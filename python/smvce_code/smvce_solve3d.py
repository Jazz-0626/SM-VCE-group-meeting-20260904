"""
SM-VCE method for solving 3-D surface deformations.
Converted from MATLAB: SMVCE_solve3D.m (the first/main version used in SMVCE_main.m)
"""
import os
import sys
import time
import numpy as np

from .lk_vec import lk_vec
from .get_design_mat import get_design_mat, get_design_mat_2D
from .smvce_vce import smvce_vce
from .smvce_smad import smvce_smad
from .utils import lonlat2sub

try:
    if os.environ.get('SMVCE_DISABLE_TORCH') == '1':
        raise ImportError('Torch disabled by SMVCE_DISABLE_TORCH')
    import torch
    import torch.nn.functional as F
except Exception:  # pragma: no cover - optional dependency at runtime
    torch = None
    F = None


def smvce_solve3d(data_dict):
    """
    Solve 3-D deformations based on SM-VCE method.

    Parameters
    ----------
    data_dict : dict - loaded DATA.mat equivalent

    Returns
    -------
    Result_smvce : dict
    """
    dem = data_dict['dem'].astype(float)
    mask = data_dict['mask'].astype(float)
    flag_smad = data_dict['flag_smad']
    flag_adpws = data_dict['flag_adpws']
    flag_interWeight = data_dict['flag_interWeight']
    fsmpara = data_dict['fsmpara']
    inc = data_dict['inc']
    azi = data_dict['azi']
    losazienu = data_dict['losazienu']
    leftorright = data_dict['leftorright']
    data = data_dict['data'].astype(float)
    flag_if_2D = data_dict['flag_if_2D']
    windowsize0 = data_dict['windowsize']
    coor = data_dict['coor']
    fault = data_dict['fault']
    flag_BOI = np.zeros(data.shape[2], dtype=int)

    row, col, data_num = data.shape

    # Create fault mask
    maskfault = np.zeros_like(mask)
    if isinstance(fault, list):
        from shapely.geometry import LineString
        y_grid, x_grid = np.mgrid[0:row, 0:col]
        for fi in fault:
            subi = lonlat2sub(fi, coor)
            try:
                line = LineString(subi[:, ::-1])  # col, row
                buf = line.buffer(windowsize0 * 2.5)
                for ri in range(row):
                    for ci in range(col):
                        from shapely.geometry import Point
                        if buf.contains(Point(ci, ri)):
                            maskfault[ri, ci] = 1
            except Exception:
                pass

    fault0 = fault

    data = data.astype(float)
    data[np.isnan(azi)] = np.nan
    data[np.isnan(data)] = 0

    # Determine unknowns
    if fsmpara == 2 and flag_if_2D == 0:
        unknum = 9
        npara_sm = 6
    elif fsmpara == 2 and flag_if_2D == 1:
        unknum = 6
        npara_sm = 4
    elif fsmpara == 3 and flag_if_2D == 0:
        unknum = 12
        npara_sm = 9
    elif fsmpara == 3 and flag_if_2D == 1:
        unknum = 8
        npara_sm = 6
    else:
        unknum = 9
        npara_sm = 6

    defo_e = np.zeros((row, col))
    defo_n = np.zeros((row, col))
    defo_u = np.zeros((row, col))
    para_sm = np.zeros((row, col, npara_sm))
    var_e = np.zeros((row, col))
    var_n = np.zeros((row, col))
    var_u = np.zeros((row, col))
    var_obs = np.zeros((row, col, data_num))
    y_grid, x_grid = np.mgrid[0:row, 0:col]
    # Convert to 1-based for compatibility with MATLAB logic
    xx = x_grid + 1
    yy = y_grid + 1

    xy_m_s = 1000
    x_m = abs(coor['post_lon']) * 108 * xy_m_s
    y_m = abs(coor['post_lat']) * 108 * xy_m_s
    sitaxy = np.zeros((row, col, data_num))
    SHPcount = np.ones((row, col, data_num))

    log_file = 'sm_vce.log'
    if os.path.exists(log_file):
        os.remove(log_file)

    Bgeo = lk_vec(azi, inc, losazienu, leftorright)

    use_cuda = bool(data_dict.get('use_cuda', True))
    cuda_rows_per_batch = max(1, int(data_dict.get('cuda_rows_per_batch', 2)))

    if _can_use_cuda_fast_path(use_cuda, fault0, flag_smad, flag_adpws, flag_interWeight, fsmpara):
        try:
            return _smvce_solve3d_cuda_fixed_window(
                data_dict=data_dict,
                Bgeo=Bgeo,
                unknum=unknum,
                npara_sm=npara_sm,
                log_file=log_file,
                cuda_rows_per_batch=cuda_rows_per_batch,
            )
        except Exception as exc:
            print(f'CUDA 快路径执行失败，回落到 CPU 路径：{exc}')

    start_time = time.time()

    for i in range(row):
        _print_progress(i + 1, row, 'SM-VCE 求解三维形变')
        for j in range(col):
            if mask[i, j] == 0:
                continue

            windowsize = windowsize0
            maxwin = 5 * windowsize
            maxdis = (maxwin / 2) * np.sqrt(x_m ** 2 + y_m ** 2) / 1000
            diswa = 1 / (1 - maxdis)
            diswb = maxdis / (maxdis - 1)

            if isinstance(fault0, list) and maskfault[i, j] > 0:
                fault = fault0
            else:
                fault = 0

            L_i0, iis, jjs, SHPcounti = _getL(data, i, j, windowsize, coor, fault, flag_smad, flag_adpws)
            SHPcount[i, j, :] = SHPcounti

            L, B, k, P = _getLBkP(L_i0, iis, jjs, Bgeo, xx, yy, x_m, y_m,
                                   fsmpara, dem, i, j, flag_interWeight, flag_BOI,
                                   diswa, diswb, flag_if_2D)

            # Check rank condition
            w2 = 5
            ii_s = max(0, i - w2)
            ii_e = min(i + w2 + 1, row)
            jj_s = max(0, j - w2)
            jj_e = min(j + w2 + 1, col)
            Bgeo_block = Bgeo[ii_s:ii_e, jj_s:jj_e, :]
            Bgeo_i = np.nanmean(np.nanmean(Bgeo_block, axis=0), axis=0).reshape(data_num, 3)

            neq0_check = np.where(k >= 0.1 * windowsize ** 2)[0]
            Bgeo_i_check = Bgeo_i[neq0_check, :]

            if flag_if_2D == 0:
                rankB = 3
            else:
                rankB = 2
                Bgeo_i_check = np.delete(Bgeo_i_check, 1, axis=1)

            if Bgeo_i_check.shape[0] == 0 or np.linalg.matrix_rank(Bgeo_i_check) < rankB:
                continue
            if np.linalg.cond(Bgeo_i_check) > 1e10:
                continue

            # VCE
            sita, P_vce, unkn, f_flag, _ = smvce_vce(L, B, P)

            sitaxy[i, j, :] = sita

            if np.sum(sita) == 0:
                with open(log_file, 'a') as flog:
                    flog.write(f'像元 row:{i + 1}, col:{j + 1} 因 VCE 不收敛未求解\n')

            # Observation variance
            p = np.zeros(data_num)
            var_obs_i = np.zeros(data_num)
            for kk in range(data_num):
                pii = np.asarray(P_vce[kk], dtype=float).ravel()
                pii_nz = pii[pii != 0]
                if len(pii_nz) == 0:
                    p[kk] = 0
                    continue
                p[kk] = np.nanmean(pii_nz)
                if np.isnan(p[kk]):
                    p[kk] = 0
                    continue
                var_obs_i[kk] = sita[kk] / (p[kk] + np.finfo(float).eps)

            var_obs[i, j, :] = np.sqrt(np.abs(var_obs_i))

            if flag_if_2D == 0:
                defo_e[i, j] = unkn[0]
                defo_n[i, j] = unkn[1]
                defo_u[i, j] = unkn[2]
                if len(unkn) > 3:
                    para_sm[i, j, :] = unkn[3:3 + npara_sm]

                # Variance of 3-D deformations
                Bgeo_var = Bgeo[i, j, :].reshape(data_num, 3)
                Bgeo_var[np.isnan(Bgeo_var)] = 0
                sita_var = sita.copy()
                sita_var[sita_var == 0] = np.nan
                d0 = np.nanmean(sita_var)
                if not np.isnan(d0):
                    try:
                        var_enu_mat = d0 * np.linalg.pinv(Bgeo_var.T @ np.diag(p) @ Bgeo_var)
                        var_e[i, j] = np.sqrt(np.abs(var_enu_mat[0, 0]))
                        var_n[i, j] = np.sqrt(np.abs(var_enu_mat[1, 1]))
                        var_u[i, j] = np.sqrt(np.abs(var_enu_mat[2, 2]))
                    except Exception:
                        pass
            else:
                defo_e[i, j] = unkn[0]
                defo_u[i, j] = unkn[1]
                if len(unkn) > 2:
                    para_sm[i, j, :] = unkn[2:2 + npara_sm]

                Bgeo_var = Bgeo[i, j, :].reshape(data_num, 3)
                Bgeo_var[np.isnan(Bgeo_var)] = 0
                sita_var = sita.copy()
                sita_var[sita_var == 0] = np.nan
                d0 = np.nanmean(sita_var)
                if not np.isnan(d0):
                    try:
                        Bgeo_var_2d = np.delete(Bgeo_var, 1, axis=1)
                        var_enu_mat = d0 * np.linalg.pinv(Bgeo_var_2d.T @ np.diag(p) @ Bgeo_var_2d)
                        var_e[i, j] = np.sqrt(np.abs(var_enu_mat[0, 0]))
                        var_u[i, j] = np.sqrt(np.abs(var_enu_mat[1, 1]))
                    except Exception:
                        pass

    print()  # newline after progress bar
    total_time = time.time() - start_time

    enu = np.stack([defo_e, defo_n, defo_u], axis=2)
    var_enu = np.stack([var_e, var_n, var_u], axis=2)

    enu[enu == 0] = np.nan
    para_sm[para_sm == 0] = np.nan
    var_enu[var_enu == 0] = np.nan
    var_obs[var_obs == 0] = np.nan
    sitaxy[sitaxy == 0] = np.nan
    SHPcount[SHPcount == 0] = np.nan

    Result_smvce = {
        'enu': enu,
        'var': {'obs': var_obs, 'enu': var_enu},
        'para_sm': para_sm,
        'sita': sitaxy,
        'coor': coor,
        'InputData': data_dict,
        'SHPcount': SHPcount,
        'total_time': total_time
    }

    print(f'SM-VCE 反演完成。总耗时：{total_time:.2f} s')

    return Result_smvce


def _can_use_cuda_fast_path(use_cuda, fault, flag_smad, flag_adpws, flag_interWeight, fsmpara):
    """CUDA path currently supports the default fixed-window configuration."""
    return (
        use_cuda
        and torch is not None
        and F is not None
        and torch.cuda.is_available()
        and not isinstance(fault, list)
        and flag_smad == 0
        and flag_adpws == 0
        and flag_interWeight == 0
        and fsmpara == 2
    )


def _smvce_solve3d_cuda_fixed_window(data_dict, Bgeo, unknum, npara_sm, log_file, cuda_rows_per_batch):
    """CUDA fast path for the common fixed-window SM-VCE configuration."""
    device = torch.device('cuda')
    data = data_dict['data'].astype(np.float32)
    mask = data_dict['mask'].astype(float)
    flag_if_2D = data_dict['flag_if_2D']
    windowsize = int(data_dict['windowsize'])
    coor = data_dict['coor']

    row, col, data_num = data.shape
    rankB = 3 if flag_if_2D == 0 else 2
    threshold_k = 0.1 * windowsize ** 2
    half_w = (windowsize - 1) // 2
    window_len = 2 * half_w + 1
    small_half_w = 5

    xy_m_s = 1000.0
    x_m = abs(coor['post_lon']) * 108 * xy_m_s
    y_m = abs(coor['post_lat']) * 108 * xy_m_s

    dx_axis = x_m * np.arange(-half_w, half_w + 1, dtype=np.float32)
    dy_axis = y_m * np.arange(-half_w, half_w + 1, dtype=np.float32)
    dx_grid, dy_grid = np.meshgrid(dx_axis, dy_axis)
    dx_flat = torch.as_tensor(dx_grid.reshape(-1), device=device, dtype=torch.float64)
    dy_flat = torch.as_tensor(dy_grid.reshape(-1), device=device, dtype=torch.float64)

    data_t = torch.as_tensor(np.moveaxis(data, 2, 0)[np.newaxis, ...], device=device)

    Bgeo_nan_to_num = np.nan_to_num(Bgeo, nan=0.0).astype(np.float32)
    Bgeo_t = torch.as_tensor(np.moveaxis(Bgeo_nan_to_num, 2, 0)[np.newaxis, ...], device=device)
    Bgeo_valid_t = torch.as_tensor(
        np.moveaxis(np.isfinite(Bgeo).astype(np.float32), 2, 0)[np.newaxis, ...],
        device=device
    )

    # Precompute the 11x11 local mean geometry used by the rank/condition screening.
    kernel_small = 2 * small_half_w + 1
    Bgeo_mean_num = F.avg_pool2d(Bgeo_t, kernel_small, stride=1, padding=small_half_w) * (kernel_small ** 2)
    Bgeo_mean_den = F.avg_pool2d(Bgeo_valid_t, kernel_small, stride=1, padding=small_half_w) * (kernel_small ** 2)
    Bgeo_mean_t = Bgeo_mean_num / torch.clamp(Bgeo_mean_den, min=1.0)

    defo_e = np.zeros((row, col))
    defo_n = np.zeros((row, col))
    defo_u = np.zeros((row, col))
    para_sm = np.zeros((row, col, npara_sm))
    var_e = np.zeros((row, col))
    var_n = np.zeros((row, col))
    var_u = np.zeros((row, col))
    var_obs = np.zeros((row, col, data_num))
    sitaxy = np.zeros((row, col, data_num))
    SHPcount = np.zeros((row, col, data_num))
    log_lines = []

    start_time = time.time()

    with torch.no_grad():
        for row_start in range(0, row, cuda_rows_per_batch):
            row_end = min(row, row_start + cuda_rows_per_batch)
            block = _solve_cuda_block(
                row_start=row_start,
                row_end=row_end,
                col=col,
                data_num=data_num,
                unknum=unknum,
                npara_sm=npara_sm,
                rankB=rankB,
                threshold_k=threshold_k,
                flag_if_2D=flag_if_2D,
                data_t=data_t,
                mask_block=mask[row_start:row_end, :],
                Bgeo_t=Bgeo_t,
                Bgeo_valid_t=Bgeo_valid_t,
                Bgeo_mean_t=Bgeo_mean_t,
                dx_flat=dx_flat,
                dy_flat=dy_flat,
                half_w=half_w,
                window_len=window_len,
                device=device,
            )

            defo_e[row_start:row_end, :] = block['defo_e']
            defo_n[row_start:row_end, :] = block['defo_n']
            defo_u[row_start:row_end, :] = block['defo_u']
            para_sm[row_start:row_end, :, :] = block['para_sm']
            var_e[row_start:row_end, :] = block['var_e']
            var_n[row_start:row_end, :] = block['var_n']
            var_u[row_start:row_end, :] = block['var_u']
            var_obs[row_start:row_end, :, :] = block['var_obs']
            sitaxy[row_start:row_end, :, :] = block['sita']
            SHPcount[row_start:row_end, :, :] = block['SHPcount']
            log_lines.extend(block['log_lines'])
            _print_progress(row_end, row, 'SM-VCE 求解三维形变 (CUDA)')

    if log_lines:
        with open(log_file, 'a') as flog:
            flog.writelines(log_lines)

    torch.cuda.synchronize()
    print()
    total_time = time.time() - start_time

    enu = np.stack([defo_e, defo_n, defo_u], axis=2)
    var_enu = np.stack([var_e, var_n, var_u], axis=2)

    enu[enu == 0] = np.nan
    para_sm[para_sm == 0] = np.nan
    var_enu[var_enu == 0] = np.nan
    var_obs[var_obs == 0] = np.nan
    sitaxy[sitaxy == 0] = np.nan
    SHPcount[SHPcount == 0] = np.nan

    Result_smvce = {
        'enu': enu,
        'var': {'obs': var_obs, 'enu': var_enu},
        'para_sm': para_sm,
        'sita': sitaxy,
        'coor': coor,
        'InputData': data_dict,
        'SHPcount': SHPcount,
        'total_time': total_time
    }

    print(f'SM-VCE 反演（CUDA 快路径）完成。总耗时：{total_time:.2f} s')
    return Result_smvce


def _solve_cuda_block(row_start, row_end, col, data_num, unknum, npara_sm, rankB, threshold_k,
                      flag_if_2D, data_t, mask_block, Bgeo_t, Bgeo_valid_t, Bgeo_mean_t, dx_flat, dy_flat,
                      half_w, window_len, device):
    """Solve a small row stripe on CUDA."""
    block_rows = row_end - row_start
    batch_pixels = block_rows * col

    result = {
        'defo_e': np.zeros((block_rows, col)),
        'defo_n': np.zeros((block_rows, col)),
        'defo_u': np.zeros((block_rows, col)),
        'para_sm': np.zeros((block_rows, col, npara_sm)),
        'var_e': np.zeros((block_rows, col)),
        'var_n': np.zeros((block_rows, col)),
        'var_u': np.zeros((block_rows, col)),
        'var_obs': np.zeros((block_rows, col, data_num)),
        'sita': np.zeros((block_rows, col, data_num)),
        'SHPcount': np.zeros((block_rows, col, data_num)),
        'log_lines': [],
    }

    mask_block_flat = torch.as_tensor(mask_block.reshape(-1) != 0, device=device)
    if not bool(mask_block_flat.any()):
        return result

    BtB_list = []
    BtL_list = []
    LtL_list = []
    k_list = []

    for kk in range(data_num):
        L_win = _unfold_block(data_t[:, kk:kk + 1, :, :], row_start, row_end, half_w, window_len)[:, 0, :].double()
        Bgeo_win = _unfold_block(data_t=Bgeo_t[:, kk * 3:(kk + 1) * 3, :, :],
                                 row_start=row_start, row_end=row_end,
                                 half_w=half_w, window_len=window_len).double()
        Bgeo_valid_win = _unfold_block(data_t=Bgeo_valid_t[:, kk * 3:(kk + 1) * 3, :, :],
                                       row_start=row_start, row_end=row_end,
                                       half_w=half_w, window_len=window_len)

        a = Bgeo_win[:, 0, :]
        b = Bgeo_win[:, 1, :]
        c = Bgeo_win[:, 2, :]
        B_i = _build_design_mat_torch(a, b, c, dx_flat, dy_flat, flag_if_2D)

        valid_geo = torch.all(Bgeo_valid_win > 0.5, dim=1)
        valid_rows = (L_win != 0) & valid_geo & torch.isfinite(B_i).all(dim=2)
        valid_rows_f = valid_rows.to(torch.float64)
        B_valid = B_i * valid_rows_f.unsqueeze(-1)
        L_valid = L_win * valid_rows_f

        BtB_list.append(torch.matmul(B_valid.transpose(1, 2), B_valid))
        BtL_list.append(torch.matmul(B_valid.transpose(1, 2), L_valid.unsqueeze(-1)).squeeze(-1))
        LtL_list.append(torch.sum(L_valid * L_valid, dim=1))
        k_list.append(valid_rows.sum(dim=1).double())

    BtB = torch.stack(BtB_list, dim=1)
    BtL = torch.stack(BtL_list, dim=1)
    LtL = torch.stack(LtL_list, dim=1)
    k_stack = torch.stack(k_list, dim=1)

    result['SHPcount'] = k_stack.cpu().numpy().reshape(block_rows, col, data_num)

    Bgeo_mean_block = Bgeo_mean_t[0, :, row_start:row_end, :].permute(1, 2, 0).reshape(batch_pixels, data_num, 3)
    strong_geom = Bgeo_mean_block.clone()
    strong_geom[k_stack < threshold_k] = 0
    if flag_if_2D == 1:
        strong_geom = torch.cat([strong_geom[:, :, :1], strong_geom[:, :, 2:]], dim=2)

    sv = torch.linalg.svdvals(strong_geom.double())
    rank = (sv > 1e-8).sum(dim=1)
    cond = torch.full((batch_pixels,), float('inf'), dtype=torch.float64, device=device)
    if sv.shape[1] > 0:
        cond_ok = sv[:, -1] > 1e-12
        cond[cond_ok] = sv[cond_ok, 0] / sv[cond_ok, -1]

    active_mask = mask_block_flat & (rank >= rankB) & torch.isfinite(cond) & (cond <= 1e10)
    active_idx = torch.where(active_mask)[0]
    if active_idx.numel() == 0:
        return result

    BtB_active = BtB[active_idx]
    BtL_active = BtL[active_idx]
    LtL_active = LtL[active_idx]
    k_active = k_stack[active_idx]

    sita, p_vec, x = _batched_vce_torch(BtB_active, BtL_active, LtL_active, k_active)

    sita_sum_zero = torch.where(torch.abs(torch.sum(sita, dim=1)) < 1e-12)[0]
    for local_idx in sita_sum_zero.tolist():
        pixel_idx = int(active_idx[local_idx].item())
        rr = row_start + pixel_idx // col
        cc = pixel_idx % col
        result['log_lines'].append(f'像元 row:{rr + 1}, col:{cc + 1} 因 VCE 不收敛未求解\n')

    x_cpu = x.cpu().numpy()
    sita_cpu = sita.cpu().numpy()
    p_cpu = p_vec.cpu().numpy()
    active_np = active_idx.cpu().numpy()

    defo_e = np.zeros(batch_pixels)
    defo_n = np.zeros(batch_pixels)
    defo_u = np.zeros(batch_pixels)
    para_sm = np.zeros((batch_pixels, npara_sm))
    var_e = np.zeros(batch_pixels)
    var_n = np.zeros(batch_pixels)
    var_u = np.zeros(batch_pixels)
    var_obs = np.zeros((batch_pixels, data_num))
    sita_out = np.zeros((batch_pixels, data_num))

    sita_out[active_np, :] = sita_cpu
    var_obs_active = np.zeros_like(sita_cpu)
    valid_p = p_cpu != 0
    var_obs_active[valid_p] = sita_cpu[valid_p] / (p_cpu[valid_p] + np.finfo(float).eps)
    var_obs[active_np, :] = np.sqrt(np.abs(var_obs_active))

    if flag_if_2D == 0:
        defo_e[active_np] = x_cpu[:, 0]
        defo_n[active_np] = x_cpu[:, 1]
        defo_u[active_np] = x_cpu[:, 2]
        if npara_sm > 0:
            para_sm[active_np, :] = x_cpu[:, 3:3 + npara_sm]
    else:
        defo_e[active_np] = x_cpu[:, 0]
        defo_u[active_np] = x_cpu[:, 1]
        if npara_sm > 0:
            para_sm[active_np, :] = x_cpu[:, 2:2 + npara_sm]

    Bgeo_center = Bgeo_t[0, :, row_start:row_end, :].permute(1, 2, 0).reshape(batch_pixels, data_num, 3)[active_idx].double()
    sita_active = sita
    sita_mask = sita_active != 0
    sita_count = torch.clamp(sita_mask.sum(dim=1), min=1)
    d0 = (sita_active * sita_mask).sum(dim=1) / sita_count
    d0 = d0.unsqueeze(-1).unsqueeze(-1)

    if flag_if_2D == 0:
        N_var = torch.matmul(Bgeo_center.transpose(1, 2), p_vec.unsqueeze(-1) * Bgeo_center)
        var_enu_mat = d0 * torch.linalg.pinv(N_var)
        diag = torch.diagonal(var_enu_mat, dim1=-2, dim2=-1).cpu().numpy()
        var_e[active_np] = np.sqrt(np.abs(diag[:, 0]))
        var_n[active_np] = np.sqrt(np.abs(diag[:, 1]))
        var_u[active_np] = np.sqrt(np.abs(diag[:, 2]))
    else:
        Bgeo_center_2d = torch.cat([Bgeo_center[:, :, :1], Bgeo_center[:, :, 2:]], dim=2)
        N_var = torch.matmul(Bgeo_center_2d.transpose(1, 2), p_vec.unsqueeze(-1) * Bgeo_center_2d)
        var_enu_mat = d0 * torch.linalg.pinv(N_var)
        diag = torch.diagonal(var_enu_mat, dim1=-2, dim2=-1).cpu().numpy()
        var_e[active_np] = np.sqrt(np.abs(diag[:, 0]))
        var_u[active_np] = np.sqrt(np.abs(diag[:, 1]))

    result['defo_e'] = defo_e.reshape(block_rows, col)
    result['defo_n'] = defo_n.reshape(block_rows, col)
    result['defo_u'] = defo_u.reshape(block_rows, col)
    result['para_sm'] = para_sm.reshape(block_rows, col, npara_sm)
    result['var_e'] = var_e.reshape(block_rows, col)
    result['var_n'] = var_n.reshape(block_rows, col)
    result['var_u'] = var_u.reshape(block_rows, col)
    result['var_obs'] = var_obs.reshape(block_rows, col, data_num)
    result['sita'] = sita_out.reshape(block_rows, col, data_num)

    return result


def _unfold_block(data_t, row_start, row_end, half_w, window_len):
    """Extract sliding windows for a row stripe using unfold."""
    row0 = max(0, row_start - half_w)
    row1 = min(data_t.shape[2], row_end + half_w)
    crop = data_t[:, :, row0:row1, :]

    pad_top = max(0, half_w - row_start)
    pad_bottom = max(0, row_end + half_w - data_t.shape[2])
    crop = F.pad(crop, (half_w, half_w, pad_top, pad_bottom))
    unfolded = F.unfold(crop, kernel_size=(window_len, window_len), stride=1)
    unfolded = unfolded.transpose(1, 2)
    channels = data_t.shape[1]
    return unfolded.reshape(-1, channels, window_len * window_len)


def _build_design_mat_torch(a, b, c, dx_flat, dy_flat, flag_if_2D):
    """Build the strain-model design matrix on CUDA for one observation type."""
    if flag_if_2D == 0:
        return torch.stack([
            a, b, c,
            a * dx_flat, a * dy_flat,
            b * dx_flat, b * dy_flat,
            c * dx_flat, c * dy_flat
        ], dim=2)

    return torch.stack([
        a, c,
        a * dx_flat, a * dy_flat,
        c * dx_flat, c * dy_flat
    ], dim=2)


def _batched_vce_torch(BtB, BtL, LtL, k_stack):
    """Batched single-iteration VCE on CUDA.

    与 MATLAB SMVCE_vce.m 和 Python CPU 路径 smvce_vce.py 行为一致：
    k=0 的观测在 VCE 中被剔除；参考类自动取每行第一个 k>0 且 sita>0 的观测，
    而非固定取第 0 个。早期实现 `ref = sita0[:, :1]` 在断层附近失相干导致
    obs 0 (常为 S1 升轨 DInSAR) 在窗口内 k=0 时整像素崩为零，被上层 NaN 化。
    """
    has_data = k_stack > 0
    # 初始权重：k=0 的观测置 0，等价于 MATLAB neq0 跳过
    weights0 = has_data.to(torch.float64)
    sita0, _, _ = _batched_sita_torch(BtB, BtL, LtL, k_stack, weights0)
    sita0 = torch.abs(sita0)

    # 参考类：每行取第一个 (has_data AND sita0>0) 的观测
    eligible = has_data & (sita0 > 1e-12)
    any_eligible = eligible.any(dim=1, keepdim=True)
    # argmax 在全 False 行返回 0，靠 any_eligible 兜底
    ref_idx = torch.argmax(eligible.to(torch.int8), dim=1, keepdim=True)
    ref = sita0.gather(1, ref_idx)

    weights1 = torch.where(
        eligible,
        ref / torch.clamp(sita0, min=1e-12),
        torch.zeros_like(sita0),
    )
    # 整行均不 eligible（极少数像素）→ 退回到等权重以仍能得到 x，sita 全 0
    weights1 = torch.where(any_eligible, weights1, weights0)

    sita1, _, x1 = _batched_sita_torch(BtB, BtL, LtL, k_stack, weights1)
    return sita1, weights1, x1


def _batched_sita_torch(BtB, BtL, LtL, k_stack, weights):
    """Solve x and variance components for a batch of pixels."""
    N = torch.sum(weights.unsqueeze(-1).unsqueeze(-1) * BtB, dim=1)
    U = torch.sum(weights.unsqueeze(-1) * BtL, dim=1)
    NN = torch.linalg.pinv(N)
    x = torch.matmul(NN, U.unsqueeze(-1)).squeeze(-1)

    A = torch.matmul(NN.unsqueeze(1), BtB)
    S = torch.einsum('bdij,beji->bde', A, A)
    traceA = torch.diagonal(A, dim1=-2, dim2=-1).sum(dim=-1)
    diag_idx = torch.arange(S.shape[1], device=S.device)
    S[:, diag_idx, diag_idx] = k_stack - 2 * traceA + S[:, diag_idx, diag_idx]

    xBtBx = torch.einsum('bi,bdij,bj->bd', x, BtB, x)
    xBtL = torch.einsum('bi,bdi->bd', x, BtL)
    residual_norm = torch.clamp(xBtBx - 2 * xBtL + LtL, min=0)
    W = weights * residual_norm

    sita = torch.matmul(torch.linalg.pinv(S), W.unsqueeze(-1)).squeeze(-1)
    return sita, NN, x


def _print_progress(current, total, prefix='Progress', bar_len=42):
    """Print a progress bar to stderr."""
    frac = current / total
    filled = int(bar_len * frac)
    bar = '█' * filled + '░' * (bar_len - filled)
    sys.stderr.write(f'\r{prefix}: |{bar}| {current}/{total} ({frac * 100:.1f}%)')
    sys.stderr.flush()


def _getLBkP(L_i0, iis, jjs, Bgeo, xx, yy, x_m, y_m, fsmpara, dem, i, j,
             flag_interWeight, flag_BOI, diswa, diswb, flag_if_2D):
    """Get observation vectors, design matrices, counts, and weight matrices."""
    data_num = len(L_i0)
    L = [None] * data_num
    B = [None] * data_num
    k = np.zeros(data_num, dtype=int)
    P = [None] * data_num

    for kk in range(data_num):
        ii_s, ii_e = int(iis[kk, 0]), int(iis[kk, 1]) + 1
        jj_s, jj_e = int(jjs[kk, 0]), int(jjs[kk, 1]) + 1

        Bgeo_i = Bgeo[ii_s:ii_e, jj_s:jj_e, kk * 3:(kk + 1) * 3]
        dx_i = x_m * (xx[ii_s:ii_e, jj_s:jj_e] - (j + 1))  # 1-based
        dy_i = y_m * (yy[ii_s:ii_e, jj_s:jj_e] - (i + 1))

        if fsmpara == 3:
            dz_i = dem[ii_s:ii_e, jj_s:jj_e] - dem[i, j]
            if flag_if_2D == 0:
                B_i = get_design_mat(Bgeo_i, dx_i, dy_i, dz_i)
            else:
                B_i = get_design_mat_2D(Bgeo_i, dx_i, dy_i, dz_i)
        else:
            if flag_if_2D == 0:
                B_i = get_design_mat(Bgeo_i, dx_i, dy_i)
            else:
                B_i = get_design_mat_2D(Bgeo_i, dx_i, dy_i)

        tem = L_i0[kk].ravel()
        tem_B = B_i

        # Find valid (non-zero, non-nan) indices
        valid_B = ~np.isnan(np.sum(B_i, axis=1))
        neq0 = np.where((tem != 0) & valid_B)[0]

        L[kk] = tem[neq0]
        B[kk] = tem_B[neq0, :]
        k[kk] = len(neq0)

        if flag_interWeight == 1 and k[kk] > 0:
            wi = np.ones(k[kk])
            # Simplified: skip iterative weighting for now
        else:
            wi = np.ones(k[kk])

        P[kk] = wi

    return L, B, k, P


def _getL(data, i, j, windowsize, coor, fault, flag_smad, flag_adpws):
    """Get observation data within a window around pixel (i, j)."""
    row, col, data_num = data.shape
    maxwin = 5 * windowsize
    minnum = 0.2 * windowsize * windowsize

    half_w = (windowsize - 1) // 2
    ii_s = max(0, i - half_w)
    ii_e = min(row - 1, i + half_w)
    jj_s = max(0, j - half_w)
    jj_e = min(col - 1, j + half_w)

    corner_loni = coor['corner_lon'] + jj_s * coor['post_lon']
    corner_lati = coor['corner_lat'] + ii_s * coor['post_lat']
    r0 = int(round(np.median(range(ii_s, ii_e + 1)) - ii_s))
    c0 = int(round(np.median(range(jj_s, jj_e + 1)) - jj_s))

    iis = np.tile([ii_s, ii_e], (data_num, 1))
    jjs = np.tile([jj_s, jj_e], (data_num, 1))

    # Fault-based homogeneous point selection
    if isinstance(fault, list):
        iffault = 1
        win_h = ii_e - ii_s + 1
        win_w = jj_e - jj_s + 1
        mask_fault = _getHomoPoints(fault, [win_h, win_w],
                                     corner_loni, corner_lati,
                                     coor['post_lon'], coor['post_lat'], r0, c0)
    else:
        iffault = 0
        win_h = ii_e - ii_s + 1
        win_w = jj_e - jj_s + 1
        mask_fault = np.ones((win_h, win_w))

    mask_fault_3d = np.repeat(mask_fault[:, :, np.newaxis], data_num, axis=2)
    L_i = data[ii_s:ii_e + 1, jj_s:jj_e + 1, :] * mask_fault_3d

    if flag_smad == 1:
        fault_has_effect = iffault * (int(np.sum(mask_fault == 1) != mask_fault.size))
        _, L_i, _ = smvce_smad(L_i, r0, c0, 1, fault_has_effect)

    mi_sum = np.sum(np.sum(L_i != 0, axis=0), axis=0)
    mi_sum00 = mi_sum.copy()

    L_i0 = [L_i[:, :, di] for di in range(data_num)]

    windowsizeij = windowsize
    mask_all_ij = (mi_sum > 0).astype(float)

    if flag_adpws == 1:
        while np.sum(mi_sum[mask_all_ij == 1] > minnum) != np.sum(mask_all_ij):
            mi_sum0 = mi_sum.copy()
            windowsizeij = windowsizeij * 1.2
            if windowsizeij > maxwin:
                break

            half_w2 = int(round((windowsizeij - 1) / 2))
            ii_s2 = max(0, i - half_w2)
            ii_e2 = min(row - 1, i + half_w2)
            jj_s2 = max(0, j - half_w2)
            jj_e2 = min(col - 1, j + half_w2)

            corner_loni2 = coor['corner_lon'] + jj_s2 * coor['post_lon']
            corner_lati2 = coor['corner_lat'] + ii_s2 * coor['post_lat']
            r0_2 = int(round(np.median(range(ii_s2, ii_e2 + 1)) - ii_s2))
            c0_2 = int(round(np.median(range(jj_s2, jj_e2 + 1)) - jj_s2))

            if isinstance(fault, list):
                win_h2 = ii_e2 - ii_s2 + 1
                win_w2 = jj_e2 - jj_s2 + 1
                mask_fault2 = _getHomoPoints(fault, [win_h2, win_w2],
                                              corner_loni2, corner_lati2,
                                              coor['post_lon'], coor['post_lat'], r0_2, c0_2)
            else:
                win_h2 = ii_e2 - ii_s2 + 1
                win_w2 = jj_e2 - jj_s2 + 1
                mask_fault2 = np.ones((win_h2, win_w2))

            mask_fault_3d2 = np.repeat(mask_fault2[:, :, np.newaxis], data_num, axis=2)
            L_i2 = data[ii_s2:ii_e2 + 1, jj_s2:jj_e2 + 1, :] * mask_fault_3d2

            if flag_smad == 1:
                fault_effect2 = iffault * (int(np.sum(mask_fault2 == 1) != mask_fault2.size))
                _, L_i2, _ = smvce_smad(L_i2, r0_2, c0_2, 1, fault_effect2)

            mi_sum = np.sum(np.sum(L_i2 != 0, axis=0), axis=0)

            idxtemp1 = np.where(mask_all_ij == 1)[0]
            cond = (mi_sum0[idxtemp1] <= minnum) & (mi_sum[idxtemp1] > minnum)
            idxtemp2 = idxtemp1[cond]

            if len(idxtemp2) > 0:
                iis[idxtemp2, :] = [ii_s2, ii_e2]
                jjs[idxtemp2, :] = [jj_s2, jj_e2]
                for di in idxtemp2:
                    L_i0[di] = L_i2[:, :, di]
                mi_sum00[idxtemp2] = mi_sum[idxtemp2]

    SHPcounti = mi_sum00
    return L_i0, iis, jjs, SHPcounti


def _getHomoPoints(faulttraces, siz, corner_lon, corner_lat, post_lon, post_lat, r0, c0):
    """
    Get homogeneous points by fault separation.

    Returns a 0/1 matrix where 1 indicates pixels on the same side as the central point.
    """
    m, n = siz
    result = np.ones((m, n))

    x_grid, y_grid = np.meshgrid(np.arange(n), np.arange(m))
    lons = corner_lon + x_grid * post_lon
    lats = corner_lat + y_grid * post_lat
    lonlim = [corner_lon, corner_lon + (n - 1) * post_lon]
    latlim = [corner_lat + (n - 1) * post_lat, corner_lat]

    for faulti in faulttraces:
        for fi in range(1, len(faulti)):
            f_prev = faulti[fi - 1]
            f_curr = faulti[fi]

            # Check if segment is outside the window
            if ((f_curr[0] > lonlim[1] and f_prev[0] > lonlim[1]) or
                    (f_curr[0] < lonlim[0] and f_prev[0] < lonlim[0]) or
                    (f_curr[1] > latlim[1] and f_prev[1] > latlim[1]) or
                    (f_curr[1] < latlim[0] and f_prev[1] < latlim[0])):
                continue

            x1 = min(f_prev[0], f_curr[0])
            x2 = max(f_prev[0], f_curr[0])
            y1 = min(f_prev[1], f_curr[1])
            y2 = max(f_prev[1], f_curr[1])

            x11 = np.minimum(lons, lons[r0, c0])
            x21 = np.maximum(lons, lons[r0, c0])
            y11 = np.minimum(lats, lats[r0, c0])
            y21 = np.maximum(lats, lats[r0, c0])

            dx_fault = f_curr[0] - f_prev[0]
            if abs(dx_fault) < 1e-15:
                dx_fault = 1e-15
            k2 = (f_curr[1] - f_prev[1]) / dx_fault

            dx_grid = lons - lons[r0, c0]
            dx_grid[dx_grid == 0] = 1e-15
            k_grid = (lats - lats[r0, c0]) / dx_grid

            b2 = f_curr[1] - k2 * f_curr[0]
            b_grid = lats - k_grid * lons

            intlon = -(b_grid - b2) / (k_grid - k2 + 1e-15)
            intlon[r0, :] = (lats[r0, :] - b2) / (k2 + 1e-15)
            intlon[:, c0] = lons[r0, c0]

            intlat = -(-b2 * k_grid + b_grid * k2) / (k_grid - k2 + 1e-15)
            intlat[:, c0] = k2 * lons[r0, c0] + b2
            intlat[r0, :] = lats[r0, c0]

            temp = ((intlon > x1) & (intlon < x2) &
                    (intlat > y1) & (intlat < y2) &
                    (intlon >= x11) & (intlon <= x21) &
                    (intlat >= y11) & (intlat <= y21))
            result[temp] = 0

    result[r0, c0] = 1
    return result
