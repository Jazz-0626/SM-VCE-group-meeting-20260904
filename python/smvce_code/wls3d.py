"""
Weighted Least Squares (WLS) method for 3-D displacement solving.
Converted from MATLAB: WLS3D.m
"""
import time
import sys
import numpy as np
from .lk_vec import lk_vec


def wls3d(data_dict):
    """
    Solve 3-D deformations based on WLS method.

    Parameters
    ----------
    data_dict : dict - contains data, inc, azi, losazienu, leftorright, mask, flag_if_2D

    Returns
    -------
    Result_wls : dict with keys enu, var, obsP, total_time
    """
    start_time = time.time()

    data = data_dict['data']
    inc = data_dict['inc']
    azi = data_dict['azi']
    losazienu = data_dict['losazienu']
    leftorright = data_dict['leftorright']
    mask = data_dict['mask'].astype(float)
    flag_if_2D = data_dict['flag_if_2D']

    row, col, data_num = data.shape

    P = np.ones((row, col, data_num))
    data_work = data.copy()
    data_work[np.isnan(data_work)] = 0

    defo_e = np.zeros((row, col))
    defo_n = np.zeros((row, col))
    defo_u = np.zeros((row, col))
    var_e = np.zeros((row, col))
    var_n = np.zeros((row, col))
    var_u = np.zeros((row, col))
    var_obs = np.zeros((row, col, data_num))

    Bgeo = lk_vec(azi, inc, losazienu, leftorright)

    for i in range(row):
        # Progress bar
        _print_progress(i + 1, row, 'WLS 求解三维形变')
        for j in range(col):
            if mask[i, j] == 0:
                continue

            L_i = data_work[i, j, :]
            w2 = 30
            ii_s = max(0, i - w2)
            ii_e = min(i + w2 + 1, row)
            jj_s = max(0, j - w2)
            jj_e = min(j + w2 + 1, col)

            Bgeo_block = Bgeo[ii_s:ii_e, jj_s:jj_e, :]
            Bgeo_i = np.nanmean(np.nanmean(Bgeo_block, axis=0), axis=0).reshape(data_num, 3)

            if flag_if_2D == 1:
                Bgeo_i = np.delete(Bgeo_i, 1, axis=1)
                rankB = 2
            else:
                rankB = 3

            neq0 = np.where(L_i != 0)[0]
            if len(neq0) == 0:
                continue

            L_sel = L_i[neq0]
            Bgeo_sel = Bgeo_i[neq0, :]
            P_i = np.diag(P[i, j, neq0])

            if np.linalg.matrix_rank(Bgeo_sel) < rankB:
                continue
            cond_val = np.linalg.cond(Bgeo_sel)
            if cond_val > 1e10:
                continue

            NN = np.linalg.inv(Bgeo_sel.T @ P_i @ Bgeo_sel)
            x = NN @ Bgeo_sel.T @ P_i @ L_sel

            if flag_if_2D == 1:
                defo_e[i, j] = x[0]
                defo_u[i, j] = x[1]
                var_e[i, j] = NN[0, 0]
                var_u[i, j] = NN[1, 1]
            else:
                defo_e[i, j] = x[0]
                defo_n[i, j] = x[1]
                defo_u[i, j] = x[2]
                var_e[i, j] = NN[0, 0]
                var_n[i, j] = NN[1, 1]
                var_u[i, j] = NN[2, 2]

    print()  # newline after progress bar
    total_time = time.time() - start_time

    enu = np.stack([defo_e, defo_n, defo_u], axis=2)
    var_enu = np.stack([var_e, var_n, var_u], axis=2)

    enu[enu == 0] = np.nan
    var_obs[var_obs == 0] = np.nan
    var_enu[var_enu == 0] = np.nan

    Result_wls = {
        'enu': enu,
        'var': {'obs': var_obs, 'enu': var_enu},
        'obsP': P,
        'total_time': total_time
    }

    print(f'WLS 反演完成。总耗时：{total_time:.2f} s')
    return Result_wls


def _print_progress(current, total, prefix='Progress', bar_len=42):
    """Print a progress bar to stderr."""
    frac = current / total
    filled = int(bar_len * frac)
    bar = '█' * filled + '░' * (bar_len - filled)
    sys.stderr.write(f'\r{prefix}: |{bar}| {current}/{total} ({frac * 100:.1f}%)')
    sys.stderr.flush()
