"""
SMAD: Strain Model dependent Adaptive Neighbourhood.
Converted from MATLAB: SMVCE_SMAD.m
"""
import numpy as np


def smvce_smad(smp, r0=None, c0=None, nom=2, iffault=0):
    """
    Select points under the stationary assumption based on SMAD.

    Parameters
    ----------
    smp : ndarray, shape (m, n, data_num)
    r0, c0 : int - central point indices (0-based)
    nom : int - L1/L2-norm (default 2)
    iffault : int - if fault trace data exists

    Returns
    -------
    mask : ndarray - same shape as smp, 1/0 inclusion mask
    smp3 : ndarray - smp * mask
    defo0 : ndarray - filtered deformation for central point
    """
    m, n, data_num = smp.shape
    remainratio = 0.05

    if r0 is None:
        r0 = m // 2
    if c0 is None:
        c0 = n // 2

    sigscale1 = 3
    sigscale2 = 6
    Lratio1 = 0.65
    Lratio2 = 1
    flag_scd = 1

    defo0 = np.zeros(data_num)
    smp0 = smp.copy()
    smp0[np.isnan(smp0)] = 0

    if np.sum(smp0 == 0) == smp.size:
        return np.zeros_like(smp), np.zeros_like(smp), defo0

    # Get directional templates
    temps = _gettemps([m, n], 32)

    mask = np.zeros_like(smp)
    smp2 = np.zeros_like(smp)
    smp3 = np.zeros_like(smp)

    # Initial point selection
    if iffault == 0:
        smp1 = _getminsum(smp, temps)
    else:
        smp1 = smp.copy()

    y_grid, x_grid = np.mgrid[0:m, 0:n]
    x_grid = x_grid - c0
    y_grid = y_grid - r0

    for datai in range(data_num):
        smpi = smp[:, :, datai]
        Length_smpi = np.sum(smpi != 0)
        smp1i = smp1[:, :, datai]

        L0 = smpi.ravel()
        L01 = smp1i.ravel()
        B0 = np.column_stack([np.ones(len(L01)), x_grid.ravel(), y_grid.ravel()])

        Lneq0 = np.where(L01 != 0)[0]
        if len(Lneq0) == 0:
            continue

        L02 = L01[Lneq0]
        B02 = B0[Lneq0, :]

        # Sort by distance from median
        L02_med = np.abs(L02 - np.median(L02))
        sort_idx = np.argsort(L02_med)
        n_select = max(1, int(round(len(L02) * Lratio1)))
        Lneq2 = sort_idx[:n_select]

        if len(Lneq2) < remainratio * Lratio1 * Length_smpi:
            continue

        L = L02[Lneq2]
        B = B02[Lneq2, :]

        unk = np.linalg.pinv(B.T @ B) @ B.T @ L
        rs1 = B @ unk - L
        sig = np.sqrt(np.mean(rs1 ** 2))

        # Calculate residuals in the whole region
        rs0 = B0 @ unk - smpi.ravel()
        maski = np.zeros((m, n))

        if np.sum(np.abs(rs0 / (sig + np.finfo(float).eps)) <= sigscale1) > remainratio * Length_smpi:
            maski_flat = np.abs(rs0 / (sig + np.finfo(float).eps)) <= sigscale1
            maski = maski_flat.reshape(m, n).astype(float)

        smp2[:, :, datai] = smpi * maski

        # Second refinement
        if flag_scd == 1 and np.sum(maski) != maski.size:
            mask_flat = maski.ravel()
            idx_mask = np.where(mask_flat == 1)[0]
            if len(idx_mask) == 0:
                defo0[datai] = unk[0]
                mask[:, :, datai] = maski
                continue

            L02_2 = L0[idx_mask]
            B02_2 = B0[idx_mask, :]

            L02_med_2 = np.abs(L02_2 - np.median(L02_2))
            sort_idx_2 = np.argsort(L02_med_2)
            n_select_2 = max(1, int(round(len(L02_2) * Lratio2)))
            Lneq2_2 = sort_idx_2[:n_select_2]

            if len(Lneq2_2) < remainratio * Lratio2 * Length_smpi:
                defo0[datai] = unk[0]
                mask[:, :, datai] = maski
                continue

            L2 = L02_2[Lneq2_2]
            B2 = B02_2[Lneq2_2, :]

            unk = np.linalg.pinv(B2.T @ B2) @ B2.T @ L2
            rs1 = B2 @ unk - L2
            sig = np.sqrt(np.mean(rs1 ** 2))

            rs0 = B0 @ unk - smpi.ravel()
            if np.sum(np.abs(rs0 / (sig + np.finfo(float).eps)) <= sigscale2) > remainratio * Length_smpi:
                maski_flat2 = np.abs(rs0 / (sig + np.finfo(float).eps)) <= sigscale2
                maski2 = maski_flat2.reshape(m, n).astype(float)
                # Union of first and second masks
                maski = np.maximum(maski, maski2)

            smp3[:, :, datai] = smpi * maski
        else:
            smp3[:, :, datai] = smp2[:, :, datai]

        defo0[datai] = unk[0]
        mask[:, :, datai] = maski

    return mask, smp3, defo0


def _gettemps(siz, ndir):
    """Generate directional templates."""
    m, n = siz
    temps = np.zeros((m, n, ndir))
    r0 = m // 2
    c0 = n // 2
    r1 = m - r0 - 1
    r2 = -r0
    c1 = n - c0 - 1
    c2 = -c0

    y_grid, x_grid = np.mgrid[r2:r1 + 1, c2:c1 + 1]

    dirs = np.arange(ndir) * (360.0 / ndir)
    for idx, d in enumerate(dirs):
        x1, y1, x2, y2 = _getintpoint(d, r1, r2, c1, c2)
        xx1 = x1 - x_grid
        yy1 = y1 - y_grid
        xx2 = x2 - x_grid
        yy2 = y2 - y_grid

        # Cross product z-component
        cro_z = xx1 * yy2 - xx2 * yy1

        tempsi = np.ones(siz)
        tempsi[cro_z < 0] = -1
        temps[:, :, idx] = tempsi

    return temps


def _getintpoint(dir_deg, r1, r2, c1, c2):
    """Get intersection points for directional template."""
    k = np.tan(np.deg2rad(dir_deg))

    if np.abs(k) > 1e10:
        k = np.sign(k) * 9999 if k != 0 else 9999

    # Right side intersection
    if (k * c1) <= r1 and (k * c1) >= r2:
        x1, y1 = c1, k * c1
    elif (k * c1) < r2:
        x1 = r2 / k if k != 0 else c1
        y1 = r2
    else:
        x1 = r1 / k if k != 0 else c1
        y1 = r1

    # Left side intersection
    if (k * c2) <= r1 and (k * c2) >= r2:
        x2, y2 = c2, k * c2
    elif (k * c2) < r2:
        x2 = r2 / k if k != 0 else c2
        y2 = r2
    else:
        x2 = r1 / k if k != 0 else c2
        y2 = r1

    if dir_deg == 90:
        x1, y1 = 0, r1
        x2, y2 = 0, r2
    if dir_deg == 270:
        x1, y1 = 0, r2
        x2, y2 = 0, r1

    # Check direction
    dir1 = np.degrees(np.arctan2(y1 - y2, x1 - x2))
    if dir1 < 0:
        dir1 += 360
    if abs(dir1 - dir_deg) > 1e-6:
        x1, y1, x2, y2 = x2, y2, x1, y1

    return x1, y1, x2, y2


def _getminsum(smp, temps):
    """Select half of the window based on directional gradient."""
    data_num = smp.shape[2]
    ndir = temps.shape[2]

    smp_c = smp.copy()
    smp_c[smp_c == 0] = np.nan

    g = np.zeros((data_num, ndir))
    for di in range(data_num):
        for ti in range(ndir):
            g[di, ti] = np.nansum(smp_c[:, :, di] * temps[:, :, ti])

    # Find dominant direction for each observation
    dir2i1 = []
    dir21i1 = []
    dir2i_complex = 0
    dir21i_complex = 0

    for di in range(data_num):
        if np.all(g[di, :] == 0):
            continue

        GC = np.argmax(np.abs(g[di, :]))
        dir1_deg = (360.0 / ndir) * GC
        dir2 = dir1_deg + 90
        if dir2 >= 360:
            dir2 -= 360
        dir21 = dir1_deg + 270
        if dir21 >= 360:
            dir21 -= 360

        if dir2i_complex == 0 and dir21i_complex == 0:
            dir2i_complex = np.exp(1j * np.deg2rad(dir2))
            dir21i_complex = np.exp(1j * np.deg2rad(dir21))
            dir2i1.append(dir2)
            dir21i1.append(dir21)
        else:
            c2 = np.exp(1j * np.deg2rad(dir2))
            c21 = np.exp(1j * np.deg2rad(dir21))
            if abs(dir2i_complex + c2) > np.sqrt(2):
                dir2i_complex += c2
                dir21i_complex += c21
                dir2i_complex /= abs(dir2i_complex)
                dir21i_complex /= abs(dir21i_complex)
                dir2i1.append(dir2)
                dir21i1.append(dir21)
            else:
                dir2i_complex += c21
                dir21i_complex += c2
                dir2i_complex /= abs(dir2i_complex)
                dir21i_complex /= abs(dir21i_complex)
                dir2i1.append(dir21)
                dir21i1.append(dir2)

    if len(dir2i1) == 0:
        return smp.copy()

    # Get median direction
    dir2 = _mediandir(np.array(dir2i1))
    if dir2 < 0:
        dir2 += 360
    dir21 = _mediandir(np.array(dir21i1))
    if dir21 < 0:
        dir21 += 360

    # Use 3x3 block median analysis to determine the correct half
    mod330 = _getmod33(smp)
    ind2 = _closeddirind(dir2)
    ind21 = _closeddirind(dir21)

    ind2r, ind2c = np.unravel_index(ind2, (3, 3))
    ind21r, ind21c = np.unravel_index(ind21, (3, 3))

    mod330_c = mod330.copy()
    mod330_c[mod330_c == 0] = np.nan
    ind22v = np.nanmedian(mod330_c[1, 1, :])
    ind2v = np.nanmedian(mod330_c[ind2r, ind2c, :])
    ind21v = np.nanmedian(mod330_c[ind21r, ind21c, :])

    if np.isnan(ind22v):
        ind22v = 0
    if np.isnan(ind2v):
        ind2v = 99999999
    if np.isnan(ind21v):
        ind21v = 99999999

    if np.sum((ind2v - ind22v) ** 2) <= np.sum((ind21v - ind22v) ** 2):
        dir11 = dir2 + 90
        if dir11 >= 360:
            dir11 -= 360
    else:
        dir11 = dir2 - 90
        if dir11 < 0:
            dir11 += 360

    GC = int(round(dir11 / 360 * ndir))
    GC = max(0, min(ndir - 1, GC))

    tempsi = temps[:, :, GC].copy()
    tempsi[tempsi == -1] = 0

    smp_out = smp.copy()
    smp_out[np.isnan(smp_out)] = 0
    smp1 = smp_out * np.repeat(tempsi[:, :, np.newaxis], data_num, axis=2)

    return smp1


def _getmod33(smp, dw=0):
    """Get 3x3 block median values."""
    smp_c = smp.copy()
    smp_c[smp_c == 0] = np.nan
    m, n, data_num = smp_c.shape
    mod33 = np.zeros((3, 3, data_num))

    for i in range(3):
        for j in range(3):
            ii_start = max(0, int(round(i / 3 * m)) - dw)
            ii_end = min(m, int(round((i + 1) / 3 * m)) + dw)
            jj_start = max(0, int(round(j / 3 * n)) - dw)
            jj_end = min(n, int(round((j + 1) / 3 * n)) + dw)

            for k in range(data_num):
                block = smp_c[ii_start:ii_end, jj_start:jj_end, k]
                mk = np.nanmedian(smp_c[:, :, k])
                mkij = np.nanmedian(block)
                if mk != 0 and not np.isnan(mk) and not np.isnan(mkij):
                    mod33[i, j, k] = mkij / mk
                else:
                    mod33[i, j, k] = 0

    return mod33


def _mediandir(dirs):
    """Compute median direction handling wraparound."""
    dirs = dirs.copy()
    maxdir = np.max(dirs)
    mindir = np.min(dirs)
    if maxdir - mindir > 180:
        dirs[maxdir - dirs > 180] += 360
    result = np.median(dirs)
    if result >= 360:
        result -= 360
    return result


def _closeddirind(dir2):
    """Find the closest 3x3 grid index for a direction."""
    mod33_temp = np.zeros((3, 3))
    # Mapping: [7,8,5,2,1,0,3,6] -> 0:45:315 (MATLAB 1-based [8,9,6,3,2,1,4,7])
    positions = [(2, 0), (2, 1), (1, 2), (0, 2), (0, 1), (0, 0), (1, 0), (2, 0)]
    angles = np.arange(0, 360, 45)
    for pos, angle in zip(positions, angles):
        mod33_temp[pos] = angle
    mod33_temp[1, 1] = 500

    if dir2 > 180:
        mod33_temp[2, 0] = 360

    flat = mod33_temp.ravel()
    diffs = np.abs(dir2 - flat)
    return np.argmin(diffs)
