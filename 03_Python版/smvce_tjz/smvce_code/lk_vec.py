"""
Calculate coefficient between InSAR observations and 3-D displacements.
Converted from MATLAB: lk_vec.m
"""
import numpy as np


def lk_vec(azi, inc, losazi, leftright):
    """
    Calculate the coefficient between InSAR and 3-D displacements.

    Parameters
    ----------
    azi : ndarray, shape (row, col, data_num) - azimuth angles in degrees
    inc : ndarray, shape (row, col, data_num) - incidence angles in degrees
    losazi : ndarray, shape (data_num,) - observation geometry flag (1-5)
    leftright : ndarray, shape (data_num,) - right(1) or left(-1) looking

    Returns
    -------
    Bgeo : ndarray, shape (row, col, 3*data_num)
    """
    row, col, data_num = azi.shape
    Bgeo = np.zeros((row, col, 3 * data_num))

    for i in range(data_num):
        inc_i = np.deg2rad(inc[:, :, i])
        azi_i = np.deg2rad(azi[:, :, i] - 270)

        if losazi[i] == 1:  # LOS
            Bgeo[:, :, 3 * i] = -leftright[i] * np.sin(inc_i) * np.sin(azi_i)
            Bgeo[:, :, 3 * i + 1] = -leftright[i] * np.sin(inc_i) * np.cos(azi_i)
            Bgeo[:, :, 3 * i + 2] = np.cos(inc_i)
        elif losazi[i] == 2:  # AZI
            Bgeo[:, :, 3 * i] = -np.cos(azi_i)
            Bgeo[:, :, 3 * i + 1] = np.sin(azi_i)
            Bgeo[:, :, 3 * i + 2] = np.zeros((row, col))
        elif losazi[i] == 3:  # E-W
            Bgeo[:, :, 3 * i] = 1
            Bgeo[:, :, 3 * i + 1] = 0
            Bgeo[:, :, 3 * i + 2] = 0
        elif losazi[i] == 4:  # N-S
            Bgeo[:, :, 3 * i] = 0
            Bgeo[:, :, 3 * i + 1] = 1
            Bgeo[:, :, 3 * i + 2] = 0
        elif losazi[i] == 5:  # Vertical
            Bgeo[:, :, 3 * i] = 0
            Bgeo[:, :, 3 * i + 1] = 0
            Bgeo[:, :, 3 * i + 2] = 1

    if row * col * data_num == 1:
        Bgeo = Bgeo.reshape(1, 3)

    return Bgeo
