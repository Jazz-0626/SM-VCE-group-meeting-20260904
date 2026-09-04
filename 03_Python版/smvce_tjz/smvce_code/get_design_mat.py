"""
Establish the design matrix based on Strain Model (SM).
Converted from MATLAB: get_design_mat.m, get_design_mat_2D.m
"""
import numpy as np


def get_design_mat(Bgeo, de, dn, du=None):
    """
    Establish the design matrix based on SM for 3D displacement.

    Parameters
    ----------
    Bgeo : ndarray, shape (m, n, data_num*3)
    de, dn : ndarray, shape (m, n) - coordinate increments
    du : ndarray, shape (m, n), optional - vertical coordinate increment

    Returns
    -------
    B : ndarray - design matrix
    """
    m, n, total = Bgeo.shape
    data_num = total // 3

    if du is not None:
        # 3-D SM
        ncols = 12
        B = np.zeros((m * n * data_num, ncols))
        for i in range(data_num):
            ai = Bgeo[:, :, 3 * i].ravel()
            bi = Bgeo[:, :, 3 * i + 1].ravel()
            ci = Bgeo[:, :, 3 * i + 2].ravel()
            de_f = de.ravel()
            dn_f = dn.ravel()
            du_f = du.ravel()

            s = m * n * i
            e = m * n * (i + 1)
            B[s:e, :] = np.column_stack([
                ai, bi, ci,
                ai * de_f, ai * dn_f, ai * du_f,
                bi * de_f, bi * dn_f, bi * du_f,
                ci * de_f, ci * dn_f, ci * du_f
            ])
    else:
        # 2-D SM
        ncols = 9
        B = np.zeros((m * n * data_num, ncols))
        for i in range(data_num):
            ai = Bgeo[:, :, 3 * i].ravel()
            bi = Bgeo[:, :, 3 * i + 1].ravel()
            ci = Bgeo[:, :, 3 * i + 2].ravel()
            de_f = de.ravel()
            dn_f = dn.ravel()

            s = m * n * i
            e = m * n * (i + 1)
            B[s:e, :] = np.column_stack([
                ai, bi, ci,
                ai * de_f, ai * dn_f,
                bi * de_f, bi * dn_f,
                ci * de_f, ci * dn_f
            ])

    return B


def get_design_mat_2D(Bgeo, de, dn, du=None):
    """
    Establish the design matrix based on SM for 2D (E-W and Vertical) displacement.

    Parameters
    ----------
    Bgeo : ndarray, shape (m, n, data_num*3)
    de, dn : ndarray, shape (m, n)
    du : ndarray, shape (m, n), optional

    Returns
    -------
    B : ndarray - design matrix
    """
    m, n, total = Bgeo.shape
    data_num = total // 3

    if du is not None:
        # 3-D SM for 2D displacement
        ncols = 8
        B = np.zeros((m * n * data_num, ncols))
        for i in range(data_num):
            ai = Bgeo[:, :, 3 * i].ravel()
            ci = Bgeo[:, :, 3 * i + 2].ravel()
            de_f = de.ravel()
            dn_f = dn.ravel()
            du_f = du.ravel()

            s = m * n * i
            e = m * n * (i + 1)
            B[s:e, :] = np.column_stack([
                ai, ci,
                ai * de_f, ai * dn_f, ai * du_f,
                ci * de_f, ci * dn_f, ci * du_f
            ])
    else:
        # 2-D SM for 2D displacement
        ncols = 6
        B = np.zeros((m * n * data_num, ncols))
        for i in range(data_num):
            ai = Bgeo[:, :, 3 * i].ravel()
            ci = Bgeo[:, :, 3 * i + 2].ravel()
            de_f = de.ravel()
            dn_f = dn.ravel()

            s = m * n * i
            e = m * n * (i + 1)
            B[s:e, :] = np.column_stack([
                ai, ci,
                ai * de_f, ai * dn_f,
                ci * de_f, ci * dn_f
            ])

    return B
