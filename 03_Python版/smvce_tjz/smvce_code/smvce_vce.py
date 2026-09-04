"""
Variance Component Estimation (VCE) for determining observation weights.
Converted from MATLAB: SMVCE_vce.m
"""
import numpy as np


def smvce_vce(L, B, P=None):
    """
    Determine the weight of observations based on Variance Component Estimation.

    Parameters
    ----------
    L : list of ndarray - observations for each type
    B : list of ndarray - design matrices
    P : list of 1-D ndarray, optional - diagonal weights for each observation type

    Returns
    -------
    sita1 : ndarray - variance of unit weight
    P_vce : list - updated weight matrices
    x : ndarray - solved unknowns
    f : int - flag (0=ok, 1=negative sita, 2=exceeded max iterations)
    iterator : int - iteration count
    """
    vce_mode = 0
    t1 = 2
    t2 = 5
    thr3 = 0.0
    data_num0 = len(L)
    k = np.zeros(data_num0, dtype=int)
    f = 0

    for i in range(data_num0):
        k[i] = len(L[i]) if L[i] is not None and len(L[i]) > 0 else 0

    if P is None:
        P = [None] * data_num0
        for i in range(data_num0):
            if k[i] > 0:
                P[i] = np.ones(k[i], dtype=float)
            else:
                P[i] = np.zeros(0, dtype=float)

    max_k = max(k) if max(k) > 0 else 1
    for i in range(data_num0):
        if k[i] / max_k < thr3:
            L[i] = np.array([])
            P[i] = np.zeros(0, dtype=float)

    P0 = [_as_weight_vector(p).copy() for p in P]
    sita, v, x = _getsita(B, P, L)

    max_k = max(k) if max(k) > 0 else 1
    k_ratio = k / max_k
    k_ratio[k_ratio < thr3] = 0
    neq0 = np.where(k_ratio > 0)[0]

    iterator = 0
    # Single iteration (matching the simplified MATLAB version)
    iterator += 1
    sita = np.abs(sita)

    _, P = _getP(sita, v, P, neq0, t1, t2, vce_mode)
    sita, v, x = _getsita(B, P, L)

    P_vce = P
    sita1 = np.zeros(data_num0)
    if len(neq0) > 0 and len(sita) == len(neq0):
        sita1[neq0] = sita

    return sita1, P_vce, x, f, iterator


def _getsita(B, P, L):
    """Compute variance components."""
    data_num0 = len(L)
    k = np.zeros(data_num0, dtype=int)
    unknum = B[0].shape[1] if B[0] is not None and len(B[0]) > 0 else 0

    if unknum == 0:
        return np.zeros(0), [], np.zeros(0)

    N = np.zeros((unknum, unknum))
    U = np.zeros(unknum)

    for i in range(data_num0):
        k[i] = len(L[i]) if L[i] is not None and len(L[i]) > 0 else 0

    neq0 = np.where(k > 0)[0]
    data_num = len(neq0)
    k_nonzero = k[k > 0]

    Ni = [None] * data_num
    v = [None] * data_num

    for idx, i in enumerate(neq0):
        Bi = B[i]
        Li = L[i].ravel()
        Pi_vec = _as_weight_vector(P[i])
        Pi_vec[np.isnan(Pi_vec)] = 0

        weighted_B = Pi_vec[:, np.newaxis] * Bi
        Ni[idx] = Bi.T @ weighted_B
        N += Ni[idx]
        U += Bi.T @ (Pi_vec * Li)

    N[np.abs(N) == np.inf] = np.nan
    if np.any(np.isnan(N)):
        return np.zeros(len(neq0)), [np.zeros(0)] * data_num, np.zeros(unknum)

    NN = np.linalg.pinv(N)
    x = NN @ U

    W = np.zeros(data_num)
    for idx, i in enumerate(neq0):
        v[idx] = B[i] @ x - L[i].ravel()
        Pi_vec = _as_weight_vector(P[i])
        W[idx] = np.dot(Pi_vec * v[idx], v[idx])

    S = np.zeros((data_num, data_num))
    for i in range(data_num):
        for j in range(data_num):
            if i == j:
                NNNi = NN @ Ni[i]
                S[i, j] = k_nonzero[i] - 2 * np.trace(NNNi) + np.trace(NNNi @ NNNi)
            else:
                S[i, j] = np.trace(NN @ Ni[i] @ (NN @ Ni[j]))

    if S.size > 0:
        sita = np.linalg.pinv(S) @ W
    else:
        sita = np.zeros(0)

    return sita, v, x


def _getP(sita, v, P, neq0, t1, t2, vce_mode):
    """Update weight matrices based on VCE results."""
    P0 = [_as_weight_vector(p).copy() for p in P]

    for idx, i in enumerate(neq0):
        if idx < len(sita) and sita[idx] != 0:
            scale = sita[0] / sita[idx]
            P[i] = scale * P0[i]

            if vce_mode == 1:
                P_vec = _as_weight_vector(P[i])
                P_vec0 = _as_weight_vector(P0[i])
                P_vec = _robust_vce(P_vec, P_vec0, v[idx], sita[idx], t1, t2)
                P[i] = P_vec

    return P0, P


def _robust_vce(P_vec, P_vec0, v, sita0, t1, t2):
    """Robust VCE weight adjustment."""
    Pi_tem = P_vec.copy()
    P_vec0_nz = P_vec0[P_vec0 != 0]
    if len(P_vec0_nz) == 0:
        return Pi_tem

    d_tem = (np.median(np.abs(v)) / 0.6745) ** 2
    b1 = np.abs(v / np.sqrt(d_tem + np.finfo(float).eps))

    Pi_tem[np.abs(b1) >= t2] = 0

    mask = (b1 < t2) & (b1 >= t1)
    t1t21 = np.where(mask)[0]
    if len(t1t21) > 0:
        Pi_tem[t1t21] = t1 * Pi_tem[t1t21] / b1[t1t21] * ((t2 - np.abs(b1[t1t21])) / (t2 - t1)) ** 2

    return Pi_tem


def _as_weight_vector(Pi):
    """Convert a diagonal weight matrix or vector into a dense 1-D weight vector."""
    if Pi is None:
        return np.zeros(0, dtype=float)
    arr = np.asarray(Pi)
    if arr.ndim == 1:
        return arr.astype(float, copy=False).ravel()
    if hasattr(Pi, 'diagonal'):
        return np.asarray(Pi.diagonal(), dtype=float).ravel()
    return arr.astype(float, copy=False).ravel()
