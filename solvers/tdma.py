import numpy as np


# import numba as nb
# from numba import guvectorize, float64, prange, njit


# @guvectorize([(float64[:],) * 5], '(n), (n), (n), (n) -> (n)', nopython=False, target_backend=True)
def run_tdma(a_p: np.ndarray, a_e: np.ndarray, a_w: np.ndarray, b: np.ndarray, result: np.ndarray):
    """TDMA algorithm using numba.

    Parameters
    ----------
    a_p: np.ndarray
        Main diagonal values.
    a_e: np.ndarray
        Upper diagonal values.
    a_w: np.ndarray
        Lower diagonal values.
    b: np.ndarray
        Right-side vector.
    result: np.ndarray
        Result vector.

    """

    n = a_p.shape[0]

    p = np.zeros_like(a_p)
    q = np.zeros_like(a_p)

    p[0] = a_e[0] / a_p[0]
    q[0] = b[0] / a_p[0]

    for i in range(1, n):
        p[i] = a_e[i] / (a_p[i] - a_w[i] * p[i - 1])
        q[i] = (b[i] + a_w[i] * q[i - 1]) / (a_p[i] - a_w[i] * p[i - 1])

    result[n - 1] = q[n - 1]

    for i in range(n - 2, -1, -1):
        result[i] = p[i] * result[i + 1] + q[i]
