import numpy as np


# import numba as nb
# from numba import guvectorize, float64, prange, njit


# @guvectorize([(float64[:],) * 5], '(n), (n), (n), (n) -> (n)', nopython=False, target_backend=True)
def run_tdma(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray):
    """TDMA algorithm using numba.

    Parameters
    ----------
    a: np.ndarray
        Main diagonal values.
    b: np.ndarray
        Upper diagonal values.
    c: np.ndarray
        Lower diagonal values.
    d: np.ndarray
        Right-side vector.
    result: np.ndarray
        Result vector.

    """

    result = np.zeros_like(a)
    n = a.shape[0]

    p = np.zeros_like(a)
    q = np.zeros_like(a)

    p[0] = b[0] / a[0]
    q[0] = d[0] / a[0]

    for i in range(1, n):
        tmp = a[i] - c[i] * p[i - 1]
        p[i] = b[i] / tmp[0]
        q[i] = (d[i] + c[i] * q[i - 1]) / tmp[0]

    result[n - 1] = q[n - 1]
    p[n - 1] = 0.0

    for i in range(n - 2, -1, -1):
        result[i] = p[i] * result[i + 1] + q[i]

    return result
