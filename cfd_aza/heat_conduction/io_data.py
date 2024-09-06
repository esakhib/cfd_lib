from dataclasses import dataclass
import numpy as np


@dataclass
class InputData:
    N: int
    length: float
    k: float
    T_left: float
    T_right: float

    # These are for thermal source linearization S = S_c + S_p * T[i]
    S_p: float
    S_c: float

    c: float    # specific heat [const]
    rho: float    # density

class OutputData:
    T_numerical = np.ndarray
    T_analytical = np.ndarray
    L = np.ndarray




