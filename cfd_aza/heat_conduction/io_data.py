from dataclasses import dataclass
import numpy as np


@dataclass
class InputData:
    N: int
    length: float
    k: float
    T_left: float
    T_right: float

class OutputData:
    T_numerical = np.ndarray
    T_analytical = np.ndarray
    L = np.ndarray




