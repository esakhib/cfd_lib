from dataclasses import dataclass
import numpy as np


@dataclass
class OutputData:
    T_num = np.ndarray
    T_an = np.ndarray
    L: float
