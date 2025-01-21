from dataclasses import dataclass
import numpy as np
from enum import Enum

# TODO: add c (specific heat) and rho (density) coefficients
# TODO: add Sp and Sc for source linearizing

@dataclass
class InputData:
    N: int
    length: float
    T_init: float
    k: float
    T_left: float
    T_right: float
    q_left: float
    q_right: float
    T_env: float
    h: float

@dataclass
class TimeData:
    all_time: float
    delta_time: float

class OutputData:
    T_solution_set: np.ndarray

