from dataclasses import dataclass
import numpy as np

# TODO: add c (specific heat) and rho (density) coefficients
# TODO: add Sp and Sc for source linearizing

@dataclass
class InputDataHeatConductivity:
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
    rho: float

@dataclass
class TimeData:
    all_time: float
    delta_time: float

@dataclass
class OutputData:
    solution_set: np.ndarray
    L: np.ndarray
    delta_time: float
    all_time: float

