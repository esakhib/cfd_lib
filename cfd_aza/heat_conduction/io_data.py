from dataclasses import dataclass
import numpy as np


@dataclass
class Dirichlet: #prescribed temperature
    N: int
    length: float
    T_left: float
    T_right: float
    k: float
    #c: float    # specific heat [const]
    #rho: float    # density
    # These are for thermal source linearization S = S_c + S_p * T[i]
    # S_p: float
    # S_c: float



@dataclass
class Neumann: #specified heat flux condition
    N: int
    length: float
    k: float
    c: float
    rho: float
    T_left: float
    T_right: float
    q: float


@dataclass
class Robin: #convection boundary condition
    N: int
    length: float
    k: float
    c: float
    rho: float
    h: float
    T_env: float
    T_flat: float



class OutputData:
    T_current_solution_numerical = np.ndarray
    L = np.ndarray




