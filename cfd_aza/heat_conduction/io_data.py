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
    T_left: float | None
    T_right: float | None
    q: float | None
    q_left: float | None
    q_right: float | None
    T_env: float | None
    h: float | None


class OutputData:
    T_current_solution_numerical: np.ndarray
    L: np.ndarray


@dataclass
class Dirichlet: #prescribed temperature
    N: int
    length: float
    T_left: float
    T_right: float
    T_init: float
    k: float



@dataclass
class Neumann: #specified heat flux condition
    N: int
    length: float
    T_init: float
    k: float
    q: float | None
    q_left: float | None
    q_right: float | None


@dataclass
class Robin: #convection boundary condition
    N: int
    length: float
    T_env: float
    T_init: float
    k: float
    h: float









