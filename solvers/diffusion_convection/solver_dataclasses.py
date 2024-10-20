from dataclasses import dataclass
from enum import Enum
from typing import Callable

from scipy.constants import g


def f_c(c: float) -> float:
    return (1.0 - c) ** 4.7


@dataclass
class GridTimeData:
    nx: int = 10  # dx = x_length / nx
    ny: int = 1  # dy = y_height / ny
    x_length: float = 0.1  # m
    y_height: float = 1.0  # m
    nt: int = 10  # dt = total_time / nt
    total_time: float = 100.0


@dataclass
class InputData:
    # TODO: add source linearization
    sc: float = 0.0  # s = sc + sp * Cp
    sp: float = 0.0  # s = sc + sp * Cp

    # Used
    const_u_sed = 0.2  # 2 / 9
    g: float = g  # physical constant
    f_c: Callable = f_c  # approximation function
    r0: float = 0.001  # m
    rho1: float = 1000.0  # kg / m^3
    rho2: float = 900.0  # kg / m^3
    mu2: float = 0.6  # Pa * sec
    d: float = 9.46E-19  # diffusion coefficient, m^2 / sec

    c_init: float = 0.01  # initial concentration (t = 0)
    c_wall_left: float = 0.07  # left boundary concentration (x = 0, t)
    c_wall_right: float = 0.07  # right boundary concentration (x = L, t)
    q_source = 0.5  # source value


class BoundaryType(Enum):
    Dirichlet = 'dirichlet'  # I рода (значение на границе)
    Neumann = 'neumann'  # II рода (поток на границе)
    Robin = 'robin'  # III рода
