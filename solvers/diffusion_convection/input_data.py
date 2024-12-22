from dataclasses import dataclass, field
from typing import Callable

from scipy.constants import g

from utils.boundary_type import BoundaryConditionsModel, BoundaryConditionsType


def f_c(c: float) -> float:
    return (1.0 - c) ** 4.7


@dataclass
class GridTimeDataDiffusionConvection:
    nx: int = 10
    ny: int = 1

    x_length: float = 0.1  # m
    y_height: float = 1.0  # m

    nt: int = 10
    total_time: float = 100.0  # sec


@dataclass
class InputDataDiffusionConvection:
    sc: float = 0.0  # s = sc + sp * Cp
    sp: float = 0.0  # s = sc + sp * Cp

    # for Used calculation
    use_velocity: bool = False
    const_u_sed: float = 0.2  # 2 / 9
    g: float = g  # physical constant
    f_c: Callable = f_c  # approximation function
    r0: float = 0.001  # m
    rho1: float = 1000.0  # kg / m^3
    rho2: float = 900.0  # kg / m^3
    mu2: float = 0.6  # Pa * sec
    d: float = 9.46E-19  # diffusion coefficient, m^2 / sec

    c_init: float = 0.01  # initial concentration (t = 0)
    c_left: float = 0.07  # left boundary concentration (x = 0, t)
    c_right: float = 0.07  # right boundary concentration (x = L, t)
    c_upper: float = 0.0  # upper boundary concentration (y = H, t)
    c_lower: float = 0.0  # lower boundary concentration (y = 0, t)

    q_source: float = 0.5  # source value

    boundary_conditions_model: BoundaryConditionsModel = field(
        default_factory=BoundaryConditionsModel(
            left_side=BoundaryConditionsType.Dirichlet,
            right_side=BoundaryConditionsType.Dirichlet,
            upper_side=None,
            lower_side=None
        )
    )
