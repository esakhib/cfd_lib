from dataclasses import dataclass

from solvers.common.boundary_conditions import BoundaryCondition, BoundaryType


@dataclass
class GridTimeDataHC:
    nx: int = 1000  # step dx = x_length / nx
    ny: int = 1  # step dy = y_height / ny

    x_length: float = 100.0  # m
    y_height: float = 1.0  # m

    nt: int = 1000
    total_time = 100.0


@dataclass
class InputDataHC:
    sc: float = 0.0  # s = sc + sp * Tp
    sp: float = 0.0  # s = sc + sp * Tp

    k: float | None = 100.0  # thermal diffusivity, m^2 / sec
    cp: float | None = None  # specific heat capacity, J / (kg * K)
    lambda_coefficient: float | None = None  # thermal conductivity, W / (m * K)
    rho: float | None = None  # density, kg / m^3

    t_init: float = 100.0  # K
    t_left: float = 100.0  # K
    t_right: float = 100.0  # K
    t_upper: float = 0.0  # K
    t_lower: float = 0.0  # K

    q_source = 0.0  # source value

    boundary_model: BoundaryCondition = BoundaryCondition(
        left_side=BoundaryType.Dirichlet,
        right_side=BoundaryType.Dirichlet,
        upper_side=None,
        lower_side=None
    )
