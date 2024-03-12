from dataclasses import dataclass


@dataclass
class InputGridData:
    nx: int = 1000  # step dx = x_length / nx
    ny: int = 1  # step dy = y_height / ny
    x_length: float = 100.0  # m
    y_height: float = 1.0  # m


@dataclass
class InputHeatConductivityData:
    # TODO: add source linearization
    sc: float = 0.0  # s = sc + sp * Tp
    sp: float = 0.0  # s = sc + sp * Tp

    k: float | None = 100.0  # thermal diffusivity, m^2 / sec
    cp: float | None = None  # specific heat capacity, J / (kg * K)
    lambda_coef: float | None = None  # thermal conductivity, W / (m * K)
    rho: float | None = None  # density, kg / m^3

    # TODO: add class for boundary and initial conditions
    t_init: float = 100.0  # K
    t_left: float = 100.0  # K
    t_right: float = 300.0  # K
