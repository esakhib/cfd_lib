import logging
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from solvers.diffusion_convection.input_data import InputDataDiffusionConvection, GridTimeDataDiffusionConvection
from solvers.diffusion_convection.models import DiffusionConvectionSolverModel


class EquationTypeEnum(Enum):
    HEAT_CONDUCTIVITY = 'Heat Conduction'
    DIFFUSION_CONVECTION = 'Diffusion-Convection'


@dataclass
class SolverOutputData:
    time_scale: np.ndarray = field(default_factory=lambda: np.array([]))  # timescale
    grid: np.ndarray = field(default_factory=lambda: np.array([]))  # domain grid
    numerical_solution: np.ndarray = field(default_factory=lambda: np.array([]))  # numerical solution
    analytical_solution: np.ndarray | None = field(default_factory=lambda: np.array([]))  # analytical solution
    total_solutions: list = field(default_factory=lambda: [])
    total_velocity: list = field(default_factory=lambda: [])


@dataclass
class SolverInputData:
    grid_time_data: type(GridTimeDataDiffusionConvection)  # domain parameters
    equation_input_data: type(InputDataDiffusionConvection)  # physical input data
    equation_solver: type(DiffusionConvectionSolverModel)  # solver type


def get_input_data_by_equation(equation_type: EquationTypeEnum) -> SolverInputData | None:
    """Get data for current equation.

    Parameters
    ----------
    equation_type : EquationTypeEnum
        Equation type.

    Returns
    -------
    SolverInputData | None
        Input data for current solver.

    """

    if equation_type == EquationTypeEnum.DIFFUSION_CONVECTION:
        return SolverInputData(
            grid_time_data=GridTimeDataDiffusionConvection,
            equation_input_data=InputDataDiffusionConvection,
            equation_solver=DiffusionConvectionSolverModel
        )

    elif equation_type == EquationTypeEnum.HEAT_CONDUCTIVITY:
        # return SolverInputData(
        #     grid_time_data=GridTimeDataHC,
        #     equation_input_data=InputDataHC,
        #     equation_solver=HeatConductivity
        # )

        logging.warning('This equation type is not implemented yet.')
        return None

    else:
        logging.warning('Unknown equation type! Please, check input data and try again.')
        return None
