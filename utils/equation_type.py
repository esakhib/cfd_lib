import logging
from dataclasses import dataclass, field
from enum import Enum
from typing import Type

import numpy as np

from solvers.diffusion_convection.solver_dataclasses import InputData as InputDataDC, GridTimeData as GridTimeDataDC
from solvers.heat_conduction.solver import HeatConductivity
from solvers.diffusion_convection.solver import DiffsuionConvection
from solvers.heat_conduction.solver_dataclasses import InputData as InputDataHC, GridTimeData as GridTimeDataHC


class EquationTypeEnum(Enum):
    HEAT_CONDUCTIVITY = 'heat_conductivity'  # 2D unsteady
    DIFFUSION_CONVECTION = 'diffusion_convection'  # 2D unsteady


@dataclass
class EquationData:
    grid_time_data: dataclass  # domain parameters
    equation_input_data: dataclass  # physical input data
    equation_output_data: dataclass  # solution output data
    equation_solver: Type[HeatConductivity | DiffsuionConvection]  # solver type


@dataclass
class OutputData:
    time_grid: np.ndarray = field(default_factory=lambda: np.array([]))  # output time grid
    grid: np.ndarray = field(default_factory=lambda: np.array([]))  # output domain grid
    numerical_solution: np.ndarray = field(default_factory=lambda: np.array([]))  # output numerical solution
    analytical_solution: np.ndarray | None = field(default_factory=lambda: np.array([]))  # output analytical solution


def get_input_data_by_equation(equation_type: EquationTypeEnum) -> EquationData | None:
    """

    Parameters
    ----------
    equation_type

    Returns
    -------

    """

    if equation_type == EquationTypeEnum.HEAT_CONDUCTIVITY:
        return EquationData(
            grid_time_data=GridTimeDataHC,
            equation_input_data=InputDataHC,
            equation_output_data=OutputData,
            equation_solver=HeatConductivity
        )

    elif equation_type == EquationTypeEnum.DIFFUSION_CONVECTION:
        return EquationData(
            grid_time_data=GridTimeDataDC,
            equation_input_data=InputDataDC,
            equation_output_data=OutputData,
            equation_solver=DiffsuionConvection
        )

    else:
        logging.warning('Unknown equation type! Please, check input data and try again.')
        return None
