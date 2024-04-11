import json
import logging
import warnings

import numpy as np

from solvers.heat_conduction_discrete_analogue import FiniteVolumeScheme
from utils.common import timer
from utils.output_dataclasses import OutputData


class SteadyStateHeatConductivity(FiniteVolumeScheme):
    @timer
    def __init__(self, input_json: str):
        """Solve one dimensional steady-state heat conductivity equation.

        Parameters
        ----------
        input_json : json
            Input json with grid and solver parameters.

        Notes
        ----------
        1D steady-state heat conductivity equation:
            d/dx (k * dT / dx) + S = 0, S = 0

        Analytical solution:
            T(x) = Tl + (Tr - Tl) / L * x

        """

        logging.info('Start initialization grid and solver data...')

        # read input data
        self._input_data: dict = json.loads(input_json)

        # parse domain parameters
        self._length = self._input_data['grid_data']['x_length']
        self._height = self._input_data['grid_data']['y_height']
        nx = self._input_data['grid_data']['nx']
        ny = self._input_data['grid_data']['ny']

        # parse initial and boundary conditions
        self._t_init: float = self._input_data['heat_conductivity_data']['t_init']
        self._t_left: float = self._input_data['heat_conductivity_data']['t_left']
        self._t_right: float = self._input_data['heat_conductivity_data']['t_right']

        # calculate parameters
        dx = self._length / (nx - 1)
        dy = self._height / ny
        self._get_k_coef()

        # initialize scheme class
        super().__init__(nx=nx, ny=ny, dx=dx, dy=dy, k=self._k,
                         initial_time_value=self._t_init,
                         left_condition_value=self._t_left,
                         right_condition_value=self._t_right)

        # result container
        self._output_data: OutputData = OutputData()

        self._output_data.grid = np.arange(start=0.0, stop=self._length, step=dx)
        self._output_data.grid = np.append(self._output_data.grid, self._length)

        self._output_data.numerical_solution = self._result

        logging.info('End initialization grid and solver data.')

    def _get_k_coef(self) -> float:
        """Get thermal diffusivity coefficient.

        Returns
        ----------
        k: float
            Thermal diffusivity.

        """

        self._k = self._input_data['heat_conductivity_data']['k']

        if self._k is None:
            lambda_coef = self._input_data['heat_conductivity_data']['lambda_coef']
            rho = self._input_data['heat_conductivity_data']['rho']
            cp = self._input_data['heat_conductivity_data']['cp']

            assert lambda_coef is not None or rho is not None or cp is not None, (
                warnings.warn('Set parameters for k calculation!'))

            self._k = cp / (lambda_coef * rho)

        return self._k

    @timer
    def solve_analytical(self):
        """Get analytical solution on grid.
        """

        logging.info('Start analytical solution...')

        self._output_data.analytical_solution = (self._t_left + (self._t_right - self._t_left) /
                                                 self._length * self.output_data.grid)

        logging.info('End analytical solution.')

    @timer
    def solve_numerical(self):
        """Return numerical solution from FiniteVolumeScheme.
        """

        logging.info('Start numerical solution...')

        self.initialize_discrete_analogue()

        self.solve_equation()

        logging.info('End numerical solution.')

    @property
    def output_data(self):
        return self._output_data
