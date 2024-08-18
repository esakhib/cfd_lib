import logging
import warnings

import numpy as np

from solvers.heat_conduction.discrete_analogue import FiniteVolumeScheme
from utils.common import timer


class HeatConductivity(FiniteVolumeScheme):
    @timer
    def __init__(self, input_data):
        """Solve 2D unsteady heat conductivity equation.

        Parameters
        ----------
        input_data :
            Input data for equation. Includes domain and  parameters.

        Notes
        ----------
        2D unsteady heat conductivity equation:

        Analytical solution:

        """

        logging.info('Start initialization grid and solver data...')

        # read and parse input/otuput data
        self._input_data = input_data
        self._grid_time_data = self._input_data.grid_time_data
        self._equation_input_data = self._input_data.equation_input_data
        self._equation_output_data = self._input_data.equation_output_data

        # parse grid and time parameters
        self._length = self._grid_time_data.x_length
        self._height = self._grid_time_data.y_height
        nx = self._grid_time_data.nx
        ny = self._grid_time_data.ny

        # parse initial and boundary conditions
        self._t_init: float = self._equation_input_data.t_init
        self._t_left: float = self._equation_input_data.t_left
        self._t_right: float = self._equation_input_data.t_right

        # calculate parameters
        dx = self._length / (nx - 1)
        dy = 1.0

        # TODO: 2D grid
        # dy = self._height / (ny - 1)

        # TODO: unsteady equation
        # total_time = self._grid_time_data.total_time
        # nt = self._grid_time_data.nt
        # dt = total_time / nt

        self._get_k_coef()

        # initialize scheme class
        super().__init__(nx=nx, ny=ny, dx=dx, dy=dy, k=self._k,
                         initial_time_value=self._t_init,
                         left_condition_value=self._t_left,
                         right_condition_value=self._t_right)

        self._equation_output_data.grid = np.arange(start=0.0, stop=self._length, step=dx)
        self._equation_output_data.grid = np.append(self._equation_output_data.grid, self._length)

        self._equation_output_data.numerical_solution = self._result

        logging.info('End initialization grid and solver data.')

    def _get_k_coef(self) -> float:
        """Get thermal diffusivity coefficient.

        Returns
        ----------
        k: float
            Thermal diffusivity.

        """

        self._k = self._equation_input_data.k

        if self._k is None:
            lambda_coef = self._equation_input_data.lambda_coef
            rho = self._equation_input_data.rho
            cp = self._equation_input_data.cp

            assert lambda_coef is not None or rho is not None or cp is not None, (
                warnings.warn('Set parameters for k calculation!'))

            self._k = cp / (lambda_coef * rho)

        return self._k

    @timer
    def solve_analytical(self):
        """Get analytical solution on grid.
        """

        logging.info('Start analytical solution...')

        self._equation_output_data.analytical_solution = (self._t_left + (self._t_right - self._t_left) /
                                                          self._length * self._equation_output_data.grid)

        logging.info('End analytical solution.')

    @timer
    def solve_numerical(self):
        """Return numerical solution from FiniteVolumeScheme.
        """

        logging.info('Start numerical solution...')

        # discrete analogue
        self.initialize_discrete_analogue()

        # solve numerical using discrete scheme
        self.solve_equation()

        logging.info('End numerical solution.')

    @property
    def output_data(self):
        return self._equation_output_data
