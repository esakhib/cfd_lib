import logging

import numpy as np

from solvers.diffusion_convection.discrete_analogue import FiniteVolumeScheme
from utils.common import timer


class DiffsuionConvection(FiniteVolumeScheme):
    @timer
    def __init__(self, input_data):
        """Solve 2D unsteady diffusion convection equation.

        Parameters
        ----------
        input_data :
            Input data for equation. Includes domain and  parameters.

        Notes
        ----------
        2D unsteady  diffusion convection equation:

        Analytical solution:

        """

        logging.info('Start initialization grid and solver data...')

        # read and parse input/output data
        self._input_data = input_data
        self._grid_time_data = self._input_data.grid_time_data
        self._equation_input_data = self._input_data.equation_input_data
        self._equation_output_data = self._input_data.equation_output_data

        # parse grid and time parameters
        self._length = self._grid_time_data.x_length
        self._height = self._grid_time_data.y_height
        nx = self._grid_time_data.nx
        ny = self._grid_time_data.ny
        nt = self._grid_time_data.nt
        self._total_time = self._grid_time_data.total_time

        # parse initial and boundary conditions
        self._c_init: float = self._equation_input_data.c_init
        self._c_left: float = self._equation_input_data.c_left
        self._c_right: float = self._equation_input_data.c_right

        # list of solutions by each time iteration
        self._solutions: list = []

        # calculate parameters
        dx = self._length / (nx - 1)
        dy = 1.0
        self._dt = self._total_time / (nt - 1)

        # TODO: 2D grid
        # dy = self._height / (ny - 1)

        self._d = self._equation_input_data.d
        self._u_sed = self._calc_u_sed(c=self._c_init)

        # initialize scheme class
        super().__init__(
            nx=nx,
            ny=ny,
            nt=nt,
            dx=dx,
            dy=dy,
            dt=self._dt,
            d=self._d,
            u_sed=self._u_sed,
            c_initial_time_value=self._c_init,
            c_left_condition_value=self._c_left,
            c_right_condition_value=self._c_right
        )

        self._equation_output_data.grid = np.arange(start=0.0, stop=self._length, step=dx)
        self._equation_output_data.grid = np.append(self._equation_output_data.grid, self._length)

        self._equation_output_data.numerical_solution = self._result
        self._equation_output_data.total_solutions = self._solutions

        logging.info('End initialization grid and solver data.')

    def _calc_u_sed(self, c) -> float:
        """Calculate U_sed by C value.

        Returns
        ----------
        C: float
            Concentration coefficient.

        """

        const = self._equation_input_data.const_u_sed
        r0 = self._equation_input_data.r0
        g = self._equation_input_data.g
        rho1 = self._equation_input_data.rho1
        rho2 = self._equation_input_data.rho2
        mu2 = self._equation_input_data.mu2
        f_c = self._equation_input_data.f_c

        return const * r0 ** 2.0 * g * (rho1 - rho2) * f_c(c) / mu2

    @timer
    def solve_numerical(self):
        """Return numerical solution from FiniteVolumeScheme.
        """

        logging.info('Start numerical solution...')

        # set initial condition
        c_old_time_step = np.full_like(self._result, self._c_init)

        # set boundary condition
        c_old_time_step[0, 0] = self._c_left
        c_old_time_step[-1, -1] = self._c_right

        # list of solutions by each time iteration
        solutions = []

        # start time
        current_time = 0.0

        # loop through time layers
        while current_time <= self._total_time:
            # discrete analogue
            self.initialize_discrete_analogue(old_time_solution=c_old_time_step)

            # solve numerical using discrete scheme
            self.solve_equation()

            # next time step
            current_time += self._dt

            # update time solution
            c_old_time_step = self._result

            # save current time solution
            self._solutions.append(self._result)

        logging.info('End numerical solution.')

    @property
    def output_data(self):
        return self._equation_output_data
