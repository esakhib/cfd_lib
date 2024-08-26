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
        self._solutions: dict = {}

        # calculate parameters
        dx = self._length / (nx - 1)
        dy = 1.0
        self._dt = self._total_time / (nt - 1)

        # TODO: 2D grid
        # dy = self._height / (ny - 1)

        self._d = self._equation_input_data.d

        # initialize scheme class
        super().__init__(
            nx=nx,
            ny=ny,
            nt=nt,
            dx=dx,
            dy=dy,
            dt=self._dt,
            d=self._d,
            c_initial_time_value=self._c_init,
            c_left_condition_value=self._c_left,
            c_right_condition_value=self._c_right
        )

        self._equation_output_data.grid = np.arange(start=0.0, stop=self._length, step=dx)
        self._equation_output_data.grid = np.append(self._equation_output_data.grid, self._length)

        self._equation_output_data.numerical_solution = self._current_solution
        self._equation_output_data.total_solutions = self._solutions

        logging.info('End initialization grid and solver data.')

    def _calc_u_sed(self, c: np.ndarray):
        """Calculate U_sed by C values.
        """

        const = self._equation_input_data.const_u_sed
        r0 = self._equation_input_data.r0
        g = self._equation_input_data.g
        rho1 = self._equation_input_data.rho1
        rho2 = self._equation_input_data.rho2
        mu2 = self._equation_input_data.mu2
        f_c = self._equation_input_data.f_c

        for i in np.arange(1, self._nx - 1):
            self._u_sed[i] = const * r0 ** 2.0 * g * (rho1 - rho2) * f_c((c[i] + c[i + 1]) / 2.0) / mu2

    @timer
    def solve_numerical(self):
        """Return numerical solution from FiniteVolumeScheme.
        """

        logging.info('Start numerical solution...')

        # задали начальное условие
        self._old_solution = np.full_like(self._current_solution, self._c_init)

        # начальное время
        current_time = 0.0

        # цикл через временные слои
        while current_time <= self._total_time:
            logging.info(f'Solving for time = {current_time}')

            # посчитаем скорость, используя концентрацию на текущем временном слое
            self._calc_u_sed(self._old_solution)
            self.update_u_sed()

            # инициализиурем дискретный аналог, используя решение на текущем временном слое
            self.initialize_discrete_analogue(
                old_time_solution=self._old_solution
            )

            # получаем решение на следующем временном слое
            self.solve_equation()

            # обновляем решение на текущем временном слое
            self._old_solution = self._current_solution

            # сохраняем решение для временного слоя current_time в словарь
            self._solutions[current_time] = self._current_solution

            # переключились на следующий временной слой
            current_time += self._dt

        logging.info('End numerical solution.')

    @property
    def output_data(self):
        return self._equation_output_data
