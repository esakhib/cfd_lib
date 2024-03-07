import numpy as np
from numba import prange

from utils.tdma import run_tdma


class FiniteVolumeScheme:
    def __init__(self, nx: int, ny: int, dx: float, dy: float, k: float, left_condition_value: float,
                 right_condition_value: float, initial_time_value: float | None):
        """Finite volume method scheme by describing discrete analogue of the equation.

        Parameters
        ----------
        nx: int
            Number of grid points by X.
        ny: int
            Number of grid points by Y.
        dx: float
            Step size in X direction.
        dy: float
            Step size in Y direction.
        k: float
            Thermal diffusion coefficient.
        left_condition_value: float
            Left boundary condition value.
        right_condition_value: float
            Right boundary condition value.
        initial_time_value: float
            Initial condition value.

        """

        self._nx: int = nx
        self._ny: int = ny

        self._dx_e: float = dx
        self._dx_w: float = dx
        self._dy_n: float = dy
        self._dy_s: float = dy

        self._k_e: float = k
        self._k_w: float = k
        self._k_s: float = k
        self._k_n: float = k

        self._left_condition_value = left_condition_value
        self._right_condition_value = right_condition_value
        self._initial_time_value = initial_time_value

        self._a_p: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_e: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_w: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_n: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_s: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._b: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)

        self._result = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)

    def initialize_discrete_analogue(self):
        """Initialize discrete analogue by scheme.
        """

        self._b[0] = self._k_w / (self._dx_w / 2.0) * self._left_condition_value
        self._a_e[0] = self._k_e / self._dx_e
        self._a_p[0] = self._a_e[0] + self._a_w[0] + self._k_w / (self._dx_w / 2.0)

        for i in range(1, self._nx - 1):
            self._a_e[i] = self._k_e / self._dx_e
            self._a_w[i] = self._k_w / self._dx_w
            self._a_p[i] = self._a_w[i] + self._a_e[i]

        self._b[self._nx - 1] = self._k_e / (self._dx_e / 2.0) * self._right_condition_value
        self._a_w[self._nx - 1] = self._k_w / self._dx_w
        self._a_p[self._nx - 1] = self._a_e[self._nx - 1] + self._a_w[self._nx - 1] + self._k_e / (self._dx_e / 2.0)

    def solve_equation(self):
        """Solve the equation by TDMA algorithm.
        """

        run_tdma(self._a_p, self._a_e, self._a_w, self._b, self._result)

        # TODO: something wrong with boundary conditions
        self._result[0] = self._left_condition_value

    @property
    def result(self):
        return self._result
