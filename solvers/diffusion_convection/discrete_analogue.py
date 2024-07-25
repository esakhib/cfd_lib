import numpy as np

from solvers.tdma import run_tdma


class FiniteVolumeScheme:
    def __init__(self, nx: int, ny: int, dx: float, dy: float, d: float, c_left_condition_value: float,
                 c_right_condition_value: float, c_initial_time_value: float,
                 u_left_condition_value: float, u_right_condition_value: float, u_initial_time_value: float):
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
        d: float
            Thermal diffusion coefficient.
        c_left_condition_value: float
            Left boundary condition value for concentration.
        c_right_condition_value: float
            Right boundary condition value for concentration.
        c_initial_time_value: float
            Initial condition value for concentration.
        u_left_condition_value: float
            Left boundary condition value for X velocity.
        u_right_condition_value: float
            Right boundary condition value for X velocity.
        u_initial_time_value: float
            Initial condition value for X velocity.

        """

        self._nx: int = nx
        self._ny: int = ny

        self._dx_e: float = dx
        self._dx_w: float = dx
        self._dy_n: float = dy
        self._dy_s: float = dy

        self._d_e: float = d
        self._d_w: float = d
        self._d_s: float = d
        self._d_n: float = d

        self._c_left_condition_value = c_left_condition_value
        self._c_right_condition_value = c_right_condition_value
        self._c_initial_time_value = c_initial_time_value

        self._u_left_condition_value = u_left_condition_value
        self._u_right_condition_value = u_right_condition_value
        self._u_initial_time_value = u_initial_time_value

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

        self._b[0] = ...
        self._a_e[0] = ...
        self._a_w[0] = ...
        self._a_p[0] = ...

        for i in range(1, self._nx - 1):
            self._a_e[i] = ...
            self._a_w[i] = ...
            self._a_p[i] = ...

        self._b[self._nx - 1] = ...
        self._a_w[self._nx - 1] = ...
        self._a_e[self._nx - 1] = ...
        self._a_p[self._nx - 1] = ...

    def solve_equation(self):
        """Solve the equation by TDMA algorithm.
        """

        run_tdma(self._a_p, self._a_e, self._a_w, self._b, self._result)

    @property
    def result(self):
        return self._result
