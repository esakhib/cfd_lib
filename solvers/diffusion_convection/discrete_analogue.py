import numpy as np

from solvers.tdma import run_tdma


class FiniteVolumeScheme:
    def __init__(self,
                 nx: int,
                 ny: int,
                 nt: int,
                 dx: float,
                 dy: float,
                 dt: float,
                 d: float,
                 u_sed: float,
                 c_left_condition_value: float,
                 c_right_condition_value: float,
                 c_initial_time_value: float):
        """Finite volume method scheme by describing discrete analogue of the equation.

        Parameters
        ----------
        nx: int
            Number of grid points by X.
        ny: int
            Number of grid points by Y.
        nt: int
            Number of time points..
        dx: float
            Step size in X direction.
        dy: float
            Step size in Y direction.
        dt: float
            Step size in time.
        d: float
            Thermal diffusion coefficient.
        u_sed: float
            Sedation velocity.
        c_left_condition_value: float
            Left boundary condition value for concentration.
        c_right_condition_value: float
            Right boundary condition value for concentration.
        c_initial_time_value: float
            Initial condition value for concentration.

        """

        self._nx: int = nx
        self._ny: int = ny
        self._nt: int = nt

        self._dx_e: float = dx
        self._dx_w: float = dx
        self._dy_n: float = dy
        self._dy_s: float = dy

        self._dx: float = dx
        self._dy: float = dy
        self._dt: float = dt

        self._d_e: float = d
        self._d_w: float = d
        self._d_s: float = d
        self._d_n: float = d

        self._u_sed_e: float = u_sed
        self._u_sed_w: float = u_sed
        self._u_sed_n: float = u_sed
        self._u_sed_s: float = u_sed

        self._c_left_condition_value = c_left_condition_value
        self._c_right_condition_value = c_right_condition_value
        self._c_initial_time_value = c_initial_time_value

        self._a_p: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_e: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_w: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_n: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._a_s: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)
        self._b: np.ndarray = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)

        self._result = np.zeros(shape=(self._nx, self._ny), dtype=np.float64)

    def initialize_discrete_analogue(self, old_time_solution: np.ndarray):
        """Initialize discrete analogue by scheme.
        """

        # left boundary
        self._b[0] = self._c_left_condition_value * (2.0 * self._d_w / self._dx_w + self._u_sed_w) + old_time_solution[0]
        self._a_e[0] = self._d_e / self._dx_e + max(-self._u_sed_e, 0.0)
        self._a_w[0] = 0.0
        self._a_p[0] = (self._a_e[0] + self._a_w[0] +
                        (self._u_sed_e - self._u_sed_w) +
                        self._dx / self._dt + 2.0 * self._d_w / self._dx_w + self._u_sed_w)

        # center
        for i in range(1, self._nx - 1):
            self._a_e[i] = self._d_e / self._dx_e + max(-self._u_sed_e, 0.0)
            self._a_w[i] = self._d_w / self._dx_w + max(self._u_sed_w, 0.0)

            self._a_p[i] = (self._dx / self._dt +
                            self._d_e / self._dx_e +
                            self._d_w / self._dx_w +
                            (self._u_sed_e - self._u_sed_w))

            self._b[i] = old_time_solution[i]

        # right boundary
        self._b[self._nx - 1] = self._c_right_condition_value * (2.0 * self._d_e / self._dx_e) + old_time_solution[self._nx - 1]
        self._a_w[self._nx - 1] = self._d_w / self._dx_w + max(self._u_sed_w, 0.0)
        self._a_e[self._nx - 1] = 0.0
        self._a_p[self._nx - 1] = ((self._a_e[self._nx - 1] + self._a_w[self._nx - 1] +
                                    (self._u_sed_e - self._u_sed_w)) +
                                   self._dx / self._dt + 2.0 * self._d_e / self._dx_e)

    def solve_equation(self):
        """Solve the equation by TDMA algorithm.
        """

        run_tdma(self._a_p, self._a_e, self._a_w, self._b, self._result)

    @property
    def result(self):
        return self._result

    @property
    def result_total(self):
        return self._result
