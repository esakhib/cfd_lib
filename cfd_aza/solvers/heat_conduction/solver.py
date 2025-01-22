import numpy as np

from cfd_aza.solvers.TDMA import tdma_algorithm
from cfd_aza.solvers.boundary_type import BoundaryType


class HeatConductivity:
    def __init__(self,
                 input_data,
                 boundary_condition,
                 time_data):

        """
          Solving the task of one-dimension unsteady heat conductivity
                 with discrete analogue - TDMA.
          --------------------------------------------------------------
          Parameters:
               N: int  - quantity of control volumes
               length: float  - length of whole thing, [m]
               k: float  - coefficient of heat conductivity [const], [m^2 / sec]
               T_left: float  - left boundary condition of temperature, [K]
               T_right: float  - right boundary condition of temperature, [K]
               q_left: float  - heat flux from left side, [W / m^2]
               q_right: float  - heat flux from right side, [W / m^2]
               T_env: float  - temperature of environment, [K]
               h: float  - heat-transfer coefficient, [W / (m^2 * K)]

               T_current_solution: np.ndarray  - array with numerical solutions for each control volume at This time
               T_solution: np.ndarray  - array with numerical solutions for each control volume at Last time
               T_solution_set: np.ndarray  - multidimensional array with all numerical solutions
               a, b, c, d: np.ndarray  - arrays with coefficients of discrete analogue
                                         for each control volume
          --------------------------------------------------------------
          Extra parameters:
               dt: float  - time step interval, [sec]
               all_time: float  - all considering time, [sec]
               time_steps: int  - quantity of time steps

               dx: float  - length of control volume, [m]
               delta: float  - extra length just for correct program working, [m]
               L: np.ndarray  - array with control volumes
          --------------------------------------------------------------

        """

        self._left_side = boundary_condition.left_side
        self._right_side = boundary_condition.right_side

        self._N_origin: int = input_data.N
        self._N: int = self._N_origin + 2
        self._length: float = input_data.length
        self._k: float = input_data.k
        self._T_init: float = input_data.T_init
        self._T_left: float = input_data.T_left
        self._T_right: float = input_data.T_right
        self._q_left: float = input_data.q_left
        self._q_right: float = input_data.q_right
        self._T_environment: float = input_data.T_env
        self._h: float = input_data.h

        # TODO: add Sp and Sc for source linearizing
        # linearize temperature source S = S_c + S_p * T[i]
        self._S_c = 0
        self._S_p = 0

        # TODO: add k coefficient calculation when it is not [const]
        # array filled with coefficient of heat conductivity for each control volume
        self._k_arr: np.ndarray = np.array([self._k] * (self._N + 1), float)

        self._a_p: np.ndarray = np.zeros(shape=self._N, dtype=float)
        self._a_e: np.ndarray = np.zeros(shape=self._N, dtype=float)
        self._a_w: np.ndarray = np.zeros(shape=self._N, dtype=float)
        self._b: np.ndarray = np.zeros(shape=self._N, dtype=float)

        self._T_solution: np.ndarray = np.zeros(shape=self._N, dtype=float)
        self._T_solution = np.full_like(self._T_solution, self._T_init)

        self._T_solution_set: np.array = np.array([], dtype=float)

        # length calculations
        self._delta: float = 0.1
        self._dx: float = self._length / self._N_origin
        self._L: np.ndarray = np.arange(
            start=(-self._dx / 2),
            stop=self._length + (self._dx / 2) + self._delta,
            step=self._dx)
        self._L[0], self._L[self._N - 1] = 0, self._L[self._N - 1] - (self._dx / 2)

        # time calculations
        self._all_time: float = time_data.all_time
        self._dt: float = time_data.delta_time
        self._time_steps: int = int(self._all_time / self._dt)

        # extra variables for correct program working
        self._a_o: float = self._k * self._dx / self._dt  # a_o = (rho * c * dx) / dt
        self._c: float = self._h / self._k

    def apply_bndry_cond(self):

        """ Get coefficients with rule of discrete analogue """

        # calculation boundary coefficients on left side
        if (self._left_side == BoundaryType.Dirichlet):
            self._a_p[0] = 1
            self._a_w[0] = 0
            self._a_e[0] = -1
            self._b[0] = 2 * self._T_left
        elif (self._left_side == BoundaryType.Neumann):
            self._a_p[0] = 1
            self._a_w[0] = 0
            self._a_e[0] = 1
            self._b[0] = self._q_left * (self._dx / self._k)
        elif (self._left_side == BoundaryType.Robin):
            self._a_p[0] = -(self._c * self._dx / 2) - 1
            self._a_w[0] = 0
            self._a_e[0] = (self._c * self._dx / 2) - 1
            self._b[0] = -self._c * self._dx * self._T_environment

        # calculation boundary coefficients on right side
        if (self._right_side == BoundaryType.Dirichlet):
            self._a_p[self._N - 1] = 1
            self._a_w[self._N - 1] = -1
            self._a_e[self._N - 1] = 0
            self._b[self._N - 1] = 2 * self._T_right
        elif (self._right_side == BoundaryType.Neumann):
            self._a_p[self._N - 1] = 1
            self._a_w[self._N - 1] = 1
            self._a_e[self._N - 1] = 0
            self._b[self._N - 1] = self._q_right * (self._dx / self._k)
        elif (self._right_side == BoundaryType.Robin):
            self._a_p[self._N - 1] = 1 - (self._c * self._dx / 2)
            self._a_w[self._N - 1] = 1 + (self._c * self._dx / 2)
            self._a_e[self._N - 1] = 0
            self._b[self._N - 1] = -self._c * self._dx * self._T_environment

        # calculation the rest coefficients with rule of discrete analogue
        for i in range(1, self._N - 1):
            self._a_w[i] = self._k / self._dx
            self._a_e[i] = self._k / self._dx
            self._a_p[i] = self._a_w[i] + self._a_e[i] + self._a_o - (self._S_p * self._dx)
            self._b[i] = self._S_c * self._dx + self._a_o * self._T_solution[i]

    def tdma_solver(self) -> np.ndarray:

        """ Get solution with TDMA """

        tdma_algorithm(self._a_p,
                       self._a_e,
                       self._a_w,
                       self._b,
                       self._N,
                       self._T_solution)

        return self._T_solution

    def time_solver(self) -> np.ndarray:

        """ Get solutions in time """

        self._time_iter: float = self._dt

        while (self._time_iter <= self._all_time):
            self.apply_bndry_cond()
            self._T_current_solution = self.tdma_solver()
            self._T_solution = np.copy(self._T_current_solution)

            # get temperature on endpoints deleting fictitious control volumes
            self._T_current_solution[0] = (self._T_current_solution[0] + self._T_current_solution[1]) / 2
            self._T_current_solution[self._N - 1] = (self._T_current_solution[self._N - 2] + self._T_current_solution[
                self._N - 1]) / 2

            # record all solutions
            self._T_solution_set = np.concatenate(
                (self._T_solution_set, self._T_current_solution))  # записываем отдельно все эти решения

            self._time_iter += self._dt

        self._T_solution_set = self._T_solution_set.reshape((self._time_steps, self._N))

        return self._T_solution_set

    @property
    def L(self) -> np.ndarray:
        return self._L
