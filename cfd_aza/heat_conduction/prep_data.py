import numpy as np

#from cfd_aza.heat_conduction.TDMA_solver import tdma_algorithm
from cfd_aza.heat_conduction.TDMA_solver_old import tdma_algorithm
from cfd_aza.heat_conduction.init import BoundaryType


class HeatConductivityData:
    def __init__(self, input_data, boundary_condition, time_data):
        """
          Solving the task of one-dimension non-stationary heat conductivity
                 with discrete analogue - TDMA.
          --------------------------------------------------------------
          Parameters:
               N: int  - quantity of control volumes
               length: float  - length of whole thing
               k: float  - coefficient of heat conductivity [const]
               T_left: float  - left boundary condition of temperature
               T_right: float  - right boundary condition of temperature

               T_current_solution_numerical: np.ndarray  - array with numerical solutions for each control volume at This time
               T_current_solution_analytical: np.ndarray  - array with analytical solutions for each control volume at This time
               T_old_solution_numerical: np.ndarray  - array with numerical solutions for each control volume at Last time
               a, b, c, d: np.ndarray  - arrays with coefficients of discrete analogue
                                         for each control volume
          --------------------------------------------------------------
          Extra parameters:
               dx: float  - length of control volume
               dt: float  - time step interval
               delta: float  - extra length just for correct program working
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
        self._T_left: float | None = input_data.T_left
        self._T_right: float | None = input_data.T_right
        self._q: float | None = input_data.q  # если не задана одна из сторон, то приравнять к соответсвующей стороне
                                             # если заданы обе стороны, то эта переменная не используется
        self._q_left: float | None = input_data.q_left
        self._q_right: float | None = input_data.q_right
        self._T_environment: float | None = input_data.T_env
        self._h: float | None = input_data._h
        self._c: float = self._h / self._k

        self._delta: float = 0.1  # m
        self._dx: float = self._length / self._N_origin  # m
        self._L: np.ndarray = np.arange(
            start=(-self._dx / 2),
            stop=self._length + (self._dx / 2) + self._delta,
            step=self._dx)
        self._L[0], self._L[self._N - 1] = 0, self._L[self._N - 1] - (self._dx / 2)

        # данные, касающиеся времени
        self._all_time: float = time_data.all_time  # все рассматриваемое время, sec
        self._dt: float = time_data.delta_time  # sec
        self._time_steps: int = int(self._all_time / self._dt)  # количество врем промежутков
        self._a_o: float = self._k * self._dx / self._dt  # a_o = (rho * c * dx) / dt

        self._T_solution: np.ndarray = np.zeros(shape=self._N, dtype=float)
        self._T_solution = np.full_like(self._T_solution, self._T_init)

        # array filled with coefficient of heat conductivity for each control volume
        self._k_arr: np.ndarray = np.array([self._k] * (self._N + 1), float)

        self._a_p: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._a_e: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._a_w: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._b: np.ndarray = np.zeros(shape = self._N, dtype = float)
        
        
        # TODO: add Sp and Sc for source linearizing
        #linearize temperature source S = S_c + S_p * T[i]
        self._S_c = 0
        self._S_p = 0

    def apply_bndry_cond(self):
        """ Get boundary coefficients with rule of discrete analogue """
        if (self._left_side == BoundaryType.Dirichlet):
            self._a_p[0], self._a_w[0], self._a_e[0], self._b[0] = 1, 0, -1, 2 * self._T_left
        elif (self._left_side == BoundaryType.Neumann):
            self._a_p[0], self._a_w[0], self._a_e[0], self._b[0] = 1, 0, 1, self._q_left * (self._dx / self._k)
        elif (self._left_side == BoundaryType.Robin):
            self._a_p[0], self._a_w[0], self._a_e[0], self._b[0] = (-(self._c * self._dx / 2) - 1), 0, ((self._c * self._dx / 2) - 1), -self._c * self._dx * self._T_environment

        if (self._right_side == BoundaryType.Dirichlet):
            self._a_p[self._N - 1], self._a_w[self._N - 1], self._a_e[self._N - 1], self._b[self._N - 1] = 1, -1, 0, 2 * self._T_right
        elif (self._right_side == BoundaryType.Neumann):
            self._a_p[self._N - 1], self._a_w[self._N - 1], self._a_e[self._N - 1], self._b[self._N - 1] = 1, 1, 0, self._q_right * (self._dx / self._k)
        elif (self._right_side == BoundaryType.Robin):
            self._a_p[self._N - 1], self._a_w[self._N - 1], self._a_e[self._N - 1], self._b[self._N - 1] = (1 - (self._c * self._dx / 2)), (1 + (self._c * self._dx / 2)), 0, -self._c * self._dx * self._T_environment



    def thomas_solution(self):
        """ Get solution with TDMA """
        tdma_algorithm(self._a_p, self._a_e, self._a_w, self._b, self._N, self._T_solution)
        return self._T_solution


    @property
    def N_time(self):
        return self._N_time

    @property
    def dt(self):
        return self._dt

    @property
    def L(self):
        return self._L



    # @property
    # def T_old_solution_numerical(self):
    #     return self._T_old_solution_numerical
    # @T_old_solution_numerical.setter
    # def T_old_solution_numerical(self, T_old_solution_numerical):
    #     self._T_old_solution_numerical = T_old_solution_numerical


