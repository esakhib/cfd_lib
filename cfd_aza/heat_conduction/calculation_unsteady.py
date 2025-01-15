import numpy as np

#from cfd_aza.heat_conduction.TDMA_solver import tdma_algorithm
from cfd_aza.heat_conduction.TDMA_solver_old import tdma_algorithm
from cfd_aza.heat_conduction.analytical_solver import analytical_formula


class Solutions:
    def __init__(self, main_data, boundary_condition):
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

        self._N: int = main_data.N
        self._length: float = main_data.length
        self._k: float = main_data.k
        self._T_init: float = main_data.T_init
        self._T_left: float | None = main_data.T_left
        self._T_right: float | None = main_data.T_right
        self._q: float | None = main_data.q  # если не задана одна из сторон, то приравнять к соответсвующей стороне
                                             # если заданы обе стороны, то эта переменная не используется
        self._q_left: float | None = main_data.q_left
        self._q_right: float | None = main_data.q_right
        self._T_env: float | None = main_data.T_env
        self._h: float | None = main_data._h


        self._time: float = 10
        self._N_time: int = 5
        self._dt: float = self._time / (self._N_time - 1)
        self._delta: float = 0.1
        self._dx: float = self._length / (self._N - 1)
        self._L: np.ndarray = np.arange(start=0, stop=self._length + self._delta, step=self._dx)
        self._a_o: float = (self._k * self._dx) / self._dt  # a_o = (rho * c * dx) / Dt


        self.T_old_solution_numerical: np.ndarray = T_old_solution_numerical

        # array filled with coefficient of heat conductivity for each control volume
        self._k_arr: np.ndarray = np.array([self._k] * (self._N + 1), float)

        self._a_p: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._a_e: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._a_w: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._b: np.ndarray = np.zeros(shape = self._N, dtype = float)

        #linearize temperature source S = S_c + S_p * T[i]
        self._S_c = 0
        self._S_p = 0


        # boundary conditions for coefficients
        self._a_p[0], self._a_e[0], self._a_w[0], self._b[0] = 1, 0, 0, self._T_left #self.T_old_solution_numerical[0]
        self._a_p[self._N - 1], self._a_e[self._N - 1], self._a_w[self._N - 1], self._b[self._N - 1] = 1, 0, 0, self._T_right #self.T_old_solution_numerical[self._N - 1]




        # filling arrays of coefficients with rule of discrete analogue
        for i in range(1, self._N - 1):
            self._a_e[i] = self._k_arr[i - 1] / self._dx
            self._a_w[i] = self._k_arr[i + 1] / self._dx
            self._a_p[i] = self._a_e[i] + self._a_w[i] + self._a_o - (self._S_p * self._dx)
            self._b[i] = self._S_c * self._dx + self._a_o * self.T_old_solution_numerical[i]




    def thomas_solution(self):
        ''' Get solution with TDMA '''
        tdma_algorithm(self._a_p, self._a_e, self._a_w, self._b, self._N, self.T_old_solution_numerical)
        return self.T_old_solution_numerical


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


