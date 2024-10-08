import numpy as np

from cfd_aza.heat_conduction.TDMA_solver import tdma_algorithm
from cfd_aza.heat_conduction.analytical_solver import analytical_formula


class Solutions:
    def __init__(self, main_data, T_old_solution_numerical):
        '''
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

        '''

        self._N: int = main_data.N
        self._length: float = main_data.length
        self._k: float = main_data.k
        self._T_left: float = main_data.T_left
        self._T_right: float = main_data.T_right
        # self._c: float = main_data.c
        # self._rho: float = main_data.rho

        self._time: float = 20
        self._N_time: int = 5
        self._dt: float = self._time / (self._N_time - 1)
        self._delta: float = 0.1
        self._dx: float = self._length / (self._N - 1)
        self._L: np.ndarray = np.arange(start=0, stop=self._length + self._delta, step=self._dx)
        self._a_o: float = (self._k * self._dx) / self._dt  # a_o = (rho * c * dx) / Dt

        #linearize temperature S = S_c + S_p * T[i]
        self._S_c = 0
        self._S_p = 0

        self.T_old_solution_numerical: np.ndarray = T_old_solution_numerical

        # array filled with coefficient of heat conductivity for each control volume
        self._k_arr: np.ndarray = np.array([self._k] * (self._N + 1), float)

        self._a: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._b: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._c: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._d: np.ndarray = np.zeros(shape = self._N, dtype = float)



        # boundary conditions for coefficients
        self._a[0], self._b[0], self._c[0], self._d[0] = 1, 0, 0, self._T_left
        self._a[self._N - 1], self._b[self._N - 1], self._c[self._N - 1], self._d[self._N - 1] = 1, 0, 0, self._T_right

        # filling arrays of coefficients with rule of discrete analogue
        for i in range(1, self._N - 1):
            self._b[i] = self._k_arr[i - 1] / self._dx
            self._c[i] = self._k_arr[i + 1] / self._dx
            self._a[i] = self._b[i] + self._c[i] + self._a_o - (self._S_p * self._dx)
            self._d[i] = self._S_c * self._dx + self._a_o * self.T_old_solution_numerical[i]




    def thomas_solution(self):
        ''' Get solution with TDMA '''
        tdma_algorithm(self._a, self._b, self._c, self._d, self._N, self.T_old_solution_numerical)
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


