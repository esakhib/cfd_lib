import numpy as np

from cfd_aza.heat_conduction.TDMA_solver import tdma_algorithm
from cfd_aza.heat_conduction.analytical_solver import analytical_formula


class Solutions:
    def __init__(self, main_data):
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
        self._c: float = main_data.c
        self._rho: float = main_data.rho
        self._S_p: float = main_data.S_p
        self._S_c: float = main_data.S_c

        self._dx: float = self._length / (self._N - 1)
        self._delta: float = 0.1
        self._dt: float = 0.1
        self._L: np.ndarray = np.arange(start=0, stop=self._length + self._delta, step=self._dx)

        self._T_old_solution_numerical: np.ndarray = np.zeros(shape = self._N, dtype = float)

        self._T_current_solution_numerical: np.ndarray = np.zeros(shape = self._N, dtype = float)
        # self._T_current_solution_analytical: np.ndarray = np.zeros(shape = self._N, dtype = float)

        # array filled with coefficient of heat conductivity for each control volume
        self._k_arr: np.ndarray = np.array([self._k] * (self._N + 1), float)

        self._a: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._b: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._c: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._d: np.ndarray = np.zeros(shape = self._N, dtype = float)


        # boundary conditions for coefficients
        self._a[0], self._b[0], self._c[0], self._d[0] = 1, 0, 0, self._T_left
        self._a[self._N - 1], self._b[self._N - 1], self._c[self._N - 1], self._d[self._N - 1] = 1, 0, 0, self._T_right

        self._A = (self._rho * self._c * self._dx) / self._dt

        # filling arrays of coefficients with rule of discrete analogue
        for i in range(1, self._N - 1):
            self._b[i] = self._k_arr[i - 1] / self._dx
            self._c[i] = self._k_arr[i + 1] / self._dx
            self._a[i] = self._b[i] + self._c[i] + self._A - (self._S_p * self._dx)
            self._d[i] = self._S_p * self._dx + self._A * self._T_old_solution_numerical[i]


    def thomas_solution(self):
        ''' Get solution with TDMA '''
        tdma_algorithm(self._a, self._b, self._c, self._d, self._N, self._T_current_solution_numerical)
        return self._T_current_solution_numerical


    # def analytical_solution(self):
    #     ''' Get solution with analytical formula '''
    #     analytical_formula(self._T_left, self._T_right, self._N, self._T_analytical)
    #     return self._T_analytical




