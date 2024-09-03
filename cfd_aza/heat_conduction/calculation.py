import numpy as np

from cfd_aza.heat_conduction.TDMA_solver import thomas
from cfd_aza.heat_conduction.analytical_solver import analytical


class Solutions:
    def __init__(self, main_data):
        '''
          Solving the task of one-dimension stationary heat conductivity
                 with discrete analogue - TDMA.
          --------------------------------------------------------------
          Parameters:
               N: int  - quantity of control volumes
               length: float  - length of whole thing
               k: float  - coefficient of heat conductivity
               T_left: float  - left boundary condition of temperature
               T_right: float  - right boundary condition of temperature
               T_numerical: np.ndarray  - array with numerical solutions for each control volume
               T_analytical: np.ndarray  - array with analytical solutions for each control volume
               a, b, c, d: np.ndarray  - arrays with coefficients of discrete analogue
                                         for each control volume
          --------------------------------------------------------------
          Extra parameters:
               dx: float  - length of control volume
               delta: float  - extra length just for correct program working
               L: np.ndarray  - array with control volumes
          --------------------------------------------------------------

        '''

        self._N: int = main_data.N
        self._length: float = main_data.length
        self._k: float = main_data.k
        self._T_left: float = main_data.T_left
        self._T_right: float = main_data.T_right

        self._dx: float = self._length / (self._N - 1)
        self._delta: float = 0.1
        self._L: np.ndarray = np.arange(start=0, stop=self._length + self._delta, step=self._dx)

        self._T_numerical: np.ndarray = np.zeros(shape = self._N, dtype = float)
        self._T_analytical: np.ndarray = np.zeros(shape = self._N, dtype = float)

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
            self._c[i] = self._k_arr[i - 1] / self._dx
            self._b[i] = self._k_arr[i + 1] / self._dx
            self._a[i] = self._c[i] + self._b[i]


    def thomas_solution(self):
        ''' Get solution with TDMA '''
        thomas(self._a, self._b, self._c, self._d, self._N, self._T_numerical)
        return self._T_numerical


    def analytical_solution(self):
        ''' Get solution with analytical formula '''
        analytical(self._T_left, self._T_right, self._N, self._T_analytical)
        return self._T_analytical




