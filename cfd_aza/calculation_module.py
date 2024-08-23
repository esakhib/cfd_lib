import numpy as np

from thomas_func import thomas
from analytical_func import analytical


class Solutions:
    def __init__(self, main_data):
        self.N = main_data.N
        self.length = main_data.length
        self.k = main_data.k
        self.T_left = main_data.T_left
        self.T_right = main_data.T_right

        self._dx: float = self.length / (self.N - 1)  # length of control volume
        self._delta: float = 0.1  # extra length just for correct program working
        self._L: float = np.arange(start=0, stop=self.length + self._delta, step=self._dx)  # array with control volumes

        self.T_num = np.zeros(self.N)
        self.T_an = np.zeros(self.N)
        self.k_arr = np.array([self.k] * (self.N + 1), float)
        self.a = np.zeros(self.N)
        self.b = np.zeros(self.N)
        self.c = np.zeros(self.N)
        self.d = np.zeros(self.N)


        self.a[0], self.b[0], self.c[0], self.d[0] = 1, 0, 0, self.T_left
        self.a[self.N - 1], self.b[self.N - 1], self.c[self.N - 1], self.d[
            self.N - 1] = 1, 0, 0, self.T_right


        for i in range(1, self.N - 1):  #
            self.c[i] = self.k_arr[i - 1] / self._dx
            self.b[i] = self.k_arr[i + 1] / self._dx
            self.a[i] = self.c[i] + self.b[i]

    def thomas_solution(self):
        thomas(self.a, self.b, self.c, self.d, self.N, self.T_num)

        return self.T_num

    def analytical_solution(self):
        analytical(self.T_left, self.T_right, self.N, self.T_an)

        return self.T_an




