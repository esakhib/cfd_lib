import numpy as np

from input_data import MainData
from output_data import OutputData

class Equations(MainData):
    def __init__(self):
        self.T_num = np.zeros(self.N)
        self.T_an = np.zeros(self.N)
        self.a = np.zeros(self.N)
        self.b = np.zeros(self.N)
        self.c = np.zeros(self.N)
        self.d = np.zeros(self.N)
        self.P = np.zeros(self.N)
        self.Q = np.zeros(self.N)

    def bound_val(self):
        self.a[0], self.b[0], self.c[0], self.d[0] = 1, 0, 0, self.T_left
        self.a[self.N - 1], self.b[self.N - 1], self.c[self.N - 1], self.d[
            self.N - 1] = 1, 0, 0, self.T_right
        
    def thomas_solution(self):
        self.P[0] = self.b[0] / self.a[0]
        self.Q[0] = self.d[0] / self.a[0]
        for i in range(1, self.N):
            self.P[i] = self.b[i] / (self.a[i] - (self.c[i] * self.P[i - 1]))
            self.Q[i] = (self.d[i] + (self.c[i] * self.Q[i - 1])) / (self.a[i] - (self.c[i] * self.P[i - 1]))

        self.T_num[self.N - 1] = self.Q[self.N - 1]
        for i in range(self.N - 1, 0, -1):
            self.T_num[i - 1] = self.P[i - 1] * self.T_num[i] + self.Q[i - 1]

    def analytical_solution(self):
        self. T_an[0] = self.T_left
        for i in range(1, self.N):
            self.T_an[i] = self.T_right - (self.T_right - self.T_left) * (self.N - i) / self.N


    def result_thomas_solution(self):
        return self.T_num

    def result_analytical_solution(self):
        return self.T_an



