from dataclasses import dataclass
import numpy as np

from input_data import MainData


@dataclass
class OutputData:
    T_num = np.zeros(MainData.N)
    T_an = np.zeros(MainData.N)
    a = np.zeros(MainData.N)
    b = np.zeros(MainData.N)
    c = np.zeros(MainData.N)
    d = np.zeros(MainData.N)
    P = np.zeros(MainData.N)
    Q = np.zeros(MainData.N)

    def __init__(self, d_o, d_N):
        self.d_o = d_o
        self.d_N = d_N

    def bound_val(self):
        self.a[0], self.b[0], self.c[0], self.d[0] = 1, 0, 0, self.d_o
        self.a[MainData.N - 1], self.b[MainData.N - 1], self.c[MainData.N - 1], self.d[MainData.N - 1] = 1, 0, 0, self.d_N
