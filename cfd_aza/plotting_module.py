import numpy as np
import matplotlib.pyplot as mp

from output_data import *


class Visual:
    def __init__(self, T_num, T_an, L):
        self.T_num = T_num
        self.T_an = T_an
        self.L = L
    def plotting(self):
        mp.plot(self.L, self.T_num, "-*m", label='T_num')
        mp.plot(self.L, self.T_an, "--b", label='T_an')
        mp.legend()
        mp.xlabel('Length')
        mp.ylabel('Temperature')
        mp.title('Numerical solution of heat coductivity')
        mp.show()



