import numpy as np
import matplotlib.pyplot as mp


class Visual:
    def __init__(self, output_data):
        self.T_num = output_data.T_num
        self.T_an = output_data.T_an
        self.L = output_data.L
    @property
    def plotting(self):
        mp.plot(self.L, self.T_num, "-*m", label='T_num')
        mp.plot(self.L, self.T_an, "--b", label='T_an')
        mp.legend()
        mp.xlabel('Length')
        mp.ylabel('Temperature')
        mp.title('Numerical solution of heat coductivity')
        mp.show()



