from dataclasses import dataclass
import numpy as np


@dataclass
class MainData:
    N: int  # quantity of control volumes
    length: float  # length of whole thing
    k: float  # coefficient of heat conductivity

    def __init__(self, N, length, k):
        self.N = N
        self.length = length
        self.k = k

    def extra_data(self):
        _dx: float = self.length / (self.N - 1)  # length of control volume
        _delta: float = 0.1  # extra lenght just for correct programm working
        _L = np.arange(start=0, stop=self.length + _delta, step=_dx)  # array with control volumes



