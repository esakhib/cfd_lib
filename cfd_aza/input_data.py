from dataclasses import dataclass
import numpy as np


@dataclass
class MainData:
    N: int # quantity of control volumes
    length: float  # length of whole thing
    k: float  # coefficient of heat conductivity
    T_left: float  # boundary condition
    T_right: float  # boundary condition

class ExtraData:
    _dx: float = MainData.length / (MainData.N - 1)  # length of control volume
    _delta: float = 0.1  # extra length just for correct program working
    _L: float = np.arange(start=0, stop=MainData.length + _delta, step=_dx)  # array with control volumes
