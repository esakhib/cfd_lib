from dataclasses import dataclass
import numpy as np


@dataclass
class MainData:
    N: int # quantity of control volumes
    length: float  # length of whole thing
    k: float  # coefficient of heat conductivity
    T_left: float  # boundary condition
    T_right: float  # boundary condition



