from dataclasses import dataclass
import numpy as np

from calculation_module import Equations
from input_data import ExtraData


class OutputData(ExtraData):
    T_num = Equations.T_num
    T_an = Equations.T_an
    a = Equations.a
    b = Equations.b
    c = Equations.c
    d = Equations.d
    P = Equations.P
    Q = Equations.Q

