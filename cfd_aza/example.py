from dataclasses import dataclass
from enum import Enum

import numpy as np


class BoundaryType(Enum):
    Dirichlet = 'dirichlet'  # I рода (значение на границе)
    Neumann = 'neumann'  # II рода (поток на границе)
    Robin = 'robin'  # III рода


@dataclass
class BoundaryCondition:
    left_side: BoundaryType
    right_side: BoundaryType


boundary_condition = BoundaryCondition(
    left_side=BoundaryType.Dirichlet,
    right_side=BoundaryType.Neumann
)


def apply_bndry_condition(a_p: np.ndarray, a_e: np.ndarray, a_w: np.ndarray, a_s: np.ndarray, b: np.ndarray,
                          t_left: float, t_right: float, t_env: float, q_left: float, q_right: float, dx: float,
                          k: float, c_coef: float, boundary_condition: BoundaryCondition) \
        -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # update the first and the last elements here for the left corner
    if boundary_condition.left_side == BoundaryType.Dirichlet:
        pass
    elif boundary_condition.left_side == BoundaryType.Neumann:
        pass
    elif boundary_condition.left_side == BoundaryType.Robin:
        pass

    # update the first and the last elements here for the left right corner
    if boundary_condition.right_side == BoundaryType.Dirichlet:
        pass
    elif boundary_condition.right_side == BoundaryType.Neumann:
        pass
    elif boundary_condition.right_side == BoundaryType.Robin:
        pass

    # можно:
    # 1) передать массивы ap, ae, aw, as, b целиком и изменить [0] и [size - 1] элементы
    # 2) не передавать массивы вообще, а вычислить ap, ae, aw, as, b для [0] и [size - 1] элементов и вернуть их
    # 3) передать только [0] и [size - 1] элементы, изменить их

    return a_p, a_e, a_w, a_s, b
