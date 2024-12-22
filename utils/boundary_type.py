from dataclasses import dataclass
from enum import Enum


class BoundaryConditionsType(Enum):
    Dirichlet = 'Dirichlet'  # I BC (value at the boundary)
    Neumann = 'Neumann'  # II BC (flow at the boundary)
    Robin = 'Robin'  # III BC


@dataclass
class BoundaryConditionsModel:
    left_side: BoundaryConditionsType
    right_side: BoundaryConditionsType
    upper_side: BoundaryConditionsType | None
    lower_side: BoundaryConditionsType | None
