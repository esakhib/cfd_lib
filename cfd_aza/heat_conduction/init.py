from dataclasses import dataclass
from enum import Enum

class BoundaryType(Enum):
    Dirichlet = 'Dirichlet'
    Neumann = 'Neumann'
    Robin = 'Robin'

@dataclass
class BoundaryCondition:
    left_side: BoundaryType
    right_side: BoundaryType
