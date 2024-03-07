from dataclasses import dataclass, field

import numpy as np


@dataclass
class OutputData:
    time_grid: np.ndarray = field(default_factory=lambda: np.array([]))
    grid: np.ndarray = field(default_factory=lambda: np.array([]))
    numerical_solution: np.ndarray = field(default_factory=lambda: np.array([]))
    analytical_solution: np.ndarray = field(default_factory=lambda: np.array([]))
