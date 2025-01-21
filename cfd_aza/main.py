from cfd_aza.solvers.heat_conduction.inout_data import InputData
from cfd_aza.solvers.heat_conduction.inout_data import TimeData
from cfd_aza.solvers.heat_conduction.inout_data import OutputData

from cfd_aza.solvers.heat_conduction.solver import HeatConductivity

from cfd_aza.solvers.heat_conduction.init import BoundaryCondition
from cfd_aza.solvers.heat_conduction.init import BoundaryType
#from cfd_aza.visual.plotting import *

import matplotlib.pyplot as mp


input_data = InputData(
    N = 5,
    length = 10,
    T_init = 30,
    T_right = 50,
    T_left = 0,
    k = 5,
    q_left = 50,
    q_right = 0,
    T_env = 0,
    h = 0,
)

time_data = TimeData(
    all_time = 10.0,
    delta_time = 2.0
)

boundary_condition = BoundaryCondition(
    left_side = BoundaryType.Neumann,
    right_side = BoundaryType.Dirichlet
)

equation = HeatConductivity(
    input_data = input_data,
    time_data = time_data,
    boundary_condition = boundary_condition
)

output_data = OutputData()


T_old_solution_set = equation.time_solver()


# отрисовка

time_iter: float = time_data.delta_time  # текущее время
i: int = 0  # номер итерации
while (time_iter <= time_data.all_time):
    fig, ax = mp.subplots()
    line = ax.plot(equation.L, T_old_solution_set[i], "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.axis('scaled')
    mp.show()
    time_iter += time_data.delta_time
    i += 1




