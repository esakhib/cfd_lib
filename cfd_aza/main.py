from cfd_aza.solvers.heat_conduction.inout_data import InputDataHeatConductivity
from cfd_aza.solvers.heat_conduction.inout_data import TimeData
from cfd_aza.solvers.heat_conduction.inout_data import OutputData

from cfd_aza.solvers.heat_conduction.solver import HeatConductivity

from cfd_aza.solvers.init import BoundaryCondition
from cfd_aza.solvers.init import BoundaryType

from cfd_aza.visual.plotting import Visual




input_data = InputDataHeatConductivity(
    N = 5,
    length = 10,
    T_init = 30,
    T_right = 500,
    T_left = 100,
    k = 5,
    q_left = 0,
    q_right = 100,
    T_env = 0,
    h = 0,
)

time_data = TimeData(
    all_time = 10.0,
    delta_time = 5.0
)

boundary_condition = BoundaryCondition(
    left_side = BoundaryType.Dirichlet,
    right_side = BoundaryType.Neumann
)

equation = HeatConductivity(
    input_data = input_data,
    time_data = time_data,
    boundary_condition = boundary_condition
)

T_old_solution_set = equation.time_solver()
L = equation.L

output_data = OutputData(
    solution_set = T_old_solution_set,
    L = L,
    delta_time = time_data.delta_time,
    all_time = time_data.all_time
)

plotting = Visual(output_data = output_data)
plotting.plotting()
