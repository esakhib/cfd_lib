from cfd_aza.heat_conduction.io_data import *
from cfd_aza.heat_conduction.calculation import *
from cfd_aza.visual.plotting import *



main_data = InputData(N = 20, length = 1000, T_right = 700, T_left = 50, k = 5)

equation = Solutions(main_data = main_data)

output_data = OutputData()
output_data.L = equation._L
output_data.T_numerical = equation.thomas_solution()
#output_data.T_analytical = equation.analytical_solution()

plot = Visual(output_data = output_data)
plot.plotting()

t=0



