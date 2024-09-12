from cfd_aza.heat_conduction.io_data import *
from cfd_aza.heat_conduction.calculation import *
from cfd_aza.visual.plotting import *



main_data = InputData(N = 20, length = 1000, T_right = 700, T_left = 50, k = 5, S_p = 10, S_c = 0, c = 10, rho = 10)

equation = Solutions(main_data = main_data)

output_data = OutputData()
output_data.L = equation._L
#output_data.T_analytical = equation.analytical_solution()



for i in range(0, equation.N_time):
    plot = Visual(output_data = output_data)
    output_data.T_current_solution_numerical = equation.thomas_solution()
    plot.plotting()
    equation.T_old_solution_numerical = output_data.T_current_solution_numerical




