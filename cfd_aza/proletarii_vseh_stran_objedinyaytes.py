from input_data import *
from calculation_module import *
from output_data import *
from plotting_module import *



main_data = MainData(N = 20, length = 1000, T_right = 700, T_left = 50, k = 5)

equation = Solutions(main_data = main_data)

output_data = OutputData(L = equation._L)
output_data.T_num = equation.thomas_solution()
output_data.T_an = equation.analytical_solution()

plot = Visual(output_data = output_data)
plot.plotting





