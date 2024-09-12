from cfd_aza.heat_conduction.io_data import *
from cfd_aza.heat_conduction.calculation import *
#from cfd_aza.visual.plotting import *

import matplotlib.pyplot as mp
import time



main_data = InputData(N = 20, length = 1000, T_right = 700, T_left = 50, k = 5, S_p = 10, S_c = 0, c = 10, rho = 10)

equation = Solutions(main_data = main_data)

output_data = OutputData()
output_data.L = equation.L
#output_data.T_analytical = equation.analytical_solution()


fig, ax = mp.subplots()
for i in range(0, equation.N_time):
    output_data.T_current_solution_numerical = equation.thomas_solution()
    line, = ax.plot(output_data.L, output_data.T_current_solution_numerical, "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.draw()
    mp.gcf().canvas.flush_events()
    time.sleep(0.02)
    equation.T_old_solution_numerical = output_data.T_current_solution_numerical
    mp.show()




