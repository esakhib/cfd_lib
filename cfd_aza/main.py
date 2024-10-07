from cfd_aza.heat_conduction.io_data import *
from cfd_aza.heat_conduction.calculation_unsteady import *
#from cfd_aza.visual.plotting import *

import matplotlib.pyplot as mp
import time



main_data = InputData(N = 20, length = 20, T_right = 100, T_left = 20, k = 5)

equation = Solutions(main_data = main_data, T_old_solution_numerical = main_data.T_left * np.ones(shape = main_data.N, dtype = float))

output_data = OutputData()
output_data.L = equation.L
#output_data.T_analytical = equation.analytical_solution()



for i in range(0, equation.N_time):
    fig, ax = mp.subplots()
    output_data.T_current_solution_numerical = equation.thomas_solution()
    line, = ax.plot(output_data.L, output_data.T_current_solution_numerical, "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.draw()
    mp.gcf().canvas.flush_events()
    time.sleep(0.02)
    mp.show()
    equation = Solutions(main_data = main_data, T_old_solution_numerical = output_data.T_current_solution_numerical)





