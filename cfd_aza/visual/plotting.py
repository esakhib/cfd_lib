import numpy as np
import time
import matplotlib.pyplot as mp
import matplotlib.animation as animation


class Visual:
    def __init__(self, output_data):
        '''
          Visualising the task
        '''
        self._T_numerical: np.ndarray = output_data.T_current_solution_numerical
        # self._T_analytical: np.ndarray = output_data.T_analytical
        self._L: np.ndarray = output_data.L


    def plotting(self, delta_time, all_time, L, T_old_solution_set):
        self._time_iter: float = self._delta_time  # текущее время
        self._i: int = 0  # номер итерации
        while (self._time_iter <= self._all_time):
            fig, ax = mp.subplots()
            line, = ax.plot(self._L, T_old_solution_set[i], "-*m", label='[T] numerical')
            mp.legend()
            mp.xlabel('Length, [mm]')
            mp.ylabel('Temperature, [°C]')
            mp.title('Numerical solution of heat conductivity')
            mp.axis('scaled')
            mp.show()
            time_iter += time_data.delta_time
            i += 1
