import numpy as np
import matplotlib.pyplot as mp


class Visual:
    def __init__(self, output_data):
        """
          Visualising the task
        """
        self._solution_set: np.ndarray = output_data.solution_set
        self._L: np.ndarray = output_data.L
        self._delta_time: float = output_data.delta_time
        self._all_time: float = output_data.all_time


    def plotting(self):
        self._time_iter: float = self._delta_time
        self._i: int = 0
        while (self._time_iter <= self._all_time):
            fig, ax = mp.subplots()
            line, = ax.plot(self._L, self._solution_set[self._i], "-*m", label='[T] numerical')
            mp.legend()
            mp.xlabel('Length, [mm]')
            mp.ylabel('Temperature, [°C]')
            mp.title('Numerical solution of heat conductivity')
            mp.axis('scaled')
            mp.show()
            self._time_iter += self._delta_time
            self._i += 1
