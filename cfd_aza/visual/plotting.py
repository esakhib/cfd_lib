import numpy as np
import matplotlib.pyplot as mp


class Visual:
    def __init__(self, output_data):
        """
          Visualising the task
          --------------------------------------------------------------
          Parameters:
               solution_set: np.ndarray  - multidimensional array with all numerical solutions
               L: np.ndarray  - array with control volumes
               delta_time: float  - time step interval, [sec]
               all_time: float  - all considering time, [sec]
          --------------------------------------------------------------

        """
        self._solution_set: np.ndarray = output_data.solution_set
        self._L: np.ndarray = output_data.L
        self._delta_time: float = output_data.delta_time
        self._all_time: float = output_data.all_time

        self._time_iter: float = 0.0
        self._i: int = 0


    def plot_data(self) -> None:
        self._time_iter: float = self._delta_time
        self._i: int = 0
        while (self._time_iter <= self._all_time):
            mp.plot(self._L, self._solution_set[self._i], "-*", label='%d сек' % (self._time_iter))
            self._time_iter += self._delta_time
            self._i += 1
        # mp.legend()
        mp.xlabel('Length, [mm]')
        mp.ylabel('Temperature, [°C]')
        mp.title('Numerical solution of heat conductivity')
        # mp.axis('scaled')
        mp.show()