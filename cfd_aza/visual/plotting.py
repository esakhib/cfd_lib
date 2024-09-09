import numpy as np
import matplotlib.pyplot as mp
import matplotlib.animation as animation


class Visual:
    def __init__(self, output_data):
        '''
          Visualising the task
        '''
        self._T_numerical: np.ndarray = output_data.T_numerical
        # self._T_analytical: np.ndarray = output_data.T_analytical
        self._L: np.ndarray = output_data.L


    def plotting(self):
        mp.plot(self._L, self._T_numerical, "-*m", label='[T] numerical')
        # mp.plot(self._L, self._T_analytical, "--b", label='[T] analytical')
        mp.legend()
        mp.xlabel('Length, [mm]')
        mp.ylabel('Temperature, [°C]')
        mp.title('Numerical solution of heat conductivity')
        mp.show()

#fig, ax = mp.subplots()

# analytical = ax.scatter(T[0], L[0], c="color", s=3, cmap = "plasma", label='T_analytical')
#numerical = ax.scatter(T[0], z[0], c="color", s=3, cmap = "plasma", label='T_numerical')
#ax.set(xlim=[ 'boundary values' ], ylim=[ 'boundary values' ], xlabel='Length, [mm]', ylabel='Temperature, [°C]')
#ax.legend()


#ani = animation.FuncAnimation(fig=fig, func= 'tdma non-stationary' , frames=40, interval=30)
#mp.show()


