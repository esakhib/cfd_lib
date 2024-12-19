import numpy as np

def analytical_formula(T_right: float,
               T_left: float,
               N: int,
               T: np.ndarray) -> None:
    '''
      Analytical solution for one-dimension stationary heat conductivity
    '''

    T[0] = T_right
    for i in range(1, N):
        T[i] = T_left - (T_left - T_right) * (N - i) / N
