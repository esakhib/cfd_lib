import numpy as np
import matplotlib.pyplot as mp


def tdma_algorithm(a, b, c, d, N, T) -> None:
    '''
      Numerical solution for one-dimension unsteady heat conductivity
         with TDMA (Thomas-algorithm)
    '''

    P = np.zeros(N)
    Q = np.zeros(N)
    P[0] = b[0] / a[0]
    Q[0] = d[0] / a[0]
    for i in range(1, N):
        P[i] = b[i] / (a[i] - (c[i] * P[i - 1]))
        Q[i] = (d[i] + (c[i] * Q[i - 1])) / (a[i] - (c[i] * P[i - 1]))

    T[N - 1] = Q[N - 1]
    for i in range(N - 1, 0, -1):
        T[i - 1] = P[i - 1] * T[i] + Q[i - 1]

def discrete_analogue(k_arr, a, b, c, d, dx, S_p, a_o, S_c, T_old_solution_numerical):
    for i in range(1, N - 1):
        b[i] = k_arr[i - 1] / dx
        c[i] = k_arr[i + 1] / dx
        a[i] = b[i] + c[i] + a_o - (S_p * dx)
        d[i] = S_c * dx + a_o * T_old_solution_numerical[i]




N: int = 5
length: float = 10
k: float = 5
T_left: float = 20
T_right: float = 100
# c: float = main_data.c

time: float = 10
N_time: int = 5
dt: float = time / (N_time - 1)
delta: float = 0.1
dx: float = length / (N - 1)
L: np.ndarray = np.arange(start=0, stop=length + delta, step=dx)
a_o: float = (k * dx) / dt  # a_o = (rho * c * dx) / Dt

#linearize temperature source S = S_c + S_p * T[i]
S_c = 0
S_p = 0

T_old_solution_numerical: np.ndarray = T_left * np.ones(shape = N, dtype = float)
T_current_solution_numerical: np.ndarray = np.zeros(shape = N, dtype = float)

# array filled with coefficient of heat conductivity for each control volume
k_arr: np.ndarray = np.array([k] * (N + 1), float)

a: np.ndarray = np.zeros(shape = N, dtype = float)
b: np.ndarray = np.zeros(shape = N, dtype = float)
c: np.ndarray = np.zeros(shape = N, dtype = float)
d: np.ndarray = np.zeros(shape = N, dtype = float)

# boundary conditions for coefficients
a[0], b[0], c[0], d[0] = 1, 0, 0, T_left #self.T_old_solution_numerical[0]
a[N - 1], b[N - 1], c[N - 1], d[N - 1] = 1, 0, 0, T_right #self.T_old_solution_numerical[N - 1]

# filling arrays of coefficients with rule of discrete analogue
discrete_analogue(k_arr, a, b, c, d, dx, S_p, a_o, S_c, T_old_solution_numerical)

for i in range(0, N_time):
    fig, ax = mp.subplots()
    T_current_solution_numerical = tdma_algorithm(a, b, c, d, N, T_old_solution_numerical)
    line, = ax.plot(L, T_current_solution_numerical, "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.draw()
    mp.gcf().canvas.flush_events()
    time.sleep(0.02)
    mp.show()
    T_old_solution_numerical = T_current_solution_numerical
























