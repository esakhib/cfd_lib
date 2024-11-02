import numpy as np
import matplotlib.pyplot as mp

import time

""" 
---------------------
|     ГУ 2 рода     |
---------------------
"""

def tdma_algorithm(
        a_p: np.ndarray,
        a_w: np.ndarray,
        a_e: np.ndarray,
        b: np.ndarray,
        N: int,
        T: np.ndarray) -> np.ndarray:
    '''
      Numerical solution for one-dimension unsteady heat conductivity
         with TDMA (Thomas-algorithm)
    '''

    P = np.zeros(N)
    Q = np.zeros(N)
    P[0] = a_w[0] / a_p[0]
    Q[0] = b[0] / a_p[0]
    for i in range(1, N):
        P[i] = a_w[i] / (a_p[i] - (a_e[i] * P[i - 1]))
        Q[i] = (b[i] + (a_e[i] * Q[i - 1])) / (a_p[i] - (a_e[i] * P[i - 1]))

    T[N - 1] = Q[N - 1]
    for i in range(N - 1, 0, -1):
        T[i - 1] = P[i - 1] * T[i] + Q[i - 1]

    # T[0] = T_left
    # T[N-1] = T_right

    return T

def discrete_analogue(   # это функция для коэф дискр аналога,
                         # она не используется, т.к. раскрыта
                         # в цикле в теле основной программы
        k_arr: np.ndarray,
        a_p: np.ndarray,
        a_w: np.ndarray,
        a_e: np.ndarray,
        b: np.ndarray,
        dx: float,
        S_p: float,
        a_o: float,
        S_c: float,
        T_old_solution_numerical: np.ndarray,
        T_left: float,
        T_right: float) -> None:

    # boundary conditions for coefficients
    # a[0], b[0], c[0], d[0] = 1, 0, 0, T_left  # self.T_old_solution_numerical[0]
    # a[N - 1], b[N - 1], c[N - 1], d[N - 1] = 1, 0, 0, T_right  # self.T_old_solution_numerical[N - 1]
    a_p[0], a_w[0], a_e[0], b[0] = 1, 0, -1, 2 * T_left   # коэф для фиктивного к.о.
    a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = 1, -1, 0, 2 * T_right
    # a[0], b[0], c[0], d[0] = 1, 0, 0, 2*T_left/3  # self.T_old_solution_numerical[0]
    # a[N - 1], b[N - 1], c[N - 1], d[N - 1] = 1, 0, 0, 2*T_right/3  # self.T_old_solution_numerical[N - 1]


    for i in range(1, N - 1):
        a_w[i] = k_arr[i - 1] / dx
        a_e[i] = k_arr[i + 1] / dx
        a_p[i] = a_w[i] + a_e[i] + a_o - (S_p * dx)
        b[i] = S_c * dx + a_o * T_old_solution_numerical[i]



# данные, касающиеся самого тела
N: int = 7
length: float = 10.0
k: float = 1000.0
T_left: float = 100.0
T_right: float = 500.0
q_left: float = 1000 # поток слева
q_right: float = 0 # поток справа
delta: float = 0.1
dx: float = length / (N - 1)
L: np.ndarray = np.arange(start=0, stop=length + delta, step=dx)
# c: float = main_data.c

# данные, касающиеся времени
all_time: float = 100.0
time_steps: int = 5
dt: float = all_time / (time_steps - 1)
a_o: float = (k * dx) / dt  # a_o = (rho * c * dx) / Dt

# линеаризация источника S = S_c + S_p * T[i]
S_c = 0
S_p = 0

T_init: float = T_left # начальная температура
T_old_solution_numerical: np.array = np.array([], dtype = float)
T_current_solution_numerical: np.ndarray = T_init * np.ones(shape = N, dtype = float)

k_arr: np.ndarray = np.array([k] * (N + 1), float) # коэф теплопроводности

a_p: np.ndarray = np.zeros(shape = N, dtype = float)
a_w: np.ndarray = np.zeros(shape = N, dtype = float)
a_e: np.ndarray = np.zeros(shape = N, dtype = float)
b: np.ndarray = np.zeros(shape = N, dtype = float)


# цикл с решением уравнений
time_iter: float = 0.0
while (time_iter <= all_time):
    a_p[0], a_w[0], a_e[0], b[0] = 1, 0, -1, (2 * T_left) # учитываем, что поток и темепратура заданы на фиктивном к.о.
    a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = 1, -1, 0, (2 * T_right) # учитываем, что поток и темепратура заданы на фиктивном к.о.
    for i in range(1, N - 1):
        a_w[i] = k / dx
        a_e[i] = k / dx
        a_p[i] = a_w[i] + a_e[i] + a_o - (S_p * dx)
        b[i] = S_c * dx + a_o * T_current_solution_numerical[i] + q_left + q_right
    T_current_solution_numerical = tdma_algorithm(a_p, a_w, a_e, b, N, T_current_solution_numerical) # получаем решение на данном временном шаге
    T_old_solution_numerical = np.concatenate((T_old_solution_numerical, T_current_solution_numerical)) # записываем отдельно все эти решения
    print(time_iter, ' sec:  ', T_current_solution_numerical)
    time_iter += dt

T_old_solution_numerical = T_old_solution_numerical.reshape((time_steps, N))


# отрисовка
time_iter: float = 0.0
i: int =  0
while (time_iter <= all_time):
    fig, ax = mp.subplots()
    line, = ax.plot(L, T_old_solution_numerical[i], "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.draw()
    mp.gcf().canvas.flush_events()
    time.sleep(0.02)
    mp.show()
    time_iter += dt
    i += 1























