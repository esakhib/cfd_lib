import matplotlib.pyplot as mp
import numpy as np
from requests import Req

""" 
---------------------
|     ГУ 2 рода     |
---------------------
"""


########################################################################################################################

# функции

def tdma_algorithm(
        a_p: np.ndarray,
        a_w: np.ndarray,
        a_e: np.ndarray,
        b: np.ndarray,
        temp_size: int,
        temp_arr: np.ndarray) -> np.ndarray:
    """
    Это функция для реализации
    алгоритма TDMA

    Parameters
    ----------
    a_p
    a_w
    a_e
    b
    temp_size
    temp_arr

    Returns
    -------

    """

    P = np.zeros(temp_size)
    Q = np.zeros(temp_size)

    P[0] = a_e[0] / a_p[0]
    Q[0] = b[0] / a_p[0]

    for i in range(1, temp_size):
        P[i] = a_e[i] / (a_p[i] - (a_w[i] * P[i - 1]))
        Q[i] = (b[i] + (a_w[i] * Q[i - 1])) / (a_p[i] - (a_w[i] * P[i - 1]))

    temp_arr[temp_size - 1] = Q[temp_size - 1]

    for i in range(temp_size - 1, 0, -1):
        temp_arr[i - 1] = P[i - 1] * temp_arr[i] + Q[i - 1]

    # temp_arr[0] = T_left
    # temp_arr[temp_size-1] = T_right

    return temp_arr


# -----------------------------------------------------------------------------------------------------------------------


########################################################################################################################
# тело программы
########################################################################################################################

# выбор режима расчета ГУ: справа, слева, оба одновременно
request = Req.left  #  right / left / both

# данные, касающиеся самого тела
N_origin: int = 5  # количество к.о.
N: int = N_origin + 2  # количество к.о. с учетом фиктивных к.о.
length: float = 10.0  # длина всего предмета, m
k: float = 100  # коэффициент температуропроводности, m^2 / сек
T_left: float = 100.0  # температура слева, K
T_right: float = 130.0  # температура справа, K
q_left: float = 1000.0  # поток слева, W / m^2
q_right: float = 700.0  # поток справа, W / m^2
q_v: np.ndarray = np.zeros(shape=N, dtype=float)
# q_v[3] = -50000.0
delta: float = 0.1  # m
dx: float = length / N_origin  # m
L: np.ndarray = np.arange(start=(-dx / 2), stop=length + (dx / 2) + delta, step=dx)
L[0], L[N - 1] = 0, L[N - 1] - (dx / 2)

# данные, касающиеся времени
all_time: float = 10.0  # все рассматриваемое время, sec
dt: float = 10.0  # sec
time_steps: int = int(all_time / dt)  # количество врем промежутков
a_o: float = k * dx / dt  # a_o = (rho * c * dx) / dt

# линеаризация источника S = S_c + S_p * T[i]
S_c: float = 0.0
S_p: float = 0.0

T_init: float = T_left  # начальная температура, K
T_old_solution = T_init * np.ones(shape=N, dtype=float)  # массив для записи решения на старом временном слое

T_old_solution_set: np.array = np.array([], dtype=float)  # массив для записи всех решений

k_arr: np.ndarray = np.array([k] * (N + 1), float)  # массив для коэф температуропровондости

# массивы для коэф дискретного аналога
a_p: np.ndarray = np.zeros(shape=N, dtype=float)
a_w: np.ndarray = np.zeros(shape=N, dtype=float)
a_e: np.ndarray = np.zeros(shape=N, dtype=float)
b: np.ndarray = np.zeros(shape=N, dtype=float)

########################################################################################################################

# цикл с решением уравнений

time_iter: float = dt  # текущее время
while (time_iter <= all_time):
    if (request == Req.left):
        a_p[0], a_w[0], a_e[0], b[0] = 1, 0, 1, q_left * (dx / k)  # учитываем фиктивный к.о.
        a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = 1, -1, 0, 2 * T_right

    if (request == Req.right):
        a_p[0], a_w[0], a_e[0], b[0] = 1, 0, -1, 2 * T_left  # учитываем фиктивный к
        a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = 1, 1, 0, q_right * (dx / k)

    if (request == Req.both):
        a_p[0], a_w[0], a_e[0], b[0] = 1, 0, 1, q_left * (dx / k)
        a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = 1, 1, 0, q_right * (dx / k)

    for i in range(1, N - 1):
        a_w[i] = k / dx
        a_e[i] = k / dx
        a_p[i] = a_w[i] + a_e[i] + a_o - (S_p * dx)
        b[i] = S_c * dx + a_o * T_old_solution[i] + (q_v[i] * dx)

    T_current_solution_numerical = tdma_algorithm(a_p, a_w, a_e, b, N,
                                                  T_old_solution)  # получаем решение на данном временном шаге
    T_old_solution = T_current_solution_numerical

    # получаем температуру на крайних точках, опуская фиктивные к.о.
    T_current_solution_numerical[0] = (T_current_solution_numerical[0] + T_current_solution_numerical[1]) / 2
    T_current_solution_numerical[N - 1] = (T_current_solution_numerical[N - 2] + T_current_solution_numerical[
        N - 1]) / 2

    T_old_solution_set = np.concatenate(
        (T_old_solution_set, T_current_solution_numerical))  # записываем отдельно все эти решения

    print(request.value, "\n")
    print(time_iter, ' sec:  ', T_current_solution_numerical)
    print('         a_p = ', a_p)
    print('         a_e = ', a_e)
    print('         a_w = ', a_w)
    print('         b = ', b)
    print('\n\n')
    time_iter += dt

T_old_solution_set = T_old_solution_set.reshape((time_steps, N))

# -----------------------------------------------------------------------------------------------------------------------

# отрисовка

time_iter: float = dt  # текущее время
i: int = 0  # номер итерации
while (time_iter <= all_time):
    fig, ax = mp.subplots()
    line, = ax.plot(L, T_old_solution_set[i], "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [mm]')
    mp.ylabel('Temperature, [°C]')
    mp.title('Numerical solution of heat conductivity')
    mp.axis('scaled')
    mp.show()
    time_iter += dt
    i += 1

########################################################################################################################
