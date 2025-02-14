import matplotlib.pyplot as mp
import numpy as np
from cfd_aza.solvers.boundary_type import *

""" 
--------------------------------------
|     черновой вариант программы     |
--------------------------------------
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

#ВЫБРАТЬ ГРАНИЧНЫЕ УСЛОВИЯ:  dirichlet | neumann | robin
top = 'dirichlet'
bottom = 'neumann'

# данные, касающиеся самой системы
N_origin: int = 100  # количество к.о.
N: int = N_origin + 2  # количество к.о. с учетом фиктивных к.о.
length: float = 0.010  # длина всего объекта, m
delta: float = 0.000001  # m
dz: float = length / N_origin  # m
L: np.ndarray = np.arange(start=(-dz / 2), stop=length + (dz / 2) + delta, step=dz)
L[0], L[N - 1] = 0, L[N - 1] - (dz / 2)

r: float = 2e-6  # радиус частицы
rho_partical: float = 1070  # плотность частицы
rho_fluid: float = 1000  # плотность жидкости
myu: float = 1e-4  # коэффициент вязкости
D: float = 1e-9 # коэффициент диффузии
g: float = 9.81  # ускорение свободного падения
v: float = -(2 / 9) * r ** 2 * g * ((rho_partical - rho_fluid) / myu)
print("v = ", v)

C_o: float = 0

C_init: float = 0.1  # начальная концентрация
C_old_solution = C_init * np.ones(shape=N, dtype=float)  # массив для записи решения на старом временном слое

C_old_solution_set: np.array = np.array([], dtype=float)  # массив для записи всех решений

# данные, касающиеся времени
all_time: float = 100.0  # все рассматриваемое время, sec
dt: float = 1  # sec
AAAA: float = 20
time_steps: int = int(all_time / dt)  # количество врем промежутков
a_o: float = dz / dt

integral = np.zeros(shape=time_steps + 2, dtype=float)
integral[0] = np.sum(C_old_solution * dz)

# массивы для коэф дискретного аналога
a_p: np.ndarray = np.zeros(shape=N, dtype=float)
a_w: np.ndarray = np.zeros(shape=N, dtype=float)
a_e: np.ndarray = np.zeros(shape=N, dtype=float)
b: np.ndarray = np.zeros(shape=N, dtype=float)

########################################################################################################################

# цикл с решением уравнений

time_iter: float = dt  # текущее время
while (time_iter <= all_time):
    rho = rho_partical * C_old_solution[N - 1] + rho_fluid * (1 - C_old_solution[N - 1])
    # C = np.sum(C_old_solution) / N
    # print("C = ", C)
    # rho = rho_partical * C + rho_fluid * (1 - C)
    f = rho_partical / rho
    a_p[0], a_w[0], a_e[0], b[0] = ((D / dz) + ((1 - f) * v)), 0, (D / dz), C_o * (1 - f) * v  # учитываем фиктивный к.о.
    a_p[N - 1], a_w[N - 1], a_e[N - 1], b[N - 1] = (D / dz), ((D / dz) + ((1 - f) * v)), 0, -C_o * (1 - f) * v


    for i in range(1, N - 1):
        rho = rho_partical * C_old_solution[i] + rho_fluid * (1 - C_old_solution[i])
        # C = np.sum(C_old_solution) / N_origin
        # rho = rho_partical * C + rho_fluid * (1 - C)
        f = rho_partical / rho
        a_w[i] = D / dz + v * (1 - f)
        a_e[i] = D / dz
        a_p[i] = a_w[i] + a_e[i] + a_o
        b[i] = a_o * C_old_solution[i]

    C_current_solution_numerical = tdma_algorithm(a_p, a_w, a_e, b, N,
                                                  C_old_solution)  # получаем решение на данном временном шаге
    C_old_solution = C_current_solution_numerical
    C_current_solution_numerical[0] = (C_current_solution_numerical[0] + C_current_solution_numerical[1]) / 2
    C_current_solution_numerical[N - 1] = (C_current_solution_numerical[N - 2] + C_current_solution_numerical[
        N - 1]) / 2

    if (time_iter % AAAA == 0):
        C_old_solution_set = np.concatenate(
            (C_old_solution_set, C_current_solution_numerical))  # записываем отдельно все эти решения


    print("time_iter = ", time_iter)
    print("time_iter//dt = ", time_iter//dt)
    integral[time_iter//dt] = np.sum(C_old_solution * dz)

    print(time_iter, ' sec:  ', C_current_solution_numerical)
    print('         a_p = ', a_p)
    print('         a_e = ', a_e)
    print('         a_w = ', a_w)
    print('         b = ', b)
    print('\n\n')
    time_iter += dt


C_old_solution_set = C_old_solution_set.reshape((time_steps // AAAA, N))
# C_old_solution_set = C_old_solution_set.reshape((time_steps, N))

print("integral = ", integral)

# -----------------------------------------------------------------------------------------------------------------------

# отрисовка

time_iter: float = 0  # текущее время
i: int = 0  # номер итерации
while (time_iter <= all_time):
    fig, ax = mp.subplots()
    line, = ax.plot(L, C_old_solution_set[i], "-*m", label='[T] numerical')
    mp.legend()
    mp.xlabel('Length, [m]')
    mp.ylabel('Concentration')
    mp.title('Numerical solution of diffusion')
    #mp.axis('scaled')
    mp.show()
    time_iter += dt
    i += 1

########################################################################################################################
