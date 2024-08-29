import os
from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import g

from solvers.diffusion_convection.solver_dataclasses import BoundaryType
from solvers.tdma import run_tdma

os.environ["XDG_SESSION_TYPE"] = "xcb"


def f_c(c: float) -> float:
    return (1.0 - c) ** 4.7


def calc_u_sed(u_sed: np.ndarray, c: np.ndarray):
    const_u_sed = 0.2  # 2 / 9
    r0: float = 0.01  # m
    rho1: float = 1000.0  # kg / m^3
    rho2: float = 900.0  # kg / m^3
    mu2: float = 0.6  # Pa * sec

    for i in np.arange(1, nx - 1):
        u_sed[i] = const_u_sed * r0 ** 2.0 * g * (rho1 - rho2) * f_c((c[i] + c[i + 1]) / 2.0) / mu2

    return u_sed


nx: int = 10  # dx = x_length / nx
ny: int = 1  # dy = y_height / ny
x_length: float = 0.1  # m
y_height: float = 1.0  # m
nt: int = 1000  # dt = total_time / nt
total_time: float = 10

dx: float = x_length / (nx - 1)
dy: float = 1.0
dt: float = total_time / (nt - 1)

c_init: float = 0.01  # НУ
c_left_wall: float = 0.07  # левое ГУ (I рода)
c_right_wall: float = 0.07  # правое ГУ (I рода)
q_source = 0.000001  # источник для ГУ II рода
d: float = 1E-5  # diffusion coefficient, m^2 / sec

# элементы КО
a_p: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)
a_e: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)
a_w: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)
b: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)

# скорость
u_sed: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)
u_sed_e: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)
u_sed_w: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)

# решение на текущем временном слое
current_solution: np.ndarray = np.zeros(shape=(nx, ny), dtype=np.float64)

# словари для вывода
solutions = {}
velocity = {}

# сетка
grid = np.arange(start=0.0, stop=x_length, step=dx)
grid = np.append(grid, x_length)

boundary_type: BoundaryType = BoundaryType.Robin  # тип ГУ
dx_e = dx_w = dx  # шаг сетки
d_e = d_w = d  # коэф-т диффузии
current_time = 0.0  # начальное время

# заполним решение, используя НУ
# решение на предыдущем временном слое
old_solution = np.full_like(current_solution, c_init)

# u_sed[:, :] = 0.0001

# цикл через временные слои
while current_time <= total_time:
    # посчитаем скорость, используя концентрацию на текущем временном слое
    # u_sed = calc_u_sed(
    #     u_sed=u_sed,
    #     c=old_solution
    # )

    # сохраняем решение для временного слоя current_time в словарь
    solutions[current_time] = old_solution
    velocity[current_time] = u_sed

    # обновим скорость
    u_sed_e = u_sed_w = u_sed

    a_e[0] = -1.0
    a_w[0] = 0.0
    a_p[0] = 1.0
    b[0] = 2.0 * c_left_wall

    a_e[nx - 1] = 0.0
    a_w[nx - 1] = 1.0 + u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    a_p[nx - 1] = 1.0 - u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    b[nx - 1] = 0.0

    # # инициализиурем дискретный аналог, используя решение на текущем временном слое
    # if boundary_type == BoundaryType.Dirichlet:
    #     a_e[0] = -1.0
    #     a_w[0] = 0.0
    #     a_p[0] = 1.0
    #     b[0] = 2.0 * c_left_wall
    #
    #     a_e[nx - 1] = 0.0
    #     a_w[nx - 1] = -1.0
    #     a_p[nx - 1] = 1.0
    #     b[nx - 1] = 2.0 * c_right_wall
    #
    # if boundary_type == BoundaryType.Neumann:
    #     a_e[0] = 1.0
    #     a_w[0] = 0.0
    #     a_p[0] = 1.0
    #     b[0] = -q_source / d_e * dx_e
    #
    #     a_e[nx - 1] = 0.0
    #     a_w[nx - 1] = 1.0
    #     a_p[nx - 1] = 1.0
    #     b[nx - 1] = q_source / d_w * dx_w
    #
    # if boundary_type == BoundaryType.Robin:
    #     a_e[0] = 1.0 - u_sed_e[0] * dx_e / (2.0 * d_e)
    #     a_w[0] = 0.0
    #     a_p[0] = 1.0 + u_sed_e[0] * dx_e / (2.0 * d_e)
    #     b[0] = 0.0
    #
    #     a_e[nx - 1] = 0.0
    #     a_w[nx - 1] = 1.0 + u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    #     a_p[nx - 1] = 1.0 - u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    #     b[nx - 1] = 0.0

    for i in range(1, nx - 1):
        a_e[i] = d_e / dx_e + max(-u_sed_e[i], 0.0)
        a_w[i] = d_w / dx_w + max(u_sed_w[i], 0.0)
        a_p[i] = dx / dt + d_e / dx_e + d_w / dx_w + (u_sed_e[i] - u_sed_w[i])
        b[i] = dx / dt * old_solution[i]

    # получаем решение на следующем временном слое
    current_solution = run_tdma(
        a=a_p,
        b=a_e,
        c=a_w,
        d=b
    )

    # обновляем решение на текущем временном слое
    old_solution = current_solution

    # переключились на следующий временной слой
    current_time += dt

print(f'L = {x_length} m')
print(f'dx = {dx} m')
print(f'nx = {nx} m')
print(f'total time = {total_time} sec')
print(f'time steps count = {nt}')
print(f'dt = {dt} sec')
print(f'D = {d} m^2 / sec')

max_c = max([np.amax(solution) for solution in list(solutions.values())]) + 1.E-2
cycol = cycle('bgrcmk')
avg_idx = int(len(solutions) / 2)

# plot results
fig = plt.figure(figsize=(20, 10), dpi=100)
plt.grid(True)
plt.title('Распределение концентрации C=C(x)', size=20)
# for current_time, current_solution in solutions.items():
# plt.plot(grid,
#          current_solution.reshape(-1),
#          marker='.',
#          c=next(cycol),
#          markersize=15,
#          label=f'Численное решение в момент времени t = {current_time} сек')

plt.plot(grid,
         list(solutions.values())[0].reshape(-1),
         marker='.',
         c=next(cycol),
         markersize=15,
         label=f'Численное решение в момент времени t = {list(solutions.keys())[0]} сек')

plt.plot(grid,
         list(solutions.values())[1].reshape(-1),
         marker='.',
         c=next(cycol),
         markersize=15,
         label=f'Численное решение в момент времени t = {list(solutions.keys())[1]} сек')

plt.plot(grid,
         list(solutions.values())[avg_idx].reshape(-1),
         marker='.',
         c=next(cycol),
         markersize=15,
         label=f'Численное решение в момент времени t = {list(solutions.keys())[avg_idx]} сек')

plt.plot(grid,
         list(solutions.values())[-1].reshape(-1),
         marker='.',
         c=next(cycol),
         markersize=15,
         label=f'Численное решение в момент времени t = {list(solutions.keys())[-1]} сек')

# plt.ylim(0.0, max_c)
plt.xlabel('Длина L, м', fontsize=20)
plt.ylabel('Концентрация C', fontsize=20)
plt.legend(loc='best', prop={'size': 20})
plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()
