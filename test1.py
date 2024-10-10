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
total_time: float = 100

dx: float = x_length / (nx - 1)
dy: float = 1.0
dt: float = total_time / (nt - 1)

c_init: float = 0.0  # НУ
c_left_wall: float = 0.5  # левое ГУ (I рода)
d: float = 1E-12  # diffusion coefficient, m^2 / sec  9.46E-19

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

# цикл через временные слои
while current_time <= total_time:
    # сохраняем решение для временного слоя current_time в словарь
    solutions[current_time] = old_solution
    velocity[current_time] = u_sed

    # для отрисовки
    solutions[current_time][0] = (solutions[current_time][0] + solutions[current_time][1]) / 2.0
    solutions[current_time][-1] = (solutions[current_time][-1] + solutions[current_time][-2]) / 2.0

    # обновим скорость
    u_sed_e = u_sed_w = u_sed

    # boundary_type == BoundaryType.Dirichlet
    a_e[0] = -1.0
    a_w[0] = 0.0
    a_p[0] = 1.0
    b[0] = 2.0 * c_left_wall

    # boundary_type == BoundaryType.Robin
    a_e[nx - 1] = 0.0
    a_w[nx - 1] = 1.0 + u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    a_p[nx - 1] = 1.0 - u_sed_w[nx - 1] * dx_w / (2.0 * d_w)
    b[nx - 1] = 0.0

    # инициализиурем дискретный аналог, используя решение на текущем временном слое
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

# для отрисовки
grid[1] = (grid[0] + grid[1]) / 2.0
grid[-2] = (grid[-1] + grid[-2]) / 2.0

# plot results
fig = plt.figure(figsize=(20, 10), dpi=100)
plt.grid(True)
plt.title('Распределение концентрации C=C(x)', size=20)
plt.ticklabel_format(useOffset=False)

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
         label=f'Численное решение в момент времени t = {list(solutions.keys())[-2]} сек')

plt.plot(grid,
         list(solutions.values())[-1].reshape(-1),
         marker='.',
         c=next(cycol),
         markersize=15,
         label=f'Численное решение в момент времени t = {list(solutions.keys())[-1]} сек')

plt.xlabel('Длина L, м', fontsize=20)
plt.ylabel('Концентрация C', fontsize=20)
plt.legend(loc='best', prop={'size': 20})
plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()
