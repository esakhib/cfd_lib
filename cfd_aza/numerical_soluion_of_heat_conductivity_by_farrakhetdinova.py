import numpy as np   # importing mathematical library
import matplotlib.pyplot as mp   # importing graphical library

import thomas_func as Tf  # importing TDMA function
import analitical_func as Af


# considering of task with constant coefficient of heat-conductivity

N: int = 20  # quantity of control volumes
length: float = 1000   # length of whole thing
dx: float = length/(N-1)   # length of control volume
delta: float = 0.1   # extra lenght just for correct programm working
L = np.arange(start = 0, stop = length + delta, step = dx)  # array with control volumes

k = np.array([8.418 * 10 ** (-5)] * (N + 1), float)  # coefficient of temperature conductivity
T_num = np.zeros(N)
T_an = np.zeros(N)
a = np.zeros(N)
b = np.zeros(N)
c = np.zeros(N)
d = np.zeros(N)
P = np.zeros(N)
Q = np.zeros(N)

# boundary values
d_o = 700
a[0], b[0], c[0], d[0] = 1, 0, 0, d_o
d_N = 50
a[N - 1], b[N - 1], c[N - 1], d[N - 1] = 1, 0, 0, d_N


# filling coefficients of discrete analogue
for i in range(1, N - 1):  #
    c[i] = k[i - 1] / dx
    b[i] = k[i + 1] / dx
    a[i] = c[i] + b[i]


T_num = Tf.thomas(a, b, c, d, P, Q, N, T_num)  # run TDMA

T_an = Af.analitical(d_o, d_N, dx, N, T_an)  # run analitical solution



# visualising the task of heat conductivity
mp.plot(L, T_num, "-*m", label = 'T_num')
mp.plot(L, T_an, "--b", label = 'T_an')
mp.legend()
mp.xlabel('Length')
mp.ylabel('Temperature')
mp.title('Numerical solution of heat coductivity')
mp.show()


