import numpy as np   # importing mathematical library
import matplotlib.pyplot as mp   # importing graphical library

import thomas_func as Tf  # importing TDMA function

# considering of task with constant coefficient of heat-conductivity

N: int = 1000   # quantity of control volumes
length: float = 1000   # length of whole thing
dx: float = length/(N-1)   # length of control volume
delta = 0.1   # extra lenght just for correct programm working
L = np.arange(start = 0, stop = length + delta, step = dx)  # array with control volumes

k = np.array([203] * (N + 1), float)  # NEED CORRECT VALUE OF TEMPERATURE CONDUCTIVITY COEF,
#                                       NOT HEAT CONDUCTIVITY!!!
T = np.zeros(N)
a = np.zeros(N)
b = np.zeros(N)
c = np.zeros(N)
d = np.zeros(N)
P = np.zeros(N)
Q = np.zeros(N)

# boundary values
T[0] = 300
a[0], b[0], c[0], d[0] = 1, 0, 0, T[0]
T[N - 1] = 50
a[N - 1], b[N - 1], c[N - 1], d[N - 1] = 1, 0, 0, T[N - 1]


# filling magic arrays
for i in range(1, N - 1):  #
    c[i] = k[i - 1] / dx
    b[i] = k[i + 1] / dx
    a[i] = c[i] + b[i]

T = Tf.thomas(a, b, c, d, P, Q, N, T)  # run TDMA (namely magic)


# visualising the task of heat conductivity
mp.plot(L, T)
mp.xlabel('Length')
mp.ylabel('Temperature')
mp.title('Numerical solution of heat coductivity')
mp.show()


