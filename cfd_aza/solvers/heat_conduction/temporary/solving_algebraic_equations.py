import numpy as np

def tdma_algorithm(a: np.ndarray,
           b: np.ndarray,
           c: np.ndarray,
           d: np.ndarray,
           N: int,
           T: np.ndarray):
    '''
      Numerical solution for one-dimension stationary heat conductivity
         with TDMA (Thomas-algorithm)
    '''

    P = np.zeros(N)
    Q = np.zeros(N)
    P[0] = (-b[0]) / (a[0])
    Q[0] = d[0] / (a[0])
    for i in range(1, N):
        P[i] = (-b[i]) / ((a[i]) - ((-c[i]) * P[i - 1]))
        Q[i] = (d[i] + ((-c[i]) * Q[i - 1])) / ((a[i]) - ((-c[i]) * P[i - 1]))

    T[N - 1] = Q[N - 1]
    for i in range(N - 1, 0, -1):
        T[i - 1] = P[i - 1] * T[i] + Q[i - 1]




N: int = 4

a: np.ndarray = np.ones(shape = N, dtype = float)
b: np.ndarray = np.zeros(shape = N, dtype = float)
c: np.ndarray = np.zeros(shape = N, dtype = float)
d: np.ndarray = np.zeros(shape = N, dtype = float)
T: np.ndarray = np.zeros(shape = N, dtype = float)


for i in range(0, N):
    c[i], a[i], b[i], d[i] = map(float, input().split())

print(a)
print(b)
print(c)
print(d)

tdma_algorithm(a, b, c, d, N, T)

print(T)



