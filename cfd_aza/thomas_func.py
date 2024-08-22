import numpy as np

def thomas(a, b, c, d, N, T):
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

