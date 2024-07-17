def analitical(T_right, T_left, dx, N, T):
    T[0] = T_right
    for i in range(1, N):
        T[i] = T_left - (T_left - T_right) * (N - i) / N

    return T