import numpy as np

from input_data import MainData
from output_data import OutputData


class NumericalSolution:

    def thomas_solution(self):
        OutputData.P[0] = OutputData.b[0] / OutputData.a[0]
        OutputData.Q[0] = OutputData.d[0] / OutputData.a[0]
        for i in range(1, MainData.N):
            OutputData.P[i] = OutputData.b[i] / (OutputData.a[i] - (OutputData.c[i] * OutputData.P[i - 1]))
            OutputData.Q[i] = (OutputData.d[i] + (OutputData.c[i] * OutputData.Q[i - 1])) / (OutputData.a[i] - (OutputData.c[i] * OutputData.P[i - 1]))

        OutputData.T_num[MainData.N - 1] = OutputData.Q[MainData.N - 1]
        for i in range(MainData.N - 1, 0, -1):
            OutputData.T_num[i - 1] = OutputData.P[i - 1] * OutputData.T_num[i] + OutputData.Q[i - 1]




