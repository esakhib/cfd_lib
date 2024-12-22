import logging

import matplotlib.pyplot as plt

from solvers.diffusion_convection.models import DiffusionConvectionSolverModel
from utils.equation_type import get_input_data_by_equation, EquationTypeEnum

logging.getLogger().setLevel(logging.INFO)

# set equation type
equation_type = EquationTypeEnum.DIFFUSION_CONVECTION

# get data for solving equation
input_data = get_input_data_by_equation(
    equation_type=equation_type
)

if input_data is not None:
    # initialize the equation
    equation_solver: DiffusionConvectionSolverModel = input_data.equation_solver(
        input_data=input_data
    )

    # solve numerical
    equation_solver.solve_numerical()

    time_keys = list(equation_solver.output_data.total_solutions.keys())
    needed_items = {0: 'k', 2: 'b', 4: 'r', 6: 'y', 8: 'm'}

    # plot results
    fig = plt.figure(figsize=(20, 10), dpi=100)
    plt.grid(True)
    plt.title('Распределение концентрации C=C(x)', size=20)

    for current_idx, current_color in needed_items.items():
        plt.plot(equation_solver.output_data.grid,
                 equation_solver.output_data.total_solutions[time_keys[current_idx]].reshape(-1),
                 f'.-{current_color}',
                 markersize=15,
                 label=f'Численное решение в момент времени t = {current_idx} сек')
    plt.xlabel('Длина L, м', fontsize=20)
    plt.ylabel('Концентрация C', fontsize=20)
    plt.legend(loc='best', prop={'size': 20})
    plt.tick_params(axis='both', which='major', labelsize=20)
    fig.tight_layout()
    plt.show()
