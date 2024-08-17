import logging

import matplotlib.pyplot as plt

from utils.equation_type import get_input_data_by_equation, EquationTypeEnum

logging.getLogger().setLevel(logging.INFO)

# set equation type
equation_type = EquationTypeEnum.DIFFUSION_CONVECTION

# gat data for solving equation
input_data = get_input_data_by_equation(equation_type=equation_type)

# initialize the equation
equation = input_data.equation_solver(input_data=input_data)

# solve analytical
# equation.solve_analytical()

# solve numerical
equation.solve_numerical()

# plot results
# _ = plot_results(results=equation.output_data, save_output_fig=True, delete_previous_results=True)

time_keys = list(equation.output_data.total_solutions.keys())
needed_items = {0: 'k', 2: 'b', 4: 'r', 6: 'y', 8: 'm', len(time_keys) - 1: 'g'}

# plot results
fig = plt.figure(figsize=(20, 10), dpi=100)
plt.grid(True)
plt.title('Распределение концентрации C=C(x)', size=20)

for current_idx, current_color in needed_items.items():
    plt.plot(equation.output_data.grid,
             equation.output_data.total_solutions[time_keys[current_idx]].reshape(-1),
             f'.-{current_color}',
             markersize=15,
             label=f'Численное решение в момент времени t = {current_idx} сек')

plt.xlabel('Длина L, м', fontsize=20)
plt.ylabel('Концентрация C', fontsize=20)
plt.legend(loc='best', prop={'size': 20})
plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()
