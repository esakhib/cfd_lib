import logging

import matplotlib.pyplot as plt

from utils.equation_type import get_input_data_by_equation, EquationTypeEnum

logging.getLogger().setLevel(logging.INFO)

# set equation type
equation_type = EquationTypeEnum.HEAT_CONDUCTIVITY

# gat data for solving equation
input_data = get_input_data_by_equation(equation_type=equation_type)

# initialize the equation
equation = input_data.equation_solver(input_data=input_data)

# solve analytical
equation.solve_analytical()

# solve numerical
equation.solve_numerical()

# plot results
# _ = plot_results(results=equation.output_data, save_output_fig=True, delete_previous_results=True)

# plot results
fig = plt.figure(figsize=(20, 10), dpi=100)
plt.grid(True)
plt.title('Распределение температуры в зависимости от длины', size=20)
plt.plot(equation.output_data.grid, equation.output_data.numerical_solution.reshape(-1), '.k', markersize=15,
         label='Численное решение')
plt.plot(equation.output_data.grid, equation.output_data.analytical_solution, '.r', markersize=10,
         label='Аналитическое решение')
plt.xlabel('Длина, м', fontsize=20)
plt.ylabel('Температура, К', fontsize=20)
plt.legend(loc='best', prop={'size': 20})
plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()
