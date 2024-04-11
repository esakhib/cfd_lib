import logging

import matplotlib.pyplot as plt

from prepare_input_data.get_input_json import prepare_input_json
from solvers.heat_conduction import SteadyStateHeatConductivity

logging.getLogger().setLevel(logging.INFO)

# parse input data
input_json = prepare_input_json(save_input_data=False, delete_previous_results=True)

# initialize the equation
equation = SteadyStateHeatConductivity(input_json=input_json)

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
