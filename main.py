import logging

from prepare_input_data.get_input_json import prepare_input_json
from solvers.heat_conduction import SteadyStateHeatConductivity
from utils.plot_data import plot_results

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
_ = plot_results(results=equation.output_data, save_output_fig=False, delete_previous_results=True)
