import datetime
import logging
import os
import warnings
from pathlib import Path

import plotly.graph_objs as go

from utils.common import timer
from utils.output_dataclasses import OutputData


@timer
def plot_results(results: OutputData, save_output_fig: bool = False, delete_previous_results: bool = True):
    """Plot plotly graphics.

    Parameters
    ----------
    results: OutputData
        Dataclass with output data.
    save_output_fig: bool
        Flag to save output figure to figures folder or not.
    delete_previous_results: bool
        Flag to delete previous figure file.

    Returns
    ----------
    output_fig:
        Output plotly figure.

    """

    logging.info('Start plotting results...')

    output_fig = go.Figure()

    output_fig.add_trace(go.Scatter(x=results.grid,
                                    y=results.analytical_solution,
                                    mode='lines+markers',
                                    name='Analytical solution',
                                    marker=dict(size=10, color='Black')))

    output_fig.add_trace(go.Scatter(x=results.grid,
                                    y=results.numerical_solution.reshape(-1),
                                    mode='markers',
                                    name='Numerical solution',
                                    marker=dict(size=5, color='Red')))

    output_fig.update_layout(legend_orientation="h",
                             legend=dict(x=.5, xanchor="center"),
                             hovermode="x",
                             margin=dict(l=0, r=0, t=0, b=0))

    output_fig.update_traces(hoverinfo="all",
                             hovertemplate="Аргумент: %{x}<br>Функция: %{y}")

    output_fig.show()

    if delete_previous_results:
        [f.unlink() for f in Path('temp_files').glob("*") if f.is_file()]

    # save output figure if needed
    if save_output_fig:
        filename = 'input_file-' + datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + '.json'
        tmp_path = os.getcwd() + '/figures/'

        try:
            os.stat(tmp_path)
        except Exception as e:
            warnings.warn(str(e))
            os.mkdir(tmp_path)

        output_fig.write_html(f'figures/{filename}.html')

    logging.info('End plotting results.')

    return output_fig
