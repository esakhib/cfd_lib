import datetime
import json
import logging
import os
import warnings
from pathlib import Path

from utils.common import timer
from utils.input_dataclasses import *


@timer
def prepare_input_json(save_input_data: bool = False, delete_previous_results: bool = True) -> str:
    """Parse data from dataclasses and convert it into json file.

    Parameters
    ----------
    save_input_data: bool
        Flag to save input data to tem_files folder or not.
    delete_previous_results: bool
        Flag to delete previous json files.

    Returns
    ----------
    input_json: str
        Input json file.

    """

    logging.info('Start preparing input data...')

    # grid data
    nx, ny, x_length, y_height = InputGridData.nx, InputGridData.ny, InputGridData.x_length, InputGridData.y_height

    # tolerance
    abs_tol, rel_tol = InputSolverData.abs_tol, InputSolverData.rel_tol

    # steady-state heat conductivity input data
    sc, sp = InputHeatConductivityData.sc, InputHeatConductivityData.sp

    k, cp, lambda_coef, rho = (InputHeatConductivityData.k,
                               InputHeatConductivityData.cp,
                               InputHeatConductivityData.lambda_coef,
                               InputHeatConductivityData.rho)

    t_init, t_left, t_right = (InputHeatConductivityData.t_init,
                               InputHeatConductivityData.t_left,
                               InputHeatConductivityData.t_right)

    # convert input data to dict
    input_dict = {'grid_data': {'nx': nx, 'ny': ny, 'x_length': x_length, 'y_height': y_height},

                  'solver_data': {'abs_tol': abs_tol, 'rel_tol': rel_tol},

                  'heat_conductivity_data': {'sc': sc, 'sp': sp, 'k': k, 'cp': cp, 'lambda_coef': lambda_coef,
                                             'rho': rho, 't_init': t_init, 't_left': t_left, 't_right': t_right}}

    # convert dict to str (json format)
    input_json = json.dumps(input_dict)

    if delete_previous_results:
        [f.unlink() for f in Path('figures').glob("*") if f.is_file()]

    # save input file if needed
    if save_input_data:
        filename = 'input_file-' + datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + '.json'
        tmp_path = os.getcwd() + '/temp_files/'

        try:
            os.stat(tmp_path)
        except Exception as e:
            warnings.warn(str(e))
            os.mkdir(tmp_path)

        with open(tmp_path + filename, 'w') as f:
            f.write(input_json)

    logging.info('End preparing input data.')

    return input_json
