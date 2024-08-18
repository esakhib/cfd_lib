import datetime
import json
import logging
import os
import warnings
from pathlib import Path

from utils.common import timer


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

    # function prototype for the future

    logging.info('Start preparing input data...')

    # convert input data to dict
    input_dict = {}

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
