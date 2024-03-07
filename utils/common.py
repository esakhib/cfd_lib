import logging
import time

import numpy as np
from numba import vectorize


def timer(func):
    def wrapper(*args, **kwargs):
        # start the timer
        start_time = time.time()

        # call the decorated function
        result = func(*args, **kwargs)

        # remeasure the time
        end_time = time.time()

        # compute the elapsed time and print it
        execution_time = end_time - start_time

        logging.info(msg=f'Execution time: {np.around(execution_time, 10)} seconds')

        # return the result of the decorated function execution
        return result

    # return reference to the wrapper function
    return wrapper


def parallel_nopython(*args, **kwargs):
    return vectorize(*args, **kwargs, target='parallel', nopython=True, fastmath=True)
