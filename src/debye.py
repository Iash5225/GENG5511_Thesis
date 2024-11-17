import numpy as np
from debye3 import Debye3


def Debye(X):
    """
    Compute Cv/R for one Debye distribution.

    Parameters:
        X (float): The input value for the Debye function.

    Returns:
        float: The computed Debye value.
    """
    # Define the numerical limit
    # Equivalent to MATLAB's log(sqrt(lim)) - 1
    lim = np.log(np.sqrt(np.finfo(float).max)) - 1

    if X < lim:
        Debye = 4 * Debye3(X) - 3 * X / (np.exp(X) - 1)
    else:
        Debye = 4 * Debye3(X) + 3 * X

    return Debye



