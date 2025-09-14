import numpy as np


def Einstein(X):
    """
    Computes Cv/R for one harmonic mode using the Einstein model.

    Parameters:
    X : float
        The dimensionless parameter (usually theta_E/T).

    Returns:
    float
        The value of Cv/R.
    """

    # Machine limit for safe exponentiation, same as in MATLAB
    lim = np.log(np.sqrt(9.99999999999999E+307)) - 1

    if X < lim:
        expX = np.exp(X)
        Einstein = X**2 * expX / (expX - 1)**2
    else:
        Einstein = 0.0

    return Einstein
