import numpy as np


def psub(T):
    """
    Calculate the sublimation pressure at temperature T.

    Parameters:
    T (float or np.ndarray): Temperature in Kelvin

    Returns:
    float or np.ndarray: Sublimation pressure in MPa
    """
    # Constants
    pt = 0.068891  # MPa
    Tt = 83.806    # K

    # Coefficients
    d4 = -10.763
    d5 = -1.526
    d6 = -0.4245

    # Ensure numpy array for element-wise operations
    T = np.asarray(T)

    term = (1 - T / Tt)
    p_sub = pt * np.exp((Tt / T) * (d4 * term + d5 * term**1.5 + d6 * term**5))

    return p_sub


def pmelt(T):
    """
    Calculate the melting pressure at temperature T.

    Parameters:
    T (float or np.ndarray): Temperature in Kelvin

    Returns:
    float or np.ndarray: Melting pressure in MPa
    """
    # Constants
    pt = 0.068891  # MPa
    Tt = 83.806    # K
    b1 = 1506.5415
    b2 = 1.73136
    b3 = 4677.1597
    b4 = 0.9849295

    T = np.asarray(T)  # Ensure array-like input
    term = (T / Tt - 1)
    p_melt = pt * (1 + b1 * term**b2 + b3 * term**b4)

    return p_melt
