import numpy as np


def poly_fit3(b, T, Tt, pt):
    """
    Polynomial Fit 3

    Parameters:
    b  : list or array of coefficients [b1, b2, b3, b4]
    T  : temperature (float or np.ndarray)
    Tt : triple point temperature
    pt : triple point pressure

    Returns:
    Melting pressure or equivalent value from polynomial expression
    """
    T = np.asarray(T)
    term = T / Tt - 1
    poly_fit3 = pt * (1 + b[0] * term ** b[1] + b[2] * term ** b[3])
    return poly_fit3


def poly_fit2(b, T, Tt, pt):
    """
    Polynomial Fit 2 (Exponential Form)

    Parameters:
    b  : list or array of coefficients [b1, b2, b3]
    T  : temperature (float or np.ndarray)
    Tt : triple point temperature
    pt : triple point pressure

    Returns:
    Sublimation pressure or equivalent value from exponential expression
    """
    T = np.asarray(T)
    term = 1 - T / Tt
    exponent = (Tt / T) * (b[0] * term ** 1 + b[1]
                           * term ** 1.5 + b[2] * term ** 5)
    poly_fit2 = pt * np.exp(exponent)
    return poly_fit2
