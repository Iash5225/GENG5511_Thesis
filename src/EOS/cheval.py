import numpy as np


def Cheval(N, a, T):
    """
    This function evaluates a Chebyshev series using the Clenshaw method
    with Reinsch modification as per the method analyzed by Oliver.

    Parameters:
    N : int
        The number of terms in the sequence.
    a : array_like
        Coefficients of the Chebyshev series. Should be of length N.
    T : float
        The value at which to evaluate the series.

    Returns:
    float
        The evaluated Chebyshev series at T.
    """

    ZERO = 0.0
    HALF = 0.5
    TEST = 0.6
    TWO = 2.0
    U1 = ZERO

    a = np.asarray(a)

    # Clenshaw method for |T| < 0.6
    if abs(T) < TEST:
        U0 = ZERO
        Tt = T + T
        for i in range(N-1, -1, -1):
            U2 = U1
            U1 = U0
            U0 = Tt * U1 + a[i] - U2
        Cheval = (U0 - U2) / TWO

    else:
        d1 = ZERO
        U1 = ZERO  # ensure it's reset before usage

        if T > ZERO:
            Tt = 2 * ((T - HALF) - HALF)
            for i in range(N-1, -1, -1):
                d2 = d1
                U2 = U1
                d1 = Tt * U2 + a[i] + d2
                U1 = d1 + U2
