from constants import ZERO, HALF,TWO

def Cheval(N: int, a: list, T:float):
    """
    This function evaluates a Chebyshev series using the Clenshaw method
    with the Reinsch modification for numerical stability.

    Args:
        N (int): The number of terms in the sequence.
        a (list or numpy array): The coefficients of the Chebyshev series (length N).
        T (float): The value at which the series is to be evaluated.

    Returns:
        float: The evaluated Chebyshev series at T.
    """    

    TEST = 0.6

    U1 = ZERO  # Initialize variables

    if abs(T) < TEST:
        # Standard Clenshaw method
        U0 = ZERO
        Tt = TWO * T
        for i in range(N - 1, -1, -1):  # Loop from N-1 to 0
            U2 = U1
            U1 = U0
            U0 = Tt * U1 + a[i] - U2
        Cheval = (U0 - U2) / TWO
    else:
        # Reinsch modification
        d1 = ZERO

        if T > ZERO:
            # T >= 0.6 code
            Tt = TWO * (T - HALF)
            for i in range(N - 1, -1, -1):
                d2 = d1
                U2 = U1
                d1 = Tt * U2 + a[i] + d2
                U1 = d1 + U2
            Cheval = (d1 + d2) / TWO
        else:
            # T <= -0.6 code
            Tt = TWO * (T + HALF)
            for i in range(N - 1, -1, -1):
                d2 = d1
                U2 = U1
                d1 = Tt * U2 + a[i] - d2
                U1 = d1 - U2
            Cheval = (d1 - d2) / TWO

    return Cheval
