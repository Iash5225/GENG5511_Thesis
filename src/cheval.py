import numpy as np


def _cheval_scalar(a, T):
    """
    Scalar Clenshaw evaluation with Reinsch modification (Oliver) for Chebyshev series.
    a: sequence of length N with coefficients a[0]..a[N-1]
    T: scalar point (expected in [-1, 1])
    Returns sum_{i=1..N} a[i-1] * T_{i-1}(T) in the same convention as your MATLAB code.
    """
    ZERO = 0.0
    HALF = 0.5
    TEST = 0.6
    TWO = 2.0

    # MATLAB uses U1 initialized to ZERO outside both branches
    U1 = ZERO

    if abs(T) < TEST:
        # Standard Clenshaw
        U0 = ZERO
        Tt = T + T  # 2T
        for i in range(len(a), 0, -1):  # i = N..1
            U2 = U1
            U1 = U0
            U0 = Tt * U1 + a[i-1] - U2
        return (U0 - U2) / TWO

    # Reinsch modification
    d1 = ZERO

    if T > ZERO:
        # T >= 0.6 path
        Tt = (T - HALF) - HALF   # T - 1
        Tt = Tt + Tt             # 2(T - 1) = 2T - 2
        for i in range(len(a), 0, -1):
            d2 = d1
            U2 = U1
            d1 = Tt * U2 + a[i-1] + d2
            U1 = d1 + U2
        return (d1 + d2) / TWO
    else:
        # T <= -0.6 path
        Tt = (T + HALF) + HALF   # T + 1
        Tt = Tt + Tt             # 2(T + 1) = 2T + 2
        for i in range(len(a), 0, -1):
            d2 = d1
            U2 = U1
            d1 = Tt * U2 + a[i-1] - d2
            U1 = d1 - U2
        return (d1 - d2) / TWO


def cheval(a, T):
    """
    Vector-friendly wrapper around _cheval_scalar.
    a : 1D array-like of Chebyshev coefficients (length N)
    T : scalar or ndarray of evaluation points in [-1, 1]
    """
    a = np.asarray(a, dtype=float)
    T = np.asarray(T, dtype=float)
    if T.ndim == 0:
        return _cheval_scalar(a, float(T))
    # vectorize over T
    vec = np.vectorize(lambda x: _cheval_scalar(a, float(x)))
    return vec(T)
