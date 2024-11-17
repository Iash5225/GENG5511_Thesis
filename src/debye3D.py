import numpy as np
from constants import ZERO, PT375, HALF, ONE, THREE, FOUR, EIGHT, DEBINF, ADEB3, SMALLEST_POSITIVE, SMALLEST_SPACING
from cheval import Cheval


def Debye3D(Xvalue):
    """
    Compute the first derivative of the Debye function of order 3.

    Parameters:
        Xvalue (float): Input value for the derivative calculation.

    Returns:
        float: The computed derivative of the Debye3 function.
    """
    # Compute Chebyshev coefficients for the first derivative
    D = [0] * len(ADEB3)
    D[-1] = 36 * ADEB3[-1]
    for i in range(len(ADEB3) - 2, -1, -1):
        D[i] = D[i + 1] + 2 * (i + 1) * ADEB3[i]

    X = Xvalue

    # Error check
    if X < ZERO:
        return ZERO

    # Machine-dependent constants
    XLIM1 = -np.log(SMALLEST_POSITIVE)
    XK = ONE / THREE
    XKI = (ONE / DEBINF) ** XK
    RK = SMALLEST_POSITIVE ** XK
    XLIM2 = XKI / RK

    T = SMALLEST_SPACING
    XLOW = np.sqrt(T * EIGHT)
    XUPPER = -np.log(T + T)
    T = T / 100.0

    # Loop to determine NTERMS and calculate Debye3D
    for NTERMS in range(len(D), 0, -1):
        if abs(D[NTERMS - 1]) > T:
            if X <= FOUR:
                if X < XLOW:
                    Debye3D = (4 * X - 15) / 40
                else:
                    T = ((X * X / EIGHT) - HALF) - HALF
                    Debye3D = 0.25 * X * Cheval(NTERMS, D, T) - PT375
            else:
                if X > XLIM2:
                    Debye3D = ZERO
                else:
                    Debye3D = -THREE / (DEBINF * X**4)
                    if X < XLIM1:
                        EXPMX = np.exp(-X)
                        if X > XUPPER:
                            SUM = ((((X + 3) * X + 9) * X + 18)
                                   * X + 18) / (X**4)
                        else:
                            SUM = ZERO
                            RK = np.floor(XLIM1 / X)
                            NEXP = int(RK)
                            XK = RK * X
                            for _ in range(NEXP, 0, -1):
                                XKI = ONE / XK
                                T = (((18 * XKI + 18) * XKI + 9)
                                     * XKI + 3) * XKI + 1
                                SUM = SUM * EXPMX + T
                                RK -= ONE
                                XK -= X
                            Debye3D += THREE * SUM * EXPMX
            break

    return Debye3D
