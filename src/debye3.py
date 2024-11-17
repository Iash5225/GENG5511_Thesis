import numpy as np
from constants import ADEB3, ZERO, HALF, ONE, THREE, FOUR, EIGHT, TWENTY, PT375, DEBINF


def Debye3(Xvalue):
    """
    Compute the Debye function of order 3.

    Parameters:
        Xvalue (float): Input value for the Debye3 function.

    Returns:
        float: Computed Debye3 value.
    """
    X = Xvalue

    # Error check
    if X < ZERO:
        return ZERO

    # Machine-dependent constants
    T = 2.2251E-308  # Smallest positive number
    XLIM1 = -np.log(T)
    XK = ONE / THREE
    XKI = (ONE / DEBINF) ** XK
    RK = T ** XK
    XLIM2 = XKI / RK

    T = 1.11E-16  # Smallest relative spacing
    XLOW = np.sqrt(T * EIGHT)
    XUPPER = -np.log(T + T)
    T = T / 100.0

    # Loop to determine NTERMS and calculate Debye3
    for NTERMS in range(19, 0, -1):
        if abs(ADEB3[NTERMS - 1]) > T:
            if X <= FOUR:
                if X < XLOW:
                    Debye3 = ((X - 7.5) * X + TWENTY) / TWENTY
                else:
                    T = ((X * X / EIGHT) - HALF) - HALF
                    Debye3 = Cheval(NTERMS, ADEB3, T) - PT375 * X
            else:
                if X > XLIM2:
                    Debye3 = ZERO
                else:
                    Debye3 = ONE / (DEBINF * X**3)
                    if X < XLIM1:
                        EXPMX = np.exp(-X)
                        if X > XUPPER:
                            SUM = (((X + THREE) * X + 6) * X + 6) / (X**3)
                        else:
                            SUM = ZERO
                            RK = np.floor(XLIM1 / X)
                            NEXP = int(RK)
                            XK = RK * X
                            for _ in range(NEXP, 0, -1):
                                XKI = ONE / XK
                                T = (((6 * XKI + 6) * XKI + 3) * XKI + 1) / RK
                                SUM = SUM * EXPMX + T
                                RK -= ONE
                                XK -= X
                            Debye3 -= THREE * SUM * EXPMX
            break

    return Debye3
