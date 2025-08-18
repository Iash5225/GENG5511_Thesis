import numpy as np

# ---- Sublimation (safe) ----


def psub(T, pt, Tt, coeffs=(-10.763, -1.526, -0.4245), *, clip_exp=700.0):
    """
    Sublimation pressure P_sub (MPa) at temperature T (K).

    Parameters
    ----------
    T : float or array-like
        Temperature in K.
    pt : float
        Triple-point pressure in MPa (units must match your model).
    Tt : float
        Triple-point temperature in K.
    coeffs : (e1, e2, e3)
        Empirical coefficients for the exponent.
    clip_exp : float
        |exponent| clipping to avoid overflow in exp (default ~ log(1e304)).

    Returns
    -------
    P_sub : float or ndarray
        Sublimation pressure in MPa. Returns NaN where T > Tt (out of domain).
    """
    T = np.asarray(T, dtype=float)
    T_safe = np.maximum(T, np.finfo(float).tiny)  # avoid divide-by-zero
    e1, e2, e3 = coeffs

    # theta valid only for T <= Tt; for T > Tt set to 0 to keep fractional powers safe
    theta = 1.0 - (T_safe / Tt)
    theta_pos = np.where(theta > 0.0, theta, 0.0)

    # exponent with safe fractional powers and clipping
    expo = (Tt / T_safe) * (e1*theta_pos + e2 *
                            np.power(theta_pos, 1.5) + e3*np.power(theta_pos, 5.0))
    expo = np.clip(expo, -clip_exp, clip_exp)

    P = pt * np.exp(expo)

    # out of physical domain for sublimation (T > Tt): mark as NaN
    P = np.where(T <= Tt, P, np.nan)

    # preserve scalar behavior
    return P.item() if np.isscalar(T) else P


# ---- Melting (safe) ----
def pmelt(T, pt, Tt, coeffs=(1506.5415, 1.73136, 4677.1597, 0.9849295)):
    """
    Melting pressure P_melt (MPa) at temperature T (K).

    P_melt = pt * [1 + b1*((T/Tt - 1)^b2) + b3*((T/Tt - 1)^b4)]
    with term clamped >= 0 so fractional powers are defined.

    Parameters
    ----------
    T : float or array-like
        Temperature in K.
    pt : float
        Triple-point pressure in MPa.
    Tt : float
        Triple-point temperature in K.
    coeffs : (b1, b2, b3, b4)
        Empirical coefficients.

    Returns
    -------
    P_melt : float or ndarray
        Melting pressure in MPa. For T < Tt, term=0 ⇒ returns pt (as per model).
    """
    T = np.asarray(T, dtype=float)
    b1, b2, b3, b4 = coeffs

    # no negative base to fractional powers
    term = np.clip((T / Tt) - 1.0, 0.0, None)
    P = pt * (1.0 + b1*np.power(term, b2) + b3*np.power(term, b4))

    return P.item() if np.isscalar(T) else P
