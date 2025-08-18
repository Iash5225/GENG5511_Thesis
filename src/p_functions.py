import numpy as np

# ---- Sublimation (safe) ----


def psub(
    T,
    pt,
    Tt,
    coeffs=(-10.763, -1.526, -0.4245),
    *,
    clip_exp=700.0,   # |exponent| cap to avoid exp overflow
    ratio_cap=1e6     # cap for Tt/T to avoid overflow warnings at tiny T
):
    """
    Sublimation pressure P_sub (MPa) at temperature T (K).
    Returns NaN where T > Tt (outside sublimation domain).
    """
    # Preserve scalar behavior
    was_scalar = np.isscalar(T)

    T = np.asarray(T, dtype=float)
    P = np.full_like(T, np.nan, dtype=float)  # default NaN out-of-domain

    # Physical domain mask (sublimation only below/at triple point)
    m = T <= Tt
    if not np.any(m):
        return P.item() if was_scalar else P

    Tm = T[m]
    tiny = np.finfo(float).tiny
    Tm_safe = np.maximum(Tm, tiny)

    # Safe, warning-free ratio = Tt/T
    ratio = np.empty_like(Tm_safe)
    np.divide(Tt, Tm_safe, out=ratio)             # no divide warnings
    np.clip(ratio, -ratio_cap, ratio_cap, out=ratio)

    e1, e2, e3 = coeffs
    theta = 1.0 - (Tm_safe / Tt)
    theta_pos = np.clip(theta, 0.0, None)         # keep fractional powers real

    expo = ratio * (e1*theta_pos
                    + e2*np.power(theta_pos, 1.5)
                    + e3*np.power(theta_pos, 5.0))
    np.clip(expo, -clip_exp, clip_exp, out=expo)  # avoid exp overflow

    P[m] = pt * np.exp(expo)
    return P.item() if was_scalar else P


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
