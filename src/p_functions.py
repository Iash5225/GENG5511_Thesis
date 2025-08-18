import numpy as np

# ---- Sublimation (safe) ----


def psub(
    T, pt, Tt, coeffs=(-10.763, -1.526, -0.4245),
    *, clip_exp=700.0, ratio_cap=1e6, T_floor=1e-12
):
    import numpy as np
    was_scalar = np.isscalar(T)
    T = np.asarray(T, dtype=float)

    # Only compute in physical domain: T <= Tt
    mask = T <= Tt
    P = np.full_like(T, np.nan, dtype=float)
    if not np.any(mask):
        return P.item() if was_scalar else P

    Tm = T[mask]
    Tm_safe = np.maximum(Tm, T_floor)

    # Build Tt/T without performing a dangerous divide
    invT = np.empty_like(Tm_safe)
    invT.fill(0.0)
    good = Tm_safe > T_floor
    invT[good] = 1.0 / Tm_safe[good]
    invT[~good] = 1.0 / T_floor

    ratio = Tt * invT
    np.clip(ratio, -ratio_cap, ratio_cap, out=ratio)

    e1, e2, e3 = coeffs
    theta = 1.0 - (Tm_safe / Tt)
    theta_pos = np.clip(theta, 0.0, None)

    expo = ratio * (e1*theta_pos + e2*np.power(theta_pos,
                    1.5) + e3*np.power(theta_pos, 5.0))
    np.clip(expo, -clip_exp, clip_exp, out=expo)

    P[mask] = pt * np.exp(expo)
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
