import numpy as np


def psub(T, pt, Tt):
    """
    Sublimation pressure (MPa) at temperature T (K).
    pt: triple-point pressure (MPa)
    Tt: triple-point temperature (K)
    """
    T = np.asarray(T, dtype=float)
    T_safe = np.where(T == 0.0, np.finfo(float).tiny, T)  # avoid div/0
    d4, d5, d6 = -10.763, -1.526, -0.4245
    term = 1.0 - T_safe / Tt
    out = pt * np.exp((Tt / T_safe) * (d4*term + d5*term**1.5 + d6*term**5))
    return out.item() if out.ndim == 0 else out


def pmelt(T, pt, Tt):
    """
    Melting pressure (MPa) at temperature T (K).
    Clamps below Tt so fractional powers stay well-defined.
    """
    T = np.asarray(T, dtype=float)
    b1, b2, b3, b4 = 1506.5415, 1.73136, 4677.1597, 0.9849295
    term = np.clip((T / Tt) - 1.0, 0.0, None)  # simple + safe
    out = pt * (1.0 + b1*term**b2 + b3*term**b4)
    return out.item() if out.ndim == 0 else out
