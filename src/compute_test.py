# eos_argon_vba_equiv.py
import math
import numpy as np

# ----- constants & guards (VBA uses MPa, cm^3/mol, K, J/mol) -----
R = 8.31451
_MAXEXP = 60.0          # guard exponentials to avoid overflow
_EPSV = 1e-12         # volume guard
_EPS_DP = 1e-12         # derivative guard

# ----- Debye helpers (match VBA Debye3 / Debye3D behavior) -----


def _safe_expm1(x: float) -> float:
    return math.expm1(max(min(x, _MAXEXP), -_MAXEXP))


def Debye3(x: float) -> float:
    if x <= 0.0:
        return 0.0
    n = 400
    a, b = 0.0, x
    h = (b - a) / n

    def f(t: float) -> float:
        if t < 1e-6:
            # t^3/(e^t-1) ~ t^2 - t^3/2 + t^4/12
            return t*t * (1.0 - 0.5*t + (t*t)/12.0)
        if t > 50.0:
            # asymptotic tail is tiny
            return (t**3) * math.exp(-min(t, _MAXEXP))
        return (t**3) / _safe_expm1(t)

    s = f(a) + f(b)
    for k in range(1, n, 2):
        s += 4.0 * f(a + k*h)
    for k in range(2, n, 2):
        s += 2.0 * f(a + k*h)
    return s * h / 3.0


def Debye3D(x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x < 1e-6:
        # derivative of the small-t series above
        return x*x * (1.0 - 0.5*x)
    if x > 50.0:
        return (x**3) * math.exp(-min(x, _MAXEXP))
    return (x**3) / _safe_expm1(x)


# ----- fixed “paper” parameters (exactly as in your VBA) -----
V00 = 22.555
A1, A2, A3 = 2656.5, 7298.0, 10.0
A = 3.0

TH1 = 86.44
G1 = 2.68
Q1 = 0.0024

AA = 0.0128
BB = 0.388
CC = 7.85
# s (entropy ref) is not used in these two routines

# ----- core evaluator (VBA Evaluate: pcalc, dpdv, Theta, Gamma, lnz) -----


def _evaluate(T: float, v: float):
    # Grüneisen & Debye temperature
    Gamma = G1 * (v / V00)**Q1
    if abs(Q1) > 1e-10:
        arg = (G1 / Q1) * (1.0 - (v / V00)**Q1)
        arg = max(min(arg, _MAXEXP), -_MAXEXP)
        Theta = TH1 * math.exp(arg)
    else:
        Theta = TH1 * (V00 / v)**G1

    z = V00 / max(v, _EPSV)
    lnz = math.log(z)

    pj0 = z * (A1*lnz + A2*lnz*lnz + A3*lnz**3)
    pj1 = A * (R * T * Gamma / max(v, _EPSV)) * Debye3(Theta / max(T, 1e-12))
    pcalc = pj0 + pj1

    dpj0 = -(z / max(v, _EPSV)) * (
        A1*(1.0 + lnz) +
        A2*lnz*(2.0 + lnz) +
        A3*(lnz**2)*(3.0 + lnz)
    )
    dpj1 = (pj1 / max(v, _EPSV)) * (Q1 - 1.0) \
        - A * R * Theta * (Gamma / max(v, _EPSV))**2 * \
        Debye3D(Theta / max(T, 1e-12))
    dpdv = dpj0 + dpj1

    # anharmonic correction (exact VBA form)
    zT = T / TH1
    fz = (zT**4) / (1.0 + BB*zT*zT)
    ev = math.exp(
        max(min(CC * (v - V00) / max(V00, _EPSV), _MAXEXP), -_MAXEXP))
    pcalc -= AA * (CC / max(V00, _EPSV)) * R * TH1 * fz * ev
    dpdv -= AA * (CC / max(V00, _EPSV))**2 * R * TH1 * fz * ev

    return pcalc, dpdv, Theta, Gamma, lnz

# ----- the VBA-like interface (exact same i mapping) -----


def props_vba(T: float, X: float, i: int) -> float:
    """
    Exact port of VBA 'Props':
      if i > 0: X = p (MPa)  -> solve v(T,p)
      if i <= 0: X = v (cm^3/mol) -> compute p(T,v)
    Return value chosen by |i| (same mapping as VBA Select Case).
    """
    # branch: solve for v or not
    if i > 0:
        p = float(X)
        eps = 1e-10 * (p if p > 1.0 else 1.0)
        v = 26.0 * math.exp(-p / 25000.0)  # VBA initial guess
        for _ in range(100):
            pcalc, dpdv, Theta, Gamma, lnz = _evaluate(T, v)
            perror = p - pcalc
            vstep = perror / max(dpdv, _EPS_DP)  # keep the sign of dpdv
            # limit step to 20% v00 if it exceeds 10% (matches the VBA code)
            if abs(vstep) > 0.1*V00:
                vstep = math.copysign(0.2*V00, vstep)
            v += vstep
            if abs(perror) <= eps:
                break
            if v <= 0.0:
                v = 100.0
                break
        # keep last Theta/Gamma/lnz from loop
    else:
        v = float(X)
        pcalc, dpdv, Theta, Gamma, lnz = _evaluate(T, v)
        p = pcalc

    # anharmonic pieces
    z = T / TH1
    fz = (z**4) / (1.0 + BB*z*z)
    fz1 = 2.0*z**3 * (2.0 + BB*z*z) / (1.0 + BB*z*z)**2
    fz2 = 2.0*z*z * (6.0 + 3.0*BB*z*z + (BB**2)*(z**4)) / (1.0 + BB*z*z)**3
    fv = math.exp(max(min(CC*(v - V00)/max(V00, _EPSV), _MAXEXP), -_MAXEXP))

    # cv, dp/dT|v, cp (VBA formulas)
    cv = A * R * (Debye3(Theta / T) - (Theta / T) * Debye3D(Theta / T))
    dpdt = Gamma * cv
    cv = cv - (AA * R * T / TH1) * fz2 * fv
    dpdt = (dpdt / v) - AA * (CC / V00) * R * fz1 * fv
    cp = cv - T * (dpdt**2) / max(dpdv, _EPS_DP)

    # compressibilities & expansivity
    KappaT = -1.0 / (v * max(dpdv, _EPS_DP))
    KappaS = -cv / (cp * v * max(dpdv, _EPS_DP))
    Alpha = -(dpdt / max(dpdv, _EPS_DP)) / v

    # energies
    lnz = math.log(V00 / max(v, _EPSV))
    U = V00 * (0.5*A1*lnz*lnz + A2*(lnz**3)/3.0 + 0.25*A3*(lnz**4))
    U += A * R * T * Debye3(Theta / T)
    U += AA * R * TH1 * (fz - z*fz1) * fv

    Xi = Theta / T
    if Xi < math.log(math.sqrt(1.7976931348623157e308)) - 1.0:  # VBA's Xmax logic
        S = A * R * ((4.0/3.0)*Debye3(Xi) - math.log(1.0 - math.exp(-Xi)))
    else:
        S = 0.0
    S -= AA * R * fz1 * fv

    H = U + (p if i > 0 else pcalc) * v  # VBA uses the input p if i>0
    G = H - T * S
    Ahelm = U - T * S
    Gruneisen = Alpha * v / (cv * KappaT)

    # return selection exactly like the VBA Select Case
    k = abs(int(i))
    if k == 1:
        return U
    elif k == 2:
        return S
    elif k == 3:
        return (v if i > 0 else p)
    elif k == 4:
        return G
    elif k == 5:
        return H
    elif k == 6:
        return cv
    elif k == 7:
        return cp
    elif k == 8:
        return Alpha
    elif k == 9:
        return KappaT
    elif k == 10:
        return KappaS
    elif k == 11:
        return Gruneisen
    elif k == 12:
        return Ahelm
    else:
        raise ValueError("Unsupported property code i")

# ----- your fitter-friendly API (12-vector) -----


def compute_thermo_props(T: float, p: float, parameters=None):
    """
    Returns [v, KappaT, KappaS, Alpha, cp, cv, Gamma, U, S, A, H, G]
    Uses the same model/units as the VBA code above (MPa, cm^3/mol, J/mol).
    'parameters' is optional; if provided and 31-long, the Argon constants
    will be read from it in the same slots you’ve documented.
    """
    # If a parameter vector is supplied, read it; else, use fixed paper values
    if parameters is not None and len(parameters) >= 31:
        v00 = float(parameters[0])
        a1, a2, a3 = map(float, parameters[1:4])
        th1 = float(parameters[9])
        g1 = float(parameters[15])
        q1 = float(parameters[21])
        aa, bb, cc = map(float, parameters[27:30])
    else:
        v00, a1, a2, a3 = V00, A1, A2, A3
        th1, g1, q1 = TH1, G1, Q1
        aa,  bb,  cc = AA,  BB,  CC

    # small local closures to use those read values
    def eval_local(T_, v_):
        Gamma = g1 * (v_ / v00)**q1
        if abs(q1) > 1e-10:
            arg = (g1 / q1) * (1.0 - (v_ / v00)**q1)
            arg = max(min(arg, _MAXEXP), -_MAXEXP)
            Theta = th1 * math.exp(arg)
        else:
            Theta = th1 * (v00 / v_)**g1
        z = v00 / max(v_, _EPSV)
        lnz = math.log(z)
        pj0 = z * (a1*lnz + a2*lnz*lnz + a3*lnz**3)
        pj1 = A * (R * T_ * Gamma / max(v_, _EPSV)) * \
            Debye3(Theta / max(T_, 1e-12))
        pcalc = pj0 + pj1
        dpj0 = -(z / max(v_, _EPSV)) * (a1*(1.0+lnz) +
                                        a2*lnz*(2.0+lnz) + a3*(lnz**2)*(3.0+lnz))
        dpj1 = (pj1 / max(v_, _EPSV)) * (q1 - 1.0) \
            - A * R * Theta * (Gamma / max(v_, _EPSV))**2 * \
            Debye3D(Theta / max(T_, 1e-12))
        dpdv = dpj0 + dpj1
        zT = T_ / th1
        fz = (zT**4) / (1.0 + bb*zT*zT)
        ev = math.exp(
            max(min(cc*(v_ - v00)/max(v00, _EPSV), _MAXEXP), -_MAXEXP))
        pcalc -= aa * (cc / max(v00, _EPSV)) * R * th1 * fz * ev
        dpdv -= aa * (cc / max(v00, _EPSV))**2 * R * th1 * fz * ev
        return pcalc, dpdv, Theta, Gamma, lnz

    # Newton solve for v at (T,p)
    eps = 1e-10 * (abs(p) if abs(p) > 1.0 else 1.0)
    v = 26.0 * math.exp(-p / 25000.0)
    for _ in range(100):
        pcalc, dpdv, Theta, Gamma, lnz = eval_local(T, v)
        perror = p - pcalc
        vstep = perror / (dpdv if abs(dpdv) >
                          _EPS_DP else math.copysign(_EPS_DP, dpdv))
        if abs(vstep) > 0.1*v00:
            vstep = math.copysign(0.2*v00, vstep)
        v += vstep
        if abs(perror) <= eps:
            break
        if v <= 0.0:
            v = 100.0
            break

    # properties (same algebra as VBA)
    z = T / th1
    fz = (z**4) / (1.0 + bb*z*z)
    fz1 = 2.0*z**3 * (2.0 + bb*z*z) / (1.0 + bb*z*z)**2
    fz2 = 2.0*z*z * (6.0 + 3.0*bb*z*z + (bb**2)*(z**4)) / (1.0 + bb*z*z)**3
    fv = math.exp(max(min(cc*(v - v00)/max(v00, _EPSV), _MAXEXP), -_MAXEXP))

    cv = A * R * (Debye3(Theta / T) - (Theta / T) * Debye3D(Theta / T))
    dpdt = Gamma * cv
    cv = cv - (aa * R * T / th1) * fz2 * fv
    dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv

    dpdv_safe = (dpdv if abs(dpdv) > _EPS_DP else math.copysign(_EPS_DP, dpdv))
    KappaT = -1.0 / (v * dpdv_safe)
    Alpha = -(dpdt / dpdv_safe) / v
    cp = cv - T * (dpdt**2) / dpdv_safe
    KappaS = -cv / (cp * v * dpdv_safe)

    lnz = math.log(v00 / max(v, _EPSV))
    U = v00 * (0.5*a1*lnz*lnz + a2*(lnz**3)/3.0 + 0.25*a3*(lnz**4))
    U += A * R * T * Debye3(Theta / T)
    U += aa * R * th1 * (fz - z*fz1) * fv

    Xi = Theta / T
    S = A * R * ((4.0/3.0)*Debye3(Xi) - math.log(1.0 -
                 math.exp(-Xi))) if Xi < 36.0 else 0.0
    S -= aa * R * fz1 * fv

    Ahelm = U - T * S
    H = U + p * v
    G = H - T * S
    Gr = Alpha * v / (cv * KappaT)

    return np.array([v, KappaT, KappaS, Alpha, cp, cv, Gr, U, S, Ahelm, H, G], dtype=float)
