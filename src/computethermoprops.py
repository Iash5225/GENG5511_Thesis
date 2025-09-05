import numpy as np
import math
# -----------------------
# Constants / small guards
# -----------------------
R = 8.31451
_TINY = 1e-6
_EPS_T = 1e-9
_EPS_V = 1e-9
_MAXEXP = 60.0

# Debye section:
_MACHEPS = 1.11e-16
_XLOW = math.sqrt(8.0*_MACHEPS)
_XUP = -math.log(2.0*_MACHEPS)
# ... etc. (use _MACHEPS here only)
# Debye section (machine eps only used for the Debye series logic)

_XLIM1 = -math.log(2.2251e-308)
_DEBINF = 5.13299112734217e-02
_XLIM2 = (1.0/_DEBINF)**(1.0/3.0) / (2.2251e-308)**(1.0/3.0)

SAFE_DPDV = 1e-6  # MPa / (cm^3/mol)
# _MACHEPS = 1.11e-16  # only for Debye series switches

# -----------------------
# Debye helpers (as you had)
# -----------------------


def _debye3_cheb(x):
    # your current Chebyshev path
    if x < _XLOW:
        return ((x - 7.5)*x + 20.0)/20.0
    t = (x*x/8.0) - 1.0              # map to [-1,1]
    return _cheval(_ADEB3, t) - 0.375*x


def _debye3_asym(x):
    if x > _XLIM2:
        return 0.0
    val = 1.0/(_DEBINF * x**3)
    if x < _XLIM1:
        ex = math.exp(-x)
        if x > _XUP:
            tail = ((x + 3.0)*x + 6.0)*x + 6.0
            val -= 3.0 * tail * ex / (x**3)
        else:
            nexp = int(math.floor(_XLIM1/x))
            rk = float(nexp)
            xk = rk*x
            summ = 0.0
            for _ in range(nexp, 0, -1):
                xki = 1.0/xk
                term = (((6.0*xki + 6.0)*xki + 3.0)*xki + 1.0)/rk
                summ = summ*ex + term
                rk -= 1.0
                xk -= x
            val -= 3.0*summ*ex
    return val


def Debye3(x):
    if x <= 0.0:
        return 0.0
    x = float(x)
    x0, x1 = 3.6, 4.4               # blend window
    if x <= x0:
        return _debye3_cheb(x)
    if x >= x1:
        return _debye3_asym(x)
    # C^1 smoothstep weight
    t = (x - x0) / (x1 - x0)
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    w = t**3 * (10 - 15*t + 6*t**2)   # C^2 smoothstep
    return (1.0 - w)*_debye3_cheb(x) + w*_debye3_asym(x)


def _safe_expm1(x):
    # accurate for small x, bounded for large x
    return np.expm1(np.clip(x, -_MAXEXP, _MAXEXP))


# ---- Chebyshev evaluator (Clenshaw) on [-1,1]

def _cheval(coeffs, t):
    b0 = b1 = 0.0
    for a in coeffs[::-1]:
        b0, b1 = a + 2.0*t*b0 - b1, b0
    return 0.5*(b0 - b1)


# ---- Normalized Debye D3(x)  and its derivative D3'(x)  ----
# Coefficients from your VBA (20-digit), for x in [0,4] mapped by t = (x^2/8 - 1/2) - 1/2
_ADEB3 = np.array([
    2.70773706832744,
    0.340068135211092,
    -1.29451501844409e-02,
    7.96375538017382e-04,
    -5.46360009590824e-05,
    3.92430195988049e-06,
    -2.89403282353860e-07,
    2.17317613962500e-08,
    -1.65420999498000e-09,
    1.27279618920000e-10,
    -9.87963459000000e-12,
    7.72507400000000e-13,
    -6.07797200000000e-14,
    4.80759000000000e-15,
    -3.82040000000000e-16,
    3.04800000000000e-17,
    -2.44000000000000e-18,
    2.00000000000000e-19,
    -2.00000000000000e-20
], dtype=float)

# Precompute derivative-series coefficients (matches VBA loop)
_D = np.zeros_like(_ADEB3)
_D[-1] = 0.0
if _ADEB3.size >= 2:
    _D[-2] = 36.0 * _ADEB3[-1]
    for k in range(_ADEB3.size-2, 0, -1):
        _D[k-1] = _D[k+1] + 2.0*k * _ADEB3[k]

# 18*zeta(4) = pi^4/5; 1/(18*zeta(4)) = this value
# _DEBINF = 5.13299112734217e-02
# # Machine-like constants (double)
# _EPS = 1.11e-16
# _XLOW = math.sqrt(8.0*_EPS)         # ~ sqrt(8*eps)
# _XUP = -math.log(2.0*_EPS)         # switch to exp tail
# _XLIM1 = -math.log(2.2251e-308)     # -log(xmin)
# # XLIM2 = (1/DEBINF)^(1/3) / (xmin)^(1/3)
# _XLIM2 = (1.0/_DEBINF)**(1.0/3.0) / (2.2251e-308)**(1.0/3.0)


# def Debye3(x):
#     """Normalized Debye D3(x) = 3/x^3 ∫ t^3/(e^t-1) dt  (0 ≤ D3 ≤ 1)."""
#     if x <= 0.0:
#         return 0.0
#     if x <= 4.0:
#         if x < _XLOW:
#             # small-x series: 1 - 3x/8 + x^2/20
#             return ((x - 7.5)*x + 20.0)/20.0
#         t = (x*x/8.0) - 0.5 - 0.5   # map to [-1,1]
#         return _cheval(_ADEB3, t) - 0.375*x
#     # x > 4:
#     if x > _XLIM2:
#         return 0.0
#     val = 1.0/(_DEBINF * x**3)     # = 18*zeta(4)/x^3
#     if x < _XLIM1:
#         ex = math.exp(-x)
#         if x > _XUP:
#             # 3*exp(-x)(x^3+3x^2+6x+6)/x^3 term
#             tail = ((x + 3.0)*x + 6.0)*x + 6.0
#             val -= 3.0 * tail * ex / (x**3)
#         else:
#             # refined Poisson summation tail
#             nexp = int(math.floor(_XLIM1/x))
#             rk = float(nexp)
#             xk = rk*x
#             summ = 0.0
#             for _ in range(nexp, 0, -1):
#                 xki = 1.0/xk
#                 term = (((6.0*xki + 6.0)*xki + 3.0)*xki + 1.0)/rk
#                 summ = summ*ex + term
#                 rk -= 1.0
#                 xk -= x
#             val -= 3.0*summ*ex
#     return val


# def Debye3D(x):
#     """Derivative D3'(x) of the normalized Debye function D3(x)."""
#     if x <= 0.0:
#         return 0.0
#     if x <= 4.0:
#         if x < _XLOW:
#             return (4.0*x - 15.0)/40.0
#         t = (x*x/8.0) - 0.5 - 0.5
#         return 0.25*x*_cheval(_D, t) - 0.375
#     # x > 4:
#     if x > _XLIM2:
#         return 0.0
#     val = -3.0/(_DEBINF * x**4)    # derivative of 18*zeta(4)/x^3
#     if x < _XLIM1:
#         ex = math.exp(-x)
#         if x > _XUP:
#             tail = ((((x + 3.0)*x + 9.0)*x + 18.0)*x + 18.0) / (x**4)
#             val += 3.0*tail*ex
#         else:
#             nexp = int(math.floor(_XLIM1/x))
#             rk = float(nexp)
#             xk = rk*x
#             summ = 0.0
#             for _ in range(nexp, 0, -1):
#                 xki = 1.0/xk
#                 term = (((18.0*xki + 18.0)*xki + 9.0)*xki + 3.0)*xki + 1.0
#                 summ = summ*ex + term
#                 rk -= 1.0
#                 xk -= x
#             val += 3.0*summ*ex
#     return val


def _debye3d_cheb(x):
    if x < _XLOW:
        return (4.0*x - 15.0)/40.0
    t = (x*x/8.0) - 1.0
    return 0.25*x*_cheval(_D, t) - 0.375


def _debye3d_asym(x):
    if x > _XLIM2:
        return 0.0
    val = -3.0/(_DEBINF * x**4)
    if x < _XLIM1:
        ex = math.exp(-x)
        if x > _XUP:
            tail = ((((x + 3.0)*x + 9.0)*x + 18.0)*x + 18.0) / (x**4)
            val += 3.0*tail*ex
        else:
            nexp = int(math.floor(_XLIM1/x))
            rk = float(nexp)
            xk = rk*x
            summ = 0.0
            for _ in range(nexp, 0, -1):
                xki = 1.0/xk
                term = (((18.0*xki + 18.0)*xki + 9.0)*xki + 3.0)*xki + 1.0
                summ = summ*ex + term
                rk -= 1.0
                xk -= x
            val += 3.0*summ*ex
    return val


def Debye3D(x):
    if x <= 0.0:
        return 0.0
    x = float(x)
    x0, x1 = 3.6, 4.4
    if x <= x0:
        return _debye3d_cheb(x)
    if x >= x1:
        return _debye3d_asym(x)
    t = (x - x0) / (x1 - x0)
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    w = t**3 * (10 - 15*t + 6*t**2)   # C^2 smoothstep
    return (1.0 - w)*_debye3d_cheb(x) + w*_debye3d_asym(x)


# -----------------------
# Core EOS pieces (generalized)
# -----------------------


def evaluate_gas(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Generalized version of MATLAB's evaluateArgon:
      inputs in K (T), cm^3/mol (v, v00), MPa (outputs)
    returns: pcalc, dpdv, Theta[0..0], Gamma[0..0], lnz
    """
    # Grüneisen and Debye temperature for the first (active) branch
    # Grüneisen and Debye temperature (safe exponent)
    Gamma1 = g[0] * (v / v00)**q[0]
    if abs(q[0]) > 1e-10:
        arg = (g[0] / q[0]) * (1.0 - (v / v00)**q[0])
        arg = np.clip(arg, -_MAXEXP, _MAXEXP)
        Theta1 = Th[0] * np.exp(arg)
    else:
        Theta1 = Th[0] * (v00 / v)**g[0]

    # elastic polynomial in ln z
    z = v00 / max(v, _EPS_V)
    lnz = np.log(z)

    pj_0 = z * (a1 * lnz + a2 * lnz**2 + a3 * lnz**3)
    pj_1 = a[0] * (R * T * Gamma1 / max(v, _EPS_V)) * \
        Debye3(Theta1 / max(T, _EPS_T))

    pcalc = pj_0 + pj_1

    # derivative wrt v
    dpj_0 = -(z / max(v, _EPS_V)) * (a1*(1.0 + lnz) +
                                     a2*lnz*(2.0 + lnz) + a3*lnz**2*(3.0 + lnz))
    dpj_1 = (pj_1 / max(v, _EPS_V)) * (q[0] - 1.0) \
        - a[0] * R * Theta1 * (Gamma1 / max(v, _EPS_V))**2 * \
        Debye3D(Theta1 / max(T, _EPS_T))
    dpdv = dpj_0 + dpj_1

    # thermal correction
    zT = T / max(Th[0], _EPS_T)

    # thermal correction (safe exponent)
    argV = cc * (v - v00) / max(v00, _EPS_V)
    argV = np.clip(argV, -_MAXEXP, _MAXEXP)
    exp_factor = np.exp(argV)
    fz = zT**4 / (1.0 + bb * zT**2)

    pcalc -= aa * (cc / max(v00, _EPS_V)) * R * Th[0] * fz * exp_factor
    dpdv -= aa * (cc / max(v00, _EPS_V))**2 * R * Th[0] * fz * exp_factor

    Theta = np.array([Theta1], dtype=float)
    Gamma = np.array([Gamma1], dtype=float)
    return pcalc, dpdv, Theta, Gamma, lnz

# -----------------------
# Public API (matches fitter)
# -----------------------


def compute_thermo_props(T: float, p: float, parameters):
    """
    EXACT signature the costing function expects.
    T in K, p in MPa. Returns numpy array of length 12:
      [v, KappaT, KappaS, Alpha, cp, cv, Gamma, U, S, A, H, G]
    """
    # ---- unpack parameters in the SAME layout your fitter uses ----
    # indices (from your earlier code):
    #  0            : v00
    #  1:4          : a1, a2, a3
    #  4:9          : (unused 5 slots in this EOS)
    #  9:15         : Th[0..5]
    # 15:21         : g[0..5]
    # 21:27         : q[0..5]
    # 27:31         : aa, bb, cc, s
    v00 = float(parameters[0])
    a1, a2, a3 = [float(x) for x in parameters[1:4]]
    a = [3.0]  # MATLAB constant
    Th = np.array(parameters[9:15],  dtype=float)
    g = np.array(parameters[15:21], dtype=float)
    q = np.array(parameters[21:27], dtype=float)
    aa, bb, cc, s = [float(x) for x in parameters[27:31]]

    # ---- Newton-like solve for v (to mirror MATLAB computeThermoPropsForResult) ----
    T = float(T)
    p = float(p)

    # choose a sane box for solid molar volume (tune if needed)
    VMIN, VMAX = 5.0, 40.0
    MAX_NEWTON = 60
    BACKTRACK_MIN = 1e-3
    eps_newton = 1e-10 * max(abs(p), 1.0)


    v = float(np.clip(26.0 * np.exp(-p / 25000.0), VMIN, VMAX))
    eps_newton = 1e-8 * max(abs(p), 1.0)  # slightly looser to avoid over-iterating


    for _ in range(MAX_NEWTON):
        pcalc, dpdv, Theta, Gamma, lnz = evaluate_gas(
            T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc
        )

        # if evaluation is bad at current v, retreat toward v00 and try again
        if not (np.isfinite(pcalc) and np.isfinite(dpdv)):
            v = float(np.clip(0.5*(v + v00), VMIN, VMAX))
            continue

        perror = p - pcalc
        if abs(perror) <= eps_newton:
            break

        # bounded slope + damped step
        dpdv_safe = np.sign(dpdv) * max(abs(dpdv), 1e-4)
        vstep = float(np.clip(perror / dpdv_safe, -0.2*v00, 0.2*v00))

        # try full step; if it fails, backtrack with halving
        t = 1.0
        accepted = False
        while t >= BACKTRACK_MIN:
            v_candidate = float(np.clip(v + t*vstep, VMIN, VMAX))
            p_try, dpdv_try, Theta_try, Gamma_try, lnz_try = evaluate_gas(
                T, v_candidate, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc
            )
            if np.isfinite(p_try) and np.isfinite(dpdv_try):
                # accept the candidate
                v = v_candidate
                pcalc, dpdv, Theta, Gamma, lnz = p_try, dpdv_try, Theta_try, Gamma_try, lnz_try
                accepted = True
                break
            t *= 0.5

        if not accepted:
            # fall back toward a safe center and continue
            v = float(np.clip(0.5*(v + v00), VMIN, VMAX))

    # ---- Thermodynamic properties (same algebra/order as your fitter) ----
    zT = T / max(Th[0], _EPS_T)
    fz = zT**4 / (1.0 + bb * zT**2)
    fz1 = 2.0 * zT**3 * (2.0 + bb * zT**2) / (1.0 + bb * zT**2)**2
    fz2 = 2.0 * zT**2 * (6.0 + 3.0*bb*zT**2 + bb**2 *
                         zT**4) / (1.0 + bb * zT**2)**3
    fv = np.exp(cc * (v - v00) / max(v00, _EPS_V))
    dpdv_den = dpdv if abs(dpdv) > SAFE_DPDV else math.copysign(
        SAFE_DPDV, dpdv if dpdv != 0 else -1.0)

    x = Theta[0] / max(T, _EPS_T)
    x = max(x, 1e-8)  # avoid log(0)
    D3 = Debye3(x)
    D3d = Debye3D(x)

    cv = a[0] * R * (D3 - x * D3d)
    dpdt = Gamma[0] * cv

    cv = cv - (aa * R * T / max(Th[0], _EPS_T)) * fz2 * fv
    
    dpdt = (dpdt / max(v, _EPS_V)) - aa * \
        (cc / max(v00, _EPS_V)) * R * fz1 * fv

    # --- thermo closures from slopes ---
    dpdv_den = dpdv if abs(dpdv) > SAFE_DPDV else math.copysign(
        SAFE_DPDV, dpdv if dpdv != 0 else -1.0
    )

    # exact cp via identity (decoupled from any kappa_T floor)
    cp = cv - T * (dpdt * dpdt) / dpdv_den

    # now compute reportable Alpha, KappaT using the same denominator
    Alpha = -(dpdt / dpdv_den) / max(v, _EPS_V)
    KappaT = -1.0 / (max(v, _EPS_V) * dpdv_den)
    # optional display-floor ONLY (do not feed back into cp)
    KappaT = float(np.clip(KappaT, 1e-9, 1e+2))

    # for KappaS we only need a safe dpdv in the denominator
    dpdv_safe = np.sign(dpdv) * max(abs(dpdv), _TINY)
    KappaS = -cv / (max(cp, _TINY) * max(v, _EPS_V) * dpdv_safe)

   
   # entropy: stable log1p(1 - e^{-x})
    expmx = np.exp(-min(x, _MAXEXP))
    expmx = min(expmx, 1.0 - 1e-16)  # avoid log(0)

    # --- robust cp via identity ---

    # cp = cv - T * (dpdt**2) / dpdv_safe

    lnz = np.log(v00 / max(v, _EPS_V))
    U = v00 * (0.5*a1*lnz**2 + (a2/3.0)*lnz**3 + 0.25*a3*lnz**4)
    U += a[0] * R * T * D3
    U += aa * R * Th[0] * (fz - zT * fz1) * fv


    S = a[0] * R * ((4.0/3.0)*D3 - np.log1p(-expmx))
    S -= aa * R * fz1 * fv

    A = U - T * S
    H = U + p * v                # MPa·cm^3/mol == J/mol
    G = H - T * S

    # Grüneisen parameter (thermal)
    Gruneisen = Alpha * v / (max(cv, _TINY) * max(KappaT, _TINY))

    props = np.zeros(12, dtype=float)
    props[0] = v
    props[1] = KappaT
    props[2] = KappaS
    props[3] = Alpha
    props[4] = cp
    props[5] = cv
    props[6] = Gruneisen
    props[7] = U
    props[8] = S
    props[9] = A
    props[10] = H
    props[11] = G
    return props
