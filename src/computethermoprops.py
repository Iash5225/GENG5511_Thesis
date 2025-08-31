import numpy as np

# -----------------------
# Constants / small guards
# -----------------------
R = 8.31451   # J/(mol·K)
_EPS = 1e-4     # generic tiny
_EPS_T = 1e-9      # K
_EPS_V = 1e-9      # cm^3/mol
_MAXEXP = 60.0  # cap exponent arguments to [-60, 60] (e^60 ~ 1e26; safe)

# -----------------------
# Debye helpers (as you had)
# -----------------------


def _safe_expm1(x):
    # accurate for small x, bounded for large x
    return np.expm1(np.clip(x, -_MAXEXP, _MAXEXP))


def Debye3(x: float) -> float:
    if x <= 0.0:
        return 0.0
    n = 400
    a, b = 0.0, x
    h = (b - a) / n

    def f(t):
        if t < 1e-6:
            # t^3/(e^t - 1) ~ t^2 - t^3/2 + t^4/12
            return t*t * (1.0 - 0.5*t + (t*t)/12.0)
        if t > 50.0:
            # asymptotic: t^3 e^{-t}
            return (t**3) * np.exp(-min(t, _MAXEXP))
        return t**3 / _safe_expm1(t)

    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4.0 * f(a + i*h)
    for i in range(2, n, 2):
        s += 2.0 * f(a + i*h)
    return s * h / 3.0


def Debye3D(x: float) -> float:
    # derivative of the integral is the integrand at x
    if x <= 0.0:
        return 0.0
    if x < 1e-6:
        # derivative of series above: t^2 - 3/2 t^3 + ...
        return x*x * (1.0 - 0.5*x)  # good enough in the tiny limit
    if x > 50.0:
        return (x**3) * np.exp(-min(x, _MAXEXP))
    return (x**3) / _safe_expm1(x)

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


    x = Theta[0] / max(T, _EPS_T)
    x = max(x, 1e-8)  # avoid log(0)
    D3 = Debye3(x)
    D3d = Debye3D(x)

    cv = a[0] * R * (D3 - x * D3d)
    dpdt = Gamma[0] * cv

    cv = cv - (aa * R * T / max(Th[0], _EPS_T)) * fz2 * fv
    dpdt = (dpdt / max(v, _EPS_V)) - aa * \
        (cc / max(v00, _EPS_V)) * R * fz1 * fv

    dpdv_safe = np.sign(dpdv) * max(abs(dpdv), _EPS)


    KappaT = -1.0 / (max(v, _EPS_V) * dpdv_safe)
    KAPPA_FLOOR = 1e-9   # MPa^-1  (tune if needed; must be > 0)
    KappaT = float(np.clip(KappaT, 1e-9, 1e+2))     # MPa^-1


    Alpha = -(dpdt / dpdv_safe) / max(v, _EPS_V)
    Alpha = float(np.clip(Alpha, -1e-2, 1e-2))     # K^-1 (solid-scale)

    cp = cv + T * (Alpha**2) * v / max(KappaT, KAPPA_FLOOR)
    cp = float(np.clip(cp, 1e-3, 1e+5))         # J/mol-K
    
    KappaS = -cv / (max(cp, _EPS) * max(v, _EPS_V) * dpdv_safe)
   
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
    Gruneisen = Alpha * v / (max(cv, _EPS) * max(KappaT, _EPS))

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
