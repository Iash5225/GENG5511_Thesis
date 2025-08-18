import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad
from math import pi

# -------------------------
# Constants / numerics
# -------------------------
R = 8.31451  # J/(mol·K)

# Debye argument clamp and EOS guardrails
X_MAX = 60.0           # clamp for Debye3/3D argument
EPS_T = 1e-9           # K
EPS_V = 1e-6           # cm^3/mol (avoid v<=0)
EPS_D = 1e-14          # generic tiny for divides
V_MIN = 1e-3           # cm^3/mol
V_MAX = 1e3            # cm^3/mol
BIG = 1e300           # fail-safe large

# -------------------------
# Debye helpers (stable)
# -------------------------


def _debye_integrand(t):
    # robust for all t >= 0
    if t < 1e-6:
        # t^3/(e^t-1) ~ t^2 - t^3/2 + t^4/12
        return t*t * (1.0 - 0.5*t + (t*t)/12.0)
    if t > 40.0:
        # e^t dominates denominator
        return t**3 * np.exp(-t)
    return t**3 / np.expm1(t)


def Debye3(x):
    """D_3(x) = ∫_0^x t^3/(e^t-1) dt; handle x<=0 and clamp large x."""
    if x <= 0:
        return 0.0
    x_eff = min(float(x), X_MAX)
    val, _ = quad(_debye_integrand, 0.0, x_eff, limit=200)
    return val


def Debye3D(x):
    """Exact derivative: d/dx Debye3(x) = x^3/(e^x - 1), stabilized."""
    x = float(x)
    if x <= 0:
        return 0.0
    x_eff = min(x, X_MAX)
    return (x_eff**3) / np.expm1(x_eff)

# -------------------------
# Safe math utils
# -------------------------


def _safe_div(a, b, eps=EPS_D):
    b = np.asarray(b, float)
    b_safe = np.where(np.abs(b) < eps, np.sign(b)*eps + (b == 0)*eps, b)
    return np.asarray(a, float) / b_safe


def _clip_positive(x, lo):
    return float(max(x, lo))


def _clip_range(x, lo, hi):
    return float(min(max(x, lo), hi))

# -------------------------
# EOS core
# -------------------------


def evaluate_argon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Compute pcalc and derivatives at (T, v). All pressures MPa, v in cm^3/mol.
    """
    # Guard inputs
    T = _clip_positive(T, EPS_T)
    v = _clip_range(v, V_MIN, V_MAX)
    v00 = _clip_positive(v00, EPS_V)

    # Grüneisen and Debye temperature
    Gamma = g[0] * (v / v00) ** q[0]
    if abs(q[0]) > 1e-12:
        Theta = Th[0] * np.exp((g[0] / q[0]) * (1.0 - (v / v00) ** q[0]))
    else:
        Theta = Th[0] * (v00 / v) ** g[0]

    # Polynomials in ln z
    z = v00 / v
    lnz = np.log(z)
    pj_0 = z * (a1*lnz + a2*lnz**2 + a3*lnz**3)

    # Debye contribution
    x = _clip_range(Theta / T, 0.0, X_MAX)
    pj_1 = a[0] * (_safe_div(R*T*Gamma, v)) * Debye3(x)

    pcalc = pj_0 + pj_1

    # dp/dv pieces (stabilized)
    dpj_0 = -(z / v) * (a1*(1 + lnz) + a2*lnz*(2 + lnz) + a3*lnz**2*(3 + lnz))
    dpj_1 = (pj_1 / _clip_positive(v, EPS_V)) * (q[0] - 1.0) \
        - a[0] * R * Theta * (_safe_div(Gamma, v)**2) * Debye3D(x)
    dpdv = dpj_0 + dpj_1

    # Thermal correction (your original form)
    z_T = T / Th[0]
    fv = np.exp(cc * (v - v00) / v00)
    fz = z_T**4 / (1.0 + bb * z_T**2)
    fz1 = 2.0 * z_T**3 * (2.0 + bb*z_T**2) / (1.0 + bb*z_T**2)**2
    fz2 = 2.0 * z_T**2 * (6.0 + 3.0*bb*z_T**2 + bb**2 *
                          z_T**4) / (1.0 + bb*z_T**2)**3

    correction = aa * (cc / v00) * R * Th[0] * fz * fv
    dcorrection = aa * (cc / v00)**2 * R * Th[0] * fz * fv

    pcalc -= correction
    dpdv -= dcorrection

    return pcalc, dpdv, Theta, Gamma, lnz

# -------------------------
# Robust bracketing
# -------------------------


def _find_bracket(f, v0, grow=2.0, max_expand=18, v_min=V_MIN, v_max=V_MAX):
    """Expand around v0 to find a sign change; fallback to log scan."""
    v0 = _clip_range(v0, v_min, v_max)
    a = max(v0 / grow, v_min)
    b = min(v0 * grow, v_max)
    fa, fb = f(a), f(b)
    k = 0
    while np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) == np.sign(fb) and k < max_expand:
        a = max(a / grow, v_min)
        b = min(b * grow, v_max)
        fa, fb = f(a), f(b)
        k += 1
    if np.isfinite(fa) and np.isfinite(fb) and np.sign(fa) != np.sign(fb):
        return a, b
    # fallback: coarse log scan
    grid = np.geomspace(v_min, v_max, 200)
    vals = np.array([f(x) for x in grid])
    good = np.isfinite(vals)
    idx = np.where(good[:-1] & good[1:] &
                   (np.sign(vals[:-1]) * np.sign(vals[1:]) < 0))[0]
    if idx.size:
        return grid[idx[0]], grid[idx[0] + 1]
    return None, None

# -------------------------
# Public API
# -------------------------


def compute_thermo_props(T, p_MPa, parameters):
    """
    Returns array of 12 properties, or NaNs if no valid root could be bracketed.
    Units:
      T in K, p in MPa, v in cm^3/mol, energies in J/mol.
    """
    # unpack (keep your indexing)
    v00 = parameters[0]
    a1, a2, a3 = parameters[1:4]
    a = [3.0]  # your constant
    Th = parameters[9:15]
    g = parameters[15:21]
    q = parameters[21:27]
    aa, bb, cc, s = parameters[27:31]

    # warm start around v00
    v_guess = _clip_range(v00, V_MIN, V_MAX)

    def pressure_residual(v):
        pcalc, *_ = evaluate_argon(T, v, v00, a, a1,
                                   a2, a3, Th, g, q, aa, bb, cc)
        return p_MPa - pcalc  # root when pcalc == p_MPa

    a_b, b_b = _find_bracket(pressure_residual, v_guess)
    if a_b is None:
        # No valid bracket; return NaNs (caller should penalize)
        return np.full(12, np.nan, dtype=float)

    try:
        sol = root_scalar(pressure_residual, method='brentq',
                          a=a_b, b=b_b, xtol=1e-10, rtol=1e-8, maxiter=200)
        if not sol.converged:
            return np.full(12, np.nan, dtype=float)
        v = _clip_range(sol.root, V_MIN, V_MAX)
    except Exception:
        return np.full(12, np.nan, dtype=float)

    # Evaluate final state
    pcalc, dpdv, Theta, Gamma, lnz = evaluate_argon(
        T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)

    # Derived thermodynamics (defensive)
    x = _clip_range(Theta / _clip_positive(T, EPS_T), 0.0, X_MAX)
    D3 = Debye3(x)
    D3d = Debye3D(x)

    # Debye cv before correction
    cv = a[0] * R * (D3 - x * D3d)
    # correction terms (reuse definitions)
    z_T = T / Th[0]
    fv = np.exp(cc * (v - v00) / v00)
    fz = z_T**4 / (1.0 + bb * z_T**2)
    fz1 = 2.0 * z_T**3 * (2.0 + bb*z_T**2) / (1.0 + bb*z_T**2)**2
    fz2 = 2.0 * z_T**2 * (6.0 + 3.0*bb*z_T**2 + bb**2 *
                          z_T**4) / (1.0 + bb*z_T**2)**3

    cv = cv - (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (Gamma * cv) / _clip_positive(v, EPS_V) - \
        aa * (cc / v00) * R * fz1 * fv

    dpdv_safe = _clip_positive(dpdv, EPS_D)
    cp = cv - T * (dpdt**2) / dpdv_safe

    # Bulk moduli and alpha
    KappaT = -1.0 / (_clip_positive(v, EPS_V) * dpdv_safe)
    cp_safe = _clip_positive(cp, EPS_D)
    KappaS = -cv / (cp_safe * _clip_positive(v, EPS_V) * dpdv_safe)
    Alpha = -dpdt / dpdv_safe / _clip_positive(v, EPS_V)

    # Energetics
    U = v00 * (0.5*a1*lnz**2 + (a2/3.0)*lnz**3 + 0.25*a3*lnz**4)
    U += a[0] * R * T * D3
    U += aa * R * Th[0] * (fz - (T/Th[0]) * fz1) * fv

    # log(1 - e^-x) with stability
    S = a[0] * R * ((4.0/3.0)*D3 - np.log1p(-np.exp(-x)))
    S -= aa * R * fz1 * fv

    A = U - T * S
    H = U + p_MPa * v         # MPa·cm^3/mol ≡ J/mol
    G = H - T * S

    # Grüneisen parameter
    KappaT_safe = _clip_positive(KappaT, EPS_D)
    cv_safe = _clip_positive(cv, EPS_D)
    Gruneisen = Alpha * v / (cv_safe * KappaT_safe)

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
