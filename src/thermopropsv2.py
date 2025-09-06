from debye import Debye
from debye3 import Debye3
from debye3D import Debye3D
from cheval import cheval

import numpy as np
# -----------------------
# Public API (matches fitter)
# -----------------------

# --- assumes you already defined cheval(), Debye3(), Debye3D() above ---

R = 8.31451  # J/(mol K) ~ MPa*cm^3/(mol K) when dividing by v [cm^3/mol]

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

    # ---- Newton-like solve for v (to mirror MATLAB computeThermoProps) ----
    T = float(T)
    p = float(p)
    max_iter = 100
    rtol_p = 1e-10

    # Newton solve for v at (T,p)
    # NOTE: your initial guess v = 26*exp(-p/25000). Units match your code.
    v = 26.0 * np.exp(-p / 25000.0)  # [cm^3/mol], same as MATLAB
    # Tolerance equals 1e-10 * p in MATLAB; if p=0, avoid zero tolerance
    abs_tol = max(rtol_p * max(1.0, abs(p)), 1e-12)

    for _ in range(max_iter):
        pcalc, dpdv, Theta, Gamma, lnz = evaluate_argon(
            T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc
        )
        perror = p - pcalc
        vstep = perror / dpdv

        # Step cap |Δv| ≤ 0.2*v00 (your code actually caps at 0.1 then forces 0.2*v00)
        # We'll mirror the final cap you enforce.
        cap = 0.2 * v00
        if abs(vstep) > cap:
            vstep = np.sign(vstep) * cap

        v = v + vstep

        if abs(perror) <= abs_tol:
            break
        if v <= 0:
            v = 100.0
            break

    # Debye/anharmonic helper factors
    z = T / Th[0]
    fz = z**4 / (1.0 + bb*z**2)
    fz1 = 2.0 * z**3 * (2.0 + bb*z**2) / (1.0 + bb*z**2)**2
    fz2 = 2.0 * z**2 * (6.0 + 3.0*bb*z**2 + bb**2*z**4) / (1.0 + bb*z**2)**3
    fv = np.exp(cc * (v - v00) / v00)

    # Debye arguments
    x = Theta[0] / T
    D3 = Debye3(x)
    D3d = Debye3D(x)

    # cv and dp/dT before corrections
    cv_lattice = a[0] * R * (D3 - x * D3d)
    dpdt_base = Gamma[0] * cv_lattice

    # Apply anharmonic corrections (as in your MATLAB)
    cv = cv_lattice - (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (dpdt_base / v) - aa * (cc / v00) * R * fz1 * fv

    # cp from identity
    cp = cv - T * (dpdt**2) / dpdv

    # Compressibilities & expansion
    KappaT = -1.0 / (v * dpdv)
    KappaS = -cv / (cp * v * dpdv)
    Alpha = -(dpdt / dpdv) / v

    # Energies
    # polynomial in ln z:
    U = v00 * (0.5*a1*lnz**2 + (a2/3.0)*lnz**3 + 0.25*a3*lnz**4)
    U += a[0] * R * T * D3
    U += aa * R * Th[0] * (fz - z*fz1) * fv

    # Entropy (use log1p for stability: log(1 - e^{-x}))
    # guard tiny exp(-x)
    log_term = np.log1p(-np.exp(-x)) if x < 40.0 else np.log1p(-0.0)
    S = a[0] * R * ((4.0/3.0)*D3 - log_term)
    S -= aa * R * fz1 * fv

    A = U - T * S
    H = U + p * v
    G = H - T * S

    Gruneisen = Alpha * v / (cv * KappaT)

    props = dict(
        v=v, KappaT=KappaT, KappaS=KappaS, Alpha=Alpha,
        cp=cp, cv=cv, Gruneisen=Gruneisen,
        U=U, S=S, A=A, H=H, G=G
    )
    props_vec = np.array([
        v, KappaT, KappaS, Alpha, cp, cv, Gruneisen, U, S, A, H, G
    ], dtype=float)
    # return ONLY the vector (order matches your MATLAB)
    return np.array([v, KappaT, KappaS, Alpha, cp, cv, Gruneisen, U, S, A, H, G], dtype=float)


    # return props, props_vec


def evaluate_argon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Python port of evaluateArgon(T, v, v00, ...)

    Returns
    -------
    pcalc : float   [MPa]
    dpdv  : float   [MPa / (cm^3/mol)]
    Theta : np.ndarray shape (1,)
    Gamma : np.ndarray shape (1,)
    lnz   : float   (ln(v00/v))
    """
    # j=1 mode only in your code
    Gamma = np.zeros(1, dtype=float)
    Theta = np.zeros(1, dtype=float)

    # Gamma(v) and Theta(v) with smooth branch for |q|≈0
    Gamma[0] = g[0] * (v / v00)**q[0]
    if abs(q[0]) > 1e-10:
        Theta[0] = Th[0] * np.exp((g[0] / q[0]) * (1.0 - (v / v00)**q[0]))
    else:
        Theta[0] = Th[0] * (v00 / v)**g[0]

    # Elastic/volumetric part
    z = v00 / v
    lnz = np.log(z)

    pj_0 = z * (a1*lnz + a2*lnz**2 + a3*lnz**3)
    # Lattice pressure piece
    x = Theta[0] / T
    pj_1 = a[0] * (R * T * Gamma[0] / v) * Debye3(x)

    pcalc = pj_0 + pj_1

    # dp/dv pieces
    dpj_0 = -(z / v) * (a1*(1.0 + lnz) + a2*lnz *
                        (2.0 + lnz) + a3*lnz**2*(3.0 + lnz))
    dpj_1 = (pj_1 / v) * (q[0] - 1.0) - a[0] * R * \
        Theta[0] * (Gamma[0] / v)**2 * Debye3D(x)

    dpdv = dpj_0 + dpj_1

    # Apply (aa,bb,cc) correction to p and dp/dv
    zT = T / Th[0]
    fT = (zT**4) / (1.0 + bb*zT**2)
    exp_factor = np.exp(cc * (v - v00) / v00)

    correction_p = aa * (cc / v00) * R * Th[0] * fT * exp_factor
    pcalc = pcalc - correction_p

    dpdv_correction = aa * (cc / v00)**2 * R * Th[0] * fT * exp_factor
    dpdv = dpdv - dpdv_correction

    return pcalc, dpdv, Theta, Gamma, lnz
