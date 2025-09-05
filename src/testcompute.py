# compute_thermo_props_for_result.py
import math
import numpy as np

# ---- constants ----
R = 8.31451  # J/(mol·K)

# ---- Debye helpers (match MATLAB usage) ----


def Debye3(x: float) -> float:
    """∫_0^x t^3/(e^t-1) dt  (Simpson integration; same role as MATLAB Debye3)."""
    if x <= 0.0:
        return 0.0
    n = 400  # even
    a, b = 0.0, x
    h = (b - a) / n

    def f(t: float) -> float:
        if t < 1e-6:
            # t^3/(e^t - 1) ≈ t^2 - t^3/2 + t^4/12
            return t * t * (1.0 - 0.5 * t + (t * t) / 12.0)
        return (t ** 3) / (math.exp(t) - 1.0)

    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4.0 * f(a + i * h)
    for i in range(2, n, 2):
        s += 2.0 * f(a + i * h)
    return s * h / 3.0


def Debye3D(x: float) -> float:
    """d/dx Debye3(x) = x^3/(e^x - 1)  (same as MATLAB Debye3D)."""
    if x <= 0.0:
        return 0.0
    return (x ** 3) / (math.exp(x) - 1.0)


# ---- core evaluator (exact MATLAB formulas) ----
def evaluate_argon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Returns: pcalc, dpdv, Theta, Gamma, lnz
    Units: T [K], v [cm^3/mol], p [MPa]
    """
    # j=1 only
    Gamma = np.zeros(1)
    Theta = np.zeros(1)

    Gamma[0] = g[0] * (v / v00) ** q[0]
    if abs(q[0]) > 1e-10:
        Theta[0] = Th[0] * math.exp((g[0] / q[0]) * (1.0 - (v / v00) ** q[0]))
    else:
        Theta[0] = Th[0] * (v00 / v) ** g[0]

    z = v00 / v
    lnz = math.log(z)

    pj_0 = z * (a1 * lnz + a2 * lnz ** 2 + a3 * lnz ** 3)
    pj_1 = a[0] * (R * T * Gamma[0] / v) * Debye3(Theta[0] / T)

    pcalc = pj_0 + pj_1

    dpj_0 = -(z / v) * (a1 * (1.0 + lnz) + a2 * lnz *
                        (2.0 + lnz) + a3 * (lnz ** 2) * (3.0 + lnz))
    dpj_1 = (pj_1 / v) * (q[0] - 1.0) - a[0] * R * \
        Theta[0] * (Gamma[0] / v) ** 2 * Debye3D(Theta[0] / T)

    dpdv = dpj_0 + dpj_1

    # corrections
    zT = T / Th[0]
    exp_factor = math.exp(cc * (v - v00) / v00)
    correction = aa * (cc / v00) * R * \
        Th[0] * (zT ** 4 / (1.0 + bb * zT ** 2)) * exp_factor

    pcalc = pcalc - correction
    dpdv = dpdv - aa * (cc / v00) ** 2 * R * \
        Th[0] * (zT ** 4 / (1.0 + bb * zT ** 2)) * exp_factor

    return pcalc, dpdv, Theta, Gamma, lnz


def compute_thermo_props_for_result(T, p):
    """
    Exact port of your MATLAB computeThermoPropsForResult(T,p).

    Inputs
      T : K
      p : MPa

    Returns np.array of length 12:
      [v, KappaT, KappaS, Alpha, cp, cv, Gruneisen, U, S, A, H, G]
    """
    # constants / parameters (exact values from your MATLAB)
    v00 = 22.555
    a1, a2, a3 = 2656.5, 7298.0, 10.0
    a = np.array([3.0])  # only a(1)=3 used

    Th = np.zeros(1)
    g = np.zeros(1)
    q = np.zeros(1)
    Th[0] = 86.44
    g[0] = 2.68
    q[0] = 0.0024

    aa, bb, cc = 0.0128, 0.388, 7.85

    # solve for v at (T,p)
    eps = 1e-10 * p
    v = 26.0 * math.exp(-p / 25000.0)  # initial guess

    for _ in range(100):
        pcalc, dpdv, Theta, Gamma, lnz = evaluate_argon(
            T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)
        perror = p - pcalc
        vstep = perror / dpdv
        if abs(vstep) > 0.1 * v00:
            vstep = math.copysign(0.2 * v00, vstep)
        v = v + vstep
        if abs(perror) <= eps:
            break
        if v <= 0.0:
            v = 100.0
            break

    # thermodynamic pieces
    z = T / Th[0]
    fz = z ** 4 / (1.0 + bb * z ** 2)
    fz1 = 2.0 * z ** 3 * (2.0 + bb * z ** 2) / (1.0 + bb * z ** 2) ** 2
    fz2 = 2.0 * z ** 2 * (6.0 + 3.0 * bb * z ** 2 + (bb ** 2)
                          * z ** 4) / (1.0 + bb * z ** 2) ** 3
    fv = math.exp(cc * (v - v00) / v00)

    cv = a[0] * R * (Debye3(Theta[0] / T) - (Theta[0] / T)
                     * Debye3D(Theta[0] / T))
    dpdt = Gamma[0] * cv
    cv = cv - (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv

    cp = cv - T * (dpdt ** 2) / dpdv
    KappaT = -1.0 / (v * dpdv)
    KappaS = -cv / (cp * v * dpdv)
    Alpha = -(dpdt / dpdv) / v

    lnz = math.log(v00 / v)
    U = v00 * (0.5 * a1 * lnz ** 2 + a2 * (lnz ** 3) /
               3.0 + 0.25 * a3 * (lnz ** 4))
    U = U + a[0] * R * T * Debye3(Theta[0] / T)
    U = U + aa * R * Th[0] * (fz - z * fz1) * fv

    S = a[0] * R * ((4.0 / 3.0) * Debye3(Theta[0] / T) -
                    math.log(1.0 - math.exp(-Theta[0] / T)))
    S = S - aa * R * fz1 * fv

    A = U - T * S
    H = U + p * v
    G = H - T * S

    Gruneisen = Alpha * v / (cv * KappaT)

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
