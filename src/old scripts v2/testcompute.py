# compute_thermo_props_for_result.py
import math
import numpy as np

# ---- constants ----
R = 8.31451  # J/(mol·K)

# ---- Debye helpers (match MATLAB usage) ----


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
_DEBINF = 5.13299112734217e-02
# Machine-like constants (double)
_EPS = 1.11e-16
_XLOW = math.sqrt(8.0*_EPS)         # ~ sqrt(8*eps)
_XUP = -math.log(2.0*_EPS)         # switch to exp tail
_XLIM1 = -math.log(2.2251e-308)     # -log(xmin)
# XLIM2 = (1/DEBINF)^(1/3) / (xmin)^(1/3)
_XLIM2 = (1.0/_DEBINF)**(1.0/3.0) / (2.2251e-308)**(1.0/3.0)


def Debye3(x):
    """Normalized Debye D3(x) = 3/x^3 ∫ t^3/(e^t-1) dt  (0 ≤ D3 ≤ 1)."""
    if x <= 0.0:
        return 0.0
    if x <= 4.0:
        if x < _XLOW:
            # small-x series: 1 - 3x/8 + x^2/20
            return ((x - 7.5)*x + 20.0)/20.0
        t = (x*x/8.0) - 0.5 - 0.5   # map to [-1,1]
        return _cheval(_ADEB3, t) - 0.375*x
    # x > 4:
    if x > _XLIM2:
        return 0.0
    val = 1.0/(_DEBINF * x**3)     # = 18*zeta(4)/x^3
    if x < _XLIM1:
        ex = math.exp(-x)
        if x > _XUP:
            # 3*exp(-x)(x^3+3x^2+6x+6)/x^3 term
            tail = ((x + 3.0)*x + 6.0)*x + 6.0
            val -= 3.0 * tail * ex / (x**3)
        else:
            # refined Poisson summation tail
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


def Debye3D(x):
    """Derivative D3'(x) of the normalized Debye function D3(x)."""
    if x <= 0.0:
        return 0.0
    if x <= 4.0:
        if x < _XLOW:
            return (4.0*x - 15.0)/40.0
        t = (x*x/8.0) - 0.5 - 0.5
        return 0.25*x*_cheval(_D, t) - 0.375
    # x > 4:
    if x > _XLIM2:
        return 0.0
    val = -3.0/(_DEBINF * x**4)    # derivative of 18*zeta(4)/x^3
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
