import numpy as np
from scipy.optimize import root_scalar
from math import pi
from scipy.integrate import quad
import numpy as np

R = 8.31451  # Universal gas constant (J K^-1 mol^-1)
DEBYE3_INF = (pi**4)/15.0

def _debye_integrand(t):
    if t < 1e-6:
        # series for t^3/(exp(t)-1) ~ t^2 - t^3/2 + t^4/12
        return t*t * (1.0 - 0.5*t + (t*t)/12.0)
    if t > 40.0:
        # denominator ~ exp(t)
        return t**3 * np.exp(-t)
    return t**3 / np.expm1(t)

def Debye3(x):
    #  ∫0^x t^3/(e^t-1) dt ; for x→∞, → π^4/15 ≈ 6.4939394
    if x <= 0:
        return 0.0
    x_eff = min(x, 60.0)  # cap to avoid crazy overflow
    integral, _ = quad(_debye_integrand, 0.0, x_eff, limit=200)
    # tail beyond 60 is negligible; skip or approximate by 0
    return integral


def Debye3D(x):
    # Numerical derivative of Debye3 using central difference
    h = 1e-5
    return (Debye3(x + h) - Debye3(x - h)) / (2 * h)


def evaluate_argon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    Gamma = g[0] * (v / v00) ** q[0]
    if abs(q[0]) > 1e-10:
        Theta = Th[0] * np.exp((g[0] / q[0]) * (1 - (v / v00) ** q[0]))
    else:
        Theta = Th[0] * (v00 / v) ** g[0]

    z = v00 / v
    lnz = np.log(z)
    pj_0 = z * (a1 * lnz + a2 * lnz**2 + a3 * lnz**3)
    pj_1 = a[0] * (R * T * Gamma / v) * Debye3(Theta / T)
    pcalc = pj_0 + pj_1

    dpj_0 = -(z / v) * (a1 * (1 + lnz) + a2 * lnz *
                        (2 + lnz) + a3 * lnz**2 * (3 + lnz))
    dpj_1 = (pj_1 / v) * (q[0] - 1) - a[0] * R * \
        Theta * (Gamma / v)**2 * Debye3D(Theta / T)
    dpdv = dpj_0 + dpj_1

    z_T = T / Th[0]
    exp_factor = np.exp(cc * (v - v00) / v00)
    correction = aa * (cc / v00) * R * \
        Th[0] * (z_T**4 / (1 + bb * z_T**2)) * exp_factor
    dcorrection = aa * (cc / v00)**2 * R * \
        Th[0] * (z_T**4 / (1 + bb * z_T**2)) * exp_factor

    pcalc -= correction
    dpdv -= dcorrection

    return pcalc, dpdv, Theta, Gamma, lnz


def compute_thermo_props(T, p, parameters):
    v00 = parameters[0]
    a1, a2, a3 = parameters[1:4]
    a = [3]  # constant coefficient array
    Th = parameters[9:15]
    g = parameters[15:21]
    q = parameters[21:27]
    aa, bb, cc, s = parameters[27:31]

    def pressure_root(v):
        pcalc, _, _, _, _ = evaluate_argon(
            T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)
        return p - pcalc

    sol = root_scalar(pressure_root, method='brentq',
                      bracket=[1e-4, 200], xtol=1e-6)
    v = sol.root if sol.converged else 26 * np.exp(-p / 25000)

    pcalc, dpdv, Theta, Gamma, lnz = evaluate_argon(
        T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)
    z = T / Th[0]
    fz = z**4 / (1 + bb * z**2)
    fz1 = 2 * z**3 * (2 + bb * z**2) / (1 + bb * z**2)**2
    fz2 = 2 * z**2 * (6 + 3 * bb * z**2 + bb**2 * z**4) / (1 + bb * z**2)**3
    fv = np.exp(cc * (v - v00) / v00)

    cv = a[0] * R * (Debye3(Theta / T) - (Theta / T) * Debye3D(Theta / T))
    dpdt = Gamma * cv
    cv -= (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv
    cp = cv - T * (dpdt**2) / dpdv
    KappaT = -1 / (v * dpdv)
    KappaS = -cv / (cp * v * dpdv)
    Alpha = -(dpdt / dpdv) / v

    U = v00 * (0.5 * a1 * lnz**2 + a2 * lnz**3 / 3 + 0.25 * a3 * lnz**4)
    U += a[0] * R * T * Debye3(Theta / T)
    U += aa * R * Th[0] * (fz - z * fz1) * fv
    S = a[0] * R * ((4 / 3) * Debye3(Theta / T) -
                    np.log(1 - np.exp(-Theta / T)))
    S -= aa * R * fz1 * fv

    A = U - T * S
    H = U + p * v
    G = H - T * S
    Gruneisen = Alpha * v / (cv * KappaT)

    props = np.zeros(12)
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
