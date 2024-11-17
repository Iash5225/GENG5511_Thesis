import numpy as np
from constants import R


def computeThermoProps(T: float, p:float, parameters:list)->np.array:
    """
    Compute thermodynamic properties based on temperature (T) and pressure (p) using given model parameters.

    Parameters:
        T (float): Temperature in Kelvin.
        p (float): Pressure in MPa.
        parameters (list or numpy array): Model parameters.

    Returns:
        props (numpy array): Array of thermodynamic properties:
            [v, KappaT, KappaS, Alpha, cp, cv, Gruneisen, U, S, A, H, G]
    """
    # Constants
    # R = 8.31451  # Universal gas constant (J K^-1 mol^-1)
    v00 = parameters[0]  # Reference molar volume (cm³/mol)
    a1 = parameters[1]
    a2 = parameters[2]
    a3 = parameters[3]
    a = [3]
    Th = parameters[9:15]  # Characteristic temperatures
    g = parameters[15:21]  # Grüneisen parameters
    q = parameters[21:27]  # Volume dependence of Gamma
    aa = parameters[27]
    bb = parameters[28]
    cc = parameters[29]
    s = parameters[30]

    # Root-finding for volume v
    eps = 1e-10 * p
    v = 26 * np.exp(-p / 25000)  # Initial guess for v
    for _ in range(100):
        pcalc, dpdv, Theta, Gamma, lnz = evaluateArgon(
            T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)
        perror = p - pcalc
        vstep = perror / dpdv
        vstep = np.sign(vstep) * min(abs(vstep), 0.2 * v00)
        v += vstep
        if abs(perror) <= eps:
            break
        if v <= 0:
            v = 100  # Handle non-physical volume
            break

    # Compute properties
    z = T / Th[0]
    fz = z**4 / (1 + bb * z**2)
    fz1 = 2 * z**3 * (2 + bb * z**2) / (1 + bb * z**2)**2
    fz2 = 2 * z**2 * (6 + 3 * bb * z**2 + bb**2 * z**4) / (1 + bb * z**2)**3
    fv = np.exp(cc * (v - v00) / v00)

    # Heat capacities
    cv = a[0] * R * (Debye3(Theta[0] / T) - (Theta[0] / T)
                     * Debye3D(Theta[0] / T))
    dpdt = Gamma[0] * cv
    cv -= (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (dpdt / v) - aa * (cc / v00) * R * fz1 * fv
    cp = cv - T * (dpdt**2) / dpdv

    # Compressibilities and expansivity
    KappaT = -1 / (v * dpdv)
    KappaS = -cv / (cp * v * dpdv)
    Alpha = -(dpdt / dpdv) / v

    # Thermodynamic energies
    U = v00 * (0.5 * a1 * lnz**2 + a2 * (lnz**3) / 3 + 0.25 * a3 * (lnz**4))
    U += a[0] * R * T * Debye3(Theta[0] / T)
    U += aa * R * Th[0] * (fz - z * fz1) * fv
    S = a[0] * R * ((4 / 3) * Debye3(Theta[0] / T) -
                    np.log(1 - np.exp(-Theta[0] / T)))
    S -= aa * R * fz1 * fv

    # Thermodynamic potentials
    A = U - T * S  # Helmholtz energy
    H = U + p * v  # Enthalpy
    G = H - T * S  # Gibbs energy

    # Grüneisen parameter
    Gruneisen = Alpha * v / (cv * KappaT)

    # Store properties in an array
    props = np.array([v, KappaT, KappaS, Alpha, cp,
                     cv, Gruneisen, U, S, A, H, G])
    return props


def evaluateArgon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Evaluate pressure and its derivative at given temperature (T) and volume (v).

    Args:
        T (_type_): _description_
        v (_type_): _description_
        v00 (_type_): _description_
        a (_type_): _description_
        a1 (_type_): _description_
        a2 (_type_): _description_
        a3 (_type_): _description_
        Th (_type_): _description_
        g (_type_): _description_
        q (_type_): _description_
        aa (_type_): _description_
        bb (_type_): _description_
        cc (_type_): _description_

    Returns:
        _type_: _description_
    """    
    
    """
    Evaluate pressure and its derivative at given temperature (T) and volume (v).

    Returns:
        pcalc (float): Calculated pressure.
        dpdv (float): Derivative of pressure with respect to volume.
        Theta (list): Characteristic temperatures.
        Gamma (list): Grüneisen parameters.
        lnz (float): Logarithmic volume ratio.
    """
    # R = 8.31451
    Gamma = np.zeros(1)
    Theta = np.zeros(1)

    # Compute Gamma and Theta for j = 1
    for j in range(1):
        Gamma[j] = g[j] * (v / v00) ** q[j]
        if abs(q[j]) > 1e-10:
            Theta[j] = Th[j] * np.exp((g[j] / q[j]) * (1 - (v / v00) ** q[j]))
        else:
            Theta[j] = Th[j] * (v00 / v) ** g[j]

    z = v00 / v
    lnz = np.log(z)

    # Pressure components
    pj_0 = z * (a1 * lnz + a2 * lnz**2 + a3 * lnz**3)
    pj_1 = a[0] * (R * T * Gamma[0] / v) * Debye3(Theta[0] / T)
    pcalc = pj_0 + pj_1

    # Derivative of pressure
    dpj_0 = -(z / v) * (a1 * (1 + lnz) + a2 * lnz *
                        (2 + lnz) + a3 * lnz**2 * (3 + lnz))
    dpj_1 = (pj_1 / v) * (q[0] - 1) - a[0] * R * \
        Theta[0] * (Gamma[0] / v)**2 * Debye3D(Theta[0] / T)
    dpdv = dpj_0 + dpj_1

    # Apply corrections
    z = T / Th[0]
    exp_factor = np.exp(cc * (v - v00) / v00)
    correction = aa * (cc / v00) * R * \
        Th[0] * (z**4 / (1 + bb * z**2)) * exp_factor
    pcalc -= correction
    dpdv -= aa * (cc / v00)**2 * R * Th[0] * \
        (z**4 / (1 + bb * z**2)) * exp_factor

    return pcalc, dpdv, Theta, Gamma, lnz


def Debye3(x):
    """Debye integral approximation (stub)."""
    # Implement the Debye integral for heat capacity
    return 3 * (x**3 / (np.exp(x) - 1))


def Debye3D(x):
    """Derivative of Debye integral approximation (stub)."""
    # Implement the derivative of the Debye integral
    return 3 * (x**2 / (np.exp(x) - 1)**2)
