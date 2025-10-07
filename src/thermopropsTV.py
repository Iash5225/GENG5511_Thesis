import numpy as np
from debye import Debye
from debye3 import Debye3
from debye3D import Debye3D

R = 8.31451  # J/(mol K)


def compute_thermo_props_TV(T: float, v: float, parameters):
    """
    Compute thermodynamic properties based on temperature (T) and volume (v),
    using given model parameters.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.
    v : float
        Volume in cm^3/mol.
    parameters : list
        Model parameters.

    Returns
    -------
    np.ndarray
        Array of computed thermodynamic properties.
    """
    # Unpack parameters
    v00 = parameters[0]  # Reference molar volume (cm^3/mol)
    a1 = parameters[1]
    a2 = parameters[2]
    a3 = parameters[3]
    a = [3.0]  # MATLAB constant
    Th = np.array(parameters[9:15], dtype=float)  # Characteristic temperatures
    g = np.array(parameters[15:21], dtype=float)  # Gruneisen parameters
    # Exponent for volume dependence of Gamma
    q = np.array(parameters[21:27], dtype=float)
    aa = parameters[27]
    bb = parameters[28]
    cc = parameters[29]
    s = parameters[30]

    # Evaluate pressure and its derivative
    pcalc, dpdv, Theta, Gamma, lnz = evaluate_argon(
        T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc)

    # Compute other properties
    z = T / Th[0]
    fz = z**4 / (1 + bb * z**2)
    fz1 = 2 * z**3 * (2 + bb * z**2) / (1 + bb * z**2)**2
    fz2 = 2 * z**2 * (6 + 3 * bb * z**2 + bb**2 * z**4) / (1 + bb * z**2)**3
    fv = np.exp(cc * (v - v00) / v00)

    # Compute specific heat capacities, compressibilities, and other properties
    cv_lattice = a[0] * R * (Debye3(Theta[0] / T) -
                             (Theta[0] / T) * Debye3D(Theta[0] / T))
    dpdt_base = Gamma[0] * cv_lattice
    cv = cv_lattice - (aa * R * T / Th[0]) * fz2 * fv
    dpdt = (dpdt_base / v) - aa * (cc / v00) * R * fz1 * fv
    cp = cv - T * (dpdt**2) / dpdv
    KappaT = -1 / (v * dpdv)
    KappaS = -cv / (cp * v * dpdv)
    Alpha = -(dpdt / dpdv) / v

    # Compute thermodynamic energies
    U = v00 * (0.5 * a1 * lnz**2 + (a2 / 3) * lnz**3 + 0.25 * a3 * lnz**4)
    U += a[0] * R * T * Debye3(Theta[0] / T)
    U += aa * R * Th[0] * (fz - z * fz1) * fv

    # Entropy
    log_term = np.log1p(-np.exp(-Theta[0] / T)
                        ) if Theta[0] / T < 40.0 else np.log1p(-0.0)
    S = a[0] * R * ((4 / 3) * Debye3(Theta[0] / T) - log_term)
    S -= aa * R * fz1 * fv

    # Thermodynamic potentials
    A = U - T * S  # Helmholtz energy
    H = U + pcalc * v  # Enthalpy
    G = H - T * S  # Gibbs energy

    # Thermal Grüneisen parameter
    Gruneisen = Alpha * v / (cv * KappaT)

    # Store all properties in a numpy array
    props = np.array([pcalc, KappaT, KappaS, Alpha, cp, cv,
                     Gruneisen, U, S, A, H, G], dtype=float)
    return props


def evaluate_argon(T, v, v00, a, a1, a2, a3, Th, g, q, aa, bb, cc):
    """
    Evaluate pressure (pcalc) and its derivative (dpdv) at a given T and v.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.
    v : float
        Volume in cm^3/mol.
    v00 : float
        Reference molar volume in cm^3/mol.
    a : list
        List of constants.
    a1, a2, a3 : float
        Coefficients for the equation of state.
    Th : np.ndarray
        Characteristic temperatures.
    g : np.ndarray
        Gruneisen parameters.
    q : np.ndarray
        Exponent for volume dependence of Gamma.
    aa, bb, cc : float
        Correction parameters.

    Returns
    -------
    pcalc : float
        Calculated pressure in MPa.
    dpdv : float
        Derivative of pressure with respect to volume.
    Theta : np.ndarray
        Characteristic temperature array.
    Gamma : np.ndarray
        Gruneisen parameter array.
    lnz : float
        Natural logarithm of z.
    """
    Gamma = np.zeros(1, dtype=float)
    Theta = np.zeros(1, dtype=float)

    # Calculate Gamma and Theta
    Gamma[0] = g[0] * (v / v00)**q[0]
    if abs(q[0]) > 1e-10:
        Theta[0] = Th[0] * np.exp((g[0] / q[0]) * (1 - (v / v00)**q[0]))
    else:
        Theta[0] = Th[0] * (v00 / v)**g[0]

    # Calculate pressure components
    z = v00 / v
    lnz = np.log(z)

    pj_0 = z * (a1 * lnz + a2 * lnz**2 + a3 * lnz**3)
    pj_1 = a[0] * (R * T * Gamma[0] / v) * Debye3(Theta[0] / T)

    pcalc = pj_0 + pj_1

    # Compute dp/dv
    dpj_0 = -(z / v) * (a1 * (1 + lnz) + a2 * lnz *
                        (2 + lnz) + a3 * lnz**2 * (3 + lnz))
    dpj_1 = (pj_1 / v) * (q[0] - 1) - a[0] * R * \
        Theta[0] * (Gamma[0] / v)**2 * Debye3D(Theta[0] / T)

    dpdv = dpj_0 + dpj_1

    # Apply corrections to pcalc and dpdv
    zT = T / Th[0]
    fT = (zT**4) / (1 + bb * zT**2)
    exp_factor = np.exp(cc * (v - v00) / v00)

    correction_p = aa * (cc / v00) * R * Th[0] * fT * exp_factor
    pcalc -= correction_p

    dpdv_correction = aa * (cc / v00)**2 * R * Th[0] * fT * exp_factor
    dpdv -= dpdv_correction

    return pcalc, dpdv, Theta, Gamma, lnz
