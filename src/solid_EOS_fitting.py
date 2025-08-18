from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from computethermoprops import *
from p_functions import pmelt,psub
from constants import PARAMS_INIT, LOWER_BOUND, UPPER_BOUND, KRYPTON_P_t, KRYPTON_T_t, KRYPTON_REFERENCE_ENTROPY, KRYPTON_REFERENCE_ENTHALPY, PARAM_LABELS, PERCENT_SCALE, GAMMA_POS_SLOPE_MULT, GAMMA_POS_SLOPE_OFFSET, GAMMA_NEG_SLOPE_MULT, T6, CP_TEMP_THRESHOLD_K, CP_WEIGHT_BELOW, CP_WEIGHT_ABOVE, PMELT_EXTRA_WEIGHT_T_K, PMELT_EXTRA_FACTOR
from fitting_helper import rms , _mean_sq

BIG = 1e8
# Constants
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t

# --- cost function that returns (total, deviations) ---
def combined_cost_function(params, *datasets):
    (T_Vm_sub, p_Vm_sub, Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt,
     T_Vm_highp, p_Vm_highp, Vm_highp,
     T_cp_sub, p_cp_sub, cp_sub,
     T_alpha_sub, p_alpha_sub, alpha_sub,
     T_BetaT_sub, p_BetaT_sub, BetaT_sub,
     T_BetaS_sub, p_BetaS_sub, BetaS_sub,
     T_sub, p_sub, Year_sub, G_fluid_sub, V_fluid_sub,
     T_melt, p_melt, G_fluid_melt, V_fluid_melt,
     T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
     T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt) = datasets

    # triple-point adjustments
    deltaS_triple = params[30] - St_REFPROP
    props_Triple = compute_thermo_props(Tt, pt, params)
    Ht_fitted = props_Triple[11] + Tt * params[30]
    deltaH_triple = Ht_fitted - Ht_REFPROP

    # ---- deviations (mirror your MATLAB) ----
    Vm_sub_dev = rms(PERCENT_SCALE * (Vm_sub - np.array([compute_thermo_props(T, p, params)[0]
                                               for T, p in zip(T_Vm_sub, p_Vm_sub)])) / Vm_sub)

    Vm_melt_dev = rms(PERCENT_SCALE * (Vm_melt - np.array([compute_thermo_props(T, p, params)[0]
                                                 for T, p in zip(T_Vm_melt, p_Vm_melt)])) / Vm_melt)

    Vm_highp_dev = rms(PERCENT_SCALE * (Vm_highp - np.array([compute_thermo_props(T, p, params)[0]
                                                   for T, p in zip(T_Vm_highp, p_Vm_highp)])) / Vm_highp)

    cp_sub_dev_terms = []
    for T, p, cp in zip(T_cp_sub, p_cp_sub, cp_sub):
        model_cp = compute_thermo_props(T, p, params)[4]
        # MATLAB: weight 100 below 12 K, 700 above/equal
        w = CP_WEIGHT_BELOW if T < CP_TEMP_THRESHOLD_K else CP_WEIGHT_ABOVE
        cp_sub_dev_terms.append(w * (cp - model_cp) / cp)
    cp_sub_dev = rms(cp_sub_dev_terms)

    alpha_sub_dev = rms(PERCENT_SCALE * (alpha_sub - np.array([compute_thermo_props(T, p, params)[3]
                                                     for T, p in zip(T_alpha_sub, p_alpha_sub)])) / alpha_sub)

    BetaT_sub_dev = rms(PERCENT_SCALE * (BetaT_sub - np.array([compute_thermo_props(T, p, params)[1]
                                                     for T, p in zip(T_BetaT_sub, p_BetaT_sub)])) / BetaT_sub)

    BetaS_sub_dev = rms(PERCENT_SCALE * (BetaS_sub - np.array([compute_thermo_props(T, p, params)[2]
                                                     for T, p in zip(T_BetaS_sub, p_BetaS_sub)])) / BetaS_sub)

    # Enthalpy (sublimation)
    H_solid_sub = (H_fluid_sub * 1000.0) - (delta_H_sub * 1000.0)  # → J/mol
    H_solid_sub_fitted = np.array([compute_thermo_props(T, p, params)[10]
                                   for T, p in zip(T_H_sub, p_H_sub)]) - deltaH_triple
    H_solid_sub_dev = rms(
        PERCENT_SCALE * (H_solid_sub - H_solid_sub_fitted) / H_solid_sub_fitted)

    # Enthalpy (melting)
    H_solid_melt = (H_fluid_melt * 1000.0) - (delta_H_melt * 1000.0)  # → J/mol
    H_solid_melt_fitted = np.array([compute_thermo_props(T, p, params)[10]
                                    for T, p in zip(T_H_melt, p_H_melt)]) - deltaH_triple
    H_solid_melt_dev = rms(
        PERCENT_SCALE * (H_solid_melt - H_solid_melt_fitted) / H_solid_melt_fitted)

    # Sublimation pressure (experimental p in MPa; solid model in MPa)
    p_fitted_sub = []
    for T, pPa, Gf, Vf in zip(T_sub, p_sub, G_fluid_sub, V_fluid_sub):
        props = compute_thermo_props(T, pPa, params)
        delta_G = Gf - props[11] + deltaH_triple - T * deltaS_triple
        pf = pPa - (delta_G / (Vf - props[0]))
        p_fitted_sub.append(pf)
    p_sub_dev = rms(PERCENT_SCALE * (np.array(p_sub) -
                    np.array(p_fitted_sub)) / p_sub)

    # Melting pressure (MPa)
    p_fitted_melt = []
    for T, pM, Gf, Vf in zip(T_melt, p_melt, G_fluid_melt, V_fluid_melt):
        props = compute_thermo_props(T, pM, params)
        delta_G = Gf - props[11] + deltaH_triple - T * deltaS_triple
        pf = pM - (delta_G / (Vf - props[0]))
        # MATLAB: extra weight for T > 500 K (implemented as scaling pf diff)
        if T > PMELT_EXTRA_WEIGHT_T_K:
            pf = pM - (delta_G / (Vf - props[0])) * PMELT_EXTRA_FACTOR
        p_fitted_melt.append(pf)
    p_melt_dev = rms(
        PERCENT_SCALE * (np.array(p_melt) - np.array(p_fitted_melt)) / p_melt)

    # Gamma-T penalty on a temperature grid
    
    # note: your psub signature uses (T, pt, Tt)
    p6 = np.array([psub(T, pt, Tt) for T in T6])
    Gamma_T6 = np.array([compute_thermo_props(T, p, params)[6]
                        for T, p in zip(T6, p6)])
    Vm_T6 = np.array([compute_thermo_props(T, p, params)[0]
                     for T, p in zip(T6, p6)])
    slopes = np.diff(Gamma_T6) / np.diff(Vm_T6)

    Gamma_dev_terms = []
    mu = Gamma_T6.mean()
    for i, s in enumerate(slopes, start=1):
        if s > 0:
            Gamma_dev_terms.append(
                GAMMA_POS_SLOPE_OFFSET + abs(GAMMA_POS_SLOPE_MULT * (Gamma_T6[i] - mu) / mu))
        else:
            Gamma_dev_terms.append(
                GAMMA_NEG_SLOPE_MULT * (Gamma_T6[i] - mu) / mu)
    Gamma_T6_dev = rms(Gamma_dev_terms)

    # Weighted total (your MATLAB weights)
    total_deviation = (
        Vm_sub_dev * 55 +
        Vm_melt_dev * 35 +
        Vm_highp_dev +
        cp_sub_dev * 30 +
        alpha_sub_dev * 50 +
        BetaT_sub_dev * 20 +
        BetaS_sub_dev * 30 +
        H_solid_sub_dev * 25 +
        H_solid_melt_dev * 2 +
        p_sub_dev * 2 +
        p_melt_dev * 5 +
        Gamma_T6_dev * 3.5
    )

    deviations = {
        "Vm_sub": Vm_sub_dev,
        "Vm_melt": Vm_melt_dev,
        "Vm_highp": Vm_highp_dev,
        "cp_sub": cp_sub_dev,
        "alpha_sub": alpha_sub_dev,
        "BetaT_sub": BetaT_sub_dev,
        "BetaS_sub": BetaS_sub_dev,
        "H_solid_sub": H_solid_sub_dev,
        "H_solid_melt": H_solid_melt_dev,
        "p_sub": p_sub_dev,
        "p_melt": p_melt_dev,
        "Gamma_T": Gamma_T6_dev,
    }
    return total_deviation, deviations

# --- run_optimization using the callback ---
def run_optimization(params_init, bounds, datasets):
    print("Lens:", len(PARAMS_INIT), len(LOWER_BOUND), len(UPPER_BOUND))
    assert len(PARAMS_INIT) == len(LOWER_BOUND) == len(UPPER_BOUND)
    assert len(params_init) == len(
        bounds), "params_init/bounds length mismatch"

    callback = _make_outfun(*datasets)

    # 0) make sure these are python lists/tuples, not a mis-shaped numpy array
    bounds = [(lo, hi) for lo, hi in zip(LOWER_BOUND, UPPER_BOUND)]

    # 1) dimensions must match
    print("len(params_init) =", len(params_init))
    print("len(bounds)      =", len(bounds))
    assert len(params_init) == len(bounds), "params_init/bounds length mismatch"

    # 2) show SciPy what it will actually receive
    print("bounds[0:5] =", bounds[:5])

    # 3) does the initial cost look sensible (non-zero, finite)?
    total0, dev0 = combined_cost_function(params_init, *datasets)
    print("initial total cost =", total0)
    if not np.isfinite(total0):
        raise RuntimeError("Initial cost is not finite")
    
    tot0, dev0 = combined_cost_function(params_init, *datasets)


    print("initial total =", tot0)
    for k, v in dev0.items():
        print(f"  {k}: {v}")

    result = minimize(
        fun=_cost_only,
        x0=np.asarray(params_init, dtype=float),
        args=datasets,                 # datasets is already a tuple → OK
        method='L-BFGS-B',
        bounds=bounds,
        callback=callback,             # expects cb(xk)
        options={
            'disp': True,
            'maxiter': 200,
            'ftol': 1e-8,              # function tolerance
            'gtol': 1e-6,              # projected gradient tol (stopping)
            'maxls': 40,               # line-search steps (stability)
        }
    )
    return result.x, result.fun

def extract_datasets(data):
    """
    Extracts thermodynamic property datasets from the given data dictionary.

    Parameters:
        data (dict): Dictionary containing sub-dictionaries of thermodynamic property data.

    Returns:
        tuple: A tuple containing arrays for all datasets in the specified order.
    """

    # Cell Volume Sublimation
    T_Vm_sub = data["cell_volume_sub"]['Temperature']
    p_Vm_sub = np.array([psub(T, pt, Tt) for T in T_Vm_sub])
    Vm_sub = data['cell_volume_sub']['Cell Volume']

    # Cell Volume Melting
    T_Vm_melt = data['cell_volume_melt']['Temperature']
    p_Vm_melt = np.array([pmelt(T, pt, Tt) for T in T_Vm_melt])
    Vm_melt = data['cell_volume_melt']['Cell Volume']

    # High Pressure Cell Volume (safe defaults if missing)
    if 'cell_volume_highp' in data:
        T_Vm_highp = data['cell_volume_highp']['Temperature']
        p_Vm_highp = data['cell_volume_highp']['Pressure']
        Vm_highp = data['cell_volume_highp']['Cell Volume']
    else:
        T_Vm_highp = np.array([])
        p_Vm_highp = np.array([])
        Vm_highp = np.array([])

    # Heat Capacity Sublimation
    T_cp_sub = data['heat_capacity']['Temperature']
    p_cp_sub = np.array([psub(T, pt, Tt) for T in T_cp_sub])
    cp_sub = data['heat_capacity']['Heat Capacity']

    # Thermal Expansion Sublimation
    T_alpha_sub = data['thermal_coeff']['Temperature']
    p_alpha_sub = np.array([psub(T, pt, Tt) for T in T_alpha_sub])
    alpha_sub = data['thermal_coeff']['Thermal Expansion Coefficient']

    # Bulk Modulus (S)
    T_BetaS_sub = data['bulk_s']['Temperature']
    p_BetaS_sub = data['bulk_s']['Pressure']
    BetaS_sub = data['bulk_s']['Beta S']

    # Bulk Modulus (T)
    T_BetaT_sub = data['bulk_t']['Temperature']
    p_BetaT_sub = data['bulk_t']['Pressure']
    BetaT_sub = data['bulk_t']['Beta T']

    # Melting
    T_melt = data['melting']['Temperature']
    p_melt = data['melting']['Pressure']
    G_fluid_melt = data['melting']['Gibbs Energy']
    V_fluid_melt = data['melting']['Volume']

    # Sublimation
    T_sub = data['sublimation']['Temperature']    
    p_sub = data['sublimation']['Pressure']
    G_fluid_sub = data['sublimation']['Gibbs Energy']
    V_fluid_sub = data['sublimation']['Volume']
    # Enthalpy Sublimation
    T_H_sub = data['heatsub']['Temperature']
    p_H_sub = data['heatsub']['Pressure']
    delta_H_sub = data['heatsub']['Change in Enthalpy']
    H_fluid_sub = data['heatsub']['Enthalpy']

    # Enthalpy Melting
    T_H_melt = data['fusion']['Temperature']
    p_H_melt = data['fusion']['Pressure']
    delta_H_melt = data['fusion']['Change in Enthalpy']
    H_fluid_melt = data['fusion']['Enthalpy']

    Year_sub = np.array([])   # <-- add this


    datasets = (
        T_Vm_sub, p_Vm_sub, Vm_sub,
        T_Vm_melt, p_Vm_melt, Vm_melt,
        T_Vm_highp, p_Vm_highp, Vm_highp,
        T_cp_sub, p_cp_sub, cp_sub,
        T_alpha_sub, p_alpha_sub, alpha_sub,
        T_BetaT_sub, p_BetaT_sub, BetaT_sub,
        T_BetaS_sub, p_BetaS_sub, BetaS_sub,
        # ---- EXACT ORDER NEEDED BY COMBINED_COST_FUNCTION ----
        T_sub, p_sub, Year_sub, G_fluid_sub, V_fluid_sub,
        T_melt, p_melt, G_fluid_melt, V_fluid_melt,
        T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
        T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt
    )

    assert len(datasets) == 38
    return datasets


def _cost_only(params, *datasets):
    try:
        total, _ = combined_cost_function(params, *datasets)
        return float(total) if np.isfinite(total) else BIG
    except Exception:
        return BIG

def _make_outfun(*datasets):
    def outfun(xk):
        total, dev = combined_cost_function(xk, *datasets)
        print(f"Current total_deviation: {total:.6g}")
        print(f"Vm_sub deviation: {dev['Vm_sub']:.6g}")
        print(f"Vm_melt deviation: {dev['Vm_melt']:.6g}")
        print(f"Vm_highp deviation: {dev['Vm_highp']:.6g}")
        print(f"cp deviation: {dev['cp_sub']:.6g}")
        print(f"alpha_sub deviation: {dev['alpha_sub']:.6g}")
        print(f"BetaT_sub deviation: {dev['BetaT_sub']:.6g}")
        print(f"BetaS_sub deviation: {dev['BetaS_sub']:.6g}")
        print(f"H_solid_sub deviation: {dev['H_solid_sub']:.6g}")
        print(f"H_solid_melt deviation: {dev['H_solid_melt']:.6g}")
        print(f"p_sub deviation: {dev['p_sub']:.6g}")
        print(f"p_melt deviation: {dev['p_melt']:.6g}")
        print(f"Gamma T deviation: {dev['Gamma_T']:.6g}")
        # return None to continue (SciPy callback can’t stop like MATLAB's `stop`)
    return outfun


def debug_datasets(datasets, Tt):
    (T_Vm_sub, p_Vm_sub, Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt,
     T_Vm_highp, p_Vm_highp, Vm_highp,
     T_cp_sub, p_cp_sub, cp_sub,
     T_alpha_sub, p_alpha_sub, alpha_sub,
     T_BetaT_sub, p_BetaT_sub, BetaT_sub,
     T_BetaS_sub, p_BetaS_sub, BetaS_sub,
     T_sub, p_sub, Year_sub, G_fluid_sub, V_fluid_sub,
     T_melt, p_melt, G_fluid_melt, V_fluid_melt,
     T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
     T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt) = datasets

    def _cnt(name, arr):
        arr = np.asarray(arr)
        print(f"{name:16s} n={arr.size}  finite={np.isfinite(arr).sum()}")

    print("=== DATASET COUNTS ===")
    _cnt("T_Vm_sub", T_Vm_sub)
    _cnt("p_Vm_sub", p_Vm_sub)
    _cnt("Vm_sub", Vm_sub)
    _cnt("T_Vm_melt", T_Vm_melt)
    _cnt("p_Vm_melt", p_Vm_melt)
    _cnt("Vm_melt", Vm_melt)
    _cnt("T_cp_sub", T_cp_sub)
    _cnt("p_cp_sub", p_cp_sub)
    _cnt("cp_sub", cp_sub)
    _cnt("T_alpha_sub", T_alpha_sub)
    _cnt("p_alpha_sub", p_alpha_sub)
    _cnt("alpha_sub", alpha_sub)
    _cnt("T_BetaT_sub", T_BetaT_sub)
    _cnt("p_BetaT_sub", p_BetaT_sub)
    _cnt("BetaT_sub", BetaT_sub)
    _cnt("T_BetaS_sub", T_BetaS_sub)
    _cnt("p_BetaS_sub", p_BetaS_sub)
    _cnt("BetaS_sub", BetaS_sub)
    _cnt("T_sub", T_sub)
    _cnt("p_sub", p_sub)
    _cnt("G_fluid_sub", G_fluid_sub)
    _cnt("V_fluid_sub", V_fluid_sub)
    _cnt("T_melt", T_melt)
    _cnt("p_melt", p_melt)
    _cnt("G_fluid_melt", G_fluid_melt)
    _cnt("V_fluid_melt", V_fluid_melt)
    _cnt("T_H_sub", T_H_sub)
    _cnt("p_H_sub", p_H_sub)
    _cnt("delta_H_sub", delta_H_sub)
    _cnt("H_fluid_sub", H_fluid_sub)
    _cnt("T_H_melt", T_H_melt)
    _cnt("p_H_melt", p_H_melt)
    _cnt("delta_H_melt", delta_H_melt)
    _cnt("H_fluid_melt", H_fluid_melt)

    # Phase coverage checks
    print("\n=== PHASE COVERAGE ===")
    print("sublimation (T<=Tt):", np.sum(
        np.asarray(T_sub) <= Tt), "of", len(T_sub))
    print("melting     (T>=Tt):", np.sum(
        np.asarray(T_melt) >= Tt), "of", len(T_melt))

def main():
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    # xenon_data = load_all_gas_data('xenon', read_from_excel=False)
    # neon_data = load_all_gas_data('neon', read_from_excel=False)

    bounds = list(zip(LOWER_BOUND, UPPER_BOUND))
    datasets = extract_datasets(krypton_data)
    debug_datasets(datasets, Tt)

    params_fit, fval = run_optimization(PARAMS_INIT, bounds, datasets)

    print("Optimized parameters:", params_fit)
    print("Final objective value:", fval)

    # If your vector is longer/shorter, we’ll truncate/pad labels to match length:
    n = len(params_fit)
    labels = PARAM_LABELS[:n] if len(PARAM_LABELS) >= n else PARAM_LABELS + [
        ("param_"+str(i+1), "") for i in range(len(PARAM_LABELS), n)]

main()
