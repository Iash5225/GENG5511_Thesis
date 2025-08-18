from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from thermal_script import plot_all_gas_properties
from computethermoprops import *
from p_functions import pmelt,psub
from constants import CUSTOMCOLORS, CUSTOMMARKERS, PARAMS_INIT, LOWER_BOUND, UPPER_BOUND, KRYPTON_P_t, KRYPTON_T_t, KRYPTON_REFERENCE_ENTROPY, KRYPTON_REFERENCE_ENTHALPY
from scipy.integrate import quad


# Constantsn 
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t


# solid_EOS_fitting.py
def rms(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    return 0.0 if n == 0 else np.sqrt(np.sum(x*x) / n)


# Optimisation Options
def combined_cost_function(params, *datasets):
    # Unpack all datasets
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
    
    deltaS_triple = params[30] - St_REFPROP
    props_Triple = compute_thermo_props(Tt, pt, params)
    Ht_fitted = props_Triple[11] + Tt * params[30]
    deltaH_triple = Ht_fitted - Ht_REFPROP

    # === Deviation functions (as helper lambdas) ===
    # def rms(x): return np.sqrt(np.sum(np.square(x)) / len(x))

    # === Calculate individual deviations ===
    Vm_sub_dev = rms(100 * (Vm_sub - np.array([compute_thermo_props(
        T, p, params)[0] for T, p in zip(T_Vm_sub, p_Vm_sub)])) / Vm_sub)
    Vm_melt_dev = rms(100 * (Vm_melt - np.array([compute_thermo_props(T, p, params)[
                      0] for T, p in zip(T_Vm_melt, p_Vm_melt)])) / Vm_melt)
    Vm_highp_dev = rms(100 * (Vm_highp - np.array([compute_thermo_props(T, p, params)[
                       0] for T, p in zip(T_Vm_highp, p_Vm_highp)])) / Vm_highp)
    
    cp_sub_dev = []
    for T, p, cp in zip(T_cp_sub, p_cp_sub, cp_sub):
        model_cp = compute_thermo_props(T, p, params)[4]
        factor = 700 if T >= 12 else 100
        cp_sub_dev.append(factor * (cp - model_cp) / cp)
    cp_sub_dev = rms(cp_sub_dev)

    alpha_sub_dev = rms(100 * (alpha_sub - np.array([compute_thermo_props(
        T, p, params)[3] for T, p in zip(T_alpha_sub, p_alpha_sub)])) / alpha_sub)
    BetaT_sub_dev = rms(100 * (BetaT_sub - np.array([compute_thermo_props(
        T, p, params)[1] for T, p in zip(T_BetaT_sub, p_BetaT_sub)])) / BetaT_sub)
    BetaS_sub_dev = rms(100 * (BetaS_sub - np.array([compute_thermo_props(
        T, p, params)[2] for T, p in zip(T_BetaS_sub, p_BetaS_sub)])) / BetaS_sub)

    # === Enthalpy of sublimation and melting ===
    H_solid_sub = H_fluid_sub - 1000 * delta_H_sub
    H_solid_sub_fitted = np.array([compute_thermo_props(T, p, params)[
                                  10] for T, p in zip(T_H_sub, p_H_sub)]) - deltaH_triple
    H_solid_sub_dev = rms(
        100 * (H_solid_sub - H_solid_sub_fitted) / H_solid_sub_fitted)
    
    H_solid_melt = H_fluid_melt - 1000 * delta_H_melt
    H_solid_melt_fitted = np.array([compute_thermo_props(T, p, params)[
                                   10] for T, p in zip(T_H_melt, p_H_melt)]) - deltaH_triple
    H_solid_melt_dev = rms(
        100 * (H_solid_melt - H_solid_melt_fitted) / H_solid_melt_fitted)

    # === Sublimation pressure deviation ===
    p_fitted_sub = []
    for T, p, Gf, Vf in zip(T_sub, p_sub, G_fluid_sub, V_fluid_sub):
        props = compute_thermo_props(T, p / 1e6, params)  # Convert Pa -> MPa
        delta_G = Gf - props[11] + deltaH_triple - T * deltaS_triple
        pf = p - delta_G / (Vf - props[0]) * 1e6  # Pa
        p_fitted_sub.append(pf)
    p_sub_dev = rms(100 * (np.array(p_sub) - np.array(p_fitted_sub)) / p_sub)

    # === Melting pressure deviation ===
    p_fitted_melt = []
    for T, p, Gf, Vf in zip(T_melt, p_melt, G_fluid_melt, V_fluid_melt):
        props = compute_thermo_props(T, p, params)
        delta_G = Gf - props[11] + deltaH_triple - T * deltaS_triple
        pf = p - delta_G / (Vf - props[0])
        if T > 500:
            pf = p - delta_G / (Vf - props[0]) * 5
        p_fitted_melt.append(pf)
    p_melt_dev = rms(
        100 * (np.array(p_melt) - np.array(p_fitted_melt)) / p_melt)

    # === Gamma slope penalty ===
    T6 = np.array([0.0001] + list(range(2, 84, 2)) + [83.806])
    p6 = np.array([psub(T, pt, Tt) for T in T6])
    Gamma_T6 = np.array([compute_thermo_props(T, p, params)[6]
                        for T, p in zip(T6, p6)])
    Vm_T6 = np.array([compute_thermo_props(T, p, params)[0]
                     for T, p in zip(T6, p6)])

    slopes = np.diff(Gamma_T6) / np.diff(Vm_T6)
    Gamma_dev = []
    for i, s in enumerate(slopes):
        if s > 0:
            Gamma_dev.append(
                45 + abs(200 * (Gamma_T6[i + 1] - Gamma_T6.mean()) / Gamma_T6.mean()))
        else:
            Gamma_dev.append(
                50 * (Gamma_T6[i + 1] - Gamma_T6.mean()) / Gamma_T6.mean())
    Gamma_T6_dev = rms(Gamma_dev)

    # === Weighted total deviation ===
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

    # Save individual deviations for tracking
    deviations = {
        'Vm_sub': Vm_sub_dev,
        'Vm_melt': Vm_melt_dev,
        'Vm_highp': Vm_highp_dev,
        'cp_sub': cp_sub_dev,
        'alpha_sub': alpha_sub_dev,
        'BetaT_sub': BetaT_sub_dev,
        'BetaS_sub': BetaS_sub_dev,
        'H_solid_sub': H_solid_sub_dev,
        'H_solid_melt': H_solid_melt_dev,
        'p_sub': p_sub_dev,
        'p_melt': p_melt_dev,
        'Gamma_T': Gamma_T6_dev
    }
    return total_deviation


# === Weights from the thesis table ===
W_VM_SUB = 55
W_CP_LOW = 30     # T <= 12 K
W_CP_HIGH = 210    # T > 12 K
W_ALPHA = 50
W_KT = 20
W_KS = 20
W_DH_SUB = 25
W_PSUB = 2

W_VM_MELT = 35
W_DH_MELT = 2
W_PMELT_LO = 5      # T <= 500 K
W_PMELT_HI = 25     # T > 500 K

W_VM_HIGHP = 1
W_LOW = 3.5    # low-T behaviour penalty


def _mean_sq(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 0.0 if x.size == 0 else np.mean(x*x)

# Optimisation Options


def combined_cost_function2(params, *datasets):
    # Unpack all datasets (unchanged)
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

    # Reference offsets at the triple point
    deltaS_triple = params[30] - St_REFPROP
    props_Triple = compute_thermo_props(Tt, pt, params)
    Ht_fitted = props_Triple[11] + Tt * params[30]
    deltaH_triple = Ht_fitted - Ht_REFPROP

    # ---------- Relative deviations (arrays) ----------
    # Vm on sublimation/melting/high-p
    Vm_sub_fit = np.array([compute_thermo_props(T, p, params)[0]
                          for T, p in zip(T_Vm_sub,  p_Vm_sub)])
    Vm_melt_fit = np.array([compute_thermo_props(T, p, params)[0]
                           for T, p in zip(T_Vm_melt, p_Vm_melt)])
    Vm_highp_fit = np.array([compute_thermo_props(T, p, params)[0]
                            for T, p in zip(T_Vm_highp, p_Vm_highp)])

    rel_Vm_sub = (Vm_sub - Vm_sub_fit) / Vm_sub
    rel_Vm_melt = (Vm_melt - Vm_melt_fit) / Vm_melt
    rel_Vm_highp = (Vm_highp - Vm_highp_fit) / Vm_highp

    # cp on sublimation, split at 12 K
    cp_fit = np.array([compute_thermo_props(T, p, params)[4]
                      for T, p in zip(T_cp_sub, p_cp_sub)])
    rel_cp = (cp_sub - cp_fit) / cp_sub
    mask_cp_low = np.asarray(T_cp_sub) <= 12.0
    mask_cp_high = ~mask_cp_low

    # alpha, BetaT (KT), BetaS (KS)
    alpha_fit = np.array([compute_thermo_props(T, p, params)[3]
                         for T, p in zip(T_alpha_sub, p_alpha_sub)])
    rel_alpha = (alpha_sub - alpha_fit) / alpha_sub

    KT_fit = np.array([compute_thermo_props(T, p, params)[1]
                      for T, p in zip(T_BetaT_sub, p_BetaT_sub)])
    rel_KT = (BetaT_sub - KT_fit) / BetaT_sub

    KS_fit = np.array([compute_thermo_props(T, p, params)[2]
                      for T, p in zip(T_BetaS_sub, p_BetaS_sub)])
    rel_KS = (BetaS_sub - KS_fit) / BetaS_sub

    # ΔH sublimation & melting (solid enthalpy inferred from fluid – ΔH)
    H_solid_sub = H_fluid_sub - 1000.0 * delta_H_sub
    H_solid_sub_fit = np.array([compute_thermo_props(T, p, params)[
                               10] for T, p in zip(T_H_sub,  p_H_sub)]) - deltaH_triple
    rel_DH_sub = (H_solid_sub - H_solid_sub_fit) / H_solid_sub_fit

    H_solid_melt = H_fluid_melt - 1000.0 * delta_H_melt
    H_solid_melt_fit = np.array([compute_thermo_props(T, p, params)[
                                10] for T, p in zip(T_H_melt, p_H_melt)]) - deltaH_triple
    rel_DH_melt = (H_solid_melt - H_solid_melt_fit) / H_solid_melt_fit

    # p_sub (in Pa in data; model comparison via ΔG closure)
    p_fitted_sub = []
    for T, pPa, Gf, Vf in zip(T_sub, p_sub, G_fluid_sub, V_fluid_sub):
        props = compute_thermo_props(T, pPa/1e6, params)   # MPa
        delta_G = Gf - props[11] + deltaH_triple - T*deltaS_triple
        pfPa = pPa - (delta_G / (Vf - props[0])) * 1e6
        p_fitted_sub.append(pfPa)
    p_fitted_sub = np.asarray(p_fitted_sub)
    rel_psub = (np.asarray(p_sub) - p_fitted_sub) / np.asarray(p_sub)

    # p_melt (split at 500 K)
    p_fitted_melt = []
    for T, pM, Gf, Vf in zip(T_melt, p_melt, G_fluid_melt, V_fluid_melt):
        props = compute_thermo_props(T, pM, params)
        delta_G = Gf - props[11] + deltaH_triple - T*deltaS_triple
        pfM = pM - (delta_G / (Vf - props[0]))
        p_fitted_melt.append(pfM)
    p_fitted_melt = np.asarray(p_fitted_melt)
    rel_pmelt = (np.asarray(p_melt) - p_fitted_melt) / np.asarray(p_melt)
    mask_pm_lo = np.asarray(T_melt) <= 500.0
    mask_pm_hi = ~mask_pm_lo

    # ---------- Low-T behaviour penalty (γ slope etc.) ----------
    T6 = np.array([0.0001] + list(range(2, 84, 2)) + [83.806])
    p6 = np.array([psub(T, pt, Tt) for T in T6])
    Gamma_T6 = np.array([compute_thermo_props(T, p, params)[6]
                        for T, p in zip(T6, p6)])
    Vm_T6 = np.array([compute_thermo_props(T, p, params)[0]
                     for T, p in zip(T6, p6)])
    slopes = np.diff(Gamma_T6) / np.diff(Vm_T6)

    # penalize positive slopes and excursions
    gam_pen = []
    gmean = np.mean(Gamma_T6)
    for i, s in enumerate(slopes):
        if s > 0:
            gam_pen.append(45.0 + abs(200.0*(Gamma_T6[i+1]-gmean)/gmean))
        else:
            gam_pen.append(50.0 * (Gamma_T6[i+1]-gmean)/gmean)
    gam_pen = np.asarray(gam_pen, float)
    lowT_penalty = _mean_sq(gam_pen)   # use mean square so it's χ²-like

    # ---------- χ² objective with thesis weights ----------
    chi2 = 0.0
    chi2 += W_VM_SUB * _mean_sq(rel_Vm_sub)
    chi2 += W_VM_MELT * _mean_sq(rel_Vm_melt)
    chi2 += W_VM_HIGHP * _mean_sq(rel_Vm_highp)

    chi2 += W_CP_LOW * _mean_sq(rel_cp[mask_cp_low])
    chi2 += W_CP_HIGH * _mean_sq(rel_cp[mask_cp_high])

    chi2 += W_ALPHA * _mean_sq(rel_alpha)
    chi2 += W_KT * _mean_sq(rel_KT)
    chi2 += W_KS * _mean_sq(rel_KS)

    chi2 += W_DH_SUB * _mean_sq(rel_DH_sub)
    chi2 += W_DH_MELT * _mean_sq(rel_DH_melt)

    chi2 += W_PSUB * _mean_sq(rel_psub)
    chi2 += W_PMELT_LO * _mean_sq(rel_pmelt[mask_pm_lo])
    chi2 += W_PMELT_HI * _mean_sq(rel_pmelt[mask_pm_hi])

    chi2 += W_LOW * lowT_penalty

    return chi2


def run_optimization2(params_init, bounds, datasets):
    print("Lens:", len(PARAMS_INIT), len(LOWER_BOUND), len(UPPER_BOUND))
    assert len(PARAMS_INIT) == len(LOWER_BOUND) == len(UPPER_BOUND)

    result = minimize(
        fun=combined_cost_function2,
        x0=np.asarray(params_init, float),
        args=datasets,
        method='L-BFGS-B',
        bounds=bounds,
        options={'disp': True, 'maxiter': 500, 'ftol': 1e-9}
    )
    return result.x, result.fun
def run_optimization(params_init, bounds, datasets):
    print("Lens:", len(PARAMS_INIT), len(LOWER_BOUND), len(UPPER_BOUND))
    assert len(PARAMS_INIT) == len(LOWER_BOUND) == len(UPPER_BOUND)

    result = minimize(
        fun=combined_cost_function,
        x0=params_init,
        args=datasets,
        method='L-BFGS-B',
        bounds=bounds,
        options={'disp': True, 'maxiter': 20, 'ftol': 1e-5}
    )
    return result.x, result.fun


def plot_fitted_vs_data(T_Vm_sub, Vm_sub, Year_Vm_sub, Author_Vm_sub,
                        T_Vm_melt, Vm_melt, Year_Vm_melt, Author_Vm_melt,
                        params_fit, mymarker, mycolor, fontsize=13.5, markersize=7, linewidth=0.9):
    """
    Plot fitted and experimental molar volume data for sublimation and melting.
    """
    # Calculate fitted properties across temperature range
    T_plot_calc = np.arange(1, 401)
    p_plot_calc = np.array([psub(T, pt, Tt) if T < 83.806 else pmelt(T, pt, Tt)
                           for T in T_plot_calc])
    fitted_props = np.array([compute_thermo_props(T, p, params_fit)
                            for T, p in zip(T_plot_calc, p_plot_calc)])

    # === Group data by Year and Author ===
    def group_by_reference(T_vals, V_vals, Year, Author):
        groups = []
        start_idx = 0
        for i in range(1, len(T_vals)):
            if Year[i] != Year[start_idx] or Author[i] != Author[start_idx]:
                groups.append((start_idx, i))
                start_idx = i
        groups.append((start_idx, len(T_vals)))
        return groups

    sub_groups = group_by_reference(
        T_Vm_sub, Vm_sub, Year_Vm_sub, Author_Vm_sub)
    melt_groups = group_by_reference(
        T_Vm_melt, Vm_melt, Year_Vm_melt, Author_Vm_melt)

    # === Plot ===
    plt.figure(figsize=(10, 6))
    for i, (start, end) in enumerate(sub_groups):
        plt.plot(T_Vm_sub[start:end], Vm_sub[start:end],
                 marker=mymarker[i % len(mymarker)],
                 markersize=markersize, linewidth=linewidth,
                 color=mycolor[i % len(mycolor)],
                 label=f"{Year_Vm_sub[start]}-{Author_Vm_sub[start]}")

    for i, (start, end) in enumerate(melt_groups):
        plt.plot(T_Vm_melt[start:end], Vm_melt[start:end],
                 marker=mymarker[(i + len(sub_groups)) % len(mymarker)],
                 markersize=markersize, linewidth=linewidth,
                 color=mycolor[(i + len(sub_groups)) % len(mycolor)],
                 label=f"{Year_Vm_melt[start]}-{Author_Vm_melt[start]}")

    # Plot fitted curve
    plt.plot(T_plot_calc, fitted_props[:, 0],
             'k--', linewidth=2.0, label='Fitted Vm')

    plt.xlabel("Temperature [K]", fontsize=fontsize)
    plt.ylabel("Molar Volume [cm³/mol]", fontsize=fontsize)
    plt.title("Experimental and Fitted Molar Volume vs Temperature",
              fontsize=fontsize + 1)
    plt.grid(True)
    plt.legend(fontsize=fontsize - 2)
    plt.tight_layout()
    plt.show()


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
    p_H_sub = np.array([psub(T, pt, Tt) for T in T_H_sub])
    delta_H_sub = data['heatsub']['Change in Enthalpy']
    H_fluid_sub = data['heatsub']['Enthalpy']

    # Enthalpy Melting #TODO
    T_H_melt = T_H_sub
    p_H_melt = p_H_sub
    delta_H_melt = delta_H_sub
    H_fluid_melt = H_fluid_sub

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


def plot_vm_vs_temp(T_Vm_sub, Vm_sub, Year_Vm_sub, Author_Vm_sub,
                    T_Vm_melt, Vm_melt, Year_Vm_melt, Author_Vm_melt,
                    params_fit, compute_thermo_props,
                    CUSTOMMARKERS, CUSTOMCOLORS, fontsize=13.5, markersize=7, linewidth=0.9):
    T_plot_calc = np.arange(1, 401)
    p_plot_calc = np.array([psub(T, pt, Tt) if T < 83.806 else pmelt(T, pt, Tt)
                           for T in T_plot_calc])
    fitted_props = np.array([compute_thermo_props(T, p, params_fit)
                            for T, p in zip(T_plot_calc, p_plot_calc)])

    plt.figure(figsize=(10, 6))
    for i in range(len(T_Vm_sub)):
        plt.plot(T_Vm_sub[i], Vm_sub[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(
                     CUSTOMCOLORS)], markersize=markersize, linestyle='None',
                 label=f"{Year_Vm_sub[i]}-{Author_Vm_sub[i]}" if i == 0 else "")
    for i in range(len(T_Vm_melt)):
        plt.plot(T_Vm_melt[i], Vm_melt[i], marker=CUSTOMMARKERS[(i + len(T_Vm_sub)) % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[(i + len(T_Vm_sub)) % len(CUSTOMCOLORS)
                                    ], markersize=markersize, linestyle='None',
                 label=f"{Year_Vm_melt[i]}-{Author_Vm_melt[i]}" if i == 0 else "")
    plt.plot(T_plot_calc, fitted_props[:, 0],
             'k--', linewidth=2.0, label='Fitted Vm')
    plt.xlabel("Temperature [K]", fontsize=fontsize)
    plt.ylabel("Molar Volume [cm³/mol]", fontsize=fontsize)
    plt.title("Experimental and Fitted Molar Volume vs Temperature",
              fontsize=fontsize + 1)
    plt.grid(True)
    plt.legend(fontsize=fontsize - 2)
    plt.tight_layout()
    plt.show()


def plot_vm_deviation(T_Vm, Vm, params_fit, p_func, compute_thermo_props, Year, Author, CUSTOMMARKERS, CUSTOMCOLORS, title):
    p_Vm = np.array([p_func(T) for T in T_Vm])
    Vm_calc = np.array([compute_thermo_props(T, p, params_fit)[0]
                       for T, p in zip(T_Vm, p_Vm)])
    deviation = 100 * (Vm - Vm_calc) / Vm
    plt.figure(figsize=(10, 6))
    for i in range(len(T_Vm)):
        plt.plot(T_Vm[i], deviation[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)], markersize=7, linestyle='None',
                 label=f"{Year[i]}-{Author[i]}" if i == 0 else "")
    plt.axhline(0, color='k', linewidth=0.7)
    plt.xlabel("Temperature [K]")
    plt.ylabel("100·(Vm_exp - Vm_calc) / Vm_calc [%]")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_vm_highp(T_Vm_highp, Vm_highp, p_Vm_highp, Year_Vm_highp, Author_Vm_highp,
                  params_fit, compute_thermo_props, CUSTOMMARKERS, CUSTOMCOLORS):
    plt.figure(figsize=(10, 6))
    for i in range(len(T_Vm_highp)):
        plt.plot(p_Vm_highp[i], Vm_highp[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)], markersize=7, linestyle='None',
                 label=f"{Year_Vm_highp[i]}-{Author_Vm_highp[i]}-{T_Vm_highp[i]:.1f}K" if i == 0 else "")
        p_range = np.linspace(min(p_Vm_highp[i]), max(p_Vm_highp[i]), 50)
        Vm_fit = [compute_thermo_props(T_Vm_highp[i], p, params_fit)[
            0] for p in p_range]
        plt.plot(p_range, Vm_fit, color=CUSTOMCOLORS[i % len(
            CUSTOMCOLORS)], linestyle='--')
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("Molar Volume [cm³/mol]")
    plt.title("High Pressure Cell Volume")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_vm_highp_deviation(T_Vm_highp, Vm_highp, p_Vm_highp, params_fit, compute_thermo_props, Year, Author, CUSTOMMARKERS, CUSTOMCOLORS, title):
    Vm_calc = np.array([compute_thermo_props(T, p, params_fit)[0]
                       for T, p in zip(T_Vm_highp, p_Vm_highp)])
    deviation = 100 * (Vm_highp - Vm_calc) / Vm_highp
    plt.figure(figsize=(10, 6))
    for i in range(len(T_Vm_highp)):
        plt.plot(p_Vm_highp[i], deviation[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)], markersize=7, linestyle='None',
                 label=f"{Year[i]}-{Author[i]}" if i == 0 else "")
    plt.axhline(0, color='k', linewidth=0.7)
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("100·(Vm_exp - Vm_calc) / Vm_calc [%]")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_property_vs_temp(T, prop_exp, prop_name, params_fit, p_func, compute_thermo_props, Year, Author, CUSTOMMARKERS, CUSTOMCOLORS, prop_index, title, ylabel):
    T_plot_calc = np.arange(1, 401)
    p_plot_calc = np.array([p_func(T) for T in T_plot_calc])
    fitted_props = np.array([compute_thermo_props(T, p, params_fit)
                            for T, p in zip(T_plot_calc, p_plot_calc)])
    plt.figure(figsize=(10, 6))
    for i in range(len(T)):
        plt.plot(T[i], prop_exp[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)], markersize=7, linestyle='None',
                 label=f"{Year[i]}-{Author[i]}" if i == 0 else "")
    plt.plot(T_plot_calc, fitted_props[:, prop_index],
             'k--', linewidth=2.0, label=f'Fitted {prop_name}')
    plt.xlabel("Temperature [K]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_property_deviation(T, prop_exp, params_fit, p_func, compute_thermo_props, prop_index, Year, Author, CUSTOMMARKERS, CUSTOMCOLORS, title, ylabel):
    p_vals = np.array([p_func(Ti) for Ti in T])
    prop_calc = np.array([compute_thermo_props(Ti, pi, params_fit)[
                         prop_index] for Ti, pi in zip(T, p_vals)])
    deviation = 100 * (prop_exp - prop_calc) / prop_exp
    plt.figure(figsize=(10, 6))
    for i in range(len(T)):
        plt.plot(T[i], deviation[i], marker=CUSTOMMARKERS[i % len(CUSTOMMARKERS)],
                 color=CUSTOMCOLORS[i % len(CUSTOMCOLORS)], markersize=7, linestyle='None',
                 label=f"{Year[i]}-{Author[i]}" if i == 0 else "")
    plt.axhline(0, color='k', linewidth=0.7)
    plt.xlabel("Temperature [K]")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()
# Read all data into 1 dataframe
# krypton_data = load_all_gas_data('krypton', read_from_excel=False)
# xenon_data = load_all_gas_data('xenon', read_from_excel=False)
neon_data = load_all_gas_data('neon', read_from_excel=False)


bounds = list(zip(LOWER_BOUND, UPPER_BOUND))
datasets = extract_datasets(neon_data)

params_fit, fval = run_optimization2(PARAMS_INIT, bounds, datasets)

print("Optimized parameters:", params_fit)
print("Final objective value:", fval)
# datasets["T"]
T_Vm_sub = neon_data["cell_volume_sub"]["Temperature"]
Vm_sub = neon_data["cell_volume_sub"]["Cell Volume"]
Year_Vm_sub = neon_data["cell_volume_sub"]["Year"]
Author_Vm_sub = neon_data["cell_volume_sub"]["Author"]

T_Vm_melt = neon_data["cell_volume_melt"]["Temperature"]
Vm_melt = neon_data["cell_volume_melt"]["Cell Volume"]
Year_Vm_melt = neon_data["cell_volume_melt"]["Year"]
Author_Vm_melt = neon_data["cell_volume_melt"]["Author"]

PARAM_LABELS = [
    ("c1", "MPa"),
    ("c2", "MPa"),
    ("c3", "MPa"),
    ("Theta_D,0", "K"),
    ("gamma_D,0", ""),
    ("q_D", ""),
    ("b1", ""),
    ("b2", ""),
    ("b3", ""),
    ("S_m(g, T_t, p_t)", "J mol-1 K-1"),
]

# If your vector is longer/shorter, we’ll truncate/pad labels to match length:
n = len(params_fit)
labels = PARAM_LABELS[:n] if len(PARAM_LABELS) >= n else PARAM_LABELS + [
    ("param_"+str(i+1), "") for i in range(len(PARAM_LABELS), n)]

df_params = pd.DataFrame({
    "Parameter": [lbl for (lbl, _) in labels],
    "Value": [float(v) for v in params_fit[:len(labels)]],
    "Unit": [unit for (_, unit) in labels]
})

# Nice formatting (no scientific notation unless needed)
with pd.option_context('display.float_format', lambda x: f"{x:.6g}"):
    print("\nOptimized parameter values for the solid EOS\n")
    print(df_params.to_string(index=False))
