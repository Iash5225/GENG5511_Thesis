from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from thermal_script import plot_all_gas_properties
from computethermoprops import *
from p_functions import pmelt,psub
from constants import CUSTOMCOLORS,CUSTOMMARKERS,PARAMS_INIT,LOWER_BOUND,UPPER_BOUND


# Constantsn 
St_REFPROP = 131.1500
Ht_REFPROP = 1689.10
Tt = 83.8058
pt = 0.068891  # MPa




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
    def rms(x): return np.sqrt(np.sum(np.square(x)) / len(x))

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
    p6 = np.array([psub(T) for T in T6])
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


def run_optimization(params_init, bounds, datasets):
    result = minimize(
        fun=combined_cost_function,
        x0=params_init,
        args=datasets,
        method='L-BFGS-B',
        bounds=bounds,
        options={'disp': True, 'maxiter': 200000, 'ftol': 1e-25}
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
    p_plot_calc = np.array([psub(T) if T < 83.806 else pmelt(T)
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
        psub (callable): Function to calculate sublimation pressure from temperature.
        pmelt (callable): Function to calculate melting pressure from temperature.

    Returns:
        tuple: A tuple containing arrays for all datasets in the specified order.
    """

    # Cell Volume Sublimation
    T_Vm_sub = data["cell_volume_sub"]['Temperature']
    p_Vm_sub = np.array([psub(T) for T in T_Vm_sub])
    Vm_sub = data['cell_volume_sub']['Cell Volume']

    # Cell Volume Melting
    T_Vm_melt = data['cell_volume_melt']['Temperature']
    p_Vm_melt = np.array([pmelt(T) for T in T_Vm_melt])
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
    p_cp_sub = np.array([psub(T) for T in T_cp_sub])
    cp_sub = data['heat_capacity']['Heat Capacity']

    # Thermal Expansion Sublimation
    T_alpha_sub = data['thermal_coeff']['Temperature']
    p_alpha_sub = np.array([psub(T) for T in T_alpha_sub])
    alpha_sub = data['thermal_coeff']['Thermal Expansion Coefficient']

    # Empty placeholders
    T_BetaT_sub = p_BetaT_sub = BetaT_sub = np.array([])
    T_BetaS_sub = p_BetaS_sub = BetaS_sub = np.array([])
    T_sub = p_sub = G_fluid_sub = V_fluid_sub = np.array([])
    T_melt = p_melt = G_fluid_melt = V_fluid_melt = np.array([])
    T_H_sub = p_H_sub = delta_H_sub = H_fluid_sub = np.array([])
    T_H_melt = p_H_melt = delta_H_melt = H_fluid_melt = np.array([])
    Year_sub = np.array([])   # <-- add this

    datasets = (
        T_Vm_sub, p_Vm_sub, Vm_sub,
        T_Vm_melt, p_Vm_melt, Vm_melt,
        T_Vm_highp, p_Vm_highp, Vm_highp,
        T_cp_sub, p_cp_sub, cp_sub,
        T_alpha_sub, p_alpha_sub, alpha_sub,
        T_BetaT_sub, p_BetaT_sub, BetaT_sub,
        T_BetaS_sub, p_BetaS_sub, BetaS_sub,
        T_sub, p_sub, G_fluid_sub, V_fluid_sub, Year_sub,
        T_melt, p_melt, G_fluid_melt, V_fluid_melt,
        T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
        T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt
    )

    return datasets


# Read all data into 1 dataframe
# krypton_data = load_all_gas_data('krypton', read_from_excel=False)
# xenon_data = load_all_gas_data('xenon', read_from_excel=False)
neon_data = load_all_gas_data('neon', read_from_excel=False)


bounds = list(zip(LOWER_BOUND, UPPER_BOUND))
datasets = extract_datasets(neon_data)

params_fit, fval = run_optimization(PARAMS_INIT, bounds, datasets)

print("Optimized parameters:", params_fit)
print("Final objective value:", fval)
# datasets["T"]
# plot_fitted_vs_data(
#     T_Vm_sub, Vm_sub, Year_Vm_sub, Author_Vm_sub,
#     T_Vm_melt, Vm_melt, Year_Vm_melt, Author_Vm_melt,
#     params_fit, CUSTOMMARKERS, CUSTOMCOLORS
# )
