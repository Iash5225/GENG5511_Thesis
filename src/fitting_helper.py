from datetime import datetime
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from thermal_script import *
from thermopropsv2 import compute_thermo_props
from constants import *
from scipy.interpolate import UnivariateSpline
import os
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)
# solid_EOS_fitting.py
BIG = 1e4
Vm_anchor = 0.8
T_MARGIN = 2.0  # K, stay below Tt by this margin
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t

def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB,  KRYPTON_E_3_SUB,  KRYPTON_T_t, KRYPTON_P_t
)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t
)


def plot_variable_deviation(
    data,
    gas_name: str,
    x_col: str,
    y_exp_col: str,
    y_model_col: str,
    y_label: str,
    title: str,
    filename: str = None,
    xlim=None,
    ylim=None,
    markersize=50,
    custom_markers=None,
    custom_colors=None,
    legend_outside=True,
    fontsize=14,
    output_folder=None,
):
    """
    Plot percent deviation between experimental and model values for any variable.

    Args:
        data (pd.DataFrame): DataFrame with columns including x_col, y_exp_col, y_model_col, 'Year', 'Author'
        gas_name (str): Gas name for title and filename
        x_col (str): Column name for x-axis (e.g., 'Temperature')
        y_exp_col (str): Experimental value column (e.g., 'y_exp')
        y_model_col (str): Model value column (e.g., 'y_model')
        y_label (str): Y-axis label (LaTeX string)
        title (str): Plot title
        filename (str): Output filename (optional, .png will be appended if missing)
        xlim, ylim: Axis limits (optional)
        markersize (int): Marker size
        custom_markers (list): List of marker styles
        custom_colors (list): List of colors
        legend_outside (bool): Place legend outside plot
        fontsize (int): Font size for labels
        output_folder (str): If given, save plot in this folder
    """
    data = data.copy()
    data['Deviation_percent'] = 100 * \
        (data[y_exp_col] - data[y_model_col]) / data[y_exp_col]
    grouped = data.groupby(['Year', 'Author'])
    if custom_markers is None:
        custom_markers = ['o', 's', '^', 'v', 'D', 'x', '*', 'P', 'h', 'X']
    if custom_colors is None:
        custom_colors = plt.cm.tab10.colors

    # Compute percent deviation


    plt.figure(figsize=(10, 6))
    plt.tick_params(direction='in', top=True, right=True)
    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group[x_col], group['Deviation_percent'],
            label=f"{year}, {author}",
            s=markersize,
            edgecolor=custom_colors[i % len(custom_colors)],
            marker=custom_markers[i % len(custom_markers)],
            facecolors='none',
        )
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel(r'$\mathit{T}$ / K', fontsize=fontsize)
    plt.ylabel(y_label, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    if xlim:
        plt.xlim(*xlim)
    if ylim:
        plt.ylim(*ylim)
    if legend_outside:
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        plt.legend(fontsize=8)
        plt.tight_layout()

    # Save if requested
    if filename:
        if not filename.lower().endswith('.png'):
            filename += '.png'
        if output_folder:
            os.makedirs(output_folder, exist_ok=True)
            filename = os.path.join(output_folder, filename)
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()


def plot_thermo_variable(
    data,
    gas_name: str,
    x_col: str,
    y_col: str,
    y_label: str,
    title: str,
    model_x=None,
    model_y=None,
    logy=False,
    xlim=None,
    ylim=None,
    filename=None,
    output_folder=None,  # <-- add this line
    legend_outside=True,
    markersize=50,
    linewidth=2,
    fontsize=14,
    custom_markers=None,
    custom_colors=None,
):
    """
    General-purpose plot for thermodynamic variables (e.g., enthalpy, heat capacity, etc.)
    grouped by author/year, with optional model curve.

    Args:
        data (pd.DataFrame): DataFrame with columns including x_col, y_col, 'Year', 'Author'
        gas_name (str): Gas name for title and filename
        x_col (str): Column name for x-axis (e.g., 'Temperature')
        y_col (str): Column name for y-axis (e.g., 'Change in Enthalpy')
        y_label (str): Y-axis label (LaTeX string)
        title (str): Plot title
        model_x (array-like): X values for model curve (optional)
        model_y (array-like): Y values for model curve (optional)
        logy (bool): Use log scale for y-axis
        xlim, ylim: Axis limits (optional)
        filename (str): If given, save plot to this file
        legend_outside (bool): Place legend outside plot
        markersize (int): Marker size
        linewidth (int): Model curve line width
        fontsize (int): Font size for labels
        custom_markers (list): List of marker styles
        custom_colors (list): List of colors
    """
    grouped = data.groupby(['Year', 'Author'])
    if custom_markers is None:
        custom_markers = ['o', 's', '^', 'v', 'D', 'x', '*', 'P', 'h', 'X']
    if custom_colors is None:
        custom_colors = plt.cm.tab10.colors

    plt.figure(figsize=(10, 6))
    for i, ((year, author), group) in enumerate(grouped):
        plt.scatter(
            group[x_col], group[y_col],
            s=markersize,
            label=f"{year}, {author}",
            edgecolors=custom_colors[i % len(custom_colors)],
            facecolors='none',
            linewidths=1.5,
            marker=custom_markers[i % len(custom_markers)],
        )

    # Plot model curve if provided
    if model_x is not None and model_y is not None:
        order = np.argsort(model_x)
        plt.plot(
            np.array(model_x)[order], np.array(model_y)[order],
            color='black', linewidth=linewidth, label='Model'
        )

    plt.xlabel(r'$\mathit{T}$ / K', fontsize=fontsize)
    plt.ylabel(y_label, fontsize=fontsize)
    plt.title(title, fontsize=fontsize)
    if logy:
        plt.yscale('log')
    if xlim:
        plt.xlim(*xlim)
    if ylim:
        plt.ylim(*ylim)

    if legend_outside:
        plt.legend(
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            fontsize=8,
            ncol=1
        )
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        plt.legend(fontsize=8)
        plt.tight_layout()

    if filename:
        if not filename.lower().endswith('.png'):
            filename += '.png'
        if output_folder:
            os.makedirs(output_folder, exist_ok=True)
            filename = os.path.join(output_folder, filename)
        plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

def safe_psub(T):
    """
    Real, positive sublimation pressure; NaN outside domain.
    """
    p = np.asarray(psub_curve(T))
    if np.iscomplexobj(p):
        near_real = np.abs(p.imag) < 1e-12
        p = np.where(near_real, p.real, np.nan)
    p = p.astype(float, copy=False)
    p[~np.isfinite(p) | (p <= 0.0)] = np.nan
    return p


def rms_log(y_true, y_pred):
    """
    Root-mean-square error in log space.

    Computes sqrt(mean((ln(y_pred) - ln(y_true))^2)) over finite, positive pairs.
    Useful when errors are multiplicative or span several orders of magnitude.

    Parameters
    ----------
    y_true, y_pred : array-like
        True and predicted values (must be > 0 to be included).

    Returns
    -------
    float
        RMS of log-differences. Returns BIG if no valid pairs remain.
    """
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    m = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true > 0) & (y_pred > 0)
    if not np.any(m):
        return BIG
    r = np.log(y_pred[m]) - np.log(y_true[m])
    return float(np.sqrt(np.mean(r*r)))

def rms_percent(y_true, y_pred, scale=PERCENT_SCALE):
    """
    Percent RMS error.

    Computes sqrt(mean((scale * (y_true - y_pred) / y_true)^2)) over finite, nonzero y_true.
    With `scale=100`, this becomes percent RMS (%).

    Parameters
    ----------
    y_true, y_pred : array-like
        True and predicted values.
    scale : float, optional
        Multiplicative factor applied to the relative error (default from constants: 100 for %).

    Returns
    -------
    float
        Percent RMS (same units as `scale`). Returns BIG if no valid pairs remain.
    """
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true != 0.0)
    if not np.any(mask):
        return BIG
    err = scale * (y_true[mask] - y_pred[mask]) / y_true[mask]
    err = err[np.isfinite(err)]
    return BIG if err.size == 0 else float(np.sqrt(np.mean(err*err)))

def rms(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    return 0.0 if n == 0 else np.sqrt(np.sum(x*x) / n)

def extract_datasets_with_meta(data):
    """
    Build both numeric arrays (like extract_datasets) AND metadata (Author, Year).
    Returns:
        datasets : tuple (length 38, same order as extract_datasets)
        meta : dict mapping each property to its Author/Year arrays
    """
    # --- reuse most of the original code ---
    T_Vm_sub = data["cell_volume_sub"]["Temperature"].to_numpy()
    p_Vm_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t
    ) for T in T_Vm_sub], dtype=float)
    Vm_sub = data["cell_volume_sub"]["Cell Volume"].to_numpy()

    T_Vm_melt = data["cell_volume_melt"]["Temperature"].to_numpy()
    p_Vm_melt = np.array([melting_pressure_equation(
        T, KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t
    ) for T in T_Vm_melt], dtype=float)
    Vm_melt = data["cell_volume_melt"]["Cell Volume"].to_numpy()

    if "cell_volume_highp" in data:
        T_Vm_highp = data["cell_volume_highp"]["Temperature"].to_numpy()
        p_Vm_highp = data["cell_volume_highp"]["Pressure"].to_numpy()
        Vm_highp = data["cell_volume_highp"]["Cell Volume"].to_numpy()
    else:
        T_Vm_highp = p_Vm_highp = Vm_highp = np.array([])

    T_cp_sub = data["heat_capacity"]["Temperature"].to_numpy()
    p_cp_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t
    ) for T in T_cp_sub], dtype=float)
    cp_sub = data["heat_capacity"]["Heat Capacity"].to_numpy()

    T_alpha_sub = data["thermal_coeff"]["Temperature"].to_numpy()
    p_alpha_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t
    ) for T in T_alpha_sub], dtype=float)
    alpha_sub = data["thermal_coeff"]["Thermal Expansion Coefficient"].to_numpy()

    T_BetaS_sub = data["bulk_s"]["Temperature"].to_numpy()
    p_BetaS_sub = safe_psub(T_BetaS_sub)
    BetaS_sub = data["bulk_s"]["Beta S"].to_numpy()

    T_BetaT_sub = data["bulk_t"]["Temperature"].to_numpy()
    p_BetaT_sub = safe_psub(T_BetaT_sub)
    BetaT_sub = data["bulk_t"]["Beta T"].to_numpy()

    T_melt = data["melting"]["Temperature"].to_numpy()
    p_melt = data["melting"]["Pressure"].to_numpy()
    G_fluid_melt = data["melting"]["Gibbs Energy"].to_numpy()
    V_fluid_melt = data["melting"]["Volume"].to_numpy()

    T_sub = data["sublimation"]["Temperature"].to_numpy()
    p_sub = data["sublimation"]["Pressure"].to_numpy()
    G_fluid_sub = data["sublimation"]["Gibbs Energy"].to_numpy()
    V_fluid_sub = data["sublimation"]["Volume"].to_numpy()
    Year_sub = data["sublimation"]["Year"].to_numpy(
    ) if "Year" in data["sublimation"].columns else np.array([])

    T_H_sub = data["heatsub"]["Temperature"].to_numpy()
    p_H_sub = data["heatsub"]["Pressure"].to_numpy()
    delta_H_sub = pd.to_numeric(data["heatsub"]["Change in Enthalpy"],
                      errors="coerce").to_numpy()
    H_fluid_sub = pd.to_numeric(data["heatsub"]["Enthalpy"], errors="coerce").to_numpy()

    T_H_melt = data["fusion"]["Temperature"].to_numpy()
    p_H_melt = data["fusion"]["Pressure"].to_numpy()
    delta_H_melt = pd.to_numeric(data["fusion"]["Change in Enthalpy"],
                      errors="coerce").to_numpy()
    H_fluid_melt = pd.to_numeric(data["fusion"]["Enthalpy"], errors="coerce").to_numpy()

    # --- original tuple (38) ---
    datasets = (
        T_Vm_sub, p_Vm_sub, Vm_sub,
        T_Vm_melt, p_Vm_melt, Vm_melt,
        T_Vm_highp, p_Vm_highp, Vm_highp,
        T_cp_sub, p_cp_sub, cp_sub,
        T_alpha_sub, p_alpha_sub, alpha_sub,
        T_BetaT_sub, p_BetaT_sub, BetaT_sub,
        T_BetaS_sub, p_BetaS_sub, BetaS_sub,
        T_sub, p_sub, Year_sub, G_fluid_sub, V_fluid_sub,
        T_melt, p_melt, G_fluid_melt, V_fluid_melt,
        T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
        T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt
    )

    # --- metadata dict ---
    meta = {
        "Vm_sub":   {"Author": data["cell_volume_sub"].get("Author", pd.Series(["Unknown"]*len(T_Vm_sub))).to_numpy(),
                     "Year":   data["cell_volume_sub"].get("Year", pd.Series([0]*len(T_Vm_sub))).to_numpy()},
        "Vm_melt":  {"Author": data["cell_volume_melt"].get("Author", pd.Series(["Unknown"]*len(T_Vm_melt))).to_numpy(),
                     "Year":   data["cell_volume_melt"].get("Year", pd.Series([0]*len(T_Vm_melt))).to_numpy()},
        "Vm_highp": None if "cell_volume_highp" not in data else {
            "Author": data["cell_volume_highp"].get("Author", pd.Series(["Unknown"]*len(T_Vm_highp))).to_numpy(),
            "Year":   data["cell_volume_highp"].get("Year", pd.Series([0]*len(T_Vm_highp))).to_numpy()},
        "cp_sub":   {"Author": data["heat_capacity"].get("Author", pd.Series(["Unknown"]*len(T_cp_sub))).to_numpy(),
                     "Year":   data["heat_capacity"].get("Year", pd.Series([0]*len(T_cp_sub))).to_numpy()},
        "alpha_sub": {"Author": data["thermal_coeff"].get("Author", pd.Series(["Unknown"]*len(T_alpha_sub))).to_numpy(),
                      "Year":   data["thermal_coeff"].get("Year", pd.Series([0]*len(T_alpha_sub))).to_numpy()},
        "BetaT_sub": {"Author": data["bulk_t"].get("Author", pd.Series(["Unknown"]*len(T_BetaT_sub))).to_numpy(),
                      "Year":   data["bulk_t"].get("Year", pd.Series([0]*len(T_BetaT_sub))).to_numpy()},
        "BetaS_sub": {"Author": data["bulk_s"].get("Author", pd.Series(["Unknown"]*len(T_BetaS_sub))).to_numpy(),
                      "Year":   data["bulk_s"].get("Year", pd.Series([0]*len(T_BetaS_sub))).to_numpy()},
        "H_solid_sub": {"Author": data["heatsub"].get("Author", pd.Series(["Unknown"]*len(T_H_sub))).to_numpy(),
                        "Year":   data["heatsub"].get("Year", pd.Series([0]*len(T_H_sub))).to_numpy()},
        "H_solid_melt": {"Author": data["fusion"].get("Author", pd.Series(["Unknown"]*len(T_H_melt))).to_numpy(),
                         "Year":   data["fusion"].get("Year", pd.Series([0]*len(T_H_melt))).to_numpy()},
        "sublimation": {
            "Author": data["sublimation"].get("Author", pd.Series(["Unknown"]*len(T_sub))).to_numpy(),
            "Year":   data["sublimation"].get("Year",   pd.Series([0]*len(T_sub))).to_numpy(),
        },
        "melting": {
            "Author": data["melting"].get("Author", pd.Series(["Unknown"]*len(T_melt))).to_numpy(),
            "Year":   data["melting"].get("Year",   pd.Series([0]*len(T_melt))).to_numpy(),
        },
    }

    return datasets, meta


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
    p_Vm_sub = safe_psub(T_Vm_sub)
    # p_Vm_sub = np.array([sublimation_pressure_equation(
    #     T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_Vm_sub])
    Vm_sub = data['cell_volume_sub']['Cell Volume']

    # Cell Volume Melting
    T_Vm_melt = data['cell_volume_melt']['Temperature']
    p_Vm_melt = pmelt_curve(T_Vm_melt)
    # p_Vm_melt = np.array([melting_pressure_equation(
    #     T, KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t) for T in T_Vm_melt])
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
    p_cp_sub = safe_psub(T_cp_sub)
    # p_cp_sub = np.array([sublimation_pressure_equation(
    #     T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_cp_sub])
    cp_sub = data['heat_capacity']['Heat Capacity']

    # Thermal Expansion Sublimation
    T_alpha_sub = data['thermal_coeff']['Temperature']
    p_alpha_sub = safe_psub(T_alpha_sub)
    # p_alpha_sub = np.array([sublimation_pressure_equation(
    #     T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_alpha_sub])
    alpha_sub = data['thermal_coeff']['Thermal Expansion Coefficient']

    # Bulk Modulus (S)
    T_BetaS_sub = data['bulk_s']['Temperature']
    p_BetaS_sub = safe_psub(T_BetaS_sub)
    # p_BetaS_sub = data['bulk_s']['Pressure']
    # p_BetaS_sub = None
    BetaS_sub = data['bulk_s']['Beta S']

    # Bulk Modulus (T)
    T_BetaT_sub = data['bulk_t']['Temperature']
    p_BetaT_sub = safe_psub(T_BetaT_sub)
    BetaT_sub = data['bulk_t']['Beta T']
    # print("T_BetaT_sub:", T_BetaT_sub)
    # print("p_BetaT_sub:", p_BetaT_sub)
    # print("BetaT_sub:", BetaT_sub)

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
    delta_H_sub = 1000.0 * data['heatsub']['Change in Enthalpy']  # kJ → J
    H_fluid_sub = 1000.0 * data['heatsub']['Enthalpy']

    # Enthalpy Melting
    T_H_melt = data['fusion']['Temperature']
    p_H_melt = data['fusion']['Pressure']
    delta_H_melt = 1000.0 * data['fusion']['Change in Enthalpy']
    H_fluid_melt = 1000.0 * data['fusion']['Enthalpy']

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

def _safe_props_vector(T_arr, p_arr, params, idx):
    """Return model property vector (same length) with NaN where compute fails."""
    out = np.full_like(np.asarray(T_arr, float), np.nan, dtype=float)
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        try:
            props = compute_thermo_props(float(T), float(p), params)
            if np.all(np.isfinite(props)):
                out[i] = props[idx]
        except Exception:
            pass
    return out


def _rmse_percent(y_exp, y_mod, floor=1e-6):
    """
    Relative RMS (%) = sqrt(mean( [100*(y_exp - y_mod)/max(|y_exp|,floor)]^2 ))
    """
    y_exp = np.asarray(y_exp, float)
    y_mod = np.asarray(y_mod, float)
    den = np.maximum(np.abs(y_exp), floor)
    rel = 100.0 * (y_exp - y_mod) / den
    return float(np.sqrt(np.mean(rel**2)))
def combined_cost_vm(params, *datasets, Tt=None, weights=None):
    """
    Vm-only objective using *percentage RMSE* (in %), no anchors/teacher/shape terms.

    total = W["SUB"] * %RMSE(Vm on sub) + W["MELT"] * %RMSE(Vm on melt)
    """
    (T_Vm_sub,  p_Vm_sub,  Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt, *_) = datasets

    Tt_val = float(Tt if Tt is not None else KRYPTON_T_t)

    # default weights (both count equally by default)
    W = dict(SUB=1.0, MELT=1.0)
    if weights:
        W.update(weights)

    _cache = {}

    # --- %RMSE on sublimation branch (T <= Tt) ---
    Vm_sub_pct = BIG
    msub = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub) &
            np.isfinite(Vm_sub) & (T_Vm_sub <= Tt_val))
    if np.any(msub):
        Vmod_sub = _safe_props_vector_cached(
            T_Vm_sub[msub], p_Vm_sub[msub], params, idx=IDX["Vm"], cache=_cache)
        Vm_sub_pct = _rmse_percent(Vm_sub[msub], Vmod_sub)

    # --- %RMSE on melting branch (T >= Tt) ---
    Vm_melt_pct = BIG
    mmelt = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt) &
             np.isfinite(Vm_melt) & (T_Vm_melt >= Tt_val))
    if np.any(mmelt):
        Vmod_melt = _safe_props_vector_cached(
            T_Vm_melt[mmelt], p_Vm_melt[mmelt], params, idx=IDX["Vm"], cache=_cache)
        Vm_melt_pct = _rmse_percent(Vm_melt[mmelt], Vmod_melt)

    total = W["SUB"] * Vm_sub_pct + W["MELT"] * Vm_melt_pct

    deviations = dict(
        Vm_sub_pct=Vm_sub_pct,
        Vm_melt_pct=Vm_melt_pct,
    )
    return total, deviations

def _cost_only_vm(params, *datasets):
    try:
        total, _ = combined_cost_vm(params, *datasets)
        return float(total) if np.isfinite(total) else BIG
    except Exception:
        return BIG

def quick_plot_vm_only(params, datasets, Tt):
    (T_Vm_sub, p_Vm_sub, Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt, *_) = datasets

    # --- convert to NumPy early ---
    x_sub = np.asarray(T_Vm_sub,  dtype=float)
    p_sub = np.asarray(p_Vm_sub,  dtype=float)
    y_sub = np.asarray(Vm_sub,    dtype=float)

    x_melt = np.asarray(T_Vm_melt, dtype=float)
    p_melt = np.asarray(p_Vm_melt, dtype=float)
    y_melt = np.asarray(Vm_melt,   dtype=float)

    # model predictions (NumPy arrays)
    y_sub_mod = np.asarray(_safe_props_vector(
        x_sub,  p_sub,  params, idx=IDX["Vm"]), float)
    y_melt_mod = np.asarray(_safe_props_vector(
        x_melt, p_melt, params, idx=IDX["Vm"]), float)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)

    # --- Sublimation (T <= Tt) ---
    msub = np.isfinite(x_sub) & (x_sub <= Tt) & np.isfinite(
        y_sub) & np.isfinite(y_sub_mod)
    ax[0].scatter(x_sub[msub], y_sub[msub], s=20, alpha=0.6, label="Exp (sub)")
    order_sub = np.argsort(x_sub[msub])
    ax[0].plot(x_sub[msub][order_sub], y_sub_mod[msub][order_sub],
               color="red", lw=2, label="Model (sub)")
    ax[0].axvline(Tt, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    # --- Melting (T >= Tt) ---
    mmelt = np.isfinite(x_melt) & (x_melt >= Tt) & np.isfinite(
        y_melt) & np.isfinite(y_melt_mod)
    ax[1].scatter(x_melt[mmelt], y_melt[mmelt], s=20,
                  alpha=0.6, label="Exp (melt)")
    order_melt = np.argsort(x_melt[mmelt])
    ax[1].plot(x_melt[mmelt][order_melt], y_melt_mod[mmelt][order_melt],
               color="red", lw=2, label="Model (melt)")
    ax[1].axvline(Tt, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()


def extract_pointwise_datasets(data, property_map=None):
    """
    Build a tidy DataFrame of all experimental points with metadata.

    Parameters
    ----------
    data : dict
        Dict of property DataFrames (same input you give to extract_datasets).
    property_map : dict, optional
        Maps dataset keys -> Property names for reporting (e.g., "cell_volume_sub" -> "Vm_sub").

    Returns
    -------
    df : pd.DataFrame
        Columns: [Property, Source, Year, T, y_exp, p, Author, used]
    """
    if property_map is None:
        property_map = {
            "cell_volume_sub": "Vm_sub",
            "cell_volume_melt": "Vm_melt",
            "cell_volume_highp": "Vm_highp",
            "heat_capacity": "cp_sub",
            "thermal_coeff": "alpha_sub",
            "bulk_t": "BetaT_sub",
            "bulk_s": "BetaS_sub",
            "heatsub": "H_solid_sub",
            "fusion": "H_solid_melt",
        }

    rows = []

    for key, prop_name in property_map.items():
        if key not in data:
            continue
        df_block = data[key].copy()

        # Ensure numeric
        T = pd.to_numeric(df_block.get(
            "Temperature", pd.Series([], dtype=float)), errors="coerce")
        y = None
        if "Cell Volume" in df_block:
            y = pd.to_numeric(df_block["Cell Volume"], errors="coerce")
        elif "Heat Capacity" in df_block:
            y = pd.to_numeric(df_block["Heat Capacity"], errors="coerce")
        elif "Thermal Expansion Coefficient" in df_block:
            y = pd.to_numeric(
                df_block["Thermal Expansion Coefficient"], errors="coerce")
        elif "Beta T" in df_block:
            y = pd.to_numeric(df_block["Beta T"], errors="coerce")
        elif "Beta S" in df_block:
            y = pd.to_numeric(df_block["Beta S"], errors="coerce")
        elif "Change in Enthalpy" in df_block:
            y = pd.to_numeric(df_block["Change in Enthalpy"], errors="coerce")

        authors = df_block.get("Author", ["Unknown"]*len(T))
        years = df_block.get("Year", [0]*len(T))

        for Ti, yi, auth, yr in zip(T, y, authors, years):
            rows.append(dict(
                Property=prop_name,
                Source=str(auth),
                Year=int(yr) if str(yr).isdigit() else str(yr),
                T=float(Ti) if np.isfinite(Ti) else np.nan,
                y_exp=float(yi) if np.isfinite(yi) else np.nan,
                used=True   # placeholder; you can set phase-condition masks later
            ))

    return pd.DataFrame(rows)


def build_master_pointwise_df(datasets, meta, params_fit):
    """
    Returns one long DataFrame with columns:
      Property, Author, Year, T, p_exp, p_sub_curve, p_melt_curve, dp_to_curve, abs_dp_over_p,
      y_exp, y_model
    For sub-path properties, dp_to_curve compares p_exp to p_sub_curve.
    For melt-path properties, dp_to_curve compares p_exp to p_melt_curve.
    For high-p properties, dp_to_curve is NaN (no canonical curve).
    """
    (
        T_Vm_sub,   p_Vm_sub,   Vm_sub,
        T_Vm_melt,  p_Vm_melt,  Vm_melt,
        T_Vm_highp, p_Vm_highp, Vm_highp,
        T_cp_sub,   p_cp_sub,   cp_sub,
        T_alpha_sub, p_alpha_sub, alpha_sub,
        T_BetaT_sub, p_BetaT_sub, BetaT_sub,
        T_BetaS_sub, p_BetaS_sub, BetaS_sub,
        T_sub, p_sub, Year_sub, G_fluid_sub, V_fluid_sub,
        T_melt, p_melt, G_fluid_melt, V_fluid_melt,
        T_H_sub, p_H_sub, delta_H_sub, H_fluid_sub,
        T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt
    ) = datasets

    rows = []

    def add_block(prop_name, phase,           # phase: "sub", "melt", or None
                  T_arr, p_arr, y_exp_arr, idx_model, meta_key,
                  y_exp_transform=None):
        if T_arr is None or len(T_arr) == 0:
            return
        authors = meta[meta_key]["Author"] if meta.get(
            meta_key) is not None else np.array(["Unknown"]*len(T_arr))
        years = meta[meta_key]["Year"] if meta.get(
            meta_key) is not None else np.array([0]*len(T_arr))

        # Precompute curves
        p_sub_curve = psub_curve(T_arr)  # sublimation curve in Pa
        p_melt_curve = pmelt_curve(T_arr)

        for Ti, Pi, yexp, auth, yr, psub_i, pmelt_i in zip(T_arr, p_arr, y_exp_arr, authors, years, p_sub_curve, p_melt_curve):
            if not (np.isfinite(Ti) and np.isfinite(Pi) and np.isfinite(yexp)):
                continue

            props = compute_thermo_props(Ti, Pi, params_fit)
            y_model = props[idx_model] if np.all(
                np.isfinite(props)) else np.nan
            if y_exp_transform is not None:
                # e.g., for enthalpy derived values
                yexp_use = y_exp_transform(Ti, Pi, yexp)
            else:
                yexp_use = yexp

            # choose the “relevant” curve for delta-p
            if phase == "sub":
                dp_to_curve = Pi - psub_i  # Pa vs MPa
                denom = psub_i
            elif phase == "melt":
                dp_to_curve = Pi - pmelt_i
                denom = pmelt_i
            else:
                dp_to_curve = np.nan
                denom = np.nan

            abs_dp_over_p = np.abs(dp_to_curve) / np.abs(denom) if np.isfinite(
                dp_to_curve) and np.isfinite(denom) and denom != 0 else np.nan

            rows.append(dict(
                Property=prop_name,
                Author=str(auth),
                Year=int(yr) if str(yr).isdigit() else yr,
                T=float(Ti),
                p_exp=float(Pi),
                p_sub_curve=float(psub_i),
                p_melt_curve=float(pmelt_i),
                dp_to_curve=float(dp_to_curve) if np.isfinite(
                    dp_to_curve) else np.nan,
                abs_dp_over_p=float(abs_dp_over_p) if np.isfinite(
                    abs_dp_over_p) else np.nan,
                y_exp=float(yexp_use),
                y_model=float(y_model)
            ))
        # ===== NEW: pure pressure blocks =====
    # Treat pressure as the "quantity": y_exp = p_exp (from dataset), y_model = curve p(T)
    # This lets you compute RMS/AAD (%) on pressure per author/year too.

    def add_pressure_block(prop_name, phase, T_arr, p_arr, meta_key):
        if T_arr is None or len(T_arr) == 0:
            return
        authors = meta.get(meta_key, {}).get(
            "Author", np.array(["Unknown"]*len(T_arr)))
        years = meta.get(meta_key, {}).get("Year",   np.array([0]*len(T_arr)))
        p_sub_curve = psub_curve(T_arr)*1e6  # sublimation curve in Pa
        p_melt_curve = pmelt_curve(T_arr)

        for Ti, Pi, auth, yr, psub_i, pmelt_i in zip(T_arr, p_arr, authors, years, p_sub_curve, p_melt_curve):
            if not (np.isfinite(Ti) and np.isfinite(Pi)):
                continue
            if phase == "sub":
                curve = psub_i
                dp_to_curve = Pi - psub_i
                denom = psub_i
            else:  # "melt"
                curve = pmelt_i
                dp_to_curve = Pi - pmelt_i
                denom = pmelt_i

            abs_dp_over_p = np.abs(dp_to_curve) / np.abs(denom) if (np.isfinite(
                dp_to_curve) and np.isfinite(denom) and denom != 0) else np.nan

            rows.append(dict(
                Property=prop_name,          # "psub" or "pmelt"
                Author=str(auth),
                Year=int(yr) if str(yr).isdigit() else yr,
                T=float(Ti),
                p_exp=float(Pi),
                p_sub_curve=float(psub_i),
                p_melt_curve=float(pmelt_i),
                dp_to_curve=float(dp_to_curve) if np.isfinite(
                    dp_to_curve) else np.nan,
                abs_dp_over_p=float(abs_dp_over_p) if np.isfinite(
                    abs_dp_over_p) else np.nan,
                y_exp=float(Pi),             # for RMS/AAD on pressure
                y_model=float(curve)         # curve value at T
            ))

    # Main properties (use your IDX map)
    add_block("Vm_sub",    "sub",  T_Vm_sub,   p_Vm_sub,
              Vm_sub,    IDX["Vm"],    "Vm_sub")
    add_block("Vm_melt",   "melt", T_Vm_melt,  p_Vm_melt,
              Vm_melt,   IDX["Vm"],    "Vm_melt")
    add_block("Vm_highp",  None,   T_Vm_highp, p_Vm_highp, Vm_highp,
              IDX["Vm"],    "Vm_highp" if meta.get("Vm_highp") is not None else "Vm_sub")
    add_block("cp_sub",    "sub",  T_cp_sub,   p_cp_sub,
              cp_sub,    IDX["cp"],    "cp_sub")
    add_block("alpha_sub", "sub",  T_alpha_sub, p_alpha_sub,
              alpha_sub,  IDX["Alpha"], "alpha_sub")
    add_block("BetaT_sub", "sub",  T_BetaT_sub, p_BetaT_sub,
              BetaT_sub,  IDX["KappaT"], "BetaT_sub")
    add_block("BetaS_sub", "sub",  T_BetaS_sub, p_BetaS_sub,
              BetaS_sub,  IDX["KappaS"], "BetaS_sub")
    add_pressure_block("psub",  "sub",  T_sub,  p_sub,  "sublimation")
    add_pressure_block("pmelt", "melt", T_melt, p_melt, "melting")

    # If you later want enthalpy sets, you can uncomment and define the exact y_exp_transform to match your thesis
    add_block("H_solid_sub", "sub",  T_H_sub,  p_H_sub,  (H_fluid_sub - delta_H_sub),  IDX["H"], "H_solid_sub")
    add_block("H_solid_melt","melt", T_H_melt, p_H_melt, (H_fluid_melt - delta_H_melt),IDX["H"], "H_solid_melt")

    return pd.DataFrame(rows)


def _relative_errors(y_exp, y_model, rel_epsilon=0.0):
    y_exp = np.asarray(y_exp, dtype=float)
    y_model = np.asarray(y_model, dtype=float)
    mask = np.isfinite(y_exp) & np.isfinite(y_model)
    y_exp = y_exp[mask]
    y_model = y_model[mask]
    if y_exp.size == 0:
        return np.array([])
    if rel_epsilon == 0.0:
        nz = np.abs(y_exp) > 0
        y_exp = y_exp[nz]
        y_model = y_model[nz]
        if y_exp.size == 0:
            return np.array([])
        return (y_model - y_exp) / np.abs(y_exp)
    else:
        return (y_model - y_exp) / (np.abs(y_exp) + rel_epsilon)


def rms_percent(y_exp, y_model, rel_epsilon=0.0):
    r = _relative_errors(y_exp, y_model, rel_epsilon)
    return np.nan if r.size == 0 else 100.0 * np.sqrt(np.mean(r**2))


def aad_percent(y_exp, y_model, rel_epsilon=0.0):
    r = _relative_errors(y_exp, y_model, rel_epsilon)
    return np.nan if r.size == 0 else 100.0 * np.mean(np.abs(r))


def summarise_by_author(master_df, rel_epsilon=0.0, float_fmt="{:.2f}"):
    rows = []
    for (prop, auth, yr), g in master_df.groupby(["Property", "Author", "Year"], dropna=False):
        if len(g) == 0:
            continue
        rms = rms_percent(g["y_exp"], g["y_model"], rel_epsilon)
        aad = aad_percent(g["y_exp"], g["y_model"], rel_epsilon)
        Tmin, Tmax = np.nanmin(g["T"]), np.nanmax(g["T"])
        rows.append({
            "Property": prop,
            "Author": auth,
            "Year": yr,
            "T_range": f"{Tmin:g}--{Tmax:g}" if Tmin != Tmax else f"{Tmin:g}",
            "N": int(len(g)),
            "RMS_%": float_fmt.format(rms) if np.isfinite(rms) else "--",
            "AAD_%": float_fmt.format(aad) if np.isfinite(aad) else "--",
        })
    return pd.DataFrame(rows).sort_values(["Property", "Year", "Author"]).reset_index(drop=True)
