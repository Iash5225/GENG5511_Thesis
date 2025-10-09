from thermopropsTV import compute_thermo_props_TV
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
from plot_eos import plot_all_overlays_grid
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)
# solid_EOS_fitting.py
BIG = 1e4
Vm_anchor = 0.8
T_MARGIN = 2.0  # K, stay below Tt by this margin
St_REFPROP = NEON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = NEON_REFERENCE_ENTHALPY
Tt = NEON_T_t
pt = NEON_P_t

def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), NEON_E_1_SUB,  NEON_E_2_SUB,  NEON_E_3_SUB,  NEON_T_t, NEON_P_t
)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), NEON_E_4, NEON_E_5, NEON_E_6, NEON_E_7, NEON_T_t, NEON_P_t
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
    """

    data = data.copy()
    data['Deviation_percent'] = 100 * \
        (data[y_exp_col] - data[y_model_col]) / data[y_exp_col]

    grouped = data.groupby(['Year', 'Author'])
    if custom_markers is None:
        custom_markers = ['o', 's', '^', 'v', 'D', 'x', '*', 'P', 'h', 'X']
    if custom_colors is None:
        custom_colors = plt.cm.tab10.colors

    fig, ax = plt.subplots(figsize=(10, 6))

    # ✅ Thicker border (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    # ✅ Ticks: inward, on all sides
    ax.tick_params(
        direction='in',
        which='both',
        top=True,
        right=True,
        length=5,
        width=1.2
    )

    # Plot data points
    for i, ((year, author), group) in enumerate(grouped):
        ax.scatter(
            group[x_col], group['Deviation_percent'],
            label=f"{year}, {author}",
            s=markersize,
            edgecolor=custom_colors[i % len(custom_colors)],
            marker=custom_markers[i % len(custom_markers)],
            facecolors='none',
            linewidths=1.2,
        )

    # Horizontal zero line
    ax.axhline(0, color='black', linewidth=1)

    # Labels and title
    ax.set_xlabel(r'$\mathit{T}$ / K', fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)

    # Axis limits
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    # Legend
    if legend_outside:
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            fontsize=8
        )
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        ax.legend(fontsize=8)
        plt.tight_layout()

    # ✅ Save if requested
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
    output_folder=None,
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
    """

    grouped = data.groupby(['Year', 'Author'])
    if custom_markers is None:
        custom_markers = ['o', 's', '^', 'v', 'D', 'x', '*', 'P', 'h', 'X']
    if custom_colors is None:
        custom_colors = plt.cm.tab10.colors

    plt.figure(figsize=(10, 6))
    ax = plt.gca()  # get current axis

    for i, ((year, author), group) in enumerate(grouped):
        ax.scatter(
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
        ax.plot(
            np.array(model_x)[order], np.array(model_y)[order],
            color='black', linewidth=linewidth, label='EOS prediction'
        )

    ax.set_xlabel(r'$\mathit{T}$ / K', fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)

    # Axis limits
    if logy:
        ax.set_yscale('log')
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    # ✅ Ticks: inward, on all sides
    ax.tick_params(
        direction='in',
        which='both',   # both major and minor ticks
        top=True,
        right=True,
        length=5,
        width=1.2
    )

    # ✅ Thicker black border
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color('black')

    # Legend
    if legend_outside:
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            fontsize=8,
            ncol=1
        )
        plt.tight_layout(rect=[0, 0, 0.85, 1])
    else:
        ax.legend(fontsize=8)
        plt.tight_layout()

    # ✅ Save figure if filename is given
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

def safe_pmelt(T):
    """
    Real, positive melting pressure; NaN outside domain.
    """
    p = np.asarray(pmelt_curve(T))
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


def _ensure_meta_columns(df: pd.DataFrame, property_name: str) -> pd.DataFrame:
    """
    Ensure df has 'Property', 'Author', 'Year' columns for grouping.
    Robust to case and stray whitespace; will also match columns that only
    contain the substring (e.g., 'measurement year').
    """
    if df is None:
        return df
    out = df.copy()

    def find_col(name: str, contains: bool = False):
        lname = name.lower()
        for c in out.columns:
            if not isinstance(c, str):
                continue
            cs = c.strip().lower()           # strip whitespace
            if (not contains and cs == lname) or (contains and lname in cs):
                return c
        return None

    # Property
    if "Property" not in out.columns:
        out["Property"] = property_name

    # Author (exact or case-insensitive/whitespace-robust)
    if "Author" not in out.columns:
        ac = find_col("Author") or find_col(
            "Authors") or find_col("author", contains=True)
        out["Author"] = out[ac] if ac is not None else "Unknown"

    # Year (exact; else substring; else leave NaN)
    if "Year" not in out.columns:
        yc = find_col("Year") or find_col("year", contains=True)
        if yc is None:
            print(
                f"[ensure-meta] '{property_name}': Year column not found; filling with NaN.")
            out["Year"] = np.nan
        else:
            # coerce to numeric years
            out["Year"] = pd.to_numeric(out[yc], errors="coerce")

    return out


def extract_datasets_with_meta(data):
    """
    Reuses extract_datasets() for numeric arrays, then builds meta dict
    (Author, Year) per property key. No duplication of numeric work.
    """
    # IMPORTANT: fix columns BEFORE we read them into meta
    def ensure(key, prop):
        if key in data and isinstance(data[key], pd.DataFrame):
            data[key] = _ensure_meta_columns(data[key], prop)
    ensure("heatsub", "H_solid_sub")     # <-- heatsub fixed here BEFORE meta
    ensure("fusion", "H_solid_melt")

    datasets = extract_datasets(data)

    def meta_block(df, n):
        if df is None:
            return None
        authors = df.get("Author", pd.Series(["Unknown"] * n)).to_numpy()
        years = df.get("Year", pd.Series([0] * n))
        # Ensure year is integer (handle non-numeric gracefully)
        years = pd.to_numeric(years, errors="coerce").fillna(
            0).astype(int).to_numpy()
        return {"Author": authors, "Year": years}

    meta = {
        "Vm_sub":       meta_block(data.get("cell_volume_sub"),   len(datasets[0])),
        "Vm_melt":      meta_block(data.get("cell_volume_melt"),  len(datasets[3])),
        "Vm_highp":     meta_block(data.get("cell_volume_highp"), len(datasets[6])) if len(datasets[6]) else None,
        "cp_sub":       meta_block(data.get("heat_capacity"),     len(datasets[9])),
        "alpha_sub":    meta_block(data.get("thermal_coeff"),     len(datasets[12])),
        "BetaT_sub":    meta_block(data.get("bulk_t"),            len(datasets[15])),
        "BetaS_sub":    meta_block(data.get("bulk_s"),            len(datasets[18])),
        "sublimation":  meta_block(data.get("sublimation"),       len(datasets[21])),
        "melting":      meta_block(data.get("melting"),           len(datasets[26])),
        "H_solid_sub":  meta_block(data.get("heatsub"),           len(datasets[31])),
        "H_solid_melt": meta_block(data.get("fusion"),            len(datasets[34])),
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
    T_Vm_sub = data["cell_volume_sub"]['Temperature']

    cv_sub_df = data["cell_volume_sub"]

    # don’t drop NaNs yet, keep them so we can inspect
    cv_sub_df["Temperature"] = pd.to_numeric(
        cv_sub_df["Temperature"], errors="coerce")
    cv_sub_df["Cell Volume"] = pd.to_numeric(
        cv_sub_df["Cell Volume"], errors="coerce")

    # Print any bad rows with author info
    bad = cv_sub_df[~np.isfinite(cv_sub_df["Temperature"])]
    if not bad.empty:
        print("⚠️ Non-finite temperatures in sublimation dataset:")
        print(bad[["Temperature", "Cell Volume", "Author", "Year"]])

    # Check before calling
    bad_T = T_Vm_sub[T_Vm_sub <= 0]
    if not bad_T.empty:
        print("⚠️ Found non-physical temperatures (<= 0 K):")
        print(bad_T.to_string(index=True))
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
    H_fluid_sub = data['heatsub']['Enthalpy']

    # Enthalpy Melting
    T_H_melt = data['fusion']['Temperature']
    p_H_melt = data['fusion']['Pressure']
    delta_H_melt = 1000.0 * data['fusion']['Change in Enthalpy']
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


def plot_init():
    neon_data = load_all_gas_data('neon', read_from_excel=False)
    datasets = extract_datasets(neon_data)
    params_init = PARAMS_INIT_NEON
    # master_df = build_master_pointwise_df(datasets, meta, params_init)
    # summarise_by_author(master_df).to_csv(os.path.join(IMG_OUTPUT_FOLDER, 'krypton_summary_by_author_init.csv'), index=False)
    plot_all_overlays_grid(params_init, datasets, Tt=Tt, pt=pt, compute_thermo_props=compute_thermo_props,
                           St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP, psub_curve=psub_curve, pmelt_curve=pmelt_curve)


def plot_deviation(variable='Vm_melt'):
    neon_data = load_all_gas_data('neon', read_from_excel=False)
    datasets, meta = extract_datasets_with_meta(neon_data)
    params_fit = PARAMS_INIT_NEON
    master_df = build_master_pointwise_df(datasets, meta, params_fit)

    # 2. Filter for the property you want to plot

    df_cell_volume_melt = master_df[master_df["Property"] == "Vm_melt"]
    df_cell_volume_sub = master_df[master_df["Property"] == "Vm_sub"]
    df_cp_sub = master_df[master_df["Property"] == "cp_sub"]
    df_alpha_sub = master_df[master_df["Property"] == "alpha_sub"]
    df_alpha_sub['y_exp'] = df_alpha_sub['y_exp'] * \
        10**4  # convert to 1e-4 K^-1
    df_alpha_sub['y_model'] = df_alpha_sub['y_model'] * \
        10**4  # convert to 1e-4 K^-1
    df_BetaT_sub = master_df[master_df["Property"] == "BetaT_sub"]
    df_BETA_T_sub = df_BetaT_sub.copy()
    df_BETA_T_sub['y_exp'] = 1 / df_BETA_T_sub['y_exp']  # Now in MPa
    df_BETA_T_sub['y_model'] = 1 / df_BETA_T_sub['y_model']  # Now in MPa
    df_BetaS_sub = master_df[master_df["Property"] == "BetaS_sub"]
    df_BETA_S_sub = df_BetaS_sub.copy()
    df_BETA_S_sub['y_exp'] = 1 / df_BETA_S_sub['y_exp']  # Now in MPa
    df_BETA_S_sub['y_model'] = 1 / df_BETA_S_sub['y_model']  # Now in MPa
    df_enthalpy_solid_melt = master_df[master_df["Property"] == "H_solid_melt"]
    df_enthalpy_solid_sub = master_df[master_df["Property"] == "H_solid_sub"]
    df_pressure_sub = master_df[master_df["Property"] == "psub"]
    df_pressure_melt = master_df[master_df["Property"] == "pmelt"]
    # df_BetaT_sub['KappaT_exp'] = 1 / df_BetaT_sub['BetaT_exp']  # Now in MPa
    if variable == 'Vm_melt':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_cell_volume_melt,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3mol^{-1}}$',
            title=None,
            model_x=df_cell_volume_melt['T'],
            model_y=df_cell_volume_melt['y_model'],
            logy=False,
            xlim=(25, 55),
            ylim=(12.25, 14.25),
            filename='neon_melt_cellvolume.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cell_volume_melt,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}}) / V_{\mathrm{m,exp}}$',
            title=None,
            filename='neon_melt_cellvolume_deviation',
            xlim=(25, 55),
            ylim=(-0.5, 1.5),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'Vm_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_cell_volume_sub,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3mol^{-1}}$',
            title=None,
            xlim=(0, 25),
            ylim=(13.2, 14.2),
            model_x=df_cell_volume_sub['T'],
            model_y=df_cell_volume_sub['y_model'],
            logy=False,
            filename='neon_sub_cellvolume.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cell_volume_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}}) / V_{\mathrm{m,exp}}$',
            title=None,
            filename='neon_sub_cellvolume_deviation',
            xlim=(0, 25),
            ylim=(-0.5, 0.1),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'cp_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_cp_sub,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$c_p\,/\,\mathrm{JK^{-1}mol^{-1}}$',
            title=None,
            xlim=(10, 25),
            ylim=(7.5, 27.5),
            model_x=df_cp_sub['T'],
            model_y=df_cp_sub['y_model'],
            logy=False,
            filename='neon_sub_heatcapacity.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cp_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (c_{p,\mathrm{exp}} - c_{p,\mathrm{calc}}) / c_{p,\mathrm{exp}}$',
            title=None,
            filename='neon_sub_heatcapacity_deviation',
            xlim=(10, 25),
            ylim=(-7.5, 12.5),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'alpha_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_alpha_sub,  # convert to 1e-4 K^-1
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\alpha \cdot 10^{4}/\,\mathrm{K^{-1}}$',
            title=None,
            xlim=(0, 25),
            ylim=(0, 55),
            model_x=df_alpha_sub['T'],
            model_y=df_alpha_sub['y_model'],  # convert to 1e-4 K^-1
            logy=False,
            filename='neon_sub_thermal_expansion.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_alpha_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (\alpha_{\mathrm{exp}} - \alpha_{\mathrm{calc}}) / \alpha_{\mathrm{exp}}$',
            title=None,
            filename='neon_sub_thermal_expansion_deviation',
            xlim=(0, 25),
            ylim=(-10, 10),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'BetaT_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_BETA_T_sub,  # convert to KappaT in MPa
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$K_T\,/\,\mathrm{MPa}$',
            title=None,
            xlim=(0, 25),
            ylim=(500, 1200),
            model_x=df_BETA_T_sub['T'],
            model_y=df_BETA_T_sub['y_model'],  # convert to KappaT in MPa
            logy=False,
            filename='neon_sub_isothermal_compressibility.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_BETA_T_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (K_{T,\mathrm{exp}} - K_{T,\mathrm{calc}}) / K_{T,\mathrm{exp}}$',
            title=None,
            filename='neon_sub_isothermal_compressibility_deviation',
            xlim=(0, 25),
            ylim=(-10, 15),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'BetaS_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_BETA_S_sub,  # convert to KappaS in MPa
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$K_S\,/\,\mathrm{MPa}$',
            title=None,
            xlim=(0, 25),
            ylim=(600, 1200),
            model_x=df_BETA_S_sub['T'],
            model_y=df_BETA_S_sub['y_model'],  # convert to KappaS in MPa
            logy=False,
            filename='neon_sub_isentropic_compressibility.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_BETA_S_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (K_{S,\mathrm{exp}} - K_{S,\mathrm{calc}}) / K_{S,\mathrm{exp}}$',
            title=None,
            filename='neon_sub_isentropic_compressibility_deviation',
            xlim=(0, 25),
            ylim=(-6, 12),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'H_solid_sub':
        # shift model enthalpy by the triple-point offset (kJ/mol)
        dHtr_kJ, _ = _triple_offsets_plot(params_fit, Tt=Tt, pt=pt,
                                          St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP)

        dfH = df_enthalpy_solid_sub.copy()
        # y_exp is H_solid_exp in J/mol from builder → convert to kJ/mol
        dfH['y_exp'] = dfH['y_exp'] / 1000.0
        # y_model is H_model (J/mol) → convert and shift to solid reference
        dfH['y_model'] = dfH['y_model'] / 1000.0 - dHtr_kJ

        # Top panel: ΔH (kJ/mol) vs T with model line
        plot_thermo_variable(
            data=dfH,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\Delta H\,/\,\mathrm{kJ\,mol^{-1}}$',
            title=None,
            model_x=dfH['T'],
            model_y=dfH['y_model'],
            # xlim=(50, 170),
            ylim=(-0.525, -0.510),
            logy=False,
            filename='neon_sub_enthalpy.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        # Bottom: percent deviation with model in denominator (as in your MATLAB figure)
        plot_variable_deviation(
            data=dfH,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(50, 170),
            ylim=(-0.1, 0.8),
            y_label=r'$100\cdot(\Delta H_{\mathrm{exp}}-\Delta H_{\mathrm{calc}})/\Delta H_{\mathrm{exp}}$',
            title=None,
            filename='neon_sub_enthalpy_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    elif variable == 'H_solid_melt':
        # shift model enthalpy by the triple-point offset (kJ/mol)
        dHtr_kJ, _ = _triple_offsets_plot(params_fit, Tt=Tt, pt=pt,
                                          St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP)

        dfH = df_enthalpy_solid_melt.copy()
        # y_exp is H_solid_exp in J/mol from builder → convert to kJ/mol
        dfH['y_exp'] = dfH['y_exp'] / 1000.0
        # y_model is H_model (J/mol) → convert and shift to solid reference
        dfH['y_model'] = dfH['y_model'] / 1000.0 - dHtr_kJ

        # Top panel: ΔH (kJ/mol) vs T with model line
        plot_thermo_variable(
            data=dfH,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\Delta H\,/\,\mathrm{kJ\,mol^{-1}}$',
            title=None,
            model_x=dfH['T'],
            model_y=dfH['y_model'],
            # xlim=(110, 220),
            ylim=(12.5, 32.5),
            logy=False,
            filename='neon_melt_enthalpy.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        # Bottom: percent deviation with model in denominator (as in your MATLAB figure)
        plot_variable_deviation(
            data=dfH,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(110, 220),
            ylim=(-0.5, 3.5),
            y_label=r'$100\cdot(\Delta H_{\mathrm{exp}}-\Delta H_{\mathrm{calc}})/\Delta H_{\mathrm{exp}}$',
            title=None,
            filename='neon_melt_enthalpy_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    elif variable == 'psub':
        plot_thermo_variable(
            data=df_pressure_sub,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$p/\,\mathrm{MPa}$',
            title=None,
            model_x=df_pressure_sub['T'],
            model_y=df_pressure_sub['y_model'],
            # xlim=(60, 120),
            # ylim=(10**1, 10**5),
            logy=True,
            filename='neon_sub_pressure.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_pressure_sub,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(60, 120),
            ylim=(-40, 100),
            y_label=r'$100\cdot(p_{\mathrm{exp}}-p_{\mathrm{calc}})/p_{\mathrm{exp}}$',
            title=None,
            filename='neon_sub_pressure_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    elif variable == 'pmelt':
        plot_thermo_variable(
            data=df_pressure_melt,
            gas_name='neon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$p/\,\mathrm{MPa}$',
            title=None,
            model_x=df_pressure_melt['T'],
            model_y=df_pressure_melt['y_model'],
            # xlim=(60, 120),
            # ylim=(10**1, 10**5),
            logy=True,
            filename='neon_melt_pressure.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_pressure_melt,
            gas_name='neon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(60, 120),
            ylim=(-15, 5),
            y_label=r'$100\cdot(p_{\mathrm{exp}}-p_{\mathrm{calc}})/p_{\mathrm{exp}}$',
            title=None,
            filename='neon_melt_pressure_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )


def _triple_offsets_plot(params, Tt=Tt, pt=pt, St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP):
    """
    Same logic as fitting_xenon._triple_offsets but local here to avoid circular imports.
    Returns (dH_triple_kJ, dS_triple).
    """
    tp = compute_thermo_props(float(Tt), float(pt), params)
    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        meta = tp[1]
        dH = float(meta.get("deltaH_triple"))
        dS = float(meta.get("deltaS_triple"))
    else:
        props = np.asarray(tp, float)
        Sstar = float(params[30])
        dS = Sstar - float(St_REFPROP)
        # use G + T*S* to align enthalpy at triple
        Ht_fitted = float(props[IDX["G"]]) + float(Tt) * Sstar
        dH = Ht_fitted - float(Ht_REFPROP)  # J/mol
    return dH/1000.0, dS  # kJ/mol, dimensionless


#

def plot_tv_property_panels(
    params,
    V_vals,
    T_list=(0.0001, 50.0, 83.806, 150.0, 300.0),
    Tt=None,
    psub_curve=None,   # function T -> MPa
    pmelt_curve=None,  # function T -> MPa
    figsize=(10, 8),
    save_path=None
):
    """
    Plot 4 panels vs Vm using the TV EOS:
      a) p (MPa)
      b) K_S (MPa)  [isentropic bulk modulus = 1/KappaS]
      c) alpha * 1e4 (1/K)
      d) beta = alpha / KappaT (MPa^-1)

    Overlays:
      - Isotherms at T in T_list with styles:
          0 K: black solid
         50 K: purple dashed
      83.806 K: dark red dash-dot
        150 K: green dotted
        300 K: bold orange solid
      - Phase curves:
        Sublimation: bold blue dashed
        Melting:     bold grey dash-dot

    Notes
    - V_vals is a 1D array of molar volumes (cm^3/mol) to sweep, e.g., np.arange(18, 25.01, 0.01).
    - psub_curve and pmelt_curve should return pressures in MPa.
    """
    V_vals = np.asarray(V_vals, float)

    # Styling maps
    styles = {
        0.0001: dict(color='k',      ls='-',  lw=1.8, label='0 K'),
        # purple
        50.0:   dict(color='#8e44ad', ls='--', lw=1.8, label='50 K'),
        # dark red
        83.806: dict(color='#8B0000', ls='-.', lw=1.8, label='83.806 K'),
        # teal/green
        150.0:  dict(color='#16a085', ls=':',  lw=1.8, label='150 K'),
        # orange, bold
        300.0:  dict(color='#e67e22', ls='-',  lw=2.6, label='300 K'),
    }
    sub_style = dict(color='#1f77b4', ls='--', lw=2.6,
                     label='sublimation')   # bold blue dashed
    melt_style = dict(color='#7f7f7f', ls='-.', lw=2.6,
                      label='melting')       # bold grey dash-dot

    # Helpers to compute properties at fixed T across V_vals
    def tv_props_along_V(T):
        p_list, Ks_list, a_list, beta_list, V_list = [], [], [], [], []
        for v in V_vals:
            try:
                props = compute_thermo_props_TV(float(T), float(v), params)
                p = float(props[0])                    # MPa
                kT = float(props[1])                    # MPa^-1
                kS = float(props[2])                    # MPa^-1
                alp = float(props[3])                    # 1/K
                if np.isfinite(p) and np.isfinite(kT) and np.isfinite(kS) and np.isfinite(alp):
                    V_list.append(v)
                    p_list.append(p)
                    # bulk moduli:
                    Ks_list.append(1.0 / kS if kS !=
                                   0 else np.nan)           # MPa
                    beta_list.append(alp / kT if kT !=
                                     0 else np.nan)         # MPa^-1
                    # plot α·1e4
                    a_list.append(alp * 1e4)
            except Exception:
                pass
        return np.asarray(V_list), np.asarray(p_list), np.asarray(Ks_list), np.asarray(a_list), np.asarray(beta_list)

    # Helpers to sample a phase curve parametrically by T, map to V, then compute properties at (T,V)
    def sample_phase_curve(T_grid, p_of_T):
        Vc, pc, KSc, ac, bc = [], [], [], [], []
        for T in T_grid:
            try:
                p = float(p_of_T(T))  # MPa
                # Get V from p-based API at (T,p)
                props_p = compute_thermo_props(float(T), float(p), params)
                V = float(props_p[0])  # Vm index=0
                # Now compute properties from TV API at (T,V)
                props_tv = compute_thermo_props_TV(float(T), float(V), params)
                pMPa = float(props_tv[0])
                kT = float(props_tv[1])
                kS = float(props_tv[2])
                alp = float(props_tv[3])
                if np.all(np.isfinite([V, pMPa, kT, kS, alp])):
                    Vc.append(V)
                    pc.append(pMPa)
                    KSc.append(1.0 / kS if kS != 0 else np.nan)
                    ac.append(alp * 1e4)
                    bc.append(alp / kT if kT != 0 else np.nan)
            except Exception:
                pass
        return (np.asarray(Vc), np.asarray(pc), np.asarray(KSc),
                np.asarray(ac), np.asarray(bc))

    fig, axs = plt.subplots(2, 2, figsize=figsize)
    ax_p, ax_Ks = axs[0]
    ax_alpha, ax_beta = axs[1]

    # Plot isotherms
    for T in T_list:
        Vc, pc, KSc, ac, bc = tv_props_along_V(T)
        if Vc.size == 0:
            continue
        st = styles.get(round(T, 3), dict(
            color='k', ls='-', lw=1.5, label=f'{T:g} K'))
        # sort by V for clean lines
        order = np.argsort(Vc)
        Vc, pc, KSc, ac, bc = Vc[order], pc[order], KSc[order], ac[order], bc[order]
        ax_p.plot(Vc, pc, **st)
        ax_Ks.plot(Vc, KSc, **st)
        ax_alpha.plot(Vc, ac, **st)
        ax_beta.plot(Vc, bc, **st)

    # Plot phase curves (if provided)
    if (psub_curve is not None) and (Tt is not None):
        # Sample safely below Tt
        T_sub = np.linspace(max(1e-3, 0.1), float(Tt), 400)
        Vc, pc, KSc, ac, bc = sample_phase_curve(T_sub, psub_curve)
        if Vc.size:
            o = np.argsort(Vc)
            ax_p.plot(Vc[o], pc[o], **sub_style)
            ax_Ks.plot(Vc[o], KSc[o], **sub_style)
            ax_alpha.plot(Vc[o], ac[o], **sub_style)
            ax_beta.plot(Vc[o], bc[o], **sub_style)

    if pmelt_curve is not None and (Tt is not None):
        # Sample above Tt (limit to a reasonable max T, e.g., last isotherm or 300 K)
        Tmax = max([t for t in T_list if np.isfinite(t)]
                   ) if T_list else (float(Tt) + 100.0)
        Tmax = max(Tmax, float(Tt) + 1.0)
        T_melt = np.linspace(float(Tt), float(Tmax), 400)
        Vc, pc, KSc, ac, bc = sample_phase_curve(T_melt, pmelt_curve)
        if Vc.size:
            o = np.argsort(Vc)
            ax_p.plot(Vc[o], pc[o], **melt_style)
            ax_Ks.plot(Vc[o], KSc[o], **melt_style)
            ax_alpha.plot(Vc[o], ac[o], **melt_style)
            ax_beta.plot(Vc[o], bc[o], **melt_style)

    # Labels, limits, legends
    ax_p.set_xlabel(r"$V_m$ / (cm$^3$ mol$^{-1}$)")
    ax_p.set_ylabel(r"$p$ / MPa")

    ax_Ks.set_xlabel(r"$V_m$ / (cm$^3$ mol$^{-1}$)")
    ax_Ks.set_ylabel(r"$K_S$ / MPa")

    ax_alpha.set_xlabel(r"$V_m$ / (cm$^3$ mol$^{-1}$)")
    ax_alpha.set_ylabel(r"$\alpha \cdot 10^{4}$ / K$^{-1}$")

    ax_beta.set_xlabel(r"$V_m$ / (cm$^3$ mol$^{-1}$)")
    ax_beta.set_ylabel(r"$\beta$ / MPa$^{-1}$")

    # Aesthetic tweaks
    for ax in (ax_p, ax_Ks, ax_alpha, ax_beta):
        ax.tick_params(direction='in', top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)
            spine.set_color('black')

    # One shared legend using handles from the first axis
    # Collect unique labels in order
    handles, labels = [], []
    for ax in (ax_p,):
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in labels:
                handles.append(hi)
                labels.append(li)
    fig.legend(handles, labels, loc='upper center',
               ncol=4, frameon=False, fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.94])

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()
    return fig, axs

