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
from thermopropsTV import compute_thermo_props_TV
from scipy.interpolate import UnivariateSpline
import os
from plot_eos import plot_all_overlays_grid
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)
# solid_EOS_fitting.py
BIG = 1e4
Vm_anchor = 0.8
T_MARGIN = 2.0  # K, stay below Tt by this margin
St_REFPROP = XENON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = XENON_REFERENCE_ENTHALPY
Tt = XENON_T_t
pt = XENON_P_t

def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), XENON_E_1_SUB,  XENON_E_2_SUB,  XENON_E_3_SUB,  XENON_T_t, XENON_P_t
)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), XENON_E_4, XENON_E_5, XENON_E_6, XENON_E_7, XENON_T_t, XENON_P_t
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


def _smooth_highp_curves(df_highp, params_fit, idx_vm, npts=250):
    """
    Return list of smooth Vm(P) curves per (Year, Author) group.
    Interpolates T(P) within each group then evaluates EOS on dense grid.
    """
    curves = []
    need_cols = {'Year', 'Author', 'p_exp', 'T'}
    if not need_cols.issubset(df_highp.columns):
        return curves
    for (yr, auth), g in df_highp.groupby(['Year', 'Author']):
        g = g[['p_exp', 'T']].dropna().sort_values('p_exp')
        if g.shape[0] < 2:
            continue
        p_raw = g['p_exp'].values
        T_raw = g['T'].values
        # remove duplicated pressures
        p_unique, idx = np.unique(p_raw, return_index=True)
        T_unique = T_raw[idx]
        p_grid = np.linspace(p_unique.min(), p_unique.max(), npts)
        T_grid = np.interp(p_grid, p_unique, T_unique)
        Vm_grid = np.full_like(p_grid, np.nan, dtype=float)
        for i, (Pi, Ti) in enumerate(zip(p_grid, T_grid)):
            try:
                props = compute_thermo_props(float(Ti), float(Pi), params_fit)
                if np.all(np.isfinite(props)):
                    Vm_grid[i] = props[idx_vm]
            except Exception:
                pass
        curves.append(dict(year=yr, author=auth, p=p_grid, Vm=Vm_grid))
    return curves


def _eval_vm_from_eos(T, p_array, params, vm_idx):
    """
    Evaluate Vm(T, p) on a grid. Works whether compute_thermo_props returns
    a dict or a sequence (tuple/list/ndarray).
    """
    vals = []
    for p in np.asarray(p_array, float):
        res = compute_thermo_props(T, p, params)
        if isinstance(res, dict):
            vm = res.get("Vm", np.nan)
        else:
            vm = res[vm_idx]
        vals.append(vm)
    return np.asarray(vals, dtype=float)

def _smooth_isotherm_curves(
    df, params, vm_idx, npts=400, T_bin=0.5, p_pad=0.0
):
    """
    Build smooth EOS curves per (Year, Author, isotherm-bin).
    Returns list of dicts with keys: year, author, T, p, Vm.
    """
    d = df.copy()
    d["Tbin"] = (d["T"] / T_bin).round() * T_bin

    curves = []
    for (year, author, Tbin), g in d.groupby(["Year", "Author", "Tbin"]):
        T0 = float(g["T"].mean())
        pmin, pmax = np.nanmin(g["p_exp"]), np.nanmax(g["p_exp"])
        if not (np.isfinite(pmin) and np.isfinite(pmax)) or pmax <= pmin:
            continue

        # Slight padding if desired
        span = pmax - pmin
        pgrid = np.linspace(max(1e-6, pmin - p_pad * span),
                            pmax + p_pad * span, npts)

        # >>> replace with your actual evaluator <<<
        Vm = _eval_vm_from_eos(T0, pgrid, params, vm_idx)


        Vm = np.asarray(Vm, dtype=float)

        m = np.isfinite(pgrid) & np.isfinite(Vm)
        if m.sum() < 5:
            continue

        curves.append({
            "year": year, "author": author, "T": T0,
            "p": pgrid[m], "Vm": Vm[m]
        })
    return curves

def _curve_pressure_from_T(T_array, Tt, psub_curve, pmelt_curve, mode='auto'):
    """
    Return pressure (MPa) from temperature using chosen curve:
      mode='sub'  -> always sublimation curve
      mode='melt' -> always melting curve
      mode='auto' -> T <= Tt -> sub, else melt
    """
    T_array = np.asarray(T_array, float)
    if mode == 'sub':
        return psub_curve(T_array)
    if mode == 'melt':
        return pmelt_curve(T_array)
    # auto
    p_sub = psub_curve(T_array)
    p_melt = pmelt_curve(T_array)
    return np.where(T_array <= Tt, p_sub, p_melt)

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
    ishighp=False, params_fit=None, vm_idx=None, npts_iso=500, T_bin=0.5, eos_linestyle="--"
):
    """
    General-purpose plot for thermodynamic variables (e.g., enthalpy, heat capacity, etc.)
    grouped by author/year, with optional model curve.
    """

    if ishighp:
        data["p_calc"] = _curve_pressure_from_T(
            data["T"], Tt, psub_curve, pmelt_curve, mode='auto')
        # x_col_local = "p_calc"
        x_label = r"$p$ / MPa"
    else:
        x_col_local = x_col
        # Temperature axis label heuristic
        x_label = (
            r'$\mathit{T}$ / K' if x_col_local.lower().startswith('t') else x_col_local)

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
    # if model_x is not None and model_y is not None:
    #     order = np.argsort(model_x)
    #     ax.plot(
    #         np.array(model_x)[order], np.array(model_y)[order],
    #         color='black', linewidth=linewidth, label='EOS prediction'
    #     )
    # --- draw EOS ---
    if not ishighp:
        # original single-curve behavior
        if model_x is not None and model_y is not None:
            order = np.argsort(model_x)
            ax.plot(np.array(model_x)[order], np.array(model_y)[order],
                    color='black', linewidth=linewidth, label='EOS prediction')
    else:
        # high-p: per (Year, Author) isotherms if params provided
        if params_fit is not None and vm_idx is not None:
            # keep color consistent with scatter: same enumeration order
            for i, ((year, author), group) in enumerate(grouped):
                col = custom_colors[i % len(custom_colors)]
                # build smooth isotherms for this dataset group
                curves = _smooth_isotherm_curves_by_group(
                    group, params_fit, vm_idx, npts=npts_iso, T_bin=T_bin
                )
                for c in curves:
                    ax.plot(
                        c["p"], c["Vm"],
                        linestyle=eos_linestyle, linewidth=linewidth,
                        color=col,  # color-match the dataset
                        # keep legend clean: no separate label for each isotherm
                    )


    ax.set_xlabel(r'$\mathit{T}$ / K', fontsize=fontsize)
    if ishighp:
        ax.set_xlabel(r"$p$ / MPa", fontsize=fontsize)
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
    # G_fluid_melt = data['melting']['Gibbs Energy']
    # V_fluid_melt = data['melting']['Volume']
    G_fluid_melt = None
    V_fluid_melt = None

    # Sublimation
    T_sub = data['sublimation']['Temperature']
    p_sub = data['sublimation']['Pressure']
    # G_fluid_sub = data['sublimation']['Gibbs Energy']
    # V_fluid_sub = data['sublimation']['Volume']
    G_fluid_sub = None
    V_fluid_sub = None
    # Enthalpy Sublimation
    T_H_sub = data['heatsub']['Temperature']
    p_H_sub = safe_psub(T_H_sub)
    # p_H_sub = data['heatsub']['Pressure']
    delta_H_sub = 1000.0 * data['heatsub']['Change in Enthalpy']  # kJ → J
    H_fluid_sub = data['heatsub']['Enthalpy']

    # Enthalpy Melting
    T_H_melt = data['fusion']['Temperature']
    p_H_melt = pmelt_curve(T_H_melt)
    # p_H_melt = data['fusion']['Pressure']
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


def _cluster_highp_isotherms(df, T_tol=0.2):
    """
    Cluster high‑pressure experimental points into quasi‑isotherms per (Year, Author).

    T_tol : float (K)
        Maximum allowed deviation from the running cluster mean temperature to keep
        adding points to the current cluster.
    Returns list of dicts:
        {year, author, T_mean, p_values (np.array), T_values (np.array)}
    """
    clusters = []
    if df is None or df.empty:
        return clusters
    required = {'Year', 'Author', 'T', 'p_exp'}
    missing = required - set(df.columns)
    if missing:
        print(f"_cluster_highp_isotherms: missing columns {missing}")
        return clusters

    # Process each (Year, Author) separately
    for (yr, auth), g in df.groupby(['Year', 'Author']):
        g = g[['T', 'p_exp']].dropna().sort_values('T')
        if g.empty:
            continue

        current_T = []
        current_p = []

        def flush_cluster():
            if not current_T:
                return
            T_arr = np.array(current_T, float)
            p_arr = np.array(current_p, float)
            clusters.append(dict(
                year=yr,
                author=auth,
                T_mean=float(T_arr.mean()),
                p_values=p_arr,
                T_values=T_arr
            ))

        for _, row in g.iterrows():
            Ti = float(row['T'])
            pi = float(row['p_exp'])
            if not current_T:
                current_T.append(Ti)
                current_p.append(pi)
                continue
            mean_T = np.mean(current_T)
            if abs(Ti - mean_T) <= T_tol:
                current_T.append(Ti)
                current_p.append(pi)
            else:
                flush_cluster()
                current_T = [Ti]
                current_p = [pi]

        flush_cluster()

    return clusters


def _eval_vm_from_eos(T, p_array, params, vm_idx):
    """Vm(T, p) on a pressure grid; works for dict or sequence returns."""
    out = []
    for p in np.asarray(p_array, float):
        res = compute_thermo_props(T, p, params)
        vm = res.get("Vm", np.nan) if isinstance(res, dict) else res[vm_idx]
        out.append(vm)
    return np.asarray(out, float)


def _smooth_isotherm_curves_by_group(df_group, params, vm_idx, npts=400, T_bin=0.5):
    """
    Build EOS curves for one (Year, Author) group.
    Returns a list of dicts: {'T': T0, 'p': pgrid, 'Vm': Vm}.
    """
    g = df_group.copy()
    g["Tbin"] = (g["T"]/T_bin).round()*T_bin

    curves = []
    for Tbin, gg in g.groupby("Tbin"):
        T0 = float(gg["T"].mean())
        pmin, pmax = np.nanmin(gg["p_exp"]), np.nanmax(gg["p_exp"])
        if not (np.isfinite(pmin) and np.isfinite(pmax)) or pmax <= pmin:
            continue
        pgrid = np.linspace(max(1e-6, pmin), pmax, npts)
        Vm = _eval_vm_from_eos(T0, pgrid, params, vm_idx)
        m = np.isfinite(pgrid) & np.isfinite(Vm)
        if m.sum() >= 5:
            curves.append({"T": T0, "p": pgrid[m], "Vm": Vm[m]})
    return curves


def _eos_isotherm_for_cluster(cluster, params, idx_vm, npts=200, p_pad_frac=0.03):
    """
    Build a smooth EOS Vm(p) line at constant T = cluster['T_mean'] over
    pressure range of the cluster (with padding).
    Returns dict with p_line, Vm_line, T_mean, year, author.
    """
    p_vals = np.asarray(cluster['p_values'], float)
    if p_vals.size < 2:
        return None
    pmin, pmax = np.nanmin(p_vals), np.nanmax(p_vals)
    if not (np.isfinite(pmin) and np.isfinite(pmax)) or pmax <= pmin:
        return None
    span = pmax - pmin
    pad = p_pad_frac * span
    p_line = np.linspace(max(1e-6, pmin - pad), pmax + pad, npts)
    Vm_line = np.full_like(p_line, np.nan, dtype=float)
    T_iso = float(cluster['T_mean'])
    for i, P in enumerate(p_line):
        try:
            props = compute_thermo_props(T_iso, float(P), params)
            if np.all(np.isfinite(props)):
                Vm_line[i] = props[idx_vm]
        except Exception:
            pass
    m = np.isfinite(Vm_line)
    if m.sum() < 5:
        return None
    return dict(
        year=cluster['year'], author=cluster['author'],
        T_mean=T_iso, p_line=p_line[m], Vm_line=Vm_line[m]
    )
def plot_deviation(variable='Vm_melt'):
    xenon_data = load_all_gas_data('xenon', read_from_excel=False)
    datasets, meta = extract_datasets_with_meta(xenon_data)
    params_fit = PARAMS_INIT_XENON
    master_df = build_master_pointwise_df(datasets, meta, params_fit)

    # 2. Filter for the property you want to plot

    df_cell_volume_melt = master_df[master_df["Property"] == "Vm_melt"]
    df_cell_volume_sub = master_df[master_df["Property"] == "Vm_sub"]
    df_cell_volume_highp = master_df[master_df["Property"] == "Vm_highp"]
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
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3mol^{-1}}$',
            title=None,
            model_x=df_cell_volume_melt['T'],
            model_y=df_cell_volume_melt['y_model'],
            logy=False,
            xlim=(160, 375),
            ylim=(32, 40),
            filename='xenon_melt_cellvolume.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cell_volume_melt,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}}) / V_{\mathrm{m,exp}}$',
            title=None,
            filename='xenon_melt_cellvolume_deviation',
            xlim=(160, 375),
            ylim=(-2, 1.5),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'Vm_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_cell_volume_sub,
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3mol^{-1}}$',
            title=None,
            xlim=(0, 170),
            ylim=(34, 40),
            model_x=df_cell_volume_sub['T'],
            model_y=df_cell_volume_sub['y_model'],
            logy=False,
            filename='xenon_sub_cellvolume.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cell_volume_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}}) / V_{\mathrm{m,exp}}$',
            title=None,
            filename='xenon_sub_cellvolume_deviation',
            xlim=(0, 170),
            ylim=(-0.5, 1.5),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'cp_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_cp_sub,
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$c_p\,/\,\mathrm{JK^{-1}mol^{-1}}$',
            title=None,
            xlim=(0, 170),
            ylim=(0, 40),
            model_x=df_cp_sub['T'],
            model_y=df_cp_sub['y_model'],
            logy=False,
            filename='xenon_sub_heatcapacity.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_cp_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (c_{p,\mathrm{exp}} - c_{p,\mathrm{calc}}) / c_{p,\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_heatcapacity_deviation',
            xlim=(0, 170),
            ylim=(-10, 20),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'alpha_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_alpha_sub,  # convert to 1e-4 K^-1
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\alpha \cdot 10^{4}/\,\mathrm{K^{-1}}$',
            title=None,
            xlim=(0, 170),
            ylim=(0, 12),
            model_x=df_alpha_sub['T'],
            model_y=df_alpha_sub['y_model'],  # convert to 1e-4 K^-1
            logy=False,
            filename='xenon_sub_thermal_expansion.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_alpha_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (\alpha_{\mathrm{exp}} - \alpha_{\mathrm{calc}}) / \alpha_{\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_thermal_expansion_deviation',
            xlim=(0, 170),
            ylim=(-50, 15),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'BetaT_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_BETA_T_sub,  # convert to KappaT in MPa
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$K_T\,/\,\mathrm{MPa}$',
            title=None,
            xlim=(0, 170),
            ylim=(1250, 3750),
            model_x=df_BETA_T_sub['T'],
            model_y=df_BETA_T_sub['y_model'],  # convert to KappaT in MPa
            logy=False,
            filename='xenon_sub_isothermal_compressibility.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_BETA_T_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (K_{T,\mathrm{exp}} - K_{T,\mathrm{calc}}) / K_{T,\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_isothermal_compressibility_deviation',
            xlim=(0, 170),
            ylim=(-15, 15),
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
    elif variable == 'BetaS_sub':
        # # 3. Plot the variable
        plot_thermo_variable(
            data=df_BETA_S_sub,  # convert to KappaS in MPa
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$K_S\,/\,\mathrm{MPa}$',
            title=None,
            xlim=(0, 170),
            ylim=(1500, 4000),
            model_x=df_BETA_S_sub['T'],
            model_y=df_BETA_S_sub['y_model'],  # convert to KappaS in MPa
            logy=False,
            filename='xenon_sub_isentropic_compressibility.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_BETA_S_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (K_{S,\mathrm{exp}} - K_{S,\mathrm{calc}}) / K_{S,\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_isentropic_compressibility_deviation',
            xlim=(0, 170),
            ylim=(-40, 10),
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
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\Delta H\,/\,\mathrm{kJ\,mol^{-1}}$',
            title=None,
            model_x=dfH['T'],
            model_y=dfH['y_model'],
            xlim=(50, 170),
            ylim=(-6, -2.5),
            logy=False,
            filename='xenon_sub_enthalpy.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        # Bottom: percent deviation with model in denominator (as in your MATLAB figure)
        plot_variable_deviation(
            data=dfH,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            xlim=(50, 170),
            ylim=(-8, 1),
            y_label=r'$100\cdot(\Delta H_{\mathrm{exp}}-\Delta H_{\mathrm{calc}})/\Delta H_{\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_enthalpy_deviation',
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
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$\Delta H\,/\,\mathrm{kJ\,mol^{-1}}$',
            title=None,
            model_x=dfH['T'],
            model_y=dfH['y_model'],
            # xlim=(110, 220),
            ylim=(12.5, 32.5),
            logy=False,
            filename='xenon_melt_enthalpy.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        # Bottom: percent deviation with model in denominator (as in your MATLAB figure)
        plot_variable_deviation(
            data=dfH,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(110, 220),
            ylim=(-0.5, 3.5),
            y_label=r'$100\cdot(\Delta H_{\mathrm{exp}}-\Delta H_{\mathrm{calc}})/\Delta H_{\mathrm{exp}}$',
            title=None,
            filename='xenon_melt_enthalpy_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    elif variable == 'psub':
        plot_thermo_variable(
            data=df_pressure_sub,
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$p/\,\mathrm{MPa}$',
            title=None,
            model_x=df_pressure_sub['T'],
            model_y=df_pressure_sub['y_model'],
            # xlim=(60, 120),
            # ylim=(10**1, 10**5),
            logy=True,
            filename='xenon_sub_pressure.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_pressure_sub,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(60, 120),
            ylim=(-40, 100),
            y_label=r'$100\cdot(p_{\mathrm{exp}}-p_{\mathrm{calc}})/p_{\mathrm{exp}}$',
            title=None,
            filename='xenon_sub_pressure_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    elif variable == 'pmelt':
        plot_thermo_variable(
            data=df_pressure_melt,
            gas_name='xenon',
            x_col='T',
            y_col='y_exp',
            y_label=r'$p/\,\mathrm{MPa}$',
            title=None,
            model_x=df_pressure_melt['T'],
            model_y=df_pressure_melt['y_model'],
            # xlim=(60, 120),
            # ylim=(10**1, 10**5),
            logy=True,
            filename='xenon_melt_pressure.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )
        plot_variable_deviation(
            data=df_pressure_melt,
            gas_name='xenon',
            x_col='T',
            y_exp_col='y_exp',
            y_model_col='y_model',
            # xlim=(60, 120),
            ylim=(-15, 25),
            y_label=r'$100\cdot(p_{\mathrm{exp}}-p_{\mathrm{calc}})/p_{\mathrm{exp}}$',
            title=None,
            filename='xenon_melt_pressure_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
        )
    # elif variable == 'Vm_highp':
    #         # If no data, exit early
    #     if df_cell_volume_highp.empty:
    #         print("No high-pressure Vm data available.")
    #         return
    #     plot_thermo_variable(
    #         data=df_cell_volume_highp,  # reuse helper
    #         gas_name='xenon',
    #         x_col='p_exp',                # use pressure on x-axis
    #         y_col='y_exp',
    #         y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3\,mol^{-1}}$',
    #         title=None,
    #         model_x=df_cell_volume_highp['p_exp'],
    #         model_y=df_cell_volume_highp['y_model'],
    #         logy=False,
    #         filename='xenon_highp_cellvolume.png',
    #         output_folder=IMG_OUTPUT_FOLDER,
    #         custom_colors=CUSTOMCOLORS,
    #         custom_markers=CUSTOMMARKERS,
    #         ishighp=True,

    #     )
    #     # Overlay smooth per-group EOS curves
    #     ax = plt.gca()
    #     curves = _smooth_highp_curves(
    #         df_cell_volume_highp, params_fit, IDX["Vm"], npts=300)
    #     for c in curves:
    #         mask = np.isfinite(c['Vm'])
    #         if mask.sum() < 5:
    #             continue
    #         ax.plot(c['p'][mask], c['Vm'][mask],
    #                 linewidth=1.3, alpha=0.9,
    #                 label=f"{c['year']}, {c['author']} model")

    #     # Deduplicate legend labels
    #     h, lab = ax.get_legend_handles_labels()
    #     seen = set()
    #     h2, lab2 = [], []
    #     for hh, ll in zip(h, lab):
    #         if ll not in seen:
    #             seen.add(ll)
    #             h2.append(hh)
    #             lab2.append(ll)
    #     ax.legend(h2, lab2, loc='upper left',
    #               bbox_to_anchor=(1.05, 1), fontsize=8)
    #     plt.tight_layout(rect=[0, 0, 0.85, 1])
    #     plt.savefig(os.path.join(IMG_OUTPUT_FOLDER, 'krypton_highp_cellvolume.png'),
    #                 dpi=300, bbox_inches='tight')
    #     plt.show()
    #     plot_variable_deviation(
    #         data=df_cell_volume_highp,
    #         gas_name='xenon',
    #         x_col='p_exp',
    #         y_exp_col='y_exp',
    #         y_model_col='y_model',
    #         y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}})/V_{\mathrm{m,exp}}$',
    #         title=None,
    #         filename='xenon_highp_cellvolume_deviation',
    #         # xlim=(0, 120000),
    #         # ylim example: (-2, 2),
    #         output_folder=IMG_OUTPUT_FOLDER,
    #         custom_colors=CUSTOMCOLORS,
    #         custom_markers=CUSTOMMARKERS
    #     )
    elif variable == 'Vm_highp':
        if df_cell_volume_highp.empty:
            print("No high-pressure Vm data available.")
            return

        # 1) Scatter with your house style (no global model line)
        plot_thermo_variable(
            data=df_cell_volume_highp,
            gas_name='xenon',
            x_col='p_exp',
            y_col='y_exp',
            y_label=r'$V_{\mathrm{m}}\,/\,\mathrm{cm^3\,mol^{-1}}$',
            title=None,
            model_x=None, model_y=None,            # <- important
            logy=False,
            filename='xenon_highp_cellvolume.png',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS,
            ishighp=True,
            params_fit=params_fit,     # << pass your EOS params
            vm_idx=IDX["Vm"],          # << position of Vm in your model result
        )
        # Optional: deviation plot stays as-is
        plot_variable_deviation(
            data=df_cell_volume_highp,
            gas_name='xenon',
            x_col='p_exp',
            y_exp_col='y_exp',
            y_model_col='y_model',
            y_label=r'$100 \cdot (V_{\mathrm{m,exp}} - V_{\mathrm{m,calc}})/V_{\mathrm{m,exp}}$',
            title=None,
            filename='xenon_highp_cellvolume_deviation',
            output_folder=IMG_OUTPUT_FOLDER,
            custom_colors=CUSTOMCOLORS,
            custom_markers=CUSTOMMARKERS
        )



def RMS_AAD():
    xenon_data = load_all_gas_data('xenon', read_from_excel=False)
    datasets, meta = extract_datasets_with_meta(xenon_data)
    params_fit = PARAMS_INIT
    master_df = build_master_pointwise_df(datasets, meta, params_fit)
    summary = summarise_by_author(master_df)
    # Save to CSV
    output_path = os.path.join(
        IMG_OUTPUT_FOLDER, 'xenon_summary_by_author.csv')
    summary.to_csv(output_path, index=False)


def plot_init():
    xenon_data = load_all_gas_data('xenon', read_from_excel=False)
    datasets = extract_datasets(xenon_data)
    params_init = PARAMS_INIT
    plot_all_overlays_grid(params_init, datasets, Tt=Tt, pt=pt, compute_thermo_props=compute_thermo_props,
                           St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP, psub_curve=psub_curve, pmelt_curve=pmelt_curve)


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
