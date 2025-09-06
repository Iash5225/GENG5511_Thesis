from datetime import datetime
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from thermal_script import *
from computethermoprops import *
from constants import *
from scipy.interpolate import UnivariateSpline

# solid_EOS_fitting.py
BIG = 1e4
W_SUB = 1.0
W_MELT = 0.0
W_ANCH = 1.2   # try 0.3–0.5
W_SLOPE = 0.8   # try 0.05–0.15
Vm_anchor = 0.8
T_MARGIN = 2.0  # K, stay below Tt by this margin

def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB,  KRYPTON_E_3_SUB,  KRYPTON_T_t, KRYPTON_P_t
)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t
)

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

# Hidden Functions
def _rmse_abs(y, yhat):
    m = np.isfinite(y) & np.isfinite(yhat)
    return BIG if not np.any(m) else float(np.sqrt(np.mean((y[m]-yhat[m])**2)))

def _rms_weighted(y_true, y_pred, w, floor=1e-12):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    w = np.asarray(w,      float)
    m = np.isfinite(y_true) & np.isfinite(y_pred) & np.isfinite(w)
    if not np.any(m):
        return BIG
    # relative error with small floor to avoid divide-by-zero
    err = (y_true[m] - y_pred[m]) / np.maximum(np.abs(y_true[m]), floor)
    w = w[m]
    w /= np.mean(w)  # normalize so scale is comparable to unweighted RMS
    return float(np.sqrt(np.mean(w * err * err)))

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
    p_Vm_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_Vm_sub])
    Vm_sub = data['cell_volume_sub']['Cell Volume']

    # Cell Volume Melting
    T_Vm_melt = data['cell_volume_melt']['Temperature']
    p_Vm_melt = np.array([melting_pressure_equation(
        T, KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t) for T in T_Vm_melt])
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
    p_cp_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_cp_sub])
    cp_sub = data['heat_capacity']['Heat Capacity']

    # Thermal Expansion Sublimation
    T_alpha_sub = data['thermal_coeff']['Temperature']
    p_alpha_sub = np.array([sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t) for T in T_alpha_sub])
    alpha_sub = data['thermal_coeff']['Thermal Expansion Coefficient']

    # Bulk Modulus (S)
    T_BetaS_sub = data['bulk_s']['Temperature']
    # p_BetaS_sub = data['bulk_s']['Pressure']
    p_BetaS_sub = None
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


def _abs_rmse(y, yhat):
    """Absolute RMSE in the same units as y (returns BIG if nothing valid)."""
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    m = np.isfinite(y) & np.isfinite(yhat)
    if not np.any(m):
        return BIG
    return float(np.sqrt(np.mean((y[m] - yhat[m])**2)))


def diagnose_subline(params, datasets, Tt=KRYPTON_T_t):
    (T_Vm_sub, p_Vm_sub, Vm_sub, *_) = datasets
    T_hi = float(min(110.0, Tt-2.0))
    win = (np.isfinite(T_Vm_sub) & np.isfinite(Vm_sub) &
           (T_Vm_sub >= 8.0) & (T_Vm_sub <= T_hi))
    Tsub = np.asarray(T_Vm_sub[win], float)
    Vexp = np.asarray(Vm_sub[win],   float)

    def med(Tc, half=5.0):
        m = np.abs(Tsub - Tc) <= half
        return float(np.nanmedian(Vexp[m])) if np.any(m) else float(np.nanmedian(Vexp))
    T_lo = 20.0
    V_lo_exp, V_hi_exp = med(T_lo), med(T_hi)
    V_lo_mod = _safe_props_vector(
        [T_lo], [safe_psub(T_lo)], params, idx=IDX["Vm"])[0]
    V_hi_mod = _safe_props_vector(
        [T_hi], [safe_psub(T_hi)], params, idx=IDX["Vm"])[0]
    Vmod = _safe_props_vector(Tsub, safe_psub(Tsub), params, idx=IDX["Vm"])
    ok = np.isfinite(Vmod)
    k_exp = np.polyfit(Tsub, Vexp, 1)[0]
    k_mod = np.polyfit(Tsub[ok], Vmod[ok], 1)[
        0] if np.sum(ok) >= 2 else float('nan')
    print(f"Δ@20K={V_lo_mod-V_lo_exp:+.3f}, Δ@{T_hi:.0f}K={V_hi_mod-V_hi_exp:+.3f}, slope_mod={k_mod:.5f}, slope_exp={k_exp:.5f}")

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

def _make_outfun_vm(*datasets):
    """Per-iteration logger for Stage A (Vm-only)."""
    def cb(xk):
        tot, dev = combined_cost_vm(xk, *datasets)
        # history["total"].append(tot)
        # for k, v in dev.items():
        #     history[k].append(v)
        print(f"[Stage A] total={tot:.6g} | Vm_sub={dev.get('Vm_sub', np.nan):.6g} | "
              f"Vm_melt={dev.get('Vm_melt', np.nan):.6g} | "
              f"anchor={dev.get('Vm_anchor', np.nan):.6g} | slope={dev.get('Vm_slope', np.nan):.6g}")
    return cb

def _safe_props_vector_cached(T_arr, p_arr, params, idx, cache):
    T_arr = np.asarray(T_arr, float)
    p_arr = np.asarray(p_arr, float)
    out = np.full(T_arr.shape, np.nan, dtype=float)
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        key = (float(T), float(p))
        props = cache.get(key)
        if props is None:
            try:
                props = compute_thermo_props(float(T), float(p), params)
            except Exception:
                props = None
            cache[key] = props
        if props is not None and np.all(np.isfinite(props)):
            out[i] = props[IDX["Vm"]]
    return out

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

