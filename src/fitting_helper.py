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
slope_pen = 0.2
# --- safe pressure: never feed invalid (T,p) to the EOS ---
T_MARGIN = 2.0  # K, stay below Tt by this margin
# St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
# Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
# Tt = KRYPTON_T_t
# pt = KRYPTON_P_t
# IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)

def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB,  KRYPTON_E_3_SUB,  KRYPTON_T_t, KRYPTON_P_t
)

# --- one-time teacher pack from experimental subline ---


T_MARGIN = 2.0  # keep 2 K below Tt


def _sort_and_dedup(T, V):
    """Return T,V sorted by T and deduplicated (average V for same T)."""
    T = np.asarray(T, float)
    V = np.asarray(V, float)
    m = np.isfinite(T) & np.isfinite(V)
    T, V = T[m], V[m]

    # sort
    order = np.argsort(T)
    T, V = T[order], V[order]

    # collapse duplicates: average Vm for identical T
    Tu, inv = np.unique(T, return_inverse=True)
    if Tu.size != T.size:
        counts = np.bincount(inv)
        sums = np.bincount(inv, weights=V)
        Vu = sums / counts
        T, V = Tu, Vu
    return T, V


def make_subline_teacher(datasets, Tt, smooth_scale=0.02, Ngrid=120):
    (T_Vm_sub, _p_Vm_sub, Vm_sub, *_) = datasets

    Tmin = 8.0
    Tmax = float(min(110.0, Tt - T_MARGIN))

    mask = (np.isfinite(T_Vm_sub) & np.isfinite(Vm_sub) &
            (T_Vm_sub >= Tmin) & (T_Vm_sub <= Tmax))
    Tsub = np.asarray(T_Vm_sub[mask], float)
    Vexp = np.asarray(Vm_sub[mask],   float)

    # sort + dedup so UnivariateSpline gets strictly increasing x
    Tsub, Vexp = _sort_and_dedup(Tsub, Vexp)

    # need at least k+1 points (k=3 for cubic)
    if Tsub.size < 4:
        raise ValueError("Not enough subline points for teacher spline.")

    # light smoothing; scale s to data size/variance to be unit-aware
    s = max(1, Tsub.size) * smooth_scale * max(1e-6, np.var(Vexp))

    spl = UnivariateSpline(Tsub, Vexp, s=s)  # now x is strictly increasing

    Tg = np.linspace(Tmin, Tmax, Ngrid)
    Vg_target = spl(Tg)

    # anchors & experimental slope for consistency
    def _med(Tc, half=5.0):
        m = np.abs(Tsub - Tc) <= half
        return float(np.nanmedian(Vexp[m])) if np.any(m) else float(np.nanmedian(Vexp))

    T_lo, T_hi = 20.0, Tmax
    V_lo_exp, V_hi_exp = _med(T_lo), _med(T_hi)
    k_exp = np.polyfit(Tsub, Vexp, 1)[0]

    return dict(Tg=Tg, Vg_target=Vg_target,
                T_lo=T_lo, T_hi=T_hi,
                V_lo_exp=V_lo_exp, V_hi_exp=V_hi_exp,
                k_exp=k_exp)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t
)

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

def _mean_sq(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 0.0 if x.size == 0 else np.mean(x*x)

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


# ---- helpers ----


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

# def combined_cost_vm(params, *datasets):
#     (T_Vm_sub,  p_Vm_sub,  Vm_sub,
#      T_Vm_melt, p_Vm_melt, Vm_melt, *_) = datasets

#     try:
#         Tt_val = float(Tt)
#     except NameError:
#         Tt_val = float(KRYPTON_T_t)

#     Wc = {
#         "SUB":   W_SUB,
#         "MELT":  W_MELT,
#         "ANCH":  W_ANCH,
#         "SLOPE": W_SLOPE,
#         "CURV":  0.05,
#         "MONO":  5.0,
#         "SLOPE_SIGN": 1.0,
#     }

#     # ---- cache for this evaluation ----
#     _cache = {}

#     # ---- RMSE (all points) ----
#     Vm_sub_dev = BIG
#     msub = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub) &
#             np.isfinite(Vm_sub) & (T_Vm_sub <= Tt_val))
#     if np.any(msub):
#         Vmod_sub = _safe_props_vector_cached(
#             T_Vm_sub[msub], p_Vm_sub[msub], params, idx=IDX["Vm"], cache=_cache)
#         Vm_sub_dev = _abs_rmse(Vm_sub[msub], Vmod_sub)

#     Vm_melt_dev = BIG
#     mmelt = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt) &
#              np.isfinite(Vm_melt) & (T_Vm_melt >= Tt_val))
#     if np.any(mmelt):
#         Vmod_melt = _safe_props_vector_cached(
#             T_Vm_melt[mmelt], p_Vm_melt[mmelt], params, idx=IDX["Vm"], cache=_cache)
#         Vm_melt_dev = _abs_rmse(Vm_melt[mmelt], Vmod_melt)

#     # ---- Subline helpers ----
#     Vm_anchor = 0.0
#     slope_pen = 0.0
#     curv_pen = 0.0
#     mono_pen = 0.0
#     slope_sign_pen = 0.0

#     try:
#         T_hi_cap = float(min(110.0, Tt_val - 2.0))
#         win = (np.isfinite(T_Vm_sub) & np.isfinite(Vm_sub) &
#                (T_Vm_sub >= 8.0) & (T_Vm_sub <= T_hi_cap))
#         if np.any(win):
#             Tsub = np.asarray(T_Vm_sub[win], float)
#             Vexp = np.asarray(Vm_sub[win],   float)

#             # anchors
#             def _median_near(Tc, half=5.0):
#                 m = np.abs(Tsub - Tc) <= half
#                 return float(np.nanmedian(Vexp[m])) if np.any(m) else float(np.nanmedian(Vexp))
#             T_lo, T_hi = 20.0, T_hi_cap
#             V_lo_exp = _median_near(T_lo)
#             V_hi_exp = _median_near(T_hi)

#             V_lo_model = _safe_props_vector_cached(
#                 [T_lo], [safe_psub(T_lo)], params, idx=IDX["Vm"], cache=_cache)[0]
#             V_hi_model = _safe_props_vector_cached(
#                 [T_hi], [safe_psub(T_hi)], params, idx=IDX["Vm"], cache=_cache)[0]
#             if np.isfinite(V_lo_model) and np.isfinite(V_hi_model):
#                 Vm_anchor = (V_lo_model - V_lo_exp)**2 + \
#                     (V_hi_model - V_hi_exp)**2

#             # slope/curv/monotonicity (reuse cache)
#             p_win = safe_psub(Tsub)
#             Vmod = _safe_props_vector_cached(
#                 Tsub, p_win, params, idx=IDX["Vm"], cache=_cache)
#             m_ok = np.isfinite(Vmod)
#             if np.sum(m_ok) >= 4:
#                 k_exp = np.polyfit(Tsub, Vexp, 1)[0]
#                 k_mod = np.polyfit(Tsub[m_ok], Vmod[m_ok], 1)[0]
#                 slope_pen = (k_mod - k_exp)**2
#                 if np.sign(k_mod) != np.sign(k_exp):
#                     slope_sign_pen = (abs(k_mod) + abs(k_exp))**2
#                 # optional curvature (tiny)
#                 if Tsub.size >= 3 and np.sum(m_ok) >= 3:
#                     c_exp = np.polyfit(Tsub, Vexp, 2)[0]
#                     c_mod = np.polyfit(Tsub[m_ok], Vmod[m_ok], 2)[0]
#                     curv_pen = (c_mod - c_exp)**2
#                 if k_mod < 0.0:
#                     mono_pen = (-k_mod)**2
#     except Exception:
#         pass

#     total = (Wc["SUB"] * Vm_sub_dev +
#              Wc["MELT"] * Vm_melt_dev +
#              Wc["ANCH"] * Vm_anchor +
#              Wc["SLOPE"] * slope_pen +
#              Wc["CURV"] * curv_pen +
#              Wc["MONO"] * mono_pen +
#              Wc["SLOPE_SIGN"] * slope_sign_pen)

#     deviations = {
#         "Vm_sub":        Vm_sub_dev,
#         "Vm_melt":       Vm_melt_dev,
#         "Vm_anchor":     Vm_anchor,
#         "Vm_slope":      slope_pen,
#         "Vm_curv":       curv_pen,
#         "Vm_mono":       mono_pen,
#         "Vm_slope_sign": slope_sign_pen,
#     }
#     return total, deviations


# def combined_cost_vm(params, *datasets, teacher=None, Tt=None,
#                      weights=None):
#     """
#     Vm-only objective with:
#       - RMSE on sub & melt
#       - anchors at ~20 K and ~T_hi
#       - slope match on sub window
#       - teacher penalty (dense, smoothed target curve)
#       - monotonicity (dV/dT >= 0) and convexity (d2V/dT2 >= 0) on sub window
#     """
#     (T_Vm_sub,  p_Vm_sub,  Vm_sub,
#      T_Vm_melt, p_Vm_melt, Vm_melt, *_) = datasets

#     Tt_val = float(Tt if Tt is not None else KRYPTON_T_t)
#     Tmax_sub = float(min(110.0, Tt_val - T_MARGIN))

#     # default weights
#     W = dict(SUB=1.0, MELT=0.0, ANCH=1.2, SLOPE=0.7,
#              TEACH=1.0, MONO=4.0, CONV=0.3)
#     if weights:
#         W.update(weights)

#     _cache = {}  # per-eval cache

#     # ---------- core RMSE terms ----------
#     Vm_sub_dev = BIG
#     msub = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub) &
#             np.isfinite(Vm_sub) & (T_Vm_sub <= Tt_val))
#     if np.any(msub):
#         Vmod_sub = _safe_props_vector_cached(
#             T_Vm_sub[msub], p_Vm_sub[msub], params, idx=IDX["Vm"], cache=_cache)
#         Vm_sub_dev = _abs_rmse(Vm_sub[msub], Vmod_sub)

#     Vm_melt_dev = BIG
#     mmelt = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt) &
#              np.isfinite(Vm_melt) & (T_Vm_melt >= Tt_val))
#     if np.any(mmelt):
#         Vmod_melt = _safe_props_vector_cached(
#             T_Vm_melt[mmelt], p_Vm_melt[mmelt], params, idx=IDX["Vm"], cache=_cache)
#         Vm_melt_dev = _abs_rmse(Vm_melt[mmelt], Vmod_melt)

#     # ---------- sub-window helpers ----------
#     Vm_anchor = 0.0
#     slope_pen = 0.0
#     teacher_pen = 0.0
#     mono_pen = 0.0
#     conv_pen = 0.0

#     # window (avoid ultra-low T noise; remain below Tt)
#     win = (np.isfinite(T_Vm_sub) & np.isfinite(Vm_sub) &
#            (T_Vm_sub >= 8.0) & (T_Vm_sub <= Tmax_sub))
#     if np.any(win):
#         Tsub = np.asarray(T_Vm_sub[win], float)
#         Vexp = np.asarray(Vm_sub[win],   float)

#         # anchors (experiment medians vs model at same T on p_sub)
#         def _median_near(Tc, half=5.0):
#             m = np.abs(Tsub - Tc) <= half
#             return float(np.nanmedian(Vexp[m])) if np.any(m) else float(np.nanmedian(Vexp))
#         T_lo = 20.0
#         T_hi = Tmax_sub
#         V_lo_exp = teacher["V_lo_exp"] if teacher else _median_near(T_lo)
#         V_hi_exp = teacher["V_hi_exp"] if teacher else _median_near(T_hi)

#         V_lo_model = _safe_props_vector_cached([T_lo], [safe_psub(T_lo)],
#                                                params, idx=IDX["Vm"], cache=_cache)[0]
#         V_hi_model = _safe_props_vector_cached([T_hi], [safe_psub(T_hi)],
#                                                params, idx=IDX["Vm"], cache=_cache)[0]
#         if np.isfinite(V_lo_model) and np.isfinite(V_hi_model):
#             Vm_anchor = (V_lo_model - V_lo_exp)**2 + (V_hi_model - V_hi_exp)**2

#         # slope match (robust)
#         Vmod_subline = _safe_props_vector_cached(
#             Tsub, safe_psub(Tsub), params, idx=IDX["Vm"], cache=_cache)
#         ok = np.isfinite(Vmod_subline) & (
#             26.0 < Vmod_subline) & (Vmod_subline < 31.5)
#         if np.sum(ok) >= 4:
#             k_exp = (teacher["k_exp"] if teacher is not None
#                      else np.polyfit(Tsub, Vexp, 1)[0])
#             k_mod = np.polyfit(Tsub[ok], Vmod_subline[ok], 1)[0]
#             slope_pen = (k_mod - k_exp)**2

#     # ---------- teacher curve + shape enforcement ----------
#     if teacher is not None:
#         Tg = teacher["Tg"]
#         Vg_target = teacher["Vg_target"]
#         Vg_model = _safe_props_vector_cached(
#             Tg, safe_psub(Tg), params, idx=IDX["Vm"], cache=_cache)

#         mask = np.isfinite(Vg_model)
#         if np.any(mask):
#             teacher_pen = np.nanmean((Vg_model[mask] - Vg_target[mask])**2)

#             # local monotonicity & convexity (forward differences)
#             Tg_m = Tg[mask]
#             Vg_m = Vg_model[mask]
#             if Vg_m.size >= 4:
#                 dT = np.diff(Tg_m)
#                 dV = np.diff(Vg_m) / dT
#                 mono_pen = np.nanmean(np.clip(-dV, 0.0, None)**2)

#                 d2V = np.diff(Vg_m, 2) / (dT[:-1]*dT[1:])
#                 conv_pen = np.nanmean(np.clip(-d2V, 0.0, None)**2)

#     # ---------- total ----------
#     total = (W["SUB"]*Vm_sub_dev +
#              W["MELT"]*Vm_melt_dev +
#              W["ANCH"]*Vm_anchor +
#              W["SLOPE"]*slope_pen +
#              W["TEACH"]*teacher_pen +
#              W["MONO"]*mono_pen +
#              W["CONV"]*conv_pen)

#     deviations = dict(
#         Vm_sub=Vm_sub_dev, Vm_melt=Vm_melt_dev,
#         Vm_anchor=Vm_anchor, Vm_slope=slope_pen,
#         Vm_teacher=teacher_pen, Vm_mono=mono_pen, Vm_conv=conv_pen
#     )
#     return total, deviations

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
def prefit_v00(params, datasets, grid=np.linspace(20, 32, 25)):
    best = None
    base = np.array(params, float)
    for v00 in grid:
        p = base.copy()
        p[0] = float(v00)
        tot, _ = combined_cost_vm(p, *datasets)
        if np.isfinite(tot) and (best is None or tot < best[0]):
            best = (tot, v00)
    print(f"[Stage A] prefit v00 -> {best[1]:.3f} (cost {best[0]:.3f})")
    base[0] = best[1]
    return base


def make_stageA_bounds(params_init,
                       lower, upper,
                       free_elastic=(0, 1, 2, 3),
                       v00_idx=None,
                       v00_range=(25,30),
                       elastic_range=(-1e6, 1e6),
                       free_theta_idx=None,
                       theta_range=(200.0, 2000.0)):
    """
    Build bounds for Stage A (Vm-only):
      - Free elastic block indices (default {0,1,2,3})
        * v00 gets v00_range
        * the rest get elastic_range
      - Clamp all other parameters to their starting values
      - Optionally free one thermal knob (e.g., θD0) with theta_range

    Parameters
    ----------
    params_init : array-like
        Starting parameter vector (used for clamping).
    lower, upper : array-like
        Global lower/upper bounds (same length as params_init).
    free_elastic : tuple[int]
        Indices to free as "elastic" parameters (e.g., (0,1,2,3)).
    v00_idx : int or None
        Index of v00 (defaults to the first element of free_elastic).
    v00_range : (float, float)
        Bounds for v00 (default (10, 40)).
    elastic_range : (float, float)
        Bounds for other elastic params (default (-1e6, 1e6)).
    free_theta_idx : int or None
        Optional index of a thermal knob to free (e.g., θD0).
    theta_range : (float, float)
        Bounds for the thermal knob if freed (default (200, 2000)).

    Returns
    -------
    list[tuple[float, float]]
        Bounds list suitable for scipy.optimize.minimize(..., bounds=...).
    """
    p0 = np.asarray(params_init, float)
    LB = np.asarray(lower, float)
    UB = np.asarray(upper, float)
    assert p0.shape == LB.shape == UB.shape, "params/bounds length mismatch"

    # Start with clamped bounds at the initial values
    B = [(float(v), float(v)) for v in p0]

    # Determine v00 index
    if v00_idx is None:
        if len(free_elastic) == 0:
            raise ValueError("free_elastic is empty; cannot infer v00_idx")
        v00_idx = free_elastic[0]

    # Free v00
    B[v00_idx] = (float(v00_range[0]), float(v00_range[1]))

    # Free other elastic parameters
    for i in free_elastic:
        if i == v00_idx:
            continue
        B[i] = (float(elastic_range[0]), float(elastic_range[1]))

    # Optionally free a thermal knob (e.g., θD0)
    if free_theta_idx is not None:
        B[free_theta_idx] = (float(theta_range[0]), float(theta_range[1]))

    return B
# Helper: cost wrapper used above (Vm-only)


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



# 1) Safe wrapper around psub


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


def target_v00(datasets):
    (T_Vm_sub, _, Vm_sub, *_) = datasets
    m = np.isfinite(T_Vm_sub) & np.isfinite(
        Vm_sub) & (np.abs(T_Vm_sub-20.0) <= 5.0)
    return float(np.nanmedian(Vm_sub[m])) if np.any(m) else float(np.nanmedian(Vm_sub))
__all__ = [
        "psub_curve",
        "pmelt_curve",
        "rms_log",
        "rms_percent",
        "rms",
        "_mean_sq",
        "_rmse_abs",
        "_rms_weighted",
        "extract_datasets",
        "_safe_props_vector",
        "_abs_rmse",
        "combined_cost_vm",
        "prefit_v00",
        "make_stageA_bounds",
        "_cost_only_vm",
        "_make_outfun_vm",
        "quick_plot_vm_only",
        "_safe_props_vector_cached",
        "diagnose_subline",
        "safe_psub",
        "make_subline_teacher",
        "target_v00",
    ]
