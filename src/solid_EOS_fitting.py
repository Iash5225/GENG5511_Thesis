from math import isfinite
from scipy.optimize import brentq
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from computethermoprops import *
from p_functions import pmelt,psub
from constants import *
from fitting_helper import rms , _mean_sq
import math
from plot_eos import plot_all_overlays_grid

BIG = 1e4
# Constants
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t

IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)

# --- keep only THIS version of _make_outfun ---
def _make_outfun(*datasets):
    def outfun(xk):
        total, dev = combined_cost_function(xk, *datasets)
        # record values
        history["total"].append(total)
        for k in dev:
            history[k].append(dev[k])
        # optional: still print
        print(f"Current total_deviation: {total:.6g}")
        for k, v in dev.items():
            print(f"  {k}: {v:.6g}")
    return outfun


def _pf_equilibrium(T, p_guess, Gf_at_pM, Vf_at_pM, params, compute_thermo_props,
                    p_lo=1e-6, p_hi=8.0e4,  # MPa range wide enough for your data
                    expand=3.0, max_expand=1e6):
    """
    Find p_f(T) such that G_model(T,p) - G_fluid(T,p) = 0.
    We start near p_guess and expand a bracket until the sign changes.
    """
    # Define ΔG(T,p) using your already-available fluid data at this T:
    # If you can evaluate fluid G,V at arbitrary p, use that. Otherwise,
    # approximate G_fluid(T,p) ~ Gf_at_pM + (p - pM)*Vf_at_pM.
    pM = float(p_guess)
    GfM = float(Gf_at_pM)
    VfM = float(Vf_at_pM)  # cm^3/mol

    def dG(p):
        pr = compute_thermo_props(T, p, params)
        if not (isfinite(pr[11]) and isfinite(pr[0])):   # G and Vm
            return np.nan
        # linearized fluid G around pM (if no EOS for the fluid):
        G_fluid_p = GfM + (p - pM) * VfM
        return pr[11] - G_fluid_p

    # Try to bracket around pM
    lo, hi = max(p_lo, 1e-9), min(p_hi, 1e9)
    a, b = max(lo, 0.5*pM), min(hi, 1.5*pM if pM > 0 else 1.0)
    fa, fb = dG(a), dG(b)

    # Expand until sign change or limits reached
    k = 0
    while (not (isfinite(fa) and isfinite(fb) and fa*fb < 0.0)) and k < 20:
        a = max(lo, a/expand)
        b = min(hi, b*expand)
        fa, fb = dG(a), dG(b)
        k += 1
        if b/a > max_expand:
            break

    if isfinite(fa) and isfinite(fb) and fa*fb < 0.0:
        try:
            return float(brentq(dG, a, b, maxiter=100))
        except Exception:
            return np.nan
    else:
        return np.nan

def rms_log(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    m = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true > 0) & (y_pred > 0)
    if not np.any(m):
        return BIG
    r = np.log(y_pred[m]) - np.log(y_true[m])
    return float(np.sqrt(np.mean(r*r)))


# ---------- High-T weighting + log-pressure RMS ----------

# smoothly up-weight points as T → Tmax (or as T → Tt on the low side)
def _temp_weights_high_end(T, high_anchor, exp=4.0, w_min=0.25, w_max=3.0):
    T = np.asarray(T, float)
    # map T into [0,1] against the chosen anchor (high end)
    span = max(abs(high_anchor - np.min(T)), 1e-9)
    x = np.clip((T - (high_anchor - span)) / span, 0.0, 1.0)
    w = w_min + (w_max - w_min) * (x ** exp)   # emphasize high T
    return w


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


def _rms_log_pressure(p_true, p_pred, w=None, p_floor=1e-3):
    p_true = np.asarray(p_true, float)
    p_pred = np.asarray(p_pred, float)
    m = np.isfinite(p_true) & np.isfinite(p_pred)
    if not np.any(m):
        return BIG
    lt = np.log(np.maximum(p_true[m], p_floor))
    lp = np.log(np.maximum(p_pred[m], p_floor))
    err = lt - lp
    if w is None:
        return float(np.sqrt(np.mean(err * err)))
    w = np.asarray(w, float)[m]
    w /= np.mean(w)
    return float(np.sqrt(np.mean(w * err * err)))


def _hiT_weights(T, Tt, Tscale=60.0, boost=8.0):
    # 1 below triple-point, grows ~exp as T rises
    x = np.maximum(0.0, (np.asarray(T, float) - Tt)/Tscale)
    return 1.0 + boost*(1.0 - np.exp(-x))**2

def plot_deviation_history():
    """Make two figures:
       1) total deviation vs iteration
       2) grid of per-metric deviations vs iteration
    """
    # ---- Figure 1: total deviation ----
    if len(history["total"]) == 0:
        print("No history recorded (total). Did the callback run?")
    else:
        fig1, ax1 = plt.subplots(figsize=(9, 5))
        ax1.plot(history["total"], label="total")
        ax1.set_yscale("log")
        ax1.set_xlabel("Iteration")
        ax1.set_ylabel("Total deviation")
        ax1.set_title("Total Deviation per Iteration")
        ax1.legend()
        plt.tight_layout()

    # ---- Figure 2: grid of per-metric deviations ----
    metric_keys = [k for k in history.keys() if k != "total"]
    # keep only metrics that actually have data
    metric_keys = [k for k in metric_keys if len(history[k]) > 0]

    if len(metric_keys) == 0:
        print("No per-metric history recorded. Did the callback append to history?")
        return

    # choose a neat grid shape
    n = len(metric_keys)
    # 3 rows x 4 cols works nicely up to 12 metrics; compute general shape:
    cols = min(4, n)
    rows = math.ceil(n / cols)

    fig2, axes = plt.subplots(rows, cols, figsize=(
        4*cols, 3.2*rows), squeeze=False)
    ax_iter = 0
    for r in range(rows):
        for c in range(cols):
            ax = axes[r][c]
            if ax_iter < n:
                key = metric_keys[ax_iter]
                ax.plot(history[key])
                ax.set_yscale("log")
                ax.set_title(key)
                ax.set_xlabel("Iteration")
                ax.set_ylabel("Deviation")
            else:
                # hide unused axes
                ax.axis("off")
            ax_iter += 1

    fig2.suptitle("Per-Metric Deviations per Iteration", y=1.02, fontsize=12)
    plt.tight_layout()
    plt.show()

def rms_percent(y_true, y_pred, scale=PERCENT_SCALE):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    mask = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true != 0.0)
    if not np.any(mask):
        return BIG
    err = scale * (y_true[mask] - y_pred[mask]) / y_true[mask]
    err = err[np.isfinite(err)]
    return BIG if err.size == 0 else float(np.sqrt(np.mean(err*err)))

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

    tp = compute_thermo_props(Tt, pt, params)

    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        triple_props, triple_meta = tp
        # pulled directly from compute_thermo_props
        deltaH_triple = float(triple_meta.get("deltaH_triple"))
        deltaS_triple = float(triple_meta.get("deltaS_triple"))
    else:
        # backward compatibility (old API returned only the 12-length props array)
        triple_props = np.asarray(tp, float)
        # S* at triple is params[30]; REFPROP values are global constants
        deltaS_triple = params[30] - St_REFPROP
        Ht_fitted = triple_props[11] + Tt * params[30]
        deltaH_triple = Ht_fitted - Ht_REFPROP

    # ====== BLOCKS ======

    # Vm (sublimation): T <= Tt
    m = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub)
         & np.isfinite(Vm_sub) & (T_Vm_sub <= Tt))
    Vm_sub_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_sub[m], p_Vm_sub[m], params, idx=0)
        Vm_sub_dev = rms_percent(Vm_sub[m], model)

    # Vm (melting): T >= Tt
    m = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt)
         & np.isfinite(Vm_melt) & (T_Vm_melt >= Tt))
    Vm_melt_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_melt[m], p_Vm_melt[m], params, idx=0)
        Vm_melt_dev = rms_percent(Vm_melt[m], model)

    # Vm (high-p): no phase constraint (already experimental p)
    m = (np.isfinite(T_Vm_highp) & np.isfinite(
        p_Vm_highp) & np.isfinite(Vm_highp))
    Vm_highp_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_highp[m], p_Vm_highp[m], params, idx=0)
        Vm_highp_dev = rms_percent(Vm_highp[m], model)

    # cp (sublimation)
    # m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub)
    #      & np.isfinite(cp_sub) & (T_cp_sub <= Tt))
    # cp_sub_dev = BIG
    # if np.any(m):
    #     model = _safe_props_vector(T_cp_sub[m], p_cp_sub[m], params, idx=4)
    #     # temperature-dependent weighting, but still penalize empties
    #     terms = []
    #     for Ti, cpi, cpm in zip(T_cp_sub[m], cp_sub[m], model):
    #         if np.isfinite(cpm) and cpi != 0:
    #             w = CP_WEIGHT_BELOW if Ti < CP_TEMP_THRESHOLD_K else CP_WEIGHT_ABOVE
    #             terms.append(w * (cpi - cpm) / cpi)
    #     cp_sub_dev = rms(np.array(terms))

    # # cp (sublimation)
    # m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub)
    #     & np.isfinite(cp_sub) & (T_cp_sub <= Tt))
    # cp_sub_dev = BIG

    # if np.any(m):
    #     model = _safe_props_vector(T_cp_sub[m], p_cp_sub[m], params, idx=4)
    #     # print("cp exp percentiles (J/mol-K):",
    #     #       np.nanpercentile(cp_sub[m], [5, 50, 95]))
    #     # print("cp model percentiles (J/mol-K):",
    #     #   np.nanpercentile(model,  [5, 50, 95]))
    #     cp_sub_dev = rms_percent(cp_sub[m], model)

        # cp (sublimation) with split weighting (low vs high T)
    m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub)
        & np.isfinite(cp_sub) & (T_cp_sub <= Tt))
    cp_sub_dev = BIG
    if np.any(m):
        Tm = np.asarray(T_cp_sub)[m]
        pm = np.asarray(p_cp_sub)[m]
        y_exp = np.asarray(cp_sub)[m]
        y_eos = _safe_props_vector(Tm, pm, params, idx=4)
        # weights: stronger below CP_SPLIT_K because cp -> 0
        w = np.where(Tm < CP_SPLIT_K, CP_W_BELOW, CP_W_ABOVE)
        # relative RMS with weights
        cp_sub_dev = _rms_weighted(y_exp, y_eos, w, floor=1e-6)


    # alpha (sublimation)
    m = (np.isfinite(T_alpha_sub) & np.isfinite(p_alpha_sub)
         & np.isfinite(alpha_sub) & (T_alpha_sub <= Tt))
    alpha_sub_dev = BIG
    if np.any(m):
        model = _safe_props_vector(
            T_alpha_sub[m], p_alpha_sub[m], params, idx=3)
        alpha_sub_dev = rms_percent(alpha_sub[m], model)

    # BetaT (sublimation)
    m = (np.isfinite(T_BetaT_sub) & np.isfinite(p_BetaT_sub)
         & np.isfinite(BetaT_sub) & (T_BetaT_sub <= Tt))
    BetaT_sub_dev = BIG
    if np.any(m):
        model = _safe_props_vector(
            T_BetaT_sub[m], p_BetaT_sub[m], params, idx=1)
        BetaT_sub_dev = rms_percent(BetaT_sub[m], model)

    # BetaS (sublimation) – note your finite count is only 25/74
    m = (np.isfinite(T_BetaS_sub) & np.isfinite(p_BetaS_sub)
         & np.isfinite(BetaS_sub) & (T_BetaS_sub <= Tt))
    BetaS_sub_dev = BIG
    if np.any(m):
        model = _safe_props_vector(
            T_BetaS_sub[m], p_BetaS_sub[m], params, idx=2)
        BetaS_sub_dev = rms_percent(BetaS_sub[m], model)

    # Enthalpy of sublimation: kJ/mol (no unit juggling)
    m = (np.isfinite(T_H_sub) & np.isfinite(p_H_sub) &
         np.isfinite(delta_H_sub) & np.isfinite(H_fluid_sub) & (T_H_sub <= Tt))
    H_solid_sub_dev = BIG
    if np.any(m):
        H_solid_sub = H_fluid_sub[m] - delta_H_sub[m] * 1.0  # both kJ/mol
        modelH = _safe_props_vector(
            T_H_sub[m], p_H_sub[m], params, idx=10) - deltaH_triple
        H_solid_sub_dev = rms_percent(H_solid_sub, modelH)

    # Enthalpy of melting: kJ/mol
    m = (np.isfinite(T_H_melt) & np.isfinite(p_H_melt) &
         np.isfinite(delta_H_melt) & np.isfinite(H_fluid_melt) & (T_H_melt >= Tt))
    H_solid_melt_dev = BIG
    if np.any(m):
        H_solid_melt = H_fluid_melt[m] - delta_H_melt[m] * 1.0
        modelH = _safe_props_vector(
            T_H_melt[m], p_H_melt[m], params, idx=10) - deltaH_triple
        H_solid_melt_dev = rms_percent(H_solid_melt, modelH)

    # # Sublimation pressure (all MPa)
    # m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
    #      np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (T_sub <= Tt))
    # p_sub_dev = BIG
    # if np.any(m):
    #     pf = np.full_like(p_sub[m], np.nan, dtype=float)
    #     for i, (T, pM, Gf, Vf) in enumerate(zip(T_sub[m], p_sub[m], G_fluid_sub[m], V_fluid_sub[m])):
    #         props = compute_thermo_props(T, pM, params)
    #         if np.all(np.isfinite(props)) and np.isfinite(Vf) and (Vf != props[0]):
    #             dG = Gf - props[11] + deltaH_triple - T * deltaS_triple
    #             pf[i] = pM - dG / (Vf - props[0])  # MPa
    #     #p_sub_dev = rms_percent(p_sub[m], pf)
    #     p_sub_dev = rms_log(p_sub[m], pf)


        # --- Sublimation pressure (all MPa, T <= Tt) ---
    # m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
    #     np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (T_sub <= Tt))
    # p_sub_dev = BIG
    # if np.any(m):
    #     T = np.asarray(T_sub)[m]
    #     pM = np.asarray(p_sub)[m]
    #     Gf = np.asarray(G_fluid_sub)[m]
    #     Vf = np.asarray(V_fluid_sub)[m]

    #     pf = np.full_like(pM, np.nan, dtype=float)
    #     # compute model-implied p_f for each experimental point
    #     for i, (Ti, pMi, Gfi, Vfi) in enumerate(zip(T, pM, Gf, Vf)):
    #         props = compute_thermo_props(Ti, pMi, params)
    #         if np.all(np.isfinite(props)) and np.isfinite(Vfi) and (Vfi != props[0]):
    #             dG = Gfi - props[11] + deltaH_triple - Ti * deltaS_triple
    #             pf[i] = pMi - dG / (Vfi - props[0])  # MPa

    #     # temperature weights that grow smoothly toward the triple point (high end for subl.)
    #     w_T = _temp_weights_high_end(
    #         T, high_anchor=Tt, exp=2.0, w_min=0.25, w_max=3.0)

    #     # use log-residuals (balances decades of p); also uses the high-T weights
    #     p_sub_dev = _rms_log_pressure(pM, pf, w=w_T, p_floor=1e-3)

    # # Another method
    # m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
    #      np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (T_sub <= Tt))
    # p_sub_dev = BIG
    # if np.any(m):
    #     Tm, pM, Gf, Vf = T_sub[m], p_sub[m], G_fluid_sub[m], V_fluid_sub[m]
    #     pf = np.array([_pf_equilibrium(Ti, pMi, Gfi, Vfi, params, compute_thermo_props)
    #                 for Ti, pMi, Gfi, Vfi in zip(Tm, pM, Gf, Vf)])
    #     p_sub_dev = _rms_log_pressure(pM, pf, w=_hiT_weights(Tm, Tt))
    # replace the "Another method" block for p_sub_dev
    m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
        np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (T_sub <= Tt))
    p_sub_dev = BIG
    if np.any(m):
        Tm, pM, Gf, Vf = T_sub[m], p_sub[m], G_fluid_sub[m], V_fluid_sub[m]
        pf = np.full_like(pM, np.nan, dtype=float)
        for i, (Ti, pMi, Gfi, Vfi) in enumerate(zip(Tm, pM, Gf, Vf)):
            props = compute_thermo_props(Ti, pMi, params)
            if np.all(np.isfinite(props)) and np.isfinite(Vfi) and (Vfi != props[IDX["Vm"]]):
                dG = Gfi - props[IDX["G"]] + (deltaH_triple - Ti*deltaS_triple)
                pf[i] = pMi - dG / (Vfi - props[IDX["Vm"]])
        # keep the log RMS (good over decades) and your high-T weights
        p_sub_dev = _rms_log_pressure(pM, pf, w=_hiT_weights(Tm, Tt))
    # Melting pressure (all MPa; extra weight above threshold via factor)
    # m = (np.isfinite(T_melt) & np.isfinite(p_melt) &
    #      np.isfinite(G_fluid_melt) & np.isfinite(V_fluid_melt) & (T_melt >= Tt))
    # p_melt_dev = BIG
    # if np.any(m):
    #     pf = np.full_like(p_melt[m], np.nan, dtype=float)
    #     for i, (T, pM, Gf, Vf) in enumerate(zip(T_melt[m], p_melt[m], G_fluid_melt[m], V_fluid_melt[m])):
    #         props = compute_thermo_props(T, pM, params)
    #         if np.all(np.isfinite(props)) and np.isfinite(Vf) and (Vf != props[0]):
    #             fac = PMELT_EXTRA_FACTOR if T > PMELT_EXTRA_WEIGHT_T_K else 1.0
    #             dG = Gf - props[11] + deltaH_triple - T * deltaS_triple
    #             pf[i] = pM - fac * (dG / (Vf - props[0]))
    #     p_melt_dev = rms_percent(p_melt[m], pf)


        # --- Melting pressure (all MPa, T >= Tt) ---
    m = (np.isfinite(T_melt) & np.isfinite(p_melt) &
        np.isfinite(G_fluid_melt) & np.isfinite(V_fluid_melt) & (T_melt >= Tt))
    p_melt_dev = BIG
    if np.any(m):
        T = np.asarray(T_melt)[m]
        pM = np.asarray(p_melt)[m]
        Gf = np.asarray(G_fluid_melt)[m]
        Vf = np.asarray(V_fluid_melt)[m]

        pf = np.full_like(pM, np.nan, dtype=float)
        for i, (Ti, pMi, Gfi, Vfi) in enumerate(zip(T, pM, Gf, Vf)):
            props = compute_thermo_props(Ti, pMi, params)
            if np.all(np.isfinite(props)) and np.isfinite(Vfi) and (Vfi != props[0]):
                dG = Gfi - props[11] + deltaH_triple - Ti * deltaS_triple
                pf[i] = pMi - dG / (Vfi - props[0])

        # high-T weighting grows from Tt up to max(T)
        w_T = _temp_weights_high_end(
            T, high_anchor=np.max(T), exp=2.0, w_min=0.25, w_max=3.0)

        p_melt_dev = _rms_log_pressure(pM, pf, w=w_T, p_floor=1e-3)

    # Gamma-T smoothness penalty along a subl. grid (T<=Tt)
    T6 = np.array([0.0001] + list(range(2, 84, 2)) + [83.806])
    T6 = T6[T6 <= Tt]
    p6 = np.full_like(T6, np.nan, dtype=float)
    G6 = np.full_like(T6, np.nan, dtype=float)
    V6 = np.full_like(T6, np.nan, dtype=float)
    for i, T in enumerate(T6):
        p6[i] = psub(T, pt, Tt)
        pr = compute_thermo_props(T, p6[i], params)
        if np.all(np.isfinite(pr)):
            G6[i] = pr[6]   # Gruneisen
            V6[i] = pr[0]
    # finite pairs only
    mk = np.isfinite(G6) & np.isfinite(V6)
    Gamma_T6_dev = BIG
    if np.sum(mk) >= 3:
        Gm, Vm = G6[mk], V6[mk]
        slopes = np.diff(Gm) / np.diff(Vm)
        mu = np.nanmean(Gm)
        terms = []
        for i, s in enumerate(slopes, start=1):
            if not np.isfinite(s):
                continue
            if s > 0:
                terms.append(GAMMA_POS_SLOPE_OFFSET +
                             abs(GAMMA_POS_SLOPE_MULT * (Gm[i] - mu) / mu))
            else:
                terms.append(GAMMA_NEG_SLOPE_MULT * (Gm[i] - mu) / mu)
        Gamma_T6_dev = rms(np.array(terms))


    total_deviation = (
        Vm_sub_dev * W_VM_SUB +
        Vm_melt_dev * W_VM_MELT +
        Vm_highp_dev * W_VM_HIGHP +
        cp_sub_dev * W_CP_SUB +
        alpha_sub_dev * W_ALPHA_SUB +
        BetaT_sub_dev * W_BETAT_SUB +
        BetaS_sub_dev * W_BETAS_SUB +
        H_solid_sub_dev * W_H_SOLID_SUB +
        H_solid_melt_dev * W_H_SOLID_MELT +
        p_sub_dev * W_P_SUB +
        p_melt_dev * W_P_MELT +
        Gamma_T6_dev * W_GAMMA_T
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
            'maxiter': MAX_ITERATIONS,
            'ftol': FUNCTION_TOL,              # function tolerance
            'gtol': GRADIENT_TOL,              # projected gradient tol (stopping)
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


def _cost_only(params, *datasets):
    try:
        total, _ = combined_cost_function(params, *datasets)
        return float(total) if np.isfinite(total) else BIG
    except Exception:
        return BIG

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


def reset_history():
    """Clear deviation history before starting a new optimization."""
    for k in history:
        history[k].clear()


def combined_cost_vm(params, *datasets):
    (T_Vm_sub, p_Vm_sub, Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt,
     *rest) = datasets  # ignore the rest

    Vm_sub_dev = BIG
    m = (np.isfinite(T_Vm_sub) & np.isfinite(
        Vm_sub) & (T_Vm_sub <= KRYPTON_T_t))
    if np.any(m):
        model = _safe_props_vector(T_Vm_sub[m], p_Vm_sub[m], params, idx=0)
        Vm_sub_dev = rms_percent(Vm_sub[m], model)

    Vm_melt_dev = BIG
    m = (np.isfinite(T_Vm_melt) & np.isfinite(
        Vm_melt) & (T_Vm_melt >= KRYPTON_T_t))
    if np.any(m):
        model = _safe_props_vector(T_Vm_melt[m], p_Vm_melt[m], params, idx=0)
        Vm_melt_dev = rms_percent(Vm_melt[m], model)

    total = Vm_sub_dev * 55 + Vm_melt_dev * 35
    return total, {"Vm_sub": Vm_sub_dev, "Vm_melt": Vm_melt_dev}


# =========================
# ======== STAGE A ========
# Fit Vm only (sublimation + melting); freeze all other params
# =========================

def _cost_only_vm(params, *datasets):
    try:
        total, _ = combined_cost_vm(params, *datasets)
        return float(total) if np.isfinite(total) else BIG
    except Exception:
        return BIG


def _make_outfun_vm(*datasets):
    def outfun(xk):
        tot, dev = combined_cost_vm(xk, *datasets)
        # optional live log
        print(
            f"[Stage A] total: {tot:.6g} | Vm_sub: {dev['Vm_sub']:.6g} | Vm_melt: {dev['Vm_melt']:.6g}")
    return outfun


def _stageA_bounds_free_vm_only(params_init, lower, upper):
    """
    Build bounds where ONLY [v00, a1, a2, a3] are free.
    Everything else is clamped at its initial value.
    Layout from your code:
      0: v00
      1-3: a1, a2, a3
      others: frozen
    """
    params_init = np.asarray(params_init, float)
    lb = list(lower)
    ub = list(upper)

    free_idx = {0, 1, 2, 3}  # v00, a1, a2, a3
    B = []
    for i, (lo, hi) in enumerate(zip(lb, ub)):
        if i in free_idx:
            # keep your original bounds for the free ones, but make sure they are wide enough
            B.append((float(lo), float(hi)))
        else:
            # clamp at initial value
            B.append((float(params_init[i]), float(params_init[i])))
    return B


def run_stage_A(params_init, datasets, lower_bound, upper_bound,
                maxiter=300, ftol=1e-9, gtol=1e-3):
    print("=== Stage A: Fitting Vm only (sub+melt) ===")
    # build bounds
    boundsA = _stageA_bounds_free_vm_only(
        params_init, lower_bound, upper_bound)

    # quick sanity on initial cost
    tot0, dev0 = combined_cost_vm(params_init, *datasets)
    print(
        f"[Stage A] initial total = {tot0:.6g}  (Vm_sub={dev0['Vm_sub']:.6g}, Vm_melt={dev0['Vm_melt']:.6g})")

    cb = _make_outfun_vm(*datasets)

    res = minimize(
        fun=_cost_only_vm,
        x0=np.asarray(params_init, dtype=float),
        args=datasets,
        method="L-BFGS-B",
        bounds=boundsA,
        callback=cb,
        options=dict(
            disp=True,
            maxiter=maxiter,
            ftol=ftol,
            gtol=gtol,
            maxls=40,
        )
    )
    print("\n[Stage A] Optimization status:", res.message)
    print("[Stage A] Final cost:", res.fun)

    # pretty print only the free parameters
    names = ["v00", "a1", "a2", "a3"]
    idxs = [0, 1, 2, 3]
    print("\n[Stage A] Fitted elastic parameters:")
    for nm, i in zip(names, idxs):
        print(f"  {nm:4s}: {res.x[i]:.6g}")

    # save params so we can seed Stage B later
    from datetime import datetime
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    np.save(f"stageA_params_{stamp}.npy", res.x)
    print(f"[Stage A] Saved params -> stageA_params_{stamp}.npy")

    return res.x, res.fun


def stage_A():
    # 1) load data
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)

    # 2) build Stage A bounds and run
    boundsA = list(zip(LOWER_BOUND, UPPER_BOUND))
    paramsA, fA = run_stage_A(PARAMS_INIT, datasets, LOWER_BOUND, UPPER_BOUND)

    # 3) plot overlays using Stage-A params (so you can visually check Vm)
    plot_all_overlays_grid(
        params=paramsA,
        datasets=datasets,
        Tt=KRYPTON_T_t,
        pt=KRYPTON_P_t,
        compute_thermo_props=compute_thermo_props,
        St_REFPROP=KRYPTON_REFERENCE_ENTROPY,
        Ht_REFPROP=KRYPTON_REFERENCE_ENTHALPY,
        ncols=3,
        figsize=(14, 10)
    )

    # 4) print compact list for your log
    formatted = ", ".join(f"{p:.6g}" for p in paramsA)
    print("\n[Stage A] Full parameter vector (with non-frees unchanged):")
    print(formatted)
    print("[Stage A] Final cost:", fA)

def main():
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)
    bounds = list(zip(LOWER_BOUND, UPPER_BOUND))
    params_fit, fval = run_optimization(PARAMS_INIT, bounds, datasets)
    plot_deviation_history()
    plot_all_overlays_grid(
        params=params_fit,
        datasets=datasets,
        Tt=KRYPTON_T_t,
        pt=KRYPTON_P_t,
        compute_thermo_props=compute_thermo_props,
        St_REFPROP=KRYPTON_REFERENCE_ENTROPY,
        Ht_REFPROP=KRYPTON_REFERENCE_ENTHALPY,
        ncols=3,              # change to 2/4 if you like
        figsize=(14, 10)
    )
    # --- pretty print parameters ---
    formatted = ", ".join(f"{p:.2f}" for p in params_fit)
    print("Fitted parameters:")
    print(formatted)
    # --- pretty print parameters with names ---
    param_names = [
        r"$c_1$ / MPa",
        r"$c_2$ / MPa",
        r"$c_3$ / MPa",
        r"$\theta_{D,0}$ / K",
        r"$\gamma_{D,0}$",
        r"$q_D$",
        r"$b_1$",
        r"$b_2$",
        r"$b_3$",
        r"$S_m(g,T_t,p_t)$ / (J mol$^{-1}$ K$^{-1}$)"
    ]

    print("\n=== Optimised Parameters ===")
    for name, val in zip(param_names, params_fit[:len(param_names)]):
        print(f"{name:25s}: {val:12.4f}")
    print("Final cost:", fval)
# =========================
# ======== STAGE B ========
# Fit thermal + phase-line parameters with full cost
# Keeps elastic (v00, a1, a2, a3) fixed by default
# =========================

# Indices (consistent with your unpacking):
# 0: v00
# 1-3: a1, a2, a3
# 9-15:  Th[0..5]  -> we'll free Th[0] (index 9)
# 15-21: g[0..5]   -> we'll free g[0]  (index 15)
# 21-27: q[0..5]   -> we'll free q[0]  (index 21)
# 27: aa, 28: bb, 29: cc, 30: S* (solid entropy at triple)


_STAGEB_DEFAULT_FREE = dict(
    T0=True,   # Th[0] @ index 9
    g0=True,   # g[0] @ index 15
    q0=True,   # q[0] @ index 21
    aa=True,   # index 27
    bb=True,   # index 28
    cc=True,   # index 29
    Sstar=True  # index 30 (lets triple-point H,S align)
)


def _stageB_bounds(params_start, lower, upper,
                   free_flags=_STAGEB_DEFAULT_FREE,
                   lock_elastic=True, relax_elastic_eps=0.0):
    """
    Build Stage-B bounds.
    - lock_elastic=True → clamp [0..3] to params_start values.
      (Set relax_elastic_eps>0 to allow tiny wiggle if you like.)
    - free_flags controls which thermal indices are free.
    All others are clamped at their starting values.
    """
    params_start = np.asarray(params_start, float)
    LB = list(lower)
    UB = list(upper)

    # start from all fixed at their current value:
    B = [(float(v), float(v)) for v in params_start]

    # optionally unlock elastic (0..3)
    if lock_elastic:
        if relax_elastic_eps > 0:
            for i in (0, 1, 2, 3):
                v = params_start[i]
                B[i] = (v - relax_elastic_eps*abs(v),
                        v + relax_elastic_eps*abs(v))
        # else keep as clamped
    else:
        for i in (0, 1, 2, 3):
            B[i] = (float(LB[i]), float(UB[i]))

    # Free thermal/phase-line knobs
    idx_map = dict(T0=9, g0=15, q0=21, aa=27, bb=28, cc=29, Sstar=30)
    for key, do_free in free_flags.items():
        if do_free:
            i = idx_map[key]
            B[i] = (float(LB[i]), float(UB[i]))

    return B


def _make_outfun_full(*datasets):
    # reuse your per-iteration logging with full cost
    return _make_outfun(*datasets)


def run_stage_B(params_from_stageA, datasets,
                lower_bound, upper_bound,
                free_flags=_STAGEB_DEFAULT_FREE,
                lock_elastic=True,
                maxiter=500, ftol=1e-9, gtol=1e-3):
    print("=== Stage B: Fitting thermal + phase-line parameters (full cost) ===")

    boundsB = _stageB_bounds(
        params_start=params_from_stageA,
        lower=lower_bound, upper=upper_bound,
        free_flags=free_flags,
        lock_elastic=lock_elastic,
        relax_elastic_eps=0.0  # set small (e.g., 1e-3) to allow tiny movement
    )

    # show what is free
    print("[Stage B] Free parameters:", [
          k for k, v in free_flags.items() if v])
    if lock_elastic:
        print("[Stage B] Elastic [v00,a1,a2,a3] locked to Stage-A values")

    # sanity check initial full cost
    tot0, dev0 = combined_cost_function(params_from_stageA, *datasets)
    print(f"[Stage B] initial total = {tot0:.6g}")
    for k, v in dev0.items():
        print(f"  {k}: {v:.6g}")

    cb = _make_outfun_full(*datasets)

    res = minimize(
        fun=_cost_only,                        # uses combined_cost_function
        x0=np.asarray(params_from_stageA, float),
        args=datasets,
        method="L-BFGS-B",
        bounds=boundsB,
        callback=cb,
        options=dict(
            disp=True,
            maxiter=maxiter,
            ftol=ftol,
            gtol=gtol,
            maxls=40,
        )
    )

    print("\n[Stage B] Optimization status:", res.message)
    print("[Stage B] Final cost:", res.fun)

    # Pretty print the parameters we changed
    idx_map = dict(T0=9, g0=15, q0=21, aa=27, bb=28, cc=29, Sstar=30)
    print("\n[Stage B] Fitted thermal/phase-line parameters:")
    for k, i in idx_map.items():
        if free_flags.get(k, False):
            print(f"  {k:6s} (idx {i:2d}): {res.x[i]: .6g}")

    # Save params
    from datetime import datetime
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    np.save(f"stageB_params_{stamp}.npy", res.x)
    print(f"[Stage B] Saved params -> stageB_params_{stamp}.npy")

    return res.x, res.fun


def stage_B(params_stageA):
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)

    paramsB, fB = run_stage_B(
        params_from_stageA=params_stageA,
        datasets=datasets,
        lower_bound=LOWER_BOUND,
        upper_bound=UPPER_BOUND,
        free_flags=_STAGEB_DEFAULT_FREE,  # change here if you want to freeze any
        lock_elastic=True,                # keep elastic from Stage A
        maxiter=500, ftol=1e-9, gtol=1e-3
    )

    # Plot overlays to inspect p_sub, p_melt, cp, alpha, etc.
    plot_all_overlays_grid(
        params=paramsB,
        datasets=datasets,
        Tt=KRYPTON_T_t,
        pt=KRYPTON_P_t,
        compute_thermo_props=compute_thermo_props,
        St_REFPROP=KRYPTON_REFERENCE_ENTROPY,
        Ht_REFPROP=KRYPTON_REFERENCE_ENTHALPY,
        ncols=3,
        figsize=(14, 10)
    )

    formatted = ", ".join(f"{p:.6g}" for p in paramsB)
    print("\n[Stage B] Full parameter vector:")
    print(formatted)
    print("[Stage B] Final cost:", fB)


# 1) Stage A
paramsA, fA = run_stage_A(PARAMS_INIT, extract_datasets(load_all_gas_data('krypton', False)),
                          LOWER_BOUND, UPPER_BOUND)

# 2) Stage B (uses Stage A params as the start point)
stage_B(paramsA)
