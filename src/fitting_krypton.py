from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read import load_all_gas_data
from constants import *
from thermopropsv2 import compute_thermo_props
from fitting_helper import _safe_props_vector, rms_percent, rms, extract_datasets, safe_psub, psub_curve, pmelt_curve, plot_deviation, plot_init
from plot_eos import plot_all_overlays_grid
import traceback
from deviation_recorder import DeviationRecorder, Metric
import math
import os
# Constants
BIG = 1e4
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)
GLOBAL_RECORDER = DeviationRecorder()

NAME_TO_METRIC = {
    "Vm_sub":       Metric.VM_SUB,
    "Vm_melt":      Metric.VM_MELT,
    "Vm_highp":     Metric.VM_HIGHP,
    "cp_sub":       Metric.CP_SUB,
    "alpha_sub":    Metric.ALPHA_SUB,
    "BetaT_sub":    Metric.BETAT_SUB,
    "BetaS_sub":    Metric.BETAS_SUB,
    "H_solid_sub":  Metric.H_SOLID_SUB,
    "H_solid_melt": Metric.H_SOLID_MELT,
    "Gamma_T":      Metric.GAMMA_T,
}

def _cost_only(params, *datasets):
    try:
        total, devs = combined_cost_function(params, *datasets)
        GLOBAL_RECORDER.record(Metric.TOTAL, total)
        for name, val in devs.items():
            if name in NAME_TO_METRIC:
                GLOBAL_RECORDER.record(NAME_TO_METRIC[name], float(val))
        return float(total) if np.isfinite(total) else BIG
    except Exception as e:
        print("Error in cost function evaluation.", repr(e), flush=True)
        traceback.print_exc()
        return BIG


def _triple_offsets(params, Tt, pt, compute_thermo_props, St_REFPROP, Ht_REFPROP):
    """
    Minimal helper: returns (dH_triple_kJ_per_mol, dS_triple) using the same
    convention as plotting. This is the only place S* (params[30]) enters.
    """
    tp = compute_thermo_props(float(Tt), float(pt), params)
    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        # if your compute_thermo_props ever returns meta
        meta = tp[1]
        dH = float(meta.get("deltaH_triple"))
        dS = float(meta.get("deltaS_triple"))
    else:
        props = np.asarray(tp, float)
        Sstar = float(params[30])
        dS = Sstar - float(St_REFPROP)
        # match plot_eos: use G (index 11), not H
        Ht_fitted = float(props[11]) + float(Tt) * Sstar   # G + T*S*
        dH = Ht_fitted - float(Ht_REFPROP)                 # J/mol
    return dH/1000.0, dS  # kJ/mol, dimensionless
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
    # print("Thermo props:", tp)
    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        triple_props, triple_meta = tp
        # pulled directly from compute_thermo_props
        #deltaH_triple = float(triple_meta.get("deltaH_triple"))
        deltaS_triple = float(triple_meta.get("deltaS_triple"))
    else:
        # backward compatibility (old API returned only the 12-length props array)
        triple_props = np.asarray(tp, float)
        # S* at triple is params[30]; REFPROP values are global constants
        deltaS_triple = params[30] - St_REFPROP
        Ht_fitted = triple_props[11] + Tt * params[30]
        deltaH_triple = Ht_fitted - Ht_REFPROP
        deltaH_triple_kJ = deltaH_triple / 1000        # kJ/mol

    # ====== BLOCKS ======
    Vm_sub_dev = BIG
    Vm_melt_dev = BIG
    Vm_highp_dev = BIG
    cp_sub_dev = BIG
    alpha_sub_dev = BIG
    BetaT_sub_dev = BIG
    BetaS_sub_dev = BIG
    H_solid_sub_dev = BIG
    H_solid_melt_dev = BIG
    Gamma_T6_dev = BIG

    dHtr_kJ, dStr = _triple_offsets(
        params, Tt, pt, compute_thermo_props, St_REFPROP, Ht_REFPROP)

    m = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub)
         & np.isfinite(Vm_sub) & (T_Vm_sub <= Tt))
    if np.any(m):
        model = _safe_props_vector(T_Vm_sub[m], p_Vm_sub[m], params, idx=0)
        Vm_sub_dev = rms_percent(Vm_sub[m], model)   # cm^3/mol
    # Vm (melting): T >= Tt
    m = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt)
         & np.isfinite(Vm_melt) & (T_Vm_melt >= Tt))
    if np.any(m):
        model = _safe_props_vector(T_Vm_melt[m], p_Vm_melt[m], params, idx=0)
        Vm_melt_dev = rms_percent(Vm_melt[m], model)  # cm^3/mol
    # Vm (high-p): no phase constraint (already experimental p)
    m = (np.isfinite(T_Vm_highp) & np.isfinite(
        p_Vm_highp) & np.isfinite(Vm_highp))

    if np.any(m):
        model = _safe_props_vector(T_Vm_highp[m], p_Vm_highp[m], params, idx=0)
        Vm_highp_dev = rms_percent(Vm_highp[m], model)
    # cp (sublimation) with split weighting (low vs high T)
    m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub)
         & np.isfinite(cp_sub) & (T_cp_sub <= Tt))

    if np.any(m):
        Tm = np.asarray(T_cp_sub)[m]
        pm = np.asarray(p_cp_sub)[m]
        cp_exp = np.asarray(cp_sub)[m]
        model = _safe_props_vector(Tm, pm, params, idx=4)
        # cp_sub_dev = rms_percent(cp_exp, model)
        # percent deviations
        dev = (cp_exp - model) / cp_exp
        # weights depending on T threshold
        weights = np.where(Tm < CP_TEMP_THRESHOLD_K,
                           CP_WEIGHT_BELOW, CP_WEIGHT_ABOVE)
        cp_sub_dev = float(np.sqrt(np.mean((weights * dev) ** 2)))

    # alpha (sublimation)
    m = (np.isfinite(T_alpha_sub) & np.isfinite(p_alpha_sub)
         & np.isfinite(alpha_sub) & (T_alpha_sub <= Tt))

    if np.any(m):
        model = _safe_props_vector(
            T_alpha_sub[m], p_alpha_sub[m], params, idx=3)
        alpha_sub_dev = rms_percent(alpha_sub[m], model)
    # BetaT (sublimation)
    m = (np.isfinite(T_BetaT_sub) & np.isfinite(p_BetaT_sub)
         & np.isfinite(BetaT_sub) & (T_BetaT_sub <= Tt))

    if np.any(m):
        model = _safe_props_vector(
            T_BetaT_sub[m], p_BetaT_sub[m], params, idx=1)
        BetaT_sub_dev = rms_percent(BetaT_sub[m], model)
    # BetaS (sublimation)
    m = (np.isfinite(T_BetaS_sub) & np.isfinite(p_BetaS_sub)
         & np.isfinite(BetaS_sub) & (T_BetaS_sub <= Tt))

    if np.any(m):
        model = _safe_props_vector(
            T_BetaS_sub[m], p_BetaS_sub[m], params, idx=2)
        BetaS_sub_dev = rms_percent(BetaS_sub[m], model)


    # # Enthalpy of sublimation: kJ/mol (no unit juggling)
    # m = (np.isfinite(T_H_sub) & np.isfinite(p_H_sub) &
    #      np.isfinite(delta_H_sub) & np.isfinite(H_fluid_sub) & (T_H_sub <= Tt))

    # if np.any(m):
    #     H_solid_sub = H_fluid_sub[m] - delta_H_sub[m] * 1.0  # both kJ/mol
    #     modelH = _safe_props_vector(
    #         T_H_sub[m], p_H_sub[m], params, idx=10) / 1000 - deltaH_triple_kJ
    #     H_solid_sub_dev = rms_percent(H_solid_sub, modelH)
    # # Enthalpy of melting: kJ/mol
    # m = (np.isfinite(T_H_melt) & np.isfinite(p_H_melt) &
    #      np.isfinite(delta_H_melt) & np.isfinite(H_fluid_melt) & (T_H_melt >= Tt))

    # if np.any(m):
    #     H_solid_melt = H_fluid_melt[m] - delta_H_melt[m] * 1.0
    #     modelH = _safe_props_vector(
    #         T_H_melt[m], p_H_melt[m], params, idx=10)/1000 - deltaH_triple_kJ
    #     H_solid_melt_dev = rms_percent(H_solid_melt, modelH)
    # Gamma-T smoothness penalty along a subl. grid (T<=Tt)
        # Enthalpy of sublimation: compare H_solid,exp vs (H_model - dHtr)
    m = (np.isfinite(T_H_sub) & np.isfinite(p_H_sub) &
         np.isfinite(delta_H_sub) & np.isfinite(H_fluid_sub) & (T_H_sub <= Tt))
    if np.any(m):
        H_solid_sub = H_fluid_sub[m] - delta_H_sub[m] * 1.0          # kJ/mol
        modelH = _safe_props_vector(
            T_H_sub[m], p_H_sub[m], params, idx=10)/1000.0 - dHtr_kJ
        H_solid_sub_dev = rms_percent(H_solid_sub, modelH)

    # Enthalpy of melting
    m = (np.isfinite(T_H_melt) & np.isfinite(p_H_melt) &
         np.isfinite(delta_H_melt) & np.isfinite(H_fluid_melt) & (T_H_melt >= Tt))
    if np.any(m):
        H_solid_melt = H_fluid_melt[m] - delta_H_melt[m] * 1.0       # kJ/mol
        modelH = _safe_props_vector(
            T_H_melt[m], p_H_melt[m], params, idx=10)/1000.0 - dHtr_kJ
        H_solid_melt_dev = rms_percent(H_solid_melt, modelH)
    T6 = np.array([0.0001] + list(range(2, math.ceil(Tt), 2)) + [Tt])
    T6 = T6[T6 <= Tt]
    p6 = np.full_like(T6, np.nan, dtype=float)
    G6 = np.full_like(T6, np.nan, dtype=float)
    V6 = np.full_like(T6, np.nan, dtype=float)
    for i, T in enumerate(T6):
        # p6[i] = psub(T, pt, Tt)
        p6[i] = safe_psub(T)
        pr = compute_thermo_props(T, p6[i], params)
        if np.all(np.isfinite(pr)):
            G6[i] = pr[6]   # Gruneisen
            V6[i] = pr[0]
    # finite pairs only
    mk = np.isfinite(G6) & np.isfinite(V6)

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
    
    # --- light triple-point anchor so S* is identifiable even with sparse H data ---
        # --- use the helper so S* affects the cost ---
    dHtr_kJ, dStr = _triple_offsets(
        params, Tt, pt, compute_thermo_props, St_REFPROP, Ht_REFPROP)
    s_ref = max(1.0, abs(St_REFPROP))
    h_ref_kJ = max(1.0, abs(Ht_REFPROP)/1000.0)
    triple_pen = 0.05 * (dStr / s_ref)**2 + 0.05 * (dHtr_kJ / h_ref_kJ)**2
    # ====== TOTAL COST ======    
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
        Gamma_T6_dev * W_GAMMA_T + 
        triple_pen
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
        "Gamma_T": Gamma_T6_dev,
        "dHtr_kJ": dHtr_kJ,
        "dStr": dStr,
    }
    return total_deviation, deviations

def main():
    # === 4. Package for scipy.optimize.minimize ===
    bounds = [(lo, hi) for lo, hi in zip(LOWER_BOUND, UPPER_BOUND)]
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)
    res = minimize(
        fun=lambda x, *a: _cost_only(x, *a),   
        x0=PARAMS_INIT,
        args=datasets,                         
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                     ftol=FUNCTION_TOL, gtol=GRADIENT_TOL, eps=EPS, maxls=MAXLS)
    )
    print("\n Fitting status:", res.message)
    params_fit = res.x
    # params_fit = PARAMS_INIT
    # --- pretty print parameters ---
    formatted = ", ".join(f"{p:.5f}" for p in params_fit)
    print("Fitted parameters:")
    print(formatted)
    plot_all_overlays_grid(params_fit, datasets, Tt=Tt, pt=pt, compute_thermo_props=compute_thermo_props,
                           St_REFPROP=St_REFPROP, Ht_REFPROP=Ht_REFPROP, psub_curve=psub_curve, pmelt_curve=pmelt_curve)
    GLOBAL_RECORDER.plot_history(ncols=5)
    

if __name__ == "__main__":
    plot_deviation(variable="BetaS_sub")
    # main()
    # plot_init()
    # RMS_AAD()
