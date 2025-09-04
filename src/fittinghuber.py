from math import isfinite
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from computethermoprops import *
from constants import *
from fitting_helper import *
from plot_eos import plot_all_overlays_grid
from thermal_script import melting_pressure_equation, sublimation_pressure_equation
from scipy.interpolate import UnivariateSpline

# --- robust scale & Huber weights (standard Huber model) ---
# Constants
BIG = 1e4
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)

def _mad(x):
    x = np.asarray(x, float)
    m = np.nanmedian(x)
    return 1.4826 * (np.nanmedian(np.abs(x - m)) + 1e-12)


def huber_weights(resid, delta=1.345, scale=None):
    """
    Standard Huber weights on *standardized* residuals r/s.
    delta=1.345 -> ~95% efficiency for Gaussian noise.
    Returns weights in [0,1] (1 = inlier), and the scale used.
    """
    r = np.asarray(resid, float)
    s = _mad(r) if scale is None else float(scale)
    z = r / (s + 1e-12)
    a = np.abs(z)
    w = np.ones_like(a)
    mask = a > delta
    # Huber psi(z) / z  ->  delta/|z|
    w[mask] = delta / a[mask]
    return w, s


def vm_along_subline(Tseq, params, cache=None, idx=IDX["Vm"]):
    Tseq = np.asarray(Tseq, float)
    pseq = safe_psub(Tseq)
    out = np.full_like(Tseq, np.nan, dtype=float)
    last_vm = None
    for i, (T, p) in enumerate(zip(Tseq, pseq)):
        if not (np.isfinite(p) and p > 0.0):
            continue
        key = (float(T), float(p))
        props = None if cache is None else cache.get(key)
        if props is None:
            try:
                # if vm0 not supported it's okay
                props = compute_thermo_props(T, p, params, vm0=last_vm)
            except TypeError:
                props = compute_thermo_props(T, p, params)
            if cache is not None:
                cache[key] = props
        if props is not None and np.all(np.isfinite(props)):
            vm = props[idx]
            if np.isfinite(vm):
                out[i] = vm
                last_vm = vm
    return out


def subline_residuals(params, datasets):
    # (subline Vm triple is first)
    (T_Vm_sub, p_Vm_sub, Vm_sub, *_) = datasets
    m = np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub) & np.isfinite(Vm_sub)
    Vhat = _safe_props_vector_cached(
        T_Vm_sub[m], p_Vm_sub[m], params, idx=IDX["Vm"], cache={})
    return Vm_sub[m] - Vhat, m


def melt_residuals(params, datasets, Tt):
    # melt Vm triple is second
    (T_Vm_sub, p_Vm_sub, Vm_sub,
     T_Vm_melt, p_Vm_melt, Vm_melt, *_) = datasets

    m = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt) &
         np.isfinite(Vm_melt) & (T_Vm_melt >= Tt))
    Vhat = _safe_props_vector_cached(
        T_Vm_melt[m], p_Vm_melt[m], params, idx=IDX["Vm"], cache={})
    return Vm_melt[m] - Vhat, m


def _safe_props_vector_cached(T_arr, p_arr, params, idx, cache):
    T_arr = np.asarray(T_arr, float)
    p_arr = np.asarray(p_arr)
    if np.iscomplexobj(p_arr):
        p_arr = np.where(np.abs(p_arr.imag) < 1e-12, p_arr.real, np.nan)
    p_arr = p_arr.astype(float, copy=False)

    out = np.full(T_arr.shape, np.nan, dtype=float)
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        if not (np.isfinite(T) and np.isfinite(p) and p > 0.0):
            continue
        key = (float(T), float(p))
        props = cache.get(key)
        if props is None:
            try:
                props = compute_thermo_props(float(T), float(p), params)
            except Exception:
                props = None
            cache[key] = props
        if props is not None and np.all(np.isfinite(props)):
            out[i] = props[idx]
    return out


def combined_cost_vm_huber(params, *datasets, Tt=None, teacher=None, weights=None,
                           w_sub=None, msub_mask=None, w_melt=None, mmelt_mask=None):
    """
    Stage A cost with Huber-IRLS weighting for data fits.
    w_sub / w_melt are per-point weights (fixed for this call).
    """
    (T_sub, p_sub, Vm_sub,
     T_melt, p_melt, Vm_melt, *_) = datasets

    Tt_val = float(Tt if Tt is not None else KRYPTON_T_t)
    Tmax_sub = float(min(110.0, Tt_val - 2.0))

    # default term weights
    W = dict(SUB=1.0, MELT=0.0, ANCH=1.0, SLOPE=0.7,
             TEACH=0.8, MONO=3.0, CONV=0.2)
    if weights:
        W.update(weights)

    cache = {}

    # ---- Huber-weighted data misfit (sub) ----
    Vm_sub_dev = BIG
    if msub_mask is not None and np.any(msub_mask):
        Vhat = _safe_props_vector_cached(
            T_sub[msub_mask], p_sub[msub_mask], params, idx=IDX["Vm"], cache=cache)
        r = Vm_sub[msub_mask] - Vhat
        Vm_sub_dev = float(np.sqrt(np.mean((w_sub * r) * r)))

    # ---- Huber-weighted data misfit (melt) ----
    Vm_melt_dev = BIG
    if mmelt_mask is not None and np.any(mmelt_mask):
        Vhat = _safe_props_vector_cached(
            T_melt[mmelt_mask], p_melt[mmelt_mask], params, idx=IDX["Vm"], cache=cache)
        r = Vm_melt[mmelt_mask] - Vhat
        Vm_melt_dev = float(np.sqrt(np.mean((w_melt * r) * r)))

    # ---- helper terms on sub-window ----
    Vm_anchor = 0.0
    slope_pen = 0.0
    teacher_pen = 0.0
    mono_pen = 0.0
    conv_pen = 0.0

    # window below Tt
    win = (np.isfinite(T_sub) & np.isfinite(Vm_sub)
           & (T_sub >= 8.0) & (T_sub <= Tmax_sub))
    if np.any(win):
        Tsub = np.asarray(T_sub[win], float)
        Vexp = np.asarray(Vm_sub[win],   float)

        # anchors: medians at ~20 K and ~Tmax_sub
        def _median_near(Tc, half=5.0):
            m = np.abs(Tsub - Tc) <= half
            return float(np.nanmedian(Vexp[m])) if np.any(m) else float(np.nanmedian(Vexp))
        T_lo, T_hi = 20.0, Tmax_sub
        V_lo_exp = _median_near(T_lo)
        V_hi_exp = _median_near(T_hi)
        V_lo_mod = _safe_props_vector_cached(
            [T_lo], [safe_psub(T_lo)], params, idx=IDX["Vm"], cache=cache)[0]
        V_hi_mod = _safe_props_vector_cached(
            [T_hi], [safe_psub(T_hi)], params, idx=IDX["Vm"], cache=cache)[0]
        if np.isfinite(V_lo_mod) and np.isfinite(V_hi_mod):
            Vm_anchor = (V_lo_mod - V_lo_exp)**2 + (V_hi_mod - V_hi_exp)**2

        # slope / shape on sorted subline with warm-start
        Tgrid = np.sort(np.unique(np.clip(Tsub, 8.0, Tmax_sub)))
        Vmod = vm_along_subline(Tgrid, params, cache=cache)
        ok = np.isfinite(Vmod) & (26.0 < Vmod) & (Vmod < 31.5)
        if np.sum(ok) >= 4:
            k_exp = np.polyfit(Tsub, Vexp, 1)[0]
            k_mod = np.polyfit(Tgrid[ok], Vmod[ok], 1)[0]
            slope_pen = (k_mod - k_exp)**2

            # monotonicity & convexity (gentle)
            dT = np.diff(Tgrid[ok])
            dV = np.diff(Vmod[ok]) / dT
            mono_pen = float(np.mean(np.clip(-dV, 0.0, None)**2))
            if dV.size >= 2:
                d2V = np.diff(Vmod[ok], 2) / (dT[:-1]*dT[1:])
                conv_pen = float(np.mean(np.clip(-d2V, 0.0, None)**2))

        # optional: teacher curve
        if teacher is not None:
            Tg = teacher["Tg"]
            Vt = teacher["Vg_target"]
            Vm_g = vm_along_subline(Tg, params, cache=cache)
            mask = np.isfinite(Vm_g)
            if np.any(mask):
                teacher_pen = float(np.mean((Vm_g[mask] - Vt[mask])**2))

    total = (W["SUB"]*Vm_sub_dev +
             W["MELT"]*Vm_melt_dev +
             W["ANCH"]*Vm_anchor +
             W["SLOPE"]*slope_pen +
             W["TEACH"]*teacher_pen +
             W["MONO"]*mono_pen +
             W["CONV"]*conv_pen)

    return total


def stageA_huber(datasets, x0, bounds, Tt, teacher=None,
                 weights=dict(SUB=1.0, MELT=0.0, ANCH=1.0,
                              SLOPE=0.7, TEACH=0.8, MONO=3.0, CONV=0.2),
                 delta_sub=1.345, delta_melt=1.345, n_outer=4):
    """
    Run Huber IRLS from the beginning. Start with SUB only; turn MELT on later.
    """
    x = np.asarray(x0, float)

    # initial residuals & weights
    r_sub, m_sub = subline_residuals(x, datasets)
    w_sub, s_sub = huber_weights(r_sub, delta=delta_sub)

    r_melt, m_melt = melt_residuals(x, datasets, Tt)
    w_melt, s_melt = huber_weights(r_melt, delta=delta_melt)

    for k in range(n_outer):
        res = minimize(
            fun=lambda p, *a: float(
                combined_cost_vm_huber(p, *a, Tt=Tt, teacher=teacher, weights=weights,
                                       w_sub=w_sub, msub_mask=m_sub, w_melt=w_melt, mmelt_mask=m_melt)
            ),
            x0=x, args=datasets, method="L-BFGS-B", bounds=bounds,
            options=dict(disp=True, maxiter=200, ftol=FUNCTION_TOL,
                         gtol=GRADIENT_TOL, maxls=60, eps=1e-3),
        )
        x = res.x

        # reweight with updated residuals (classic Huber IRLS)
        r_sub, _ = subline_residuals(x, datasets)
        w_sub, s_sub = huber_weights(r_sub, delta=delta_sub, scale=s_sub)

        r_melt, _ = melt_residuals(x, datasets, Tt)
        w_melt, s_melt = huber_weights(r_melt, delta=delta_melt, scale=s_melt)

    return x


print("\n=== Stage A (Huber IRLS, subline first) ===")
krypton_data = load_all_gas_data('krypton', read_from_excel=False)
datasets = extract_datasets(krypton_data)

# Optional teacher (smoothing spline); or pass teacher=None
teacher = make_subline_teacher(datasets, Tt=KRYPTON_T_t)

# Bounds: tighten v00 around 20 K median (baseline), free a1..a3 & thermal knob(s)
v00_star = target_v00(datasets)         # your helper from earlier
boundsA = make_stageA_bounds(PARAMS_INIT, LOWER_BOUND, UPPER_BOUND,
                             free_elastic=(0, 1, 2, 3), v00_idx=0,
                             v00_range=(v00_star-0.3, v00_star+0.3),
                             elastic_range=(-1e4, 1e4))
B = list(boundsA)
B[27] = (-6, 6)
B[15] = (-2, 2)
boundsA = B

# Start vector (prefit v00)
x0 = prefit_v00(PARAMS_INIT, datasets)
x0[0] = v00_star

# Subline-only Huber fit (MELT=0.0)
w_stageA = dict(SUB=1.0, MELT=0.0, ANCH=1.0, SLOPE=0.7,
                TEACH=0.8, MONO=3.0, CONV=0.2)
paramsA = stageA_huber(datasets, x0, boundsA, Tt=KRYPTON_T_t, teacher=teacher, weights=w_stageA,
                       delta_sub=1.345, delta_melt=1.345, n_outer=4)

# (Optional) short refine with melt on
# print("\n=== Short refine with MELT on (still Huber) ===")
# w_refine = dict(SUB=1.0, MELT=1.0, ANCH=0.8, SLOPE=0.4,
#                 TEACH=0.5, MONO=2.0, CONV=0.2)
# paramsA2 = stageA_huber(datasets, paramsA, boundsA, Tt=KRYPTON_T_t, teacher=teacher, weights=w_refine,
#                         delta_sub=1.345, delta_melt=1.345, n_outer=2)

# Plot & diagnostics
quick_plot_vm_only(paramsA, datasets, Tt=KRYPTON_T_t)
diagnose_subline(paramsA, datasets, Tt=KRYPTON_T_t)
print(f"Stage A (Huber) params: {paramsA}")