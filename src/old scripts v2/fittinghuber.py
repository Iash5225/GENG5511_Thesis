from scipy.optimize import minimize_scalar
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


# --- scale so a1..a3 are O(1) in the optimizer ---
S = np.ones_like(PARAMS_INIT, float)
S[1:4] = 1e-3   # a1,a2,a3 ~ 1e4 → 1e1


def to_params(z): return np.asarray(z, float) * S      # z → real params
def from_params(p): return np.asarray(p, float) / S      # real → z


def nudge_inside_vec(v, bounds, eps=1e-4):
    v = np.array(v, float)
    for i, (lb, ub) in enumerate(bounds):
        if lb < ub:
            v[i] = min(max(v[i], lb+eps), ub-eps)
        else:
            v[i] = lb  # fixed
    return v


def safe_psub(T):
    """Sublimation pressure with complex-safety; returns NaN for non-physical points."""
    T = np.asarray(T, float)
    p = sublimation_pressure_equation(
        T, KRYPTON_E_1_SUB, KRYPTON_E_2_SUB, KRYPTON_E_3_SUB, KRYPTON_T_t, KRYPTON_P_t
    )
    if np.iscomplexobj(p):
        p = np.where(np.abs(p.imag) < 1e-12, p.real, np.nan)
    return np.asarray(p, float)



def subline_residuals(params, datasets, Tt=KRYPTON_T_t):
    (T_sub, _p_dummy, Vm_sub, *_) = datasets
    T_cap = min(110.0, Tt - 2.0)
    T = np.asarray(T_sub, float)
    p = safe_psub(T)   # recompute!
    m = (np.isfinite(T) & np.isfinite(p) & np.isfinite(Vm_sub) &
         (T >= 8.0) & (T <= T_cap) & (p > 0.0))
    Vhat = _safe_props_vector_cached(
        T[m], p[m], params, idx=IDX["Vm"], cache={})
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


def slope_exp_from_data(T, V, Tt):
    m = (T >= 8.0) & (T <= min(110.0, Tt-2.0))
    return np.polyfit(T[m], V[m], 1)[0]


def slope_penalty_for(p27, base, datasets, Tt):
    p = base.copy()
    p[27] = p27
    Tgrid = np.linspace(8.0, min(110.0, Tt-2.0), 120)
    Vmod = vm_along_subline(Tgrid, p)
    ok = np.isfinite(Vmod)
    if np.sum(ok) < 6:
        return 1e9
    k_mod = np.polyfit(Tgrid[ok], Vmod[ok], 1)[0]
    (T_sub, _, V_sub, *_) = datasets
    k_exp = slope_exp_from_data(np.asarray(
        T_sub, float), np.asarray(V_sub, float), Tt)
    return (k_mod - k_exp)**2


def stageA_huber(datasets, x0, bounds, Tt, teacher=None,
                 weights=dict(SUB=1.0, MELT=0.0, ANCH=1.0,
                              SLOPE=0.7, TEACH=0.8, MONO=3.0, CONV=0.2),
                 delta_sub=1.345, delta_melt=1.345, n_outer=4):
    x = np.asarray(x0, float)  # real param space

    # initial residuals/weights (real space)
    r_sub, m_sub = subline_residuals(x, datasets)
    w_sub, s_sub = huber_weights(r_sub, delta=delta_sub)
    r_melt, m_melt = melt_residuals(x, datasets, Tt)
    w_melt, s_melt = huber_weights(r_melt, delta=delta_melt)

    for k in range(n_outer):
        # --- build scaled problem for this outer iteration ---
        bounds_z = [(lb/si, ub/si) for (lb, ub), si in zip(bounds, S)]
        z0 = from_params(x)
        z0 = nudge_inside_vec(z0, bounds_z)  # avoid "X0 at bounds"

        # optimize in z; map to real params inside objective
        res = minimize(
            fun=lambda z, *a: float(
                combined_cost_vm_huber(
                    to_params(z), *a, Tt=Tt, teacher=teacher, weights=weights,
                    w_sub=w_sub, msub_mask=m_sub, w_melt=w_melt, mmelt_mask=m_melt
                )
            ),
            x0=z0, args=datasets, method="L-BFGS-B", bounds=bounds_z,
            options=dict(disp=True, maxiter=MAX_ITERATIONS,
                         ftol=FUNCTION_TOL, gtol=GRADIENT_TOL, maxls=60),
        )
        x = to_params(res.x)  # back to real space

        # reweight for next IRLS outer loop (still in real space)
        r_sub, _ = subline_residuals(x, datasets)
        w_sub, s_sub = huber_weights(r_sub, delta=delta_sub, scale=s_sub)
        r_melt, _ = melt_residuals(x, datasets, Tt)
        w_melt, s_melt = huber_weights(r_melt, delta=delta_melt, scale=s_melt)

    return x


print("\n=== Stage A (Huber IRLS, subline first) ===")
krypton_data = load_all_gas_data('krypton', read_from_excel=False)
datasets = extract_datasets(krypton_data)
print(f"Datasets Extracted")
teacher = make_subline_teacher(datasets, Tt=KRYPTON_T_t)
print(f"Teacher Curve Prepared")
# Bounds: baseline around 20 K median + free a1..a3 + thermal knobs
v00_star = target_v00(datasets)
print(f"Target v00*: {v00_star:.4f} cm3/mol")

boundsA = make_stageA_bounds(PARAMS_INIT, LOWER_BOUND, UPPER_BOUND,
                             free_elastic=(0, 1, 2, 3), v00_idx=0,
                             v00_range=(v00_star-0.3, v00_star+0.3),
                             elastic_range=(-1e4, 1e4))

PARAMS_INIT[27] = -0.8202
boundsA = list(boundsA)
boundsA[27] = (max(-6.0, -0.8202 - 0.5), min(6.0, -0.8202 + 0.5))
boundsA = tuple(boundsA)

print(f"Stage A bounds: {boundsA}")
# 1-D pre-tune of slope knob (p27)
# res1d = minimize_scalar(slope_penalty_for, bounds=(-6, 6), method='bounded',
#                         args=(PARAMS_INIT, datasets, KRYPTON_T_t))
# PARAMS_INIT[27] = float(res1d.x)
# print(f"1-D slope knob pre-tune: p27 = {PARAMS_INIT[27]:.4f}, penalty = {res1d.fun:.4f}")
# Start vector
x0 = prefit_v00(PARAMS_INIT, datasets)
x0[0] = v00_star

# Subline-only Huber fit
w_stageA = dict(SUB=1.0, MELT=0.0, ANCH=1.0, SLOPE=0.7,
                TEACH=0.8, MONO=3.0, CONV=0.2)
paramsA = stageA_huber(datasets, x0, boundsA, Tt=KRYPTON_T_t,
                       teacher=teacher, weights=w_stageA,
                       delta_sub=1.345, delta_melt=1.345, n_outer=N_OUTER)

quick_plot_vm_only(paramsA, datasets, Tt=KRYPTON_T_t)
diagnose_subline(paramsA, datasets, Tt=KRYPTON_T_t)
print(f"Stage A (Huber) params: {paramsA}")
