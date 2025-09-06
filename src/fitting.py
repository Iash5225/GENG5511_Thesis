from math import isfinite
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
# from computethermoprops import *
from constants import *
from thermopropsv2 import compute_thermo_props
from fitting_helper import *
from plot_eos import plot_all_overlays_grid
from thermal_script import melting_pressure_equation, sublimation_pressure_equation
from scipy.interpolate import UnivariateSpline

# Constants
BIG = 1e4
St_REFPROP = KRYPTON_REFERENCE_ENTROPY  # Reference Entropy
Ht_REFPROP = KRYPTON_REFERENCE_ENTHALPY
Tt = KRYPTON_T_t
pt = KRYPTON_P_t
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)

# Auxiliary functions
def psub_curve(T): return sublimation_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_1_SUB,  KRYPTON_E_2_SUB,  KRYPTON_E_3_SUB,  KRYPTON_T_t, KRYPTON_P_t
)

def pmelt_curve(T): return melting_pressure_equation(
    np.asarray(
        T, float), KRYPTON_E_4, KRYPTON_E_5, KRYPTON_E_6, KRYPTON_E_7, KRYPTON_T_t, KRYPTON_P_t
)

# Hidden Functions
def _cost_only(params, *datasets):
    try:
        total, _ = combined_cost_function(params, *datasets)
        return float(total) if np.isfinite(total) else BIG
    except Exception:
        return BIG

def _compute_triple_offsets(params):
    """Return (deltaH_triple, deltaS_triple) using meta if provided by compute_thermo_props."""
    tp = compute_thermo_props(Tt, pt, params)
    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        _, meta = tp
        deltaH_triple = float(meta.get("deltaH_triple"))
        deltaS_triple = float(meta.get("deltaS_triple"))
    else:
        triple_props = np.asarray(tp, float)
        deltaS_triple = params[30] - St_REFPROP
        Ht_fitted = triple_props[IDX["G"]] + Tt * params[30]  # H = G + T*S*
        deltaH_triple = Ht_fitted - Ht_REFPROP
    return deltaH_triple, deltaS_triple
# Combined Cost
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

    m = (np.isfinite(T_Vm_sub) & np.isfinite(p_Vm_sub)
         & np.isfinite(Vm_sub) & (T_Vm_sub <= Tt))
    Vm_sub_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_sub[m], p_Vm_sub[m], params, idx=0)
        Vm_sub_dev = _rmse_abs(Vm_sub[m], model)   # cm^3/mol

    # Vm (melting): T >= Tt  —— use ABS RMSE in cm^3/mol
    m = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt)
         & np.isfinite(Vm_melt) & (T_Vm_melt >= Tt))
    Vm_melt_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_melt[m], p_Vm_melt[m], params, idx=0)
        Vm_melt_dev = _rmse_abs(Vm_melt[m], model)  # cm^3/mol
    # Vm (high-p): no phase constraint (already experimental p)
    m = (np.isfinite(T_Vm_highp) & np.isfinite(
        p_Vm_highp) & np.isfinite(Vm_highp))
    Vm_highp_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_highp[m], p_Vm_highp[m], params, idx=0)
        Vm_highp_dev = rms_percent(Vm_highp[m], model)

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

    # Gamma-T smoothness penalty along a subl. grid (T<=Tt)
    T6 = np.array([0.0001] + list(range(2, 84, 2)) + [83.806])
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
        "Gamma_T": Gamma_T6_dev,
    }
    return total_deviation, deviations

def export_eos_to_excel(params_fit,
                                     out_path="krypton_EOS_export.xlsx",
                                     read_from_excel=False):
    """
    Export experimental datasets and EOS model values to an Excel workbook.
    One sheet per dataset, no residual columns.

    Sheets:
      - Vm_sub, Vm_melt, Vm_highp
      - cp_sub, alpha_sub, kappaT_sub, kappaS_sub
      - H_sub, H_melt
      - parameters, readme
    """
    # 1) Load and unpack datasets
    data = load_all_gas_data('krypton', read_from_excel=read_from_excel)
    datasets = extract_datasets(data)

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

    # 2) Triple-point offsets (for enthalpy alignment)
    deltaH_triple, _ = _compute_triple_offsets(params_fit)

    sheets = {}

    # Vm: sublimation
    if np.size(T_Vm_sub):
        Vm_model = _safe_props_vector(
            T_Vm_sub, p_Vm_sub, params_fit, IDX["Vm"])
        sheets["Vm_sub"] = pd.DataFrame({
            "T [K]": T_Vm_sub, "p [MPa]": p_Vm_sub,
            "Vm_exp [cm^3/mol]": Vm_sub,
            "Vm_model [cm^3/mol]": Vm_model,
        })

    # Vm: melting
    if np.size(T_Vm_melt):
        Vm_model = _safe_props_vector(
            T_Vm_melt, p_Vm_melt, params_fit, IDX["Vm"])
        sheets["Vm_melt"] = pd.DataFrame({
            "T [K]": T_Vm_melt, "p [MPa]": p_Vm_melt,
            "Vm_exp [cm^3/mol]": Vm_melt,
            "Vm_model [cm^3/mol]": Vm_model,
        })

    # Vm: high pressure
    if np.size(T_Vm_highp):
        Vm_model = _safe_props_vector(
            T_Vm_highp, p_Vm_highp, params_fit, IDX["Vm"])
        sheets["Vm_highp"] = pd.DataFrame({
            "T [K]": T_Vm_highp, "p [MPa]": p_Vm_highp,
            "Vm_exp [cm^3/mol]": Vm_highp,
            "Vm_model [cm^3/mol]": Vm_model,
        })

    # cp (sublimation)
    if np.size(T_cp_sub):
        cp_model = _safe_props_vector(
            T_cp_sub, p_cp_sub, params_fit, IDX["cp"])
        sheets["cp_sub"] = pd.DataFrame({
            "T [K]": T_cp_sub, "p [MPa]": p_cp_sub,
            "cp_exp [J/mol/K]": cp_sub,
            "cp_model [J/mol/K]": cp_model,
        })

    # alpha (sublimation)
    if np.size(T_alpha_sub):
        a_model = _safe_props_vector(
            T_alpha_sub, p_alpha_sub, params_fit, IDX["Alpha"])
        sheets["alpha_sub"] = pd.DataFrame({
            "T [K]": T_alpha_sub, "p [MPa]": p_alpha_sub,
            "alpha_exp [1/K]": alpha_sub,
            "alpha_model [1/K]": a_model,
        })

    # kappa_T (from BetaT_sub dataset)
    if np.size(T_BetaT_sub):
        kT_model = _safe_props_vector(
            T_BetaT_sub, p_BetaT_sub, params_fit, IDX["KappaT"])
        sheets["kappaT_sub"] = pd.DataFrame({
            "T [K]": T_BetaT_sub, "p [MPa]": p_BetaT_sub,
            "kappaT_exp [1/MPa]": BetaT_sub,
            "kappaT_model [1/MPa]": kT_model,
        })

    # kappa_S (from BetaS_sub dataset)
    if np.size(T_BetaS_sub):
        kS_model = _safe_props_vector(
            T_BetaS_sub, p_BetaS_sub, params_fit, IDX["KappaS"])
        sheets["kappaS_sub"] = pd.DataFrame({
            "T [K]": T_BetaS_sub, "p [MPa]": p_BetaS_sub,
            "kappaS_exp [1/MPa]": BetaS_sub,
            "kappaS_model [1/MPa]": kS_model,
        })

    # Enthalpy of sublimation: solid = fluid − ΔH  (ensure units are J/mol)
    if np.size(T_H_sub):
        H_solid_exp = np.asarray(H_fluid_sub, float) - \
            np.asarray(delta_H_sub, float)
        H_model = _safe_props_vector(T_H_sub, p_H_sub, params_fit, IDX["H"])
        H_solid_model = H_model - deltaH_triple
        sheets["H_sub"] = pd.DataFrame({
            "T [K]": T_H_sub, "p [MPa]": p_H_sub,
            "H_solid_exp [J/mol]": H_solid_exp,
            "H_solid_model [J/mol]": H_solid_model,
        })

    # Enthalpy of fusion: solid = fluid − ΔH_melt
    if np.size(T_H_melt):
        H_solid_exp = np.asarray(H_fluid_melt, float) - \
            np.asarray(delta_H_melt, float)
        H_model = _safe_props_vector(T_H_melt, p_H_melt, params_fit, IDX["H"])
        H_solid_model = H_model - deltaH_triple
        sheets["H_melt"] = pd.DataFrame({
            "T [K]": T_H_melt, "p [MPa]": p_H_melt,
            "H_solid_exp [J/mol]": H_solid_exp,
            "H_solid_model [J/mol]": H_solid_model,
        })

    # Parameters (first N named, rest raw if longer)
    param_names = [
        r"$c_1$ / MPa", r"$c_2$ / MPa", r"$c_3$ / MPa",
        r"$\theta_{D,0}$ / K", r"$\gamma_{D,0}$", r"$q_D$",
        r"$b_1$", r"$b_2$", r"$b_3$",
        r"$S_m(g,T_t,p_t)$ / (J mol$^{-1}$ K$^{-1}$)"
    ]
    n_show = min(len(param_names), len(params_fit))
    sheets["parameters"] = pd.DataFrame({
        "name": param_names[:n_show],
        "value": [float(v) for v in np.asarray(params_fit)[:n_show]],
    })

    # Readme
    sheets["readme"] = pd.DataFrame({
        "sheet": list(sheets.keys()),
        "note": [
            "All columns are experimental vs EOS model values; no residuals."
            for _ in sheets
        ],
    })

    # 3) Write to Excel
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        for name, df in sheets.items():
            sheet = name[:31]
            df.to_excel(writer, index=False, sheet_name=sheet)
            ws = writer.sheets[sheet]
            ws.freeze_panes(1, 0)
            ws.autofilter(0, 0, max(0, len(df)), max(0, len(df.columns)-1))

    print(f"[export] Wrote {len(sheets)} sheets to {out_path}")

def stageA():


    print("\n=== Stage A: Fit elastic + optionally θ_D ===")
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)
    # 1) Build bounds: elastic free, p27/p15 with small symmetric ranges around 0
    boundsA = make_stageA_bounds(
        PARAMS_INIT, LOWER_BOUND, UPPER_BOUND,
        free_elastic=(0, 1, 2, 3),   # v00, a1, a2, a3
        v00_idx=0,
        v00_range=(26.0, 28.0),      # keep near physical range
        elastic_range=(-1e4, 1e4),
        free_theta_idx=None          # <-- don't use θ-like range for p27
    )
    B = list(boundsA)
    # was (200, 2000) – likely incorrect for a coeff that starts at 0
    B[27] = (-6.0, 6.0)
    B[15] = (-2.0, 2.0)   # second-most influential; optional but helpful
    boundsA = B

    # print("bounds[27] =", boundsA[27])        # should NOT be (p0,p0)
    # print("init p27  =", PARAMS_INIT[27])

    # 4) prefit v00 to give the solver a head start
    x0 = prefit_v00(PARAMS_INIT, datasets)
    x0[0] = 27.3        # pull baseline down toward the 20 K median
    x0[27] = 0.0         # neutral start for curvature
    x0[15] = 0.0         # get off the bound
    # x0 = sweep_p27(x0, datasets, lo=boundsA[27][0], hi=boundsA[27][1], n=37)
    # 5) optimise
    cb = _make_outfun_vm(*datasets)
    res = minimize(
        # uses your combined_cost_vm internally
        fun=lambda x, *a: _cost_only_vm(x, *a),
        x0=np.asarray(x0, dtype=float),
        args=datasets,
        method="L-BFGS-B",
        bounds=boundsA,
        callback=cb,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                     ftol=FUNCTION_TOL, gtol=GRADIENT_TOL, maxls=60),
    )
    print("\n[Stage A] status:", res.message)
    print("[Stage A] final cost:", res.fun)
    print("[Stage A] fitted elastic:",
          f"v00={res.x[0]:.6g}, a1={res.x[1]:.6g}, a2={res.x[2]:.6g}, a3={res.x[3]:.6g}")
    print(f"[Stage A] fitted p27={res.x[27]:.6g}, p15={res.x[15]:.6g}")
    diagnose_subline(res.x, datasets,  Tt=KRYPTON_T_t)
    # Quick visual checks for Vm only
    quick_plot_vm_only(res.x, datasets, Tt=KRYPTON_T_t)


def huber_rmse_rel(y_exp, y_mod, w=None, delta=2.0, floor=1e-3):
    """
    Huber RMSE on *relative* residuals: r = (y_mod - y_exp)/max(|y_exp|, floor).
    'delta' in units of sigma (relative). Returns a scalar.
    """
    y_exp = np.asarray(y_exp, float)
    y_mod = np.asarray(y_mod, float)
    den = np.maximum(np.abs(y_exp), floor)
    r = (y_mod - y_exp) / den
    if w is None:
        w = np.ones_like(r)
    w = np.asarray(w, float)

    a = np.abs(r)
    quad = a <= delta
    # classic Huber loss
    L = np.where(quad, 0.5 * r**2, delta * (a - 0.5 * delta))
    # weighted mean, then sqrt to look like RMSE
    return float(np.sqrt(np.average(L, weights=w)))


def make_stageCP_bounds(xref):
    b = [(v, v) for v in np.asarray(xref, float)]  # freeze all by default
    def setb(i, lo, hi): b[i] = (float(lo), float(hi))
    # Debye / Grüneisen
    setb(9,  40.0, 120.0)   # theta0
    setb(15, 1.0,  4.0)     # gamma0
    setb(21, 0.0,  0.02)    # q small & >=0
    # anharmonic
    setb(27, -0.08, 0.02)   # allow aa negative
    setb(28,  0.1,  2.0)
    setb(29,  5.0,  15.0)
    return b


def make_stageCP_bounds(xref):
    """
    Freeze all params at xref except thermal ones:
    [9] theta0, [15] gamma0, [21] q, [27:30] aa,bb,cc.
    """
    b = [(float(v), float(v)) for v in np.asarray(xref, float)]
    def setb(i, lo, hi): b[i] = (float(lo), float(hi))

    # Debye / Grüneisen
    setb(9,  40.0, 120.0)   # theta0 [K]
    setb(15, 1.0,   4.0)    # gamma0 [-]
    setb(21, 0.0,  0.02)    # q [-] (keep small & >=0)

    # Anharmonic correction (allow aa negative to raise cp)
    setb(27, -0.08, 0.02)   # aa
    setb(28,  0.10,  2.0)   # bb
    setb(29,  5.0,  15.0)   # cc
    return b


def stageCP(x_start=None):
    """
    Fit only c_p along sublimation (T <= Tt) by adjusting thermal parameters:
    p[9], p[15], p[21], p[27:30]. Elastic (p[0:4]) stays fixed.
    """
    # --- data ---
    kr = load_all_gas_data('krypton', read_from_excel=False)
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
     T_H_melt, p_H_melt, delta_H_melt, H_fluid_melt) = extract_datasets(kr)

    m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub) &
         np.isfinite(cp_sub) & (T_cp_sub <= KRYPTON_T_t))
    T = np.asarray(T_cp_sub, float)[m]
    P = np.asarray(p_cp_sub,  float)[m]
    Y = np.asarray(cp_sub,    float)[m]   # J/(mol·K)

    # --- starting guess ---
    if x_start is None:
        x_start = np.zeros(31, float)
        # keep your current elastic (or seed roughly if unknown)
        x_start[0] = 27.3
        x_start[1:4] = [2600., 7000., 10.0]

        # >>> Krypton thermal guess (better than argon’s) <<<
        x_start[9] = 65.0    # theta0 [K]  (lower lifts low-T cp)
        x_start[15] = 3.0     # gamma0
        x_start[21] = 0.003   # q
        x_start[27] = -0.03   # aa  (negative raises cp)
        x_start[28] = 0.9     # bb
        x_start[29] = 10.0    # cc
        x_start[30] = KRYPTON_REFERENCE_ENTROPY  # doesn’t affect cp

    # --- objective (robust, weighted) ---
    W = np.where(T < CP_SPLIT_K, CP_W_BELOW, CP_W_ABOVE)

    def obj(x):
        ymod = _safe_props_vector(T, P, x, idx=IDX["cp"])
        # robust relative RMSE with a small floor for stability
        return huber_rmse_rel(Y, ymod, w=W, delta=2.0, floor=1e-3)

    # --- quick diagnostic around ~20 K with current start ---
    def debug_terms(Ts, x):
        for TT in Ts:
            pp = psub_curve(TT)
            pr = compute_thermo_props(TT, pp, x)
            Cv, Cp = pr[5], pr[4]
            print(f"T={TT:5.1f} K  Cv={Cv:8.4f}  Cp={Cp:8.4f}")

    print("\n[cp] probe with start:")
    debug_terms((15, 18, 20, 22, 25, 30), x_start)

    # --- optimize ---
    bnds = make_stageCP_bounds(x_start)
    res = minimize(
        fun=lambda x: float(obj(x)),
        x0=np.asarray(x_start, float),
        method="L-BFGS-B",
        bounds=bnds,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                     ftol=FUNCTION_TOL, gtol=GRADIENT_TOL)
    )

    print("\n[Stage CP] status:", res.message)
    print("[Stage CP] loss:", res.fun)
    print("[Stage CP] fitted:",
          f"theta0={res.x[9]:.4g}, gamma0={res.x[15]:.4g}, q={res.x[21]:.4g}, "
          f"aa={res.x[27]:.4g}, bb={res.x[28]:.4g}, cc={res.x[29]:.4g}")

    print("\n[cp] probe with fitted params:")
    debug_terms((15, 18, 20, 22, 25, 30), res.x)

    # --- quick plot ---
    yfit = _safe_props_vector(T, P, res.x, idx=IDX["cp"])
    order = np.argsort(T)
    plt.figure(figsize=(6.8, 4.4))
    plt.scatter(T, Y, s=16, alpha=0.7, label="cp exp (sub)")
    plt.plot(T[order], yfit[order], lw=2, label="cp model")
    plt.xlabel("T [K]")
    plt.ylabel(r"$c_p$ [J mol$^{-1}$ K$^{-1}$]")
    plt.title("Krypton cp along sublimation")
    plt.legend()
    plt.tight_layout()
    plt.show()

    return res.x


def sweep_p27(params, datasets, lo=-3.0, hi=3.0, n=31):
    grid = np.linspace(lo, hi, n)
    base = np.array(params, float)
    best_cost, best_v = np.inf, None
    for v in grid:
        test = base.copy()
        test[27] = v
        cost, _ = combined_cost_vm(test, *datasets)
        if np.isfinite(cost) and cost < best_cost:
            best_cost, best_v = cost, v
    print(f"[prefit] p27* = {best_v:.4g}  cost = {best_cost:.6g}")
    base[27] = best_v
    return base
def vm_sub_sensitivities(params, idxs_to_test=None, rel_step=1e-2):
    Ts = np.linspace(20.0, min(110.0, KRYPTON_T_t-2.0), 30)
    ps = safe_psub(Ts)
    base = _safe_props_vector(Ts, ps, params, idx=IDX["Vm"])
    ok = np.isfinite(base)
    Ts0, V0 = Ts[ok], base[ok]
    k0 = np.polyfit(Ts0, V0, 1)[0]                  # slope dV/dT
    c0 = np.polyfit(Ts0, V0, 2)[0] if Ts0.size >= 3 else 0.0  # curvature

    n = len(params)
    if idxs_to_test is None:
        idxs_to_test = range(n)
    out = []
    for i in idxs_to_test:
        p = np.array(params, float)
        step = rel_step * (abs(p[i]) if p[i] != 0 else 1.0)
        p[i] += step
        Vi = _safe_props_vector(Ts, ps, p, idx=IDX["Vm"])
        ok = np.isfinite(Vi)
        if ok.sum() >= 5:
            ki = np.polyfit(Ts[ok], Vi[ok], 1)[0]
            ci = np.polyfit(Ts[ok], Vi[ok], 2)[0] if ok.sum() >= 3 else 0.0
            out.append((i, (ki-k0)/max(1e-12, step), (ci-c0)/max(1e-12, step)))
    # sort by |curvature sensitivity| then |slope sensitivity|
    out.sort(key=lambda t: (abs(t[2]), abs(t[1])), reverse=True)
    return out  # list of (param_index, dSlope/dp, dCurv/dp)

def stageA2():

        # --- Argon-style initial guess for Krypton (31 params, your EOS layout) ---
        # indices:
        # [0]=v00, [1:4]=a1..a3, [4:9]=unused, [9:15]=Theta[0..5],
        # [15:21]=gamma[0..5], [21:27]=q[0..5], [27:30]=aa,bb,cc, [30]=S*(g,Tt,pt)
    PARAMS_INIT = np.zeros(31, dtype=float)

    # elastic polynomial in ln z
    PARAMS_INIT[0] = 22.555           # v00  [cm^3/mol]  (argon guess)
    PARAMS_INIT[1:4] = [2656.5, 7298.0, 10.0]  # a1,a2,a3 [MPa]

    # Debye/Grüneisen branch 0
    PARAMS_INIT[9] = 86.44            # theta_D0 [K]
    PARAMS_INIT[15] = 2.68             # gamma0   [-]
    PARAMS_INIT[21] = 0.0024           # q        [-]

    # anharmonic thermal correction
    PARAMS_INIT[27:30] = [0.0128, 0.388, 7.85]   # aa, bb, cc

    # start entropy reference aligned to REFPROP to make ΔS_triple ~ 0 initially
    PARAMS_INIT[30] = St_REFPROP         # J/(mol·K)

    # Build once (right after you load/extract datasets)
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)
    teacher = make_subline_teacher(datasets, Tt=KRYPTON_T_t)

    # Subline-focused weights first; re-enable MELT later for a short refine
    stageA_weights = dict(SUB=1.0, MELT=0.0, ANCH=1.2,
                        SLOPE=0.8, TEACH=0.8, MONO=3.0, CONV=0.2)

    # Prefit v00 + give p27 room, p15 modest
    x0 = prefit_v00(PARAMS_INIT, datasets)
    B = list(make_stageA_bounds(PARAMS_INIT, LOWER_BOUND, UPPER_BOUND,
                                free_elastic=(0, 1, 2, 3), v00_idx=0,
                                v00_range=(target_v00(datasets)-0.4,
                                        target_v00(datasets)+0.4),
                                elastic_range=(-1e4, 1e4)))
    B[27] = (-6.0, 6.0)
    B[15] = (-2.0, 2.0)
    boundsA = B

    # Optimize (note: no heavy callback)
    res = minimize(
        fun=lambda x, *a: float(combined_cost_vm(x, *a, teacher=teacher, Tt=KRYPTON_T_t,
                                                weights=stageA_weights)[0]),
        x0=np.asarray(x0, float),
        args=datasets, method="L-BFGS-B", bounds=boundsA,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                    ftol=FUNCTION_TOL, gtol=GRADIENT_TOL)  
    )

    # Diagnose/plot
    diagnose_subline(res.x, datasets, Tt=KRYPTON_T_t)
    quick_plot_vm_only(res.x, datasets, Tt=KRYPTON_T_t)

    # # Then a short refine with melt on:
    # refine_weights = dict(SUB=1.0, MELT=1.0, ANCH=0.8,
    #                     SLOPE=0.4, TEACH=0.5, MONO=2.0, CONV=0.2)
    # res2 = minimize(
    #     fun=lambda x, *a: float(combined_cost_vm(x, *a, teacher=teacher, Tt=KRYPTON_T_t,
    #                                             weights=refine_weights)[0]),
    #     x0=res.x, args=datasets, method="L-BFGS-B", bounds=boundsA,
    #     options=dict(disp=True, maxiter=200, eps=1e-3)
    # )


# 1) Elastic & subline shape
xA = stageA2()       # (if you return res.x inside stageA2, else capture from print/log)

# 2) cp-only thermal fit (start from xA if you have it; else None uses argon guess)
xCP = stageCP(x_start=xA)

# 3) optional joint tighten
# xFinal = stage_joint_refine(xCP)
