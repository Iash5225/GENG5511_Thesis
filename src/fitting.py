from math import isfinite
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from constants import *
from thermopropsv2 import compute_thermo_props
from fitting_helper import *
from plot_eos import plot_all_overlays_grid
from thermal_script import melting_pressure_equation, sublimation_pressure_equation


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


def plot_alpha_sublimation(params, datasets, *, Tt=KRYPTON_T_t, ax=None, label="α model", save_path=None):
    """
    Plot alpha(T) (thermal expansion coefficient) along sublimation branch.
    """
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

    # mask valid points
    m = (np.isfinite(T_alpha_sub) & np.isfinite(p_alpha_sub) &
         np.isfinite(alpha_sub) & (T_alpha_sub <= Tt))
    if not np.any(m):
        raise ValueError(
            "No finite alpha data on sublimation branch (T <= Tt).")

    T = np.asarray(T_alpha_sub, float)[m]
    P = np.asarray(p_alpha_sub,  float)[m]
    Y = np.asarray(alpha_sub,    float)     # [1/K]

    # model
    yfit = (_safe_props_vector(T, P, params, idx=IDX["Alpha"])
            if "_safe_props_vector" in globals()
            else np.array([compute_thermo_props(float(t), float(p), params)[IDX["Alpha"]] for t, p in zip(T, P)], float))

    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 4.4))
    order = np.argsort(T)
    ax.scatter(T, Y, s=16, alpha=0.75, label="α exp (sub)")
    ax.plot(T[order], yfit[order], lw=2, label=label)
    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$\alpha$ [K$^{-1}$]")
    ax.set_title("Krypton: thermal expansion along sublimation")
    ax.legend()
    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=200)
    return ax


def plot_bulk_moduli_sublimation(params, datasets, *, Tt=KRYPTON_T_t, ax=None, save_path=None):
    """
    Plot isothermal (K_T) and adiabatic (K_S) bulk moduli along the sublimation branch.

    Experimental inputs are compressibilities (BetaT_sub, BetaS_sub) in [1/MPa].
    We invert them to bulk moduli: K = 1 / Beta (MPa).

    Model values are computed from compute_thermo_props:
      K_T(model) = 1 / KappaT,  K_S(model) = 1 / KappaS

    Parameters
    ----------
    params : array-like
        EOS parameter vector (len 31).
    datasets : tuple
        Output of extract_datasets(...).
    Tt : float
        Triple-point temperature [K].
    ax : matplotlib.axes.Axes or None
        If provided, plot into this axes; else create a new figure.
    save_path : str or Path or None
        If given, save figure to this path.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
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

    create_fig = ax is None
    if create_fig:
        fig, ax = plt.subplots(figsize=(7.2, 4.6))

    # ---------- K_T (isothermal bulk modulus) ----------
    mT = (np.isfinite(T_BetaT_sub) & np.isfinite(
        BetaT_sub) & (T_BetaT_sub <= Tt))
    if np.any(mT):
        Tt_exp = np.asarray(T_BetaT_sub, float)[mT]
        Pt_exp = np.asarray(p_BetaT_sub,  float)[mT]
        kappaT_exp = np.asarray(BetaT_sub, float)[mT]  # [1/MPa]
        # guard against 0 or negative compressibilities
        valid = np.isfinite(kappaT_exp) & (kappaT_exp > 0)
        KT_exp = np.full_like(kappaT_exp, np.nan)
        KT_exp[valid] = 1.0 / kappaT_exp[valid]  # MPa

        # model at same (T,P)
        KT_mod = _safe_props_vector(Tt_exp, Pt_exp, params, idx=IDX["KappaT"])
        with np.errstate(divide='ignore', invalid='ignore'):
            KT_mod = np.where((KT_mod > 0) & np.isfinite(
                KT_mod), 1.0 / KT_mod, np.nan)

        order = np.argsort(Tt_exp)
        ax.scatter(Tt_exp, KT_exp, s=18, alpha=0.75, label=r"$K_T$ exp")
        ax.plot(Tt_exp[order], KT_mod[order], lw=2.0, label=r"$K_T$ model")

    # ---------- K_S (adiabatic bulk modulus) ----------
    mS = (np.isfinite(T_BetaS_sub) & np.isfinite(
        BetaS_sub) & (T_BetaS_sub <= Tt))
    if np.any(mS):
        Ts_exp = np.asarray(T_BetaS_sub, float)[mS]
        # p_BetaS_sub may be None in your loader; if so, we pass NaN (model path only needs T,P pair)
        if p_BetaS_sub is not None:
            Ps_exp = np.asarray(p_BetaS_sub, float)[mS]
        else:
            Ps_exp = np.full_like(Ts_exp, np.nan)
        kappaS_exp = np.asarray(BetaS_sub, float)[mS]  # [1/MPa]
        validS = np.isfinite(kappaS_exp) & (kappaS_exp > 0)
        KS_exp = np.full_like(kappaS_exp, np.nan)
        KS_exp[validS] = 1.0 / kappaS_exp[validS]  # MPa

        # model at same (T,P)
        KS_mod = _safe_props_vector(Ts_exp, Ps_exp, params, idx=IDX["KappaS"])
        with np.errstate(divide='ignore', invalid='ignore'):
            KS_mod = np.where((KS_mod > 0) & np.isfinite(
                KS_mod), 1.0 / KS_mod, np.nan)

        order = np.argsort(Ts_exp)
        ax.scatter(Ts_exp, KS_exp, s=18, alpha=0.75, label=r"$K_S$ exp")
        ax.plot(Ts_exp[order], KS_mod[order], lw=2.0, label=r"$K_S$ model")

    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$K$ [MPa]")
    ax.set_title("Krypton: bulk moduli along sublimation")
    ax.legend()
    ax.grid(False)
    plt.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=200)

    return ax

def plot_kappa_sublimation(params, datasets, *, Tt=KRYPTON_T_t, ax=None, save_path=None):
    """
    Plot isothermal (κ_T) and adiabatic (κ_S) compressibility along sublimation.
    """
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

    fig, ax = (plt.subplots(figsize=(7, 4.5)) if ax is None else (None, ax))

    # --- κ_T ---
    mT = (np.isfinite(T_BetaT_sub) & np.isfinite(
        BetaT_sub) & (T_BetaT_sub <= Tt))
    if np.any(mT):
        T = np.asarray(T_BetaT_sub, float)[mT]
        P = np.asarray(p_BetaT_sub,  float)[mT]
        Y = np.asarray(BetaT_sub, float)    # [1/MPa]
        yfit = (_safe_props_vector(T, P, params, idx=IDX["KappaT"])
                if "_safe_props_vector" in globals()
                else np.array([compute_thermo_props(float(t), float(p), params)[IDX["KappaT"]] for t, p in zip(T, P)], float))
        order = np.argsort(T)
        ax.scatter(T, Y, s=16, alpha=0.7, label=r"$\kappa_T$ exp")
        ax.plot(T[order], yfit[order], lw=2, label=r"$\kappa_T$ model")

    # --- κ_S ---
    mS = (np.isfinite(T_BetaS_sub) & np.isfinite(
        BetaS_sub) & (T_BetaS_sub <= Tt))
    if np.any(mS):
        T = np.asarray(T_BetaS_sub, float)[mS]
        P = np.asarray(p_BetaS_sub,  float)[
            mS] if p_BetaS_sub is not None else np.full_like(T, np.nan)
        Y = np.asarray(BetaS_sub, float)    # [1/MPa]
        yfit = (_safe_props_vector(T, P, params, idx=IDX["KappaS"])
                if "_safe_props_vector" in globals()
                else np.array([compute_thermo_props(float(t), float(p), params)[IDX["KappaS"]] for t, p in zip(T, P)], float))
        order = np.argsort(T)
        ax.scatter(T, Y, s=16, alpha=0.7, label=r"$\kappa_S$ exp")
        ax.plot(T[order], yfit[order], lw=2, label=r"$\kappa_S$ model")

    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$\kappa$ [MPa$^{-1}$]")
    ax.set_title("Krypton: compressibility along sublimation")
    ax.legend()
    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=200)
    return ax

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
        Vm_sub_dev = rms_percent(Vm_sub[m], model)   # cm^3/mol

    # Vm (melting): T >= Tt  —— use ABS RMSE in cm^3/mol
    m = (np.isfinite(T_Vm_melt) & np.isfinite(p_Vm_melt)
         & np.isfinite(Vm_melt) & (T_Vm_melt >= Tt))
    Vm_melt_dev = BIG
    if np.any(m):
        model = _safe_props_vector(T_Vm_melt[m], p_Vm_melt[m], params, idx=0)
        Vm_melt_dev = rms_percent(Vm_melt[m], model)  # cm^3/mol
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
        cp_exp = np.asarray(cp_sub)[m]
        model = _safe_props_vector(Tm, pm, params, idx=4)
        # percent deviations
        dev = (cp_exp - model) / cp_exp
        # weights depending on T threshold
        weights = np.where(Tm < CP_TEMP_THRESHOLD_K,
                       CP_WEIGHT_BELOW, CP_WEIGHT_ABOVE)
        cp_sub_dev = float(np.sqrt(np.mean((weights * dev) ** 2)))


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

def plot_cp_sublimation(params, datasets, *, Tt=KRYPTON_T_t, ax=None, label="cp model", save_path=None):
    """
    Plot cp(T) along the sublimation branch (T <= Tt) against experimental data.

    Parameters
    ----------
    params : array-like
        Current EOS parameter vector (length 31).
    datasets : tuple
        Output of extract_datasets(...). Must contain cp_sub arrays.
    Tt : float
        Triple-point temperature [K]. Defaults to KRYPTON_T_t.
    ax : matplotlib.axes.Axes or None
        If provided, plot into this axes; else create a new figure.
    label : str
        Legend label for the model curve.
    save_path : str or Path or None
        If given, save the figure to this path.

    Returns
    -------
    ax : matplotlib.axes.Axes
        Axes with the scatter and line plotted.
    """

    # unpack datasets (matches your combined_cost_function layout)
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

    # finite points on sublimation branch
    m = (np.isfinite(T_cp_sub) & np.isfinite(p_cp_sub) &
         np.isfinite(cp_sub) & (T_cp_sub <= Tt))
    if not np.any(m):
        raise ValueError(
            "No finite cp data on the sublimation branch (T <= Tt).")

    T = np.asarray(T_cp_sub, float)[m]
    P = np.asarray(p_cp_sub,  float)[m]
    Y = np.asarray(cp_sub,    float)[m]     # J/mol/K

    # model cp at the same (T,P) points
    if '_safe_props_vector' in globals():
        yfit = _safe_props_vector(T, P, params, idx=IDX["cp"])
    else:
        yfit = np.array(
            [compute_thermo_props(float(t), float(p), params)[IDX["cp"]]
             for t, p in zip(T, P)],
            dtype=float
        )

    # make plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.8, 4.4))

    order = np.argsort(T)
    ax.scatter(T, Y, s=16, alpha=0.75, label="cp exp (sublimation)")
    ax.plot(T[order], yfit[order], lw=2.0, label=label)

    ax.set_xlabel("T [K]")
    ax.set_ylabel(r"$c_p$ [J mol$^{-1}$ K$^{-1}$]")
    ax.set_title("Krypton: heat capacity along sublimation")
    ax.legend()
    ax.grid(False)
    plt.tight_layout()

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=200)

    return ax


    return res.x

def stageA2():
    """
    Vm-only fit using percent-RMSE objective (no teacher/anchors/slope terms).
    Starts from the Argon-style initial guess and uses your LOWER/UPPER bounds.
    """
    # --- Load datasets once ---
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)

    # --- Bounds: use your global LOWER_BOUND/UPPER_BOUND arrays directly ---
    # (If you want to tweak a few entries, modify LOWER_BOUND/UPPER_BOUND before this call.)
    bounds = [(float(lo), float(hi))
              for lo, hi in zip(LOWER_BOUND, UPPER_BOUND)]

    # --- Optimize Vm-only with percent-RMSE objective (no teacher) ---
    res = minimize(
        fun=lambda x, *a: float(combined_cost_vm(x, *a, Tt=KRYPTON_T_t)[0]),
        x0=np.asarray(PARAMS_INIT, float),
        args=datasets,
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                     ftol=FUNCTION_TOL, gtol=GRADIENT_TOL)
    )

    # Quick visual check for Vm
    quick_plot_vm_only(res.x, datasets, Tt=KRYPTON_T_t)

    print("\n[Stage A2] status:", res.message)
    print("final parameters:", res.x.tolist())
    return res.x

def main():
    # === 4. Package for scipy.optimize.minimize ===
    bounds = [(lo, hi) for lo, hi in zip(LOWER_BOUND, UPPER_BOUND)]
    krypton_data = load_all_gas_data('krypton', read_from_excel=False)
    datasets = extract_datasets(krypton_data)
    res = minimize(
        fun=lambda x, *a: _cost_only(x, *a),   # your combined cost wrapper
        x0=PARAMS_INIT,
        args=datasets,                         # from extract_datasets
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(disp=True, maxiter=MAX_ITERATIONS,
                     ftol=FUNCTION_TOL, gtol=GRADIENT_TOL)
    )
    print("\n[Stage A2] status:", res.message)
    print("final parameters:", res.x.tolist())
    #plot_all_overlays_grid(res.x, datasets, Tt=KRYPTON_T_t, pt=KRYPTON_P_t, compute_thermo_props=compute_thermo_props,
    #                     St_REFPROP=KRYPTON_REFERENCE_ENTROPY, Ht_REFPROP=KRYPTON_REFERENCE_ENTHALPY, psub_curve=psub_curve, pmelt_curve=pmelt_curve)
    quick_plot_vm_only(res.x, datasets, Tt=KRYPTON_T_t)
    # cp plot
    plot_cp_sublimation(res.x, datasets, Tt=KRYPTON_T_t, save_path=None)
    # alpha plot
    plot_alpha_sublimation(res.x, datasets, Tt=KRYPTON_T_t, save_path=None) 
    # kappa plot
    # plot_bulk_moduli_sublimation(res.x, datasets, Tt=KRYPTON_T_t, save_path=None)
    plt.show()


# 1) Elastic & subline shape
# xA = stageA2()       # (if you return res.x inside stageA2, else capture from print/log)

# # 2) cp-only thermal fit (start from xA if you have it; else None uses argon guess)
# xCP = stageCP(x_start=xA

main()


