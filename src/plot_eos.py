import numpy as np
import matplotlib.pyplot as plt

# indices in compute_thermo_props
IDX = dict(Vm=0, KappaT=1, KappaS=2, Alpha=3, cp=4, H=10, G=11)


def _safe_props_vector(T, p, params, idx, compute_thermo_props):
    T = np.asarray(T, float)
    p = np.asarray(p, float)
    out = np.full_like(T, np.nan, dtype=float)
    for i, (Ti, pi) in enumerate(zip(T, p)):
        try:
            pr = compute_thermo_props(float(Ti), float(pi), params)
            if np.all(np.isfinite(pr)):
                out[i] = float(pr[idx])
        except Exception:
            pass
    return out


def _calc_deltaH_S_at_triple(params, Tt, pt, compute_thermo_props,
                             St_REFPROP, Ht_REFPROP):
    tp = compute_thermo_props(Tt, pt, params)
    if isinstance(tp, tuple) and len(tp) == 2 and isinstance(tp[1], dict):
        _, meta = tp
        dH = float(meta.get("deltaH_triple"))
        dS = float(meta.get("deltaS_triple"))
    else:
        props = np.asarray(tp, float)
        dS = params[30] - St_REFPROP
        Ht_fitted = props[11] + Tt * params[30]
        dH = Ht_fitted - Ht_REFPROP
    return dH, dS


def _overlay_ax(ax, T, y_exp, y_eos, title, ylab, xlabel="T / K"):
    if T is None or len(T) == 0:
        ax.set_axis_off()
        return
    order = np.argsort(T)
    ax.scatter(T, y_exp, s=12, label="exp")
    ax.plot(np.asarray(T)[order], np.asarray(
        y_eos)[order], lw=1.6, label="EOS")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylab)
    ax.legend(frameon=False, fontsize=8)


def plot_all_overlays_grid(params, datasets, Tt, pt,
                           compute_thermo_props,
                           St_REFPROP, Ht_REFPROP,
                           ncols=3, figsize=(13, 9), sharex=False):
    """
    Makes a grid of subplots:
      Vm_sub, Vm_melt, Vm_highp (if any), cp_sub, Alpha_sub, KappaT_sub, KappaS_sub,
      H_solid (sub & melt), p_sub (model-implied vs exp), p_melt (model-implied vs exp).
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

    dHtr, dStr = _calc_deltaH_S_at_triple(params, Tt, pt, compute_thermo_props,
                                          St_REFPROP, Ht_REFPROP)

    panels = []

    # Vm (sublimation)
    m = (np.asarray(T_Vm_sub) <= Tt)
    if np.any(m):
        T = np.asarray(T_Vm_sub)[m]
        p = np.asarray(p_Vm_sub)[m]
        y_exp = np.asarray(Vm_sub)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["Vm"], compute_thermo_props)
        panels.append(("Vm — sublimation", "T / K",
                      r"$V_m$ (cm$^3$ mol$^{-1}$)", T, y_exp, y_eos))

    # Vm (melting)
    m = (np.asarray(T_Vm_melt) >= Tt)
    if np.any(m):
        T = np.asarray(T_Vm_melt)[m]
        p = np.asarray(p_Vm_melt)[m]
        y_exp = np.asarray(Vm_melt)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["Vm"], compute_thermo_props)
        panels.append(("Vm — melting", "T / K",
                      r"$V_m$ (cm$^3$ mol$^{-1}$)", T, y_exp, y_eos))

    # Vm (high-p path)
    if len(T_Vm_highp) > 0:
        T = np.asarray(T_Vm_highp)
        p = np.asarray(p_Vm_highp)
        y_exp = np.asarray(Vm_highp)
        y_eos = _safe_props_vector(
            T, p, params, IDX["Vm"], compute_thermo_props)
        panels.append(
            ("Vm — high p", "T / K", r"$V_m$ (cm$^3$ mol$^{-1}$)", T, y_exp, y_eos))

    # cp (sublimation)
    m = (np.asarray(T_cp_sub) <= Tt)
    if np.any(m):
        T = np.asarray(T_cp_sub)[m]
        p = np.asarray(p_cp_sub)[m]
        y_exp = np.asarray(cp_sub)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["cp"], compute_thermo_props)
        panels.append(("cp — sublimation", "T / K",
                      r"$c_p$ (J mol$^{-1}$ K$^{-1}$)", T, y_exp, y_eos))

    # Alpha (sublimation)
    m = (np.asarray(T_alpha_sub) <= Tt)
    if np.any(m):
        T = np.asarray(T_alpha_sub)[m]
        p = np.asarray(p_alpha_sub)[m]
        y_exp = np.asarray(alpha_sub)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["Alpha"], compute_thermo_props)
        panels.append((r"$\alpha$ — sublimation", "T / K",
                      r"$\alpha$ (K$^{-1}$)", T, y_exp, y_eos))

    # KappaT (from BetaT_sub)
    m = (np.asarray(T_BetaT_sub) <= Tt)
    if np.any(m):
        T = np.asarray(T_BetaT_sub)[m]
        p = np.asarray(p_BetaT_sub)[m]
        y_exp = np.asarray(BetaT_sub)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["KappaT"], compute_thermo_props)
        panels.append((r"$\kappa_T$ — sublimation", "T / K",
                      r"$\kappa_T$ (MPa$^{-1}$)", T, y_exp, y_eos))

    # KappaS (from BetaS_sub)
    m = (np.asarray(T_BetaS_sub) <= Tt)
    if np.any(m):
        T = np.asarray(T_BetaS_sub)[m]
        p = np.asarray(p_BetaS_sub)[m]
        y_exp = np.asarray(BetaS_sub)[m]
        y_eos = _safe_props_vector(
            T, p, params, IDX["KappaS"], compute_thermo_props)
        panels.append((r"$\kappa_S$ — sublimation", "T / K",
                      r"$\kappa_S$ (MPa$^{-1}$)", T, y_exp, y_eos))

    # H_solid (sublimation branch)
    m = (np.isfinite(T_H_sub) & np.isfinite(delta_H_sub) &
         np.isfinite(H_fluid_sub) & (np.asarray(T_H_sub) <= Tt))
    if np.any(m):
        T = np.asarray(T_H_sub)[m]
        p = np.asarray(p_H_sub)[m]
        H_s_exp = np.asarray(H_fluid_sub)[m] - np.asarray(delta_H_sub)[m]
        H_model = _safe_props_vector(
            T, p, params, IDX["H"], compute_thermo_props) - dHtr
        panels.append((r"$H_\mathrm{solid}$ — sublimation",
                      "T / K", r"$H$ (kJ mol$^{-1}$)", T, H_s_exp, H_model))

    # H_solid (melting branch)
    m = (np.isfinite(T_H_melt) & np.isfinite(delta_H_melt) &
         np.isfinite(H_fluid_melt) & (np.asarray(T_H_melt) >= Tt))
    if np.any(m):
        T = np.asarray(T_H_melt)[m]
        p = np.asarray(p_H_melt)[m]
        H_s_exp = np.asarray(H_fluid_melt)[m] - np.asarray(delta_H_melt)[m]
        H_model = _safe_props_vector(
            T, p, params, IDX["H"], compute_thermo_props) - dHtr
        panels.append((r"$H_\mathrm{solid}$ — melting", "T / K",
                      r"$H$ (kJ mol$^{-1}$)", T, H_s_exp, H_model))

    # Sublimation pressure: model-implied p_f vs exp p_M
    m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
         np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (np.asarray(T_sub) <= Tt))
    if np.any(m):
        T = np.asarray(T_sub)[m]
        pM = np.asarray(p_sub)[m]
        Gf = np.asarray(G_fluid_sub)[m]
        Vf = np.asarray(V_fluid_sub)[m]
        pf = np.full_like(pM, np.nan, dtype=float)
        for i, (Ti, pMi, Gfi, Vfi) in enumerate(zip(T, pM, Gf, Vf)):
            pr = compute_thermo_props(Ti, pMi, params)
            if np.all(np.isfinite(pr)) and np.isfinite(Vfi) and (Vfi != pr[IDX["Vm"]]):
                dG = Gfi - pr[IDX["G"]] + dHtr - Ti*(_calc_deltaH_S_at_triple(
                    params, Tt, pt, compute_thermo_props, St_REFPROP, Ht_REFPROP)[1])
                pf[i] = pMi - dG / (Vfi - pr[IDX["Vm"]])
        panels.append(("Sublimation pressure", "T / K", "p (MPa)", T, pM, pf))

    # Melting pressure: model-implied p_f vs exp p_M
    m = (np.isfinite(T_melt) & np.isfinite(p_melt) &
         np.isfinite(G_fluid_melt) & np.isfinite(V_fluid_melt) & (np.asarray(T_melt) >= Tt))
    if np.any(m):
        T = np.asarray(T_melt)[m]
        pM = np.asarray(p_melt)[m]
        Gf = np.asarray(G_fluid_melt)[m]
        Vf = np.asarray(V_fluid_melt)[m]
        pf = np.full_like(pM, np.nan, dtype=float)
        dHtr2, dStr2 = dHtr, _calc_deltaH_S_at_triple(
            params, Tt, pt, compute_thermo_props, St_REFPROP, Ht_REFPROP)[1]
        for i, (Ti, pMi, Gfi, Vfi) in enumerate(zip(T, pM, Gf, Vf)):
            pr = compute_thermo_props(Ti, pMi, params)
            if np.all(np.isfinite(pr)) and np.isfinite(Vfi) and (Vfi != pr[IDX["Vm"]]):
                dG = Gfi - pr[IDX["G"]] + dHtr2 - Ti*dStr2
                pf[i] = pMi - dG / (Vfi - pr[IDX["Vm"]])
        panels.append(("Melting pressure", "T / K", "p (MPa)", T, pM, pf))

    # ------- draw grid -------
    n = len(panels)
    if n == 0:
        print("No panels to plot.")
        return
    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=figsize, squeeze=False, sharex=sharex)
    k = 0
    for r in range(nrows):
        for c in range(ncols):
            ax = axes[r, c]
            if k < n:
                title, xlab, ylab, T, y_exp, y_eos = panels[k]
                _overlay_ax(ax, T, y_exp, y_eos, title, ylab, xlabel=xlab)
            else:
                ax.set_axis_off()
            k += 1
    fig.suptitle("Experimental vs EOS overlays", y=1.02, fontsize=12)
    fig.tight_layout()
    plt.show()
