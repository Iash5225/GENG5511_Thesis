import numpy as np
import matplotlib.pyplot as plt
import datetime
from constants import *
import os
from thermal_script import melting_pressure_equation,sublimation_pressure_equation

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
    ax.plot(np.asarray(T)[order],
            np.asarray(y_eos)[order],
            lw=2.0, color="red", label="EOS")   # <-- bright red line
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylab)
    ax.legend(frameon=False, fontsize=8)


def plot_all_overlays_grid(params, datasets, Tt, pt,
                           compute_thermo_props,
                           St_REFPROP, Ht_REFPROP,
                           psub_curve=None, pmelt_curve=None,   # <— NEW
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
    

    if psub_curve is None:
        psub_curve = lambda T: np.full_like(
            np.asarray(T, float), np.nan, dtype=float)
    if pmelt_curve is None:
        pmelt_curve = lambda T: np.full_like(
            np.asarray(T, float), np.nan, dtype=float)

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

    # # KappaT (from BetaT_sub)
    # m = (np.asarray(T_BetaT_sub) <= Tt)
    # if np.any(m):
    #     T = np.asarray(T_BetaT_sub)[m]
    #     p = np.asarray(p_BetaT_sub)[m]
    #     y_exp = np.asarray(BetaT_sub)[m]
    #     y_eos = _safe_props_vector(
    #         T, p, params, IDX["KappaT"], compute_thermo_props)
    #     panels.append((r"$\kappa_T$ — sublimation", "T / K",
    #                   r"$\kappa_T$ (MPa$^{-1}$)", T, y_exp, y_eos))

    # # KappaS (from BetaS_sub)
    # m = (np.asarray(T_BetaS_sub) <= Tt)
    # if np.any(m):
    #     T = np.asarray(T_BetaS_sub)[m]
    #     p = np.asarray(p_BetaS_sub)[m]
    #     y_exp = np.asarray(BetaS_sub)[m]
    #     y_eos = _safe_props_vector(
    #         T, p, params, IDX["KappaS"], compute_thermo_props)
    #     panels.append((r"$\kappa_S$ — sublimation", "T / K",
    #                   r"$\kappa_S$ (MPa$^{-1}$)", T, y_exp, y_eos))

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

        # --- Sublimation pressure (aux curve vs experimental, log y) ---
    m = (np.isfinite(T_sub) & np.isfinite(p_sub) &
        np.isfinite(G_fluid_sub) & np.isfinite(V_fluid_sub) & (np.asarray(T_sub) <= Tt))
    if np.any(m):
        T_pts = np.asarray(T_sub)[m]
        p_pts = np.asarray(p_sub)[m]
        # dense curve from min(T_pts) up to Tt
        T_curve = np.linspace(np.nanmin(T_pts), Tt, 300)
        p_curve = np.asarray(psub_curve(T_curve), float)*1e6 #convert bar to MPA 
        # guard: log axis needs positive pressures
        ok = np.isfinite(p_curve) & (p_curve > 0)
        T_curve, p_curve = T_curve[ok], p_curve[ok]
        panels.append({
            "type": "pressure",
            "title": "Sublimation pressure",
            "xlabel": "T / K",
            "ylabel": "p (MPa)",
            "T_pts": T_pts, "p_pts": p_pts,
            "T_curve": T_curve, "p_curve": p_curve
        })

    # --- Melting pressure (aux curve vs experimental, log y) ---
    m = (np.isfinite(T_melt) & np.isfinite(p_melt) &
        np.isfinite(G_fluid_melt) & np.isfinite(V_fluid_melt) & (np.asarray(T_melt) >= Tt))
    if np.any(m):
        T_pts = np.asarray(T_melt)[m]
        p_pts = np.asarray(p_melt)[m]
        # dense curve from max(Tt, min exp T) to max exp T
        t_lo = max(Tt, float(np.nanmin(T_pts)))
        t_hi = float(np.nanmax(T_pts))
        if t_hi > t_lo:
            T_curve = np.linspace(t_lo, t_hi, 300)
            p_curve = np.asarray(pmelt_curve(T_curve), float)
            ok = np.isfinite(p_curve) & (p_curve > 0)
            T_curve, p_curve = T_curve[ok], p_curve[ok]
            panels.append({
                "type": "pressure",
                "title": "Melting pressure",
                "xlabel": "T / K",
                "ylabel": "p (MPa)",
                "T_pts": T_pts, "p_pts": p_pts,
                "T_curve": T_curve, "p_curve": p_curve
            })

        n = len(panels)
        if n == 0:
                print("No panels to plot.")
                return

        ncols = max(1, int(ncols))
        nrows = (n + ncols - 1) // ncols


        fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                                squeeze=False, sharex=sharex)
        k = 0
        for r in range(nrows):
            for c in range(ncols):
                ax = axes[r, c]
                if k < n:
                    panel = panels[k]
                    if isinstance(panel, dict) and panel.get("type") == "pressure":
                        _overlay_pressure_ax(
                            ax,
                            panel["T_pts"], panel["p_pts"],
                            panel["T_curve"], panel["p_curve"],
                            panel["title"], panel["ylabel"],
                            xlabel=panel["xlabel"]
                        )
                    else:
                        # your original tuple-style panels
                        title, xlab, ylab, T, y_exp, y_eos = panel
                        _overlay_ax(ax, T, y_exp, y_eos, title, ylab, xlabel=xlab)
                else:
                    ax.set_axis_off()
                k += 1


        fig.suptitle("Experimental vs EOS overlays", y=1.02, fontsize=12)
        fig.tight_layout()

        # --- save with timestamp in EOS_FILEPATH ---
        os.makedirs(EOS_FILEPATH, exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = os.path.join(
            EOS_FILEPATH, f"experimental_vs_eos_{timestamp}.png")
        fig.savefig(filename, dpi=300)
        print(f"Figure saved to {filename}")

        plt.close(fig)  # free memory


def _make_phase_curve(T_data, T_lo, T_hi, n=200):
    """Dense, sorted temperature grid spanning [T_lo, T_hi] ∩ [min(T_data), max(T_data)]."""
    if len(T_data) == 0:
        return np.array([]), np.array([])
    tmin = max(T_lo, float(np.nanmin(T_data)))
    tmax = min(T_hi, float(np.nanmax(T_data)))
    if not np.isfinite(tmin) or not np.isfinite(tmax) or tmax <= tmin:
        return np.array([]), np.array([])
    Tgrid = np.linspace(tmin, tmax, n)
    return Tgrid, None  # second value filled by caller


def _overlay_pressure_ax(ax, T_pts, p_pts, T_curve, p_curve, title, ylabel, xlabel="T / K"):
    """Scatter experimental points and draw phase-curve; set log y-scale and a neat legend."""
    ax.scatter(T_pts, p_pts, s=28, facecolors='none',
               edgecolors='C0', linewidths=1.2, label="exp.")
    ax.plot(T_curve, p_curve, 'k-', lw=1.8, label="aux curve")
    ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
