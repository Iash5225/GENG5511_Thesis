# === single_metric_fit.py (put this near your other helpers) ===
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# which index in compute_thermo_props output corresponds to which property
PROP_IDX = {
    "Vm": 0, "KappaT": 1, "KappaS": 2, "Alpha": 3,
    "cp": 4, "cv": 5, "Gamma": 6, "U": 7, "S": 8, "A": 9, "H": 10, "G": 11,
}

BIG = 1e4


def _safe_model_vector(T, p, params, prop_idx, compute_thermo_props):
    y = np.full_like(np.asarray(T, float), np.nan, dtype=float)
    for i, (Ti, pi) in enumerate(zip(T, p)):
        try:
            props = compute_thermo_props(float(Ti), float(pi), params)
            if np.all(np.isfinite(props)):
                y[i] = float(props[prop_idx])
        except Exception:
            pass
    return y


def _rms_percent(y_true, y_pred, scale=100.0):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    m = np.isfinite(y_true) & np.isfinite(y_pred) & (y_true != 0.0)
    if not np.any(m):
        return BIG
    err = scale * (y_true[m] - y_pred[m]) / y_true[m]
    return float(np.sqrt(np.mean(err*err))) if err.size else BIG


def make_single_metric_view(metric_key, datasets, Tt):
    """
    Returns T, p, y_exp, prop_idx for the chosen dataset.
    metric_key in {"Vm_sub","Vm_melt","cp_sub","Alpha_sub","KappaT_sub","KappaS_sub"}
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

    if metric_key == "Vm_sub":
        m = (np.asarray(T_Vm_sub) <= Tt)
        T = np.asarray(T_Vm_sub)[m]
        p = np.asarray(p_Vm_sub)[m]  # exact experimental path
        y_exp = np.asarray(Vm_sub)[m]
        print("T[min,max] =", T.min(), T.max(),
              "   p[min,max] =", p.min(), p.max())

        return T, p, y_exp, PROP_IDX["Vm"]
    if metric_key == "Vm_melt":
        m = (np.asarray(T_Vm_melt) >= Tt)
        return np.asarray(T_Vm_melt)[m], np.asarray(p_Vm_melt)[m], np.asarray(Vm_melt)[m], PROP_IDX["Vm"]
    if metric_key == "cp_sub":
        m = (np.asarray(T_cp_sub) <= Tt)
        return np.asarray(T_cp_sub)[m], np.asarray(p_cp_sub)[m], np.asarray(cp_sub)[m], PROP_IDX["cp"]
    if metric_key == "Alpha_sub":
        m = (np.asarray(T_alpha_sub) <= Tt)
        return np.asarray(T_alpha_sub)[m], np.asarray(p_alpha_sub)[m], np.asarray(alpha_sub)[m], PROP_IDX["Alpha"]
    if metric_key == "KappaT_sub":
        m = (np.asarray(T_BetaT_sub) <= Tt)
        return np.asarray(T_BetaT_sub)[m], np.asarray(p_BetaT_sub)[m], np.asarray(BetaT_sub)[m], PROP_IDX["KappaT"]
    if metric_key == "KappaS_sub":
        m = (np.asarray(T_BetaS_sub) <= Tt)
        return np.asarray(T_BetaS_sub)[m], np.asarray(p_BetaS_sub)[m], np.asarray(BetaS_sub)[m], PROP_IDX["KappaS"]
    raise ValueError(f"Unknown metric_key: {metric_key}")


def make_single_metric_cost(metric_key, datasets, Tt, compute_thermo_props):
    T, p, y_exp, prop_idx = make_single_metric_view(metric_key, datasets, Tt)

    def cost_only(params):
        y_calc = _safe_model_vector(T, p, params, prop_idx, compute_thermo_props)
        res = np.diff(y_calc)
        penalty = np.sum(res < 0) * 100.0  # penalize non-monotone Vm(T)
        return _rms_percent(y_exp, y_calc) + penalty


    # Lightweight callback for convergence trace
    it = {"i": 0}

    def cb(xk):
        it["i"] += 1
        c = cost_only(xk)
        print(f"[{metric_key}] iter={it['i']:3d}  RMS%={c:8.4f}")

    # Simple plot util to sanity-check fit afterwards
    def quick_plots(params, title_suffix=""):
        # overlay
        y_calc = _safe_model_vector(
            T, p, params, prop_idx, compute_thermo_props)
        idx = np.argsort(T)
        plt.figure(figsize=(6, 4))
        plt.scatter(T, y_exp, s=16, label="exp")
        plt.plot(T[idx], y_calc[idx], lw=1.6, label="EOS")
        plt.xlabel("T / K")
        plt.ylabel(metric_key)
        plt.title(f"{metric_key} — overlay {title_suffix}")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.show()
        # residuals
        res = 100.0*(y_exp - y_calc)/y_exp
        m = np.isfinite(res)
        plt.figure(figsize=(6, 4))
        plt.axhline(0, c='k', lw=1)
        plt.scatter(T[m], res[m], s=16)
        plt.xlabel("T / K")
        plt.ylabel("% error")
        plt.title(f"{metric_key} — residuals {title_suffix}")
        plt.tight_layout()
        plt.show()
        

    return cost_only, cb, quick_plots
# === end single_metric_fit.py ===
