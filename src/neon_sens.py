import numpy as np
from thermopropsv2 import compute_thermo_props
from constants import PARAMS_INIT_NEON, NEON_E_1_SUB, NEON_E_2_SUB, NEON_E_3_SUB, NEON_T_t, NEON_P_t
from thermal_script import sublimation_pressure_equation
from fitting_helper_neon import _safe_props_vector, rms_percent, rms, extract_datasets, safe_psub, extract_datasets_with_meta, psub_curve, pmelt_curve, build_master_pointwise_df, _relative_errors, summarise_by_author, plot_thermo_variable, plot_variable_deviation

# experimental table (T in K, BetaT in MPa)
T_exp = np.array([2, 3, 4, 4, 4.2, 4.7, 5, 6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12,
                 13, 13.5, 14, 14, 15, 16, 16, 17, 18, 18, 19, 20, 21, 22, 23, 24], dtype=float)
BetaT_exp = np.array([0.000885, 0.000886, 0.0009, 0.000888, 0.000911, 0.00089, 0.000891, 0.00091, 0.000896, 0.000903, 0.00092, 0.000913, 0.000926, 0.00094, 0.000942, 0.000962,
                     0.00096, 0.000987, 0.00102, 0.000974, 0.001, 0.00105, 0.0011, 0.00107, 0.00115, 0.0012, 0.00124, 0.00126, 0.00132, 0.00139, 0.00146, 0.00154, 0.00162, 0.00171], dtype=float)

KappaT_exp = 1.0 / BetaT_exp  # MPa^-1


def p_sub(T):
    # returns MPa
    return sublimation_pressure_equation(T, NEON_E_1_SUB, NEON_E_2_SUB, NEON_E_3_SUB, NEON_T_t, NEON_P_t)


def rms_percent(K_model, K_exp):
    ok = np.isfinite(K_model) & np.isfinite(K_exp)
    if not np.any(ok):
        return np.inf
    rel = (K_model[ok] - K_exp[ok]) / K_exp[ok]
    return 100.0 * np.sqrt(np.mean(rel**2))


def score_params(params):
    Kmods = []
    for T in T_exp:
        p = safe_psub(T)
        props = compute_thermo_props(float(T), p, params)
        Kmods.append(float(props[1]))  # index 1 = KappaT
    return rms_percent(np.array(Kmods), KappaT_exp)


# quick grid over g0 (index 15) and aa (index 27); keep bb=1.2, cc=1 as current seed
base = PARAMS_INIT_NEON.copy()
best = (np.inf, None)
g_vals = np.linspace(base[15]-0.2,
                     base[15]+0.2, 21)  # e.g. 1.51 .. 3.51
aa_vals = np.linspace(0.0, 1, 21)
for g in g_vals:
    for aa in aa_vals:
        p = base.copy()
        p[15] = g
        p[27] = aa
        # leave bb,cc at base (1.2,1)
        s = score_params(p)
        if s < best[0]:
            best = (s, (g, aa))
print("Best RMS%:", best[0], "params (g0, aa):", best[1])

# show a few top candidates (optional thoroughness)
candidates = []
for g in g_vals:
    for aa in aa_vals:
        p = base.copy()
        p[15] = g
        p[27] = aa
        candidates.append((score_params(p), g, aa))
candidates = sorted(candidates, key=lambda x: x[0])[:6]
print("Top candidates (RMS%, g0, aa):")
for c in candidates:
    print(c)
