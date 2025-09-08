# pip install SALib numpy

import numpy as np
from SALib.sample import saltelli
from SALib.analyze import sobol
from thermopropsv2 import compute_thermo_props

# ---------- your compute_thermo_props must be in scope ----------
# Output order (from your function):
OUTS = ["Vm", "KappaT", "KappaS", "Alpha", "cp",
        "cv", "Gruneisen", "U", "S", "A", "H", "G"]

# ----- PARAMS (from your message) -----
PARAMS_INIT = np.array([
    27.21, 2739.48, 7328.48, 122.62,
    0, 0, 0, 0, 0,
    88.97, 0, 0, 0, 0, 0,
    4.05, 0, 0, 0, 0, 0,
    -3.02, 0, 0, 0, 0, 0,
    0.0, 1.92, 7.94,
    119.23
], dtype=float)

LOWER = np.array([
    25, 0, 0, 0, 0, 0, 0, 0, 0,
    20, 0, 0, 0, 0, 0,
    -5, 0, 0, 0, 0, 0,
    -10, 0, 0, 0, 0, 0,
    0, 0, 0, 0
], dtype=float)

UPPER = np.array([
    40, 10000, 10000, 10000, 0, 0, 0, 0, 0,
    300, 0, 0, 0, 0, 0,
    5, 0, 0, 0, 0, 0,
    10, 0, 0, 0, 0, 0,
    100, 100, 100, 1000
], dtype=float)

# Parameters actually used by your function:
# v00, a1,a2,a3, Th0, g0, q0, aa,bb,cc
CANDIDATE_ACTIVE = [0, 1, 2, 3, 9, 15, 21, 27, 28, 29]
PARAM_NAMES = {0: "v00", 1: "a1", 2: "a2", 3: "a3", 9: "Th0",
               15: "g0", 21: "q0", 27: "aa", 28: "bb", 29: "cc"}

ACTIVE = [i for i in CANDIDATE_ACTIVE if UPPER[i] > LOWER[i]]
names = [PARAM_NAMES[i] for i in ACTIVE]
bounds = [[float(LOWER[i]), float(UPPER[i])] for i in ACTIVE]

problem = {"num_vars": len(ACTIVE), "names": names, "bounds": bounds}

# ----- choose a small but representative (T,p) grid -----
# Fill these with temps/pressures relevant to your dataset/branches.
grid_Tp = [
    (90.0, 0.1),    # K, MPa
    (115.8, 0.073),  # near triple-point-ish
    (140.0, 0.1),
    (170.0, 0.1),
]

# ----- sampling -----
base_N = 512  # increase for tighter CIs
X = saltelli.sample(problem, base_N, calc_second_order=False)  # shape (N, P)


def run_model(theta_full, T, p):
    y = compute_thermo_props(T, p, theta_full)
    return y  # length 12

# Pre-build full-theta array for all samples (stitch ACTIVE dims into full vector)


def theta_from_row(row):
    theta = PARAMS_INIT.copy()
    for k, idx in enumerate(ACTIVE):
        theta[idx] = row[k]
    return theta


# Evaluate once per grid point (reuse X)
# Collect ST indices per (grid point, property, parameter)
ST_all = []  # list of arrays with shape (len(OUTS), len(ACTIVE))
S1_all = []

for (T_eval, p_eval) in grid_Tp:
    # Evaluate model for all samples at this (T,p)
    Ymat = np.empty((X.shape[0], len(OUTS)), dtype=float)
    for n in range(X.shape[0]):
        theta = theta_from_row(X[n])
        y = run_model(theta, T_eval, p_eval)
        Ymat[n, :] = y

    # Analyze per output
    ST_this = np.zeros((len(OUTS), len(ACTIVE)))
    S1_this = np.zeros_like(ST_this)

    for j_out in range(len(OUTS)):
        Y = Ymat[:, j_out]
        # guard against NaNs/Infs
        if not np.all(np.isfinite(Y)):
            # if this output misbehaves at this point, set zeros
            Si = {"S1": np.zeros(len(ACTIVE)), "ST": np.zeros(len(ACTIVE))}
        else:
            Si = sobol.analyze(
                problem, Y, calc_second_order=False, print_to_console=False)
        S1_this[j_out, :] = np.nan_to_num(Si["S1"], nan=0.0)
        ST_this[j_out, :] = np.nan_to_num(Si["ST"], nan=0.0)

    S1_all.append(S1_this)
    ST_all.append(ST_this)

S1_all = np.stack(S1_all, axis=0)  # (n_grid, n_out, n_par)
ST_all = np.stack(ST_all, axis=0)

# ----- aggregate across grid (median absolute ST works well) -----
agg_ST = np.median(np.abs(ST_all), axis=0)  # (n_out, n_par)
agg_S1 = np.median(np.abs(S1_all), axis=0)

# ----- pretty print rankings per property -----


def print_top_k(agg, k=5, label="Total-order (ST)"):
    print(f"\n=== Parameter importance by property ({label}) ===")
    for j_out, outname in enumerate(OUTS):
        scores = agg[j_out, :]
        if np.all(scores == 0):
            continue
        order = np.argsort(scores)[::-1]
        top = [(names[i], scores[i]) for i in order[:k]]
        txt = ", ".join([f"{nm}: {val:.3f}" for nm, val in top])
        print(f"{outname:>10s} -> {txt}")


print_top_k(agg_ST, k=6, label="Total-order (ST, median over grid)")
# If you want to see first-order too:
print_top_k(agg_S1, k=6, label="First-order (S1, median over grid)")
