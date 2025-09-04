# fit_argon_vm_quick.py
# Minimal, fast sanity-check of Argon solid EOS against Vm(T) data
# - Reads Vm_sublimation_all_data.txt and Vm_melting_all_data.txt
# - Uses your compute_thermo_props(T, p, params)
# - Pressures along sub/melt from Argon p(T)
# - Starts from Argon paper's optimized params
# - Quick fit + overlay plot

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- your codebase imports ---
from computethermoprops import compute_thermo_props        # you already have this
from thermal_script import melting_pressure_equation
# If you already define Argon constants in constants.py, import them here instead
# from constants import ARGON_T_t, ARGON_P_t, ARGON_E_1_SUB, ARGON_E_2_SUB, ARGON_E_3_SUB, ARGON_E_4, ARGON_E_5, ARGON_E_6, ARGON_E_7

# -------------------------
# Argon triple & p(T) set-up
# -------------------------
ARGON_T_t = 83.8058   # K  (triple-point temperature)
ARGON_P_t = 0.0689    # MPa (triple-point pressure)

# If your thermal_script already has the right Argon coefficients wired, you can leave these as-is,
# otherwise put the correct ones here (numbers below are placeholders – fill from your Argon set).
ARGON_E_1_SUB = 10.763
ARGON_E_2_SUB = -1.526
ARGON_E_3_SUB = -0.4245

ARGON_E_4 = 1506.54
ARGON_E_5 = 1.731
ARGON_E_6 = 4677.16
ARGON_E_7 = 0.98493


def sublimation_pressure_equation(T, e1, e2, e3, T_t, P_t):
    """
    Argon solid–gas line (Huber-style form):
      ln(p/P_t) = -(T_t/T) * [ e1*θ + e2*θ^(3/2) + e3*θ^5 ],  θ = 1 - T/T_t
    Returns MPa. NaN for T<=0 or T>T_t (outside subline domain).
    """
    T = np.asarray(T, dtype=float)
    out = np.full_like(T, np.nan, dtype=float)

    valid = np.isfinite(T) & (T > 0.0) & (T <= T_t)  # allow exactly T_t
    if not np.any(valid):
        return out if out.ndim else np.nan

    theta = 1.0 - (T[valid] / T_t)          # θ in [0,1]
    poly = (e1*theta
            + e2*np.power(theta, 1.5)
            + e3*np.power(theta, 5.0))
    expo = -(T_t / T[valid]) * poly        # <-- critical minus sign
    out[valid] = P_t * np.exp(expo)

    # scalar-friendly return
    return out if out.ndim else float(out)


def psub_curve(T):
    """Argon solid–gas (sublimation) pressure in MPa."""
    T = np.asarray(T, float)
    return sublimation_pressure_equation(T, ARGON_E_1_SUB, ARGON_E_2_SUB, ARGON_E_3_SUB, ARGON_T_t, ARGON_P_t)


def pmelt_curve(T):
    """Argon solid–liquid (melting) pressure in MPa."""
    T = np.asarray(T, float)
    return melting_pressure_equation(T, ARGON_E_4, ARGON_E_5, ARGON_E_6, ARGON_E_7, ARGON_T_t, ARGON_P_t)


# --------------
# Data locations
# --------------
HERE = Path(__file__).resolve().parent
FILE_SUB = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"

# -------------------------
# Load Vm(T) experimental
# -------------------------


def _load_vm_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")  # whitespace-separated
    # normalize columns
    colmap = {c.lower(): c for c in df.columns}
    # Expect columns with these names (case-insensitive): Year, Author, Temperature, CellVolume
    T = df[colmap.get("temperature", "Temperature")].astype(float)
    Vm = df[colmap.get("cellvolume", "CellVolume")].astype(float)
    out = pd.DataFrame({"T": T, "Vm_exp": Vm})
    out = out[np.isfinite(out["T"]) & np.isfinite(out["Vm_exp"])]
    out = out.sort_values("T").reset_index(drop=True)
    return out


sub_df = _load_vm_table(FILE_SUB)   # T <= Tt typically (but we’ll mask below)
melt_df = _load_vm_table(FILE_MELT)  # T >= Tt

# Build pressure for each point
sub_df["p"] = psub_curve(sub_df["T"].values)
melt_df["p"] = pmelt_curve(melt_df["T"].values)

# keep proper sides of Tt
sub_df = sub_df[(sub_df["T"] <= ARGON_T_t)].copy()
melt_df = melt_df[(melt_df["T"] >= ARGON_T_t)].copy()

# -------------------------
# Parameters (Argon paper)
# -------------------------
# Layout follows your existing vector convention:
#   0=v00, 1..3 elastic (c1,c2,c3) [MPa], then thermal blocks…
# If your compute_thermo_props expects a fixed length (e.g. 31),
# we’ll build that and fill what we know; the rest set to 0.
PARAMS_INIT = np.zeros(31, dtype=float)

# Paper values you posted
v00 = 22.555            # cm^3/mol
c1 = 2656.5            # MPa
c2 = 7298.0            # MPa
c3 = 10.0              # MPa
thetaD0 = 86.44         # K   (Debye θ_0)
gamma0 = 2.68          # Grüneisen γ_0
qD = 0.0024        # q
b1, b2, b3 = 0.0128, 0.388, 7.85
# Entropy reference at triple (if used downstream)
S_ref = 130.37          # J/(mol K)

# Map into your parameter vector (same slots as your Krypton script)
PARAMS_INIT[0] = v00
PARAMS_INIT[1] = c1
PARAMS_INIT[2] = c2
PARAMS_INIT[3] = c3

# A common choice in your code was: slot 15 for a thermal “knob”
# Put Debye θ_0 somewhere sensible if your model uses it (adjust if your code expects elsewhere)
PARAMS_INIT[15] = thetaD0
# If your model uses γ0, qD, anharmonics b1–b3 at known slots, set them here too.
# (Otherwise they can stay 0 and your model will ignore them.)
# Example (adjust indices to your model):
# PARAMS_INIT[16] = gamma0
# PARAMS_INIT[17] = qD
# PARAMS_INIT[21] = b1
# PARAMS_INIT[22] = b2
# PARAMS_INIT[23] = b3

# -------------------------
# Helpers
# -------------------------
IDX = dict(Vm=0)  # your compute_thermo_props returns Vm at index 0


def _safe_props_vm(T_arr, p_arr, params):
    T_arr = np.asarray(T_arr, float)
    p_arr = np.asarray(p_arr, float)
    out = np.full_like(T_arr, np.nan, dtype=float)
    last_vm = None
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        if not (np.isfinite(T) and np.isfinite(p) and p > 0):
            continue
        try:
            vm = compute_thermo_props(T, p, params)[IDX["Vm"]]
            if np.isfinite(vm):
                out[i] = vm
                last_vm = vm
        except Exception:
            pass
    return out


def rmse(y, yhat):
    m = np.isfinite(y) & np.isfinite(yhat)
    if not np.any(m):
        return 1e9
    d = y[m] - yhat[m]
    return float(np.sqrt(np.mean(d*d)))

# -------------------------
# Cost (Vm only, simple)
# -------------------------


def cost_vm_only(params):
    vm_sub_hat = _safe_props_vm(
        sub_df["T"].values,  sub_df["p"].values,  params)
    vm_melt_hat = _safe_props_vm(
        melt_df["T"].values, melt_df["p"].values, params)
    c_sub = rmse(sub_df["Vm_exp"].values,  vm_sub_hat)
    c_melt = rmse(melt_df["Vm_exp"].values, vm_melt_hat)
    # small anchors to discourage drift (median near 20 K and near Tt-2 K)
    anch = 0.0
    if len(sub_df) > 3:
        T_lo, T_hi = 20.0, min(110.0, ARGON_T_t-2.0)

        def _median_near(df, Tc, half=5.0):
            m = np.abs(df["T"].values - Tc) <= half
            return float(np.nanmedian(df.loc[m, "Vm_exp"])) if np.any(m) else float(np.nanmedian(df["Vm_exp"]))
        V_lo_exp = _median_near(sub_df, T_lo)
        V_hi_exp = _median_near(sub_df, T_hi)
        V_lo_mod = _safe_props_vm([T_lo], [psub_curve(T_lo)], params)[0]
        V_hi_mod = _safe_props_vm([T_hi], [psub_curve(T_hi)], params)[0]
        if np.isfinite(V_lo_mod) and np.isfinite(V_hi_mod):
            anch = (V_lo_mod - V_lo_exp)**2 + (V_hi_mod - V_hi_exp)**2
    # Weights tuned for a quick check
    return 1.0*c_sub + 0.8*c_melt + 0.3*anch


# -------------------------
# Bounds & quick optimize
# -------------------------
LB = PARAMS_INIT.copy()
UB = PARAMS_INIT.copy()
# Free v00 and the 3 elastic constants (others clamped)
wiggle_v = 0.8               # +/- cm^3/mol around paper v00 (tight)
LB[0], UB[0] = v00 - wiggle_v, v00 + wiggle_v
elastic_span = 1.0e4
for i in (1, 2, 3):
    LB[i], UB[i] = PARAMS_INIT[i] - elastic_span, PARAMS_INIT[i] + elastic_span
# Optionally free thetaD0 a bit if your model uses it
LB[15], UB[15] = max(5.0, thetaD0-30.0), thetaD0+30.0

bounds = [(float(lo), float(hi)) for lo, hi in zip(LB, UB)]
for T in [10, 20, 40, 60, 80]:
    print(T, psub_curve(T))

print("Starting cost (paper params):", cost_vm_only(PARAMS_INIT))
p = PARAMS_INIT.copy()
p[0] = 30.0  # try huge change
vm_probe = compute_thermo_props(20.0, psub_curve(20.0), p)[0]
print("Vm@20K with v00=30 ->", vm_probe)

res = minimize(
    fun=lambda x: float(cost_vm_only(x)),
    x0=PARAMS_INIT,
    method="L-BFGS-B",
    bounds=bounds,
    options=dict(maxiter=60, ftol=1e-5, gtol=1e-6, maxls=40, disp=True),
)

# 1) Are we actually using ARGON parameters & correct indices?
print("param vector length & a few slots:")
for i, (name) in enumerate(["v00", "c1", "c2", "c3", "theta0", "gamma", "q", "b1", "b2", "b3"]):
    print(i, name, params[i])   # <-- replace indices with your actual layout

# 2) Units sanity (if compute_thermo_props expects Pa):
c1_MP a, c2_MPa, c3_MPa = 2656.5, 7298.0, 10.0
c1_Pa, c2_Pa, c3_Pa = [x*1e6 for x in (c1_MPa, c2_MPa, c3_MPa)]
print(c1_Pa, c2_Pa, c3_Pa)

# 3) Thermal response along the subline (should not be flat):
Ts = np.array([10, 20, 40, 60, 80.0])
# your fixed sublimation curve (MPa); convert to Pa if needed by the EOS
ps = psub(Ts)
Vm, Alpha = [], []
for T, p in zip(Ts, ps):
    props = compute_thermo_props(T, p, params)  # or (T, p*1e6, params) if Pa
    Vm.append(props[IDX["Vm"]])
    Alpha.append(props[IDX["Alpha"]])
print("Vm(T) subline:", Vm)
print("Alpha(T) subline:", Alpha)  # expect ~1e-4–1e-3 K^-1 as T→Tt


print("\nStatus:", res.message)
print("Final cost:", res.fun)
print("v00, c1, c2, c3, thetaD0:",
      res.x[0], res.x[1], res.x[2], res.x[3], res.x[15])

# -------------------------
# Plot overlays
# -------------------------


def plot_overlays(params):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)

    # sub
    vm_sub_hat = _safe_props_vm(sub_df["T"].values, sub_df["p"].values, params)
    ax[0].scatter(sub_df["T"], sub_df["Vm_exp"],
                  s=18, label="Exp (sub)", alpha=0.75)
    ax[0].scatter(sub_df["T"], vm_sub_hat, s=18,
                  label="Model (sub)", alpha=0.9, color="crimson")
    ax[0].axvline(ARGON_T_t, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    # melt
    vm_melt_hat = _safe_props_vm(
        melt_df["T"].values, melt_df["p"].values, params)
    ax[1].scatter(melt_df["T"], melt_df["Vm_exp"],
                  s=18, label="Exp (melt)", alpha=0.75)
    ax[1].scatter(melt_df["T"], vm_melt_hat, s=18,
                  label="Model (melt)", alpha=0.9, color="crimson")
    ax[1].axvline(ARGON_T_t, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()


plot_overlays(res.x)
