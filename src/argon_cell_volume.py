# fit_argon_vm_quick.py
# Minimal sanity-check of Argon solid EOS against Vm(T) data.
# - Reads Vm_sublimation_all_data.txt and Vm_melting_all_data.txt
# - Uses your compute_thermo_props(T, p, params)
# - p(T) along sub/melt from Huber-style formulas
# - Starts from Argon paper parameters; fits a small subset
# - Plots overlays

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- your codebase ---
from computethermoprops import compute_thermo_props
# expects (T, e4, e5, e6, e7, Tt, Pt)
from thermal_script import melting_pressure_equation

IDX = dict(Vm=0)  # compute_thermo_props must return Vm at index 0

# -----------------------------------------------------------------------------
# Argon triple point & p(T) coefficients (Huber)
# -----------------------------------------------------------------------------
ARGON_T_t = 83.8058    # K
ARGON_P_t = 0.0689     # MPa

# Sublimation: ln(p/Pt) = (Tt/T) * [ e1 θ + e2 θ^(3/2) + e3 θ^5 ], θ = 1 - T/Tt
# (e1,e2,e3 are NEGATIVE in Huber ⇒ exponent < 0 ⇒ p drops towards 0 at low T)
ARGON_E_1_SUB = -10.763
ARGON_E_2_SUB = -1.526
ARGON_E_3_SUB = -0.4245

# Melting: p = Pt * [ 1 + e4 (T/Tt - 1)^e5 + e6 (T/Tt - 1)^e7 ]
ARGON_E_4 = 1506.54
ARGON_E_5 = 1.731
ARGON_E_6 = 4677.16
ARGON_E_7 = 0.98493

# If your EOS takes MPa, leave 1.0; if it wants Pa, set 1e6
P_TO_EOS = 1.0


def sublimation_pressure_equation(T, e1, e2, e3, T_t, P_t):
    """
    p_sub = P_t * exp( (T_t/T) * [ e1*θ + e2*θ^(3/2) + e3*θ^5 ] ),
    θ = 1 - T/T_t. Returns MPa. NaN for T<=0 or T>T_t.
    """
    T = np.asarray(T, dtype=float)
    out = np.full_like(T, np.nan, dtype=float)

    valid = np.isfinite(T) & (T > 0.0) & (T <= T_t)
    if np.any(valid):
        theta = 1.0 - (T[valid] / T_t)
        poly = (e1 * theta
                + e2 * np.power(theta, 1.5)
                + e3 * np.power(theta, 5.0))
        expo = (T_t / T[valid]) * poly
        out[valid] = P_t * np.exp(expo)
    return out if out.ndim else float(out)


def psub_curve(T):
    T = np.asarray(T, float)
    return sublimation_pressure_equation(
        T, ARGON_E_1_SUB, ARGON_E_2_SUB, ARGON_E_3_SUB, ARGON_T_t, ARGON_P_t
    )


def pmelt_curve(T):
    T = np.asarray(T, float)
    return melting_pressure_equation(
        T, ARGON_E_4, ARGON_E_5, ARGON_E_6, ARGON_E_7, ARGON_T_t, ARGON_P_t
    )


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
FILE_SUB = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"


def _load_vm_table(path: Path) -> pd.DataFrame:
    """Load table with columns Temperature, CellVolume (case-insensitive)."""
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python", comment="#")
    except Exception:
        df = pd.read_csv(path, comment="#")

    print(f"Loaded {len(df)} rows from {path.name} | cols: {list(df.columns)}")

    # normalize column names to locate fields
    cmap = {c.lower().replace(" ", "").replace("_", ""): c for c in df.columns}
    t_col = cmap.get("temperature", "Temperature")
    v_col = cmap.get("cellvolume",  "CellVolume")

    # safer numeric coercion (won’t crash if a stray string/unit sneaks in)
    T = pd.to_numeric(df[t_col], errors="coerce")
    Vm = pd.to_numeric(df[v_col], errors="coerce")

    out = pd.DataFrame({"T": T, "Vm_exp": Vm})
    # ✅ use column names, not .T (transpose)
    m = np.isfinite(out["T"].to_numpy()) & np.isfinite(
        out["Vm_exp"].to_numpy())
    out = out.loc[m].sort_values("T").reset_index(drop=True)

    if len(out):
        print(f"  kept {len(out)} rows | T: {out['T'].min()}..{out['T'].max()} K | "
              f"Vm: {out['Vm_exp'].min()}..{out['Vm_exp'].max()} cm^3/mol")
    else:
        print("  kept 0 rows — check column mapping and numeric coercion")
    return out


sub_df = _load_vm_table(FILE_SUB)
melt_df = _load_vm_table(FILE_MELT)
print(f"Sublimation data: T {sub_df['T'].min()}..{sub_df['T'].max()} K")
print(f"Melting data:    T {melt_df['T'].min()}..{melt_df['T'].max()} K")

# Assign pressures & keep the correct side of Tt
sub_df["p"] = psub_curve(sub_df["T"].values)
melt_df["p"] = pmelt_curve(melt_df["T"].values)
sub_df = sub_df[(sub_df["T"] <= ARGON_T_t)].copy()
melt_df = melt_df[(melt_df["T"] >= ARGON_T_t)].copy()

print(f"Using {len(sub_df)} sublimation and {len(melt_df)} melting points.")

# -----------------------------------------------------------------------------
# Parameter vector (exact indices must match compute_thermo_props)
# -----------------------------------------------------------------------------
# Vector length 31; unpack in your EOS:
# 0=v00, 1..3 = c1,c2,c3; 9..14=Th[0..5]; 15..20=g[0..5]; 21..26=q[0..5]; 27..30=aa,bb,cc,s
PARAMS_INIT = np.zeros(31, dtype=float)

# Elastic
PARAMS_INIT[0] = 22.555   # v00 [cm^3/mol]
PARAMS_INIT[1] = 2656.5   # c1  [MPa]
PARAMS_INIT[2] = 7298.0   # c2  [MPa]
PARAMS_INIT[3] = 10.0     # c3  [MPa]

# Thermal ON (this is what drives the subline slope)
thetaD0, gamma0, qD = 86.44, 2.68, 0.0024
b1, b2, b3 = 0.0128, 0.388, 7.85
PARAMS_INIT[9] = thetaD0  # Th[0]
PARAMS_INIT[15] = gamma0   # g[0]
PARAMS_INIT[21] = qD       # q[0]
PARAMS_INIT[27] = b1       # aa
PARAMS_INIT[28] = b2       # bb
PARAMS_INIT[29] = b3       # cc
# PARAMS_INIT[30] = S_ref if your EOS uses an explicit reference entropy

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _safe_props_vm(T_arr, p_arr, params):
    """Call EOS with unit conversion; return Vm vector, NaN on failure."""
    T_arr = np.asarray(T_arr, float)
    p_arr = np.asarray(p_arr, float) * P_TO_EOS
    out = np.full_like(T_arr, np.nan, dtype=float)
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        if np.isfinite(T) and np.isfinite(p) and p > 0:
            try:
                vm = compute_thermo_props(T, p, params)[IDX["Vm"]]
                if np.isfinite(vm):
                    out[i] = vm
            except Exception:
                pass
    return out


def rmse(y, yhat):
    m = np.isfinite(y) & np.isfinite(yhat)
    if not np.any(m):
        return 1e9
    d = y[m] - yhat[m]
    return float(np.sqrt(np.mean(d*d)))


def finite_polyfit(x, y, deg):
    m = np.isfinite(x) & np.isfinite(y)
    if np.sum(m) < (deg+1):
        return 0.0
    return float(np.polyfit(np.asarray(x)[m], np.asarray(y)[m], deg)[0])

# -----------------------------------------------------------------------------
# Cost (Vm-only, with gentle anchors and robust slope term on subline)
# -----------------------------------------------------------------------------


def cost_vm_only(params):
    vm_sub_hat = _safe_props_vm(
        sub_df["T"].values,  sub_df["p"].values,  params)
    vm_melt_hat = _safe_props_vm(
        melt_df["T"].values, melt_df["p"].values, params)

    c_sub = rmse(sub_df["Vm_exp"].values,  vm_sub_hat)
    c_melt = rmse(melt_df["Vm_exp"].values, vm_melt_hat)

    # Slope on 8..min(110, Tt-2)
    t_hi = max(8.0, min(110.0, ARGON_T_t - 2.0))
    mwin = sub_df["T"].between(8.0, t_hi)
    k_exp = finite_polyfit(
        sub_df.loc[mwin, "T"], sub_df.loc[mwin, "Vm_exp"], 1)
    Tgrid = np.linspace(8.0, t_hi, 80)
    Vmod = _safe_props_vm(Tgrid, psub_curve(Tgrid), params)
    k_mod = finite_polyfit(Tgrid, Vmod, 1)
    slope_pen = (k_mod - k_exp)**2

    # Anchors near ~20 K and ~t_hi to prevent drift
    anch = 0.0
    if len(sub_df) > 3:
        def median_at(Tc):
            m = np.abs(sub_df["T"].values - Tc) <= 5.0
            vals = sub_df.loc[m, "Vm_exp"].values if np.any(
                m) else sub_df["Vm_exp"].values
            return float(np.nanmedian(vals))
        T_lo, T_hi = 20.0, t_hi
        V_lo_exp, V_hi_exp = median_at(T_lo), median_at(T_hi)
        V_lo_mod = _safe_props_vm([T_lo], [psub_curve(T_lo)], params)[0]
        V_hi_mod = _safe_props_vm([T_hi], [psub_curve(T_hi)], params)[0]
        if np.isfinite(V_lo_mod) and np.isfinite(V_hi_mod):
            anch = (V_lo_mod - V_lo_exp)**2 + (V_hi_mod - V_hi_exp)**2

    # Weights tuned for quick go/no-go
    return 1.0*c_sub + 0.8*c_melt + 0.3*anch + 0.5*slope_pen

# -----------------------------------------------------------------------------
# Bounds & optimize
# -----------------------------------------------------------------------------


def make_bounds(p0):
    LB, UB = p0.copy(), p0.copy()

    # Free v00 and elastics modestly
    LB[0],  UB[0] = p0[0] - 0.8, p0[0] + 0.8
    for i in (1, 2):  # c1,c2
        LB[i], UB[i] = 0.7*p0[i], 1.3*p0[i]
    LB[3],  UB[3] = max(0.0, p0[3] - 50.0), p0[3] + 250.0  # c3

    # Allow slight retune of thermal (kept realistic)
    LB[9],  UB[9] = 60.0, 120.0   # θD,0
    LB[15], UB[15] = 1.6,  3.6     # γ0
    LB[21], UB[21] = 0.0,  0.01    # q
    LB[27], UB[27] = 0.0,  0.04    # b1
    LB[28], UB[28] = 0.05, 1.0     # b2
    LB[29], UB[29] = 1.0,  20.0    # b3

    return [(float(lo), float(hi)) for lo, hi in zip(LB, UB)]


def plot_overlays(params):
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8), sharey=True)

    vm_sub_hat = _safe_props_vm(
        sub_df["T"].values,  sub_df["p"].values,  params)
    ax[0].scatter(sub_df["T"], sub_df["Vm_exp"],
                  s=18, label="Exp (sub)", alpha=0.75)
    ax[0].scatter(sub_df["T"], vm_sub_hat,       s=18,
                  label="Model (sub)", alpha=0.9, color="crimson")
    ax[0].axvline(ARGON_T_t, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    vm_melt_hat = _safe_props_vm(
        melt_df["T"].values, melt_df["p"].values, params)
    ax[1].scatter(melt_df["T"], melt_df["Vm_exp"],
                  s=18, label="Exp (melt)", alpha=0.75)
    ax[1].scatter(melt_df["T"], vm_melt_hat,       s=18,
                  label="Model (melt)", alpha=0.9, color="crimson")
    ax[1].axvline(ARGON_T_t, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()


def main():
    # quick sanity prints
    for T in (10, 20, 40, 60, 80):
        print(f"T={T:>3.0f} K  p_sub={psub_curve(T):.6g} MPa")

    print("Starting cost (paper params):", cost_vm_only(PARAMS_INIT))

    bounds = make_bounds(PARAMS_INIT)
    res = minimize(
        fun=lambda x: float(cost_vm_only(x)),
        x0=PARAMS_INIT,
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(maxiter=80, ftol=1e-6, gtol=1e-6, maxls=40, disp=True),
    )

    print("\nStatus:", res.message)
    print("Final cost:", res.fun)
    print("v00, c1, c2, c3:", res.x[0], res.x[1], res.x[2], res.x[3])
    print("thetaD0, gamma0, q0, b1..b3:",
          res.x[9], res.x[15], res.x[21], res.x[27], res.x[28], res.x[29])

    plot_overlays(res.x)


if __name__ == "__main__":
    main()
