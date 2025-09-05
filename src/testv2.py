from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# --- your codebase ---
from computethermoprops import compute_thermo_props  # returns Vm at index 0

# =============================================================================
# Argon triple point and p(T)
# =============================================================================
ARGON_TT = 83.8058      # K
ARGON_PT = 0.0689       # MPa

# Huber paper (Argon) coefficients
E1_SUB, E2_SUB, E3_SUB = -10.763, -1.526, -0.4245
E4, E5, E6, E7 = 1506.54, 1.731, 4677.16, 0.98493

# If compute_thermo_props uses MPa (as in your code), keep 1.0.
# If it needs Pa, set 1e6.
P_TO_EOS = 1.0


def psub(T):
    """Sublimation pressure (MPa). Valid for 0 < T <= Tt."""
    T = np.asarray(T, float)
    out = np.full_like(T, np.nan, dtype=float)
    m = np.isfinite(T) & (T > 0.0) & (T <= ARGON_TT)
    if not np.any(m):
        return out if out.ndim else np.nan
    th = 1.0 - T[m]/ARGON_TT
    expo = (ARGON_TT/T[m])*(E1_SUB*th + E2_SUB *
                            np.power(th, 1.5) + E3_SUB*np.power(th, 5.0))
    out[m] = ARGON_PT*np.exp(expo)
    return out if out.ndim else float(out)


def pmelt(T):
    """Melting pressure (MPa). Valid for T >= Tt."""
    T = np.asarray(T, float)
    return ARGON_PT*(1.0 + E4*((T/ARGON_TT) - 1.0)**E5 + E6*((T/ARGON_TT) - 1.0)**E7)


# =============================================================================
# Data loading (robust to header/spacing/case)
# =============================================================================
HERE = Path(__file__).resolve().parent
FILE_SUB = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"


def load_vm_table(path: Path) -> pd.DataFrame:
    """Expect columns like Temperature / CellVolume (any case, spaces/underscores ok)."""
    try:
        df = pd.read_csv(path, sep=r"\s+", engine="python")
    except Exception:
        df = pd.read_csv(path)

    # normalize names
    key = {c.lower().replace(" ", "").replace("_", ""): c for c in df.columns}
    # accept a few aliases seen in different dumps
    Tcol = key.get("temperature") or key.get("t") or list(df.columns)[-2]
    Vcol = key.get("cellvolume") or key.get(
        "vm") or key.get("v_m") or list(df.columns)[-1]

    T = pd.to_numeric(df[Tcol], errors="coerce")
    Vm = pd.to_numeric(df[Vcol], errors="coerce")
    print(f"Loaded {len(df)} rows from {path.name}, "
          f"keeping {np.isfinite(T).sum()} with valid T and {np.isfinite(Vm).sum()} with valid Vm.")
    out = pd.DataFrame({"T": T, "Vm_exp": Vm})
    out = out[np.isfinite(out.T) & np.isfinite(
        out.Vm_exp)].sort_values("T").reset_index(drop=True)
    return out


def load_vm_table_simple(path: Path) -> pd.DataFrame:
    """
    File format (exact headers):
      Year  Author  Temperature  CellVolume
    Tolerates repeated header lines later in the file.
    Returns columns:
      Year, Author, Temperature, CellVolume, Vm_exp, T
    """
    # Read whitespace-separated, using the first row as header
    df = pd.read_csv(path, sep=r"\s+", engine="python", header=0)
    # Clean possible whitespace in column names
    df.columns = df.columns.str.strip()

    # Coerce numeric cols; repeated headers will turn into NaN and be dropped
    df["Temperature"] = pd.to_numeric(df["Temperature"], errors="coerce")
    df["CellVolume"] = pd.to_numeric(df["CellVolume"],  errors="coerce")

    # Keep only valid numeric rows
    df = df[np.isfinite(df["Temperature"]) &
            np.isfinite(df["CellVolume"])].copy()

    # Sort and add aliases expected by downstream code
    df = df.sort_values("Temperature").reset_index(drop=True)
    df["Vm_exp"] = df["CellVolume"]
    df["T"] = df["Temperature"]

    return df

# =============================================================================
# Parameters
# =============================================================================
# 31-long vector layout you’ve been using:
# 0=v00, 1..3 = c1,c2,c3 (MPa), 9=thetaD0, 15=gamma0, 21=q,
# 27..29 = b1,b2,b3 (aa,bb,cc), 30=s (entropy ref; unused here)
PARAMS_INIT = np.zeros(31, float)

# Paper values
v00 = 22.555
c1 = 2656.5
c2 = 7298.0
c3 = 10.0
theta = 86.44
gamma = 2.68
qD = 0.0024
b1, b2, b3 = 0.0128, 0.388, 7.85

PARAMS_INIT[0] = v00
PARAMS_INIT[1] = c1
PARAMS_INIT[2] = c2
PARAMS_INIT[3] = c3
PARAMS_INIT[9] = theta
PARAMS_INIT[15] = gamma
PARAMS_INIT[21] = qD
PARAMS_INIT[27] = b1
PARAMS_INIT[28] = b2
PARAMS_INIT[29] = b3

# Fit only v00 and c1..c3 at first (others fixed)


def make_bounds(p0):
    LB, UB = p0.copy(), p0.copy()

    # v00: allow vertical shift
    LB[0], UB[0] = p0[0] - 1.5, p0[0] + 1.5

    # elastic c1..c3: keep reasonable room
    for i in (1, 2, 3):
        LB[i], UB[i] = p0[i] - 1.0e4, p0[i] + 1.0e4

    # thermal block – these actually set the subline slope/shape
    LB[9],  UB[9] = 60.0, 120.0     # θD,0
    LB[15], UB[15] = 2.0,  3.6       # γ0
    LB[21], UB[21] = 0.0,  0.010     # q
    LB[27], UB[27] = 0.0,  0.040     # b1 (aa)
    LB[28], UB[28] = 0.05, 1.20      # b2 (bb)
    LB[29], UB[29] = 2.0,  12.0      # b3 (cc)

    return [(float(lo), float(hi)) for lo, hi in zip(LB, UB)]

# =============================================================================
# EOS wrapper (Vm only)
# =============================================================================


def vm_hat(T_arr, p_arr, params):
    T_arr = np.asarray(T_arr, float)
    p_arr = np.asarray(p_arr, float)*P_TO_EOS
    out = np.full_like(T_arr, np.nan, dtype=float)
    for i, (T, p) in enumerate(zip(T_arr, p_arr)):
        if not (np.isfinite(T) and np.isfinite(p)):
            continue
        try:
            out[i] = compute_thermo_props(float(T), float(p), params)[0]
        except Exception:
            pass
    return out


def rmse(y, yhat):
    y = np.asarray(y, float)
    yhat = np.asarray(yhat, float)
    m = np.isfinite(y) & np.isfinite(yhat)
    if not np.any(m):
        return 1e9
    d = y[m] - yhat[m]
    return float(np.sqrt(np.mean(d*d)))

# =============================================================================
# Cost: Vm-only (sub + melt)
# =============================================================================
def cost_vm_only(params, sub_df, melt_df):
    v_sub = vm_hat(sub_df["T"].values,  sub_df["p"].values,  params)
    v_melt = vm_hat(melt_df["T"].values, melt_df["p"].values, params)
    # equal-weight RMSE; or weight by counts if you prefer
    return rmse(sub_df["Vm_exp"].values, v_sub) + rmse(melt_df["Vm_exp"].values, v_melt)

# =============================================================================
# Main
# =============================================================================
def plot_overlays(params, sub_df, melt_df):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)

    # ===== SUBLIMATION =====
    # experimental (scatter)
    ax[0].scatter(sub_df["T"], sub_df["Vm_exp"], s=18, alpha=0.75, label="Exp (sub)")

    # model (smooth line)
    Tsub_line = np.linspace(sub_df["T"].min(), sub_df["T"].max(), 400)
    psub_line = psub(Tsub_line)
    Vsub_line = vm_hat(Tsub_line, psub_line, params)
    ax[0].plot(Tsub_line, Vsub_line, lw=2.0, label="Model (sub)")

    ax[0].axvline(ARGON_TT, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    # ===== MELTING =====
    # experimental (scatter)
    ax[1].scatter(melt_df["T"], melt_df["Vm_exp"], s=18, alpha=0.75, label="Exp (melt)")

    # model (smooth line)
    Tmelt_line = np.linspace(melt_df["T"].min(), melt_df["T"].max(), 400)
    pmelt_line = pmelt(Tmelt_line)
    Vmelt_line = vm_hat(Tmelt_line, pmelt_line, params)
    ax[1].plot(Tmelt_line, Vmelt_line, lw=2.0, label="Model (melt)")

    ax[1].axvline(ARGON_TT, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()


def main():
    sub_df = load_vm_table_simple(FILE_SUB)
    melt_df = load_vm_table_simple(FILE_MELT)
    print(f"Sublimation data: T {sub_df['T'].min()}..{sub_df['T'].max()} K")
    print(f"Melting data:    T {melt_df['T'].min()}..{melt_df['T'].max()} K")
    # build pressures from T via Argon correlations
    sub_df["p"] = psub(sub_df["T"].values)
    melt_df["p"] = pmelt(melt_df["T"].values)
    # keep correct sides of Tt and drop any p NaNs just in case
    sub_df = sub_df[(sub_df["T"] > 0) & (sub_df["T"] <= ARGON_TT)
                    & np.isfinite(sub_df["p"])].copy()
    melt_df = melt_df[(melt_df["T"] >= ARGON_TT) &
                      np.isfinite(melt_df["p"])].copy()

    print(f"Sublimation points kept: {len(sub_df)}  "
          f"(T range {sub_df['T'].min():.2f}–{sub_df['T'].max():.2f} K)")
    print(f"Melting points kept:     {len(melt_df)} "
          f"(T range {melt_df['T'].min():.2f}–{melt_df['T'].max():.2f} K)")

    bounds = make_bounds(PARAMS_INIT)
    start = cost_vm_only(PARAMS_INIT, sub_df, melt_df)
    print(f"Start cost (paper params): {start:.6g}")

    res = minimize(
        fun=lambda x: cost_vm_only(x, sub_df, melt_df),
        x0=PARAMS_INIT,
        method="L-BFGS-B",
        bounds=bounds,
        options=dict(maxiter=120, ftol=1e-6, gtol=1e-6, maxls=40, disp=True),
    )

    print("\nStatus:", res.message)
    print("Final cost:", res.fun)
    print("v00, c1, c2, c3:", res.x[0], res.x[1], res.x[2], res.x[3])

    plot_overlays(res.x, sub_df, melt_df)


if __name__ == "__main__":
    main()
