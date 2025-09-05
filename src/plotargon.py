# plot_argon_vm_fixed.py
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# your EOS: must return an array where index 0 is Vm [cm^3/mol]
from computethermoprops import compute_thermo_props  # (T [K], p [MPa], params)
from testcompute import compute_thermo_props_for_result

# ----------------- Argon triple point & p(T) -----------------
T_t = 83.8058   # K
P_t = 0.0689    # MPa

# Huber-style coefficients (from the paper)
E1_sub, E2_sub, E3_sub = -10.763, -1.526, -0.4245
E4, E5, E6, E7 = 1506.54, 1.731, 4677.16, 0.98493


def p_sub(T):
    T = np.asarray(T, float)
    theta = 1.0 - T / T_t
    expo = (T_t / T) * (E1_sub*theta + E2_sub*theta**1.5 + E3_sub*theta**5)
    out = P_t * np.exp(expo)
    out[(T <= 0) | (T > T_t)] = np.nan
    return out


def p_melt(T):
    T = np.asarray(T, float)
    out = P_t * (1.0 + E4*(T/T_t - 1.0)**E5 + E6*(T/T_t - 1.0)**E7)
    out[T < T_t] = np.nan
    return out


# ----------------- Fixed Argon parameters (paper) -----------------
# Layout matches your fitter:
# [0]=v00, [1:4]=c1..c3, [9]=theta0, [15]=gamma0, [21]=q0, [27:30]=aa,bb,cc
PARAMS_ARGON = np.zeros(31, float)
PARAMS_ARGON[0] = 22.555   # v00 [cm^3/mol]
PARAMS_ARGON[1] = 2656.5   # c1 [MPa]
PARAMS_ARGON[2] = 7298.0   # c2 [MPa]
PARAMS_ARGON[3] = 10.0     # c3 [MPa]
PARAMS_ARGON[9] = 86.44    # theta_D0 [K]
PARAMS_ARGON[15] = 2.68     # gamma0
PARAMS_ARGON[21] = 0.0024   # q
PARAMS_ARGON[27] = 0.0128   # aa
PARAMS_ARGON[28] = 0.388    # bb
PARAMS_ARGON[29] = 7.85     # cc
# PARAMS_ARGON[30] could be entropy ref; not used for Vm

# Set to 1e6 if your EOS wants Pa instead of MPa
P_TO_EOS = 1.0

# ----------------- Data I/O -----------------
HERE = Path(__file__).resolve().parent
FILE_SUB = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"


def load_vm_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    # Files have: Year, Author, Temperature, CellVolume
    T = pd.to_numeric(df["Temperature"], errors="coerce")
    Vm = pd.to_numeric(df["CellVolume"],   errors="coerce")
    out = pd.DataFrame({"T": T, "Vm": Vm}).dropna(
    ).sort_values("T").reset_index(drop=True)
    return out

# ----------------- Model line helper -----------------


def vm_model_line(Ts, branch: str):
    if branch == "sub":
        ps = p_sub(Ts)
    else:
        ps = p_melt(Ts)
    Vm = np.full_like(Ts, np.nan, float)
    for i, (T, p) in enumerate(zip(Ts, ps)):
        if not np.isfinite(p):
            continue
        Vm[i] = compute_thermo_props(float(T), float(
            p)*P_TO_EOS, PARAMS_ARGON)[0]  # Vm at index 0
    return Vm

# ----------------- Main -----------------


def main():
    sub = load_vm_table(FILE_SUB)
    melt = load_vm_table(FILE_MELT)

# before (buggy: sub.T is the transpose)
# Ts_sub  = np.linspace(max(1.0, sub.T.min()),  min(T_t, sub.T.max()),  400)
# Ts_melt = np.linspace(max(T_t, melt.T.min()), melt.T.max(),          400)


    # after (correct: use the "T" column explicitly and cast to float)
    tmin_sub = float(sub["T"].min())
    tmax_sub = float(sub["T"].max())
    tmin_melt = float(melt["T"].min())
    tmax_melt = float(melt["T"].max())
    print(f"Sublimation T range: {tmin_sub}..{tmax_sub} K")
    print(f"Melting T range:     {tmin_melt}..{tmax_melt} K")

    Ts_sub = np.linspace(max(1.0, tmin_sub),  min(T_t, tmax_sub),  400)
    Ts_melt = np.linspace(max(T_t, tmin_melt), tmax_melt,           400)

    print(f"Computing model Vm for {len(Ts_sub)} sublimation and {len(Ts_melt)} melting points...")
    Vm_sub_line = vm_model_line(Ts_sub,  "sub")
    Vm_melt_line = vm_model_line(Ts_melt, "melt")
    print(f"Computed model Vm for {np.isfinite(Vm_sub_line).sum()} sublimation and {np.isfinite(Vm_melt_line).sum()} melting points.")


    # Plot
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    x_sub = sub["T"].to_numpy()
    y_sub = sub["Vm"].to_numpy()

    x_melt = melt["T"].to_numpy()
    y_melt = melt["Vm"].to_numpy()

    ax[0].scatter(x_sub, y_sub, s=16, alpha=0.7, label="Exp (sub)")
    ax[0].plot(Ts_sub, Vm_sub_line, lw=2, label="Model (sub)")
    ax[0].axvline(T_t, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    ax[1].scatter(x_melt, y_melt, s=16, alpha=0.7, label="Exp (melt)")
    ax[1].plot(Ts_melt, Vm_melt_line, lw=2, label="Model (melt)")
    ax[1].axvline(T_t, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
