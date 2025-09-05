# plot_argon_vm_and_cp.py
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from computethermoprops import compute_thermo_props
# from testcompute import compute_thermo_props_for_result  # optional

# ---------- Argon triple point & p(T) ----------
T_t = 83.8058   # K
P_t = 0.0689    # MPa

E1_sub, E2_sub, E3_sub = -10.763, -1.526, -0.4245
E4, E5, E6, E7 = 1506.54, 1.731, 4677.16, 0.98493


def p_sub(T):
    T = np.asarray(T, float)
    theta = 1.0 - T / T_t
    expo = (T_t / T) * (E1_sub*theta + E2_sub*theta**1.5 + E3_sub*theta**5)
    base = P_t * np.exp(expo)
    # works for scalar or array:
    return np.where((T <= 0) | (T > T_t), np.nan, base)


def probe_near_20K(params, Tlist, p_of_T):
    from computethermoprops import compute_thermo_props
    # or copy your evaluate_gas into scope
    from computethermoprops import evaluate_gas
    print("\n[T probe]  T    p(MPa)      v     dpdv     KappaT      cp     cv")
    for T in Tlist:
        p = float(p_of_T(T))
        props = compute_thermo_props(T, p, params)
        v = float(props[0])
        cp = float(props[4])
        cv = float(props[5])
        # re-evaluate dpdv at that (T,v)
        v00, a1, a2, a3 = params[0], params[1], params[2], params[3]
        Th = np.array(params[9:15])
        g = np.array(params[15:21])
        q = np.array(params[21:27])
        aa, bb, cc = params[27], params[28], params[29]
        a = [3.0]
        _, dpdv, *_ = evaluate_gas(T, v, v00, a, a1,
                                   a2, a3, Th, g, q, aa, bb, cc)
        KappaT = -1.0/(max(v, 1e-9)*dpdv) if dpdv != 0 else np.inf
        print(
            f"{T:8.3f} {p:9.4f} {v:8.4f} {dpdv:9.4f} {KappaT:10.4g} {cp:7.3f} {cv:7.3f}")


def p_melt(T):
    T = np.asarray(T, float)
    base = P_t * (1.0 + E4*(T/T_t - 1.0)**E5 + E6*(T/T_t - 1.0)**E7)
    return np.where(T < T_t, np.nan, base)


# ---------- Argon params (paper layout your EOS expects) ----------
PARAMS_ARGON = np.zeros(31, float)
PARAMS_ARGON[0] = 22.555
PARAMS_ARGON[1:4] = [2656.5, 7298.0, 10.0]
PARAMS_ARGON[9]  = 86.44
PARAMS_ARGON[15] = 2.68
PARAMS_ARGON[21] = 0.0024
PARAMS_ARGON[27:30] = [0.0128, 0.388, 7.85]
# PARAMS_ARGON[30] not needed for Vm/cp

P_TO_EOS = 1.0  # set to 1e6 if your EOS expects Pa

# ---------- Data I/O ----------
HERE = Path(__file__).resolve().parent
FILE_SUB  = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"

FILE_CP_SUB  = HERE / "cp_sublimation_all_data.txt"
FILE_CP_MELT = HERE / "cp_melting_all_data.txt"

def load_vm_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    T  = pd.to_numeric(df["Temperature"], errors="coerce")
    Vm = pd.to_numeric(df["CellVolume"],  errors="coerce")
    return (pd.DataFrame({"T": T, "Vm": Vm})
            .dropna().sort_values("T").reset_index(drop=True))

def load_cp_table_if_exists(path: Path):
    if not path.exists():
        print(f"[cp] file not found: {path}")
        return None
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    T = pd.to_numeric(df.get("Temperature"), errors="coerce")
    cp_col = None
    for name in ("cp", "Cp", "CP", "HeatCapacity", "HeatCapacity_J_per_mol_K"):
        if name in df.columns:
            cp_col = name
            break
    if cp_col is None:
        print(f"[cp] no cp column in: {path}")
        return None
    cp = pd.to_numeric(df[cp_col], errors="coerce")
    out = pd.DataFrame({"T": T, "cp": cp}).dropna().sort_values("T").reset_index(drop=True)
    return out

# ---------- Model lines ----------
def prop_model_line(Ts, branch: str, idx: int):
    ps = p_sub(Ts) if branch == "sub" else p_melt(Ts)
    vals = np.full_like(Ts, np.nan, float)
    for i, (T, p) in enumerate(zip(Ts, ps)):
        if not np.isfinite(p):
            continue
        vals[i] = float(compute_thermo_props(float(T), float(p)*P_TO_EOS, PARAMS_ARGON)[idx])
    return vals

def vm_model_line(Ts, branch): return prop_model_line(Ts, branch, idx=0)
def cp_model_line(Ts, branch): return prop_model_line(Ts, branch, idx=4)

# ---------- cp identity quick check ----------
def check_cp_identity(params, p_of_T, Tlist, label=""):
    print(f"\n[cp identity check] {label}")
    for T in Tlist:
        p = float(p_of_T(T))
        props = compute_thermo_props(T, p, params)
        v, kT, aL, cp, cv = props[0], props[1], props[3], props[4], props[5]
        lhs = cp - cv
        rhs = T * v * (aL**2) / kT   # MPa·cm^3 = J, so units match J/mol/K
        print(f"T={T:6.2f} K  cp-cv={lhs:9.5f}   T*V*α²/κT={rhs:9.5f}   Δ={lhs-rhs:+.2e}")

# ---------- Main ----------
def main():
    sub  = load_vm_table(FILE_SUB)
    melt = load_vm_table(FILE_MELT)
    cp_sub_df  = load_cp_table_if_exists(FILE_CP_SUB)
    cp_melt_df = load_cp_table_if_exists(FILE_CP_MELT)

    tmin_sub,  tmax_sub  = float(sub["T"].min()),  float(sub["T"].max())
    tmin_melt, tmax_melt = float(melt["T"].min()), float(melt["T"].max())
    Ts_sub  = np.linspace(max(1.0, tmin_sub),  min(T_t, tmax_sub),  500)
    Ts_melt = np.linspace(max(T_t, tmin_melt), tmax_melt,           500)

    Vm_sub_line  = vm_model_line(Ts_sub,  "sub")
    Vm_melt_line = vm_model_line(Ts_melt, "melt")
    cp_sub_line  = cp_model_line(Ts_sub,  "sub")
    cp_melt_line = cp_model_line(Ts_melt, "melt")

    # 2x2 grid
    fig, ax = plt.subplots(2, 2, figsize=(13, 9), sharex="col")
    probe_near_20K(PARAMS_ARGON, [15, 18, 19, 20, 21, 22, 24, 26], p_sub)

    # Vm — sub
    ax[0, 0].scatter(sub["T"].to_numpy(), sub["Vm"].to_numpy(), s=16, alpha=0.7, label="Exp (sub)")
    ax[0, 0].plot(Ts_sub, Vm_sub_line, lw=2, label="Model (sub)")
    ax[0, 0].axvline(T_t, ls="--", lw=1, color="gray")
    ax[0, 0].set_title("Vm along sublimation")
    ax[0, 0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0, 0].legend()

    # Vm — melt
    ax[0, 1].scatter(melt["T"].to_numpy(), melt["Vm"].to_numpy(), s=16, alpha=0.7, label="Exp (melt)")
    ax[0, 1].plot(Ts_melt, Vm_melt_line, lw=2, label="Model (melt)")
    ax[0, 1].axvline(T_t, ls="--", lw=1, color="gray")
    ax[0, 1].set_title("Vm along melting")
    ax[0, 1].legend()

    # cp — sub (model + optional exp)
    # cp — sub  (use scatter for the model instead of a line)
    if cp_sub_df is not None:
        ax[1, 0].scatter(cp_sub_df["T"].to_numpy(), cp_sub_df["cp"].to_numpy(),
                        s=16, alpha=0.7, label="cp exp (sub)")
    ax[1, 0].scatter(Ts_sub, cp_sub_line, s=14, alpha=0.8,
                    marker="o", label="cp model (sub)")
    # ax[1, 0].plot(Ts_sub, cp_sub_line, lw=2, label="cp model (sub)")  # ← remove

    ax[1, 0].axvline(T_t, ls="--", lw=1, color="gray")
    ax[1, 0].set_title(r"$c_p$ along sublimation")
    ax[1, 0].set_xlabel("T [K]")
    ax[1, 0].set_ylabel(r"$c_p$ [J mol$^{-1}$ K$^{-1}$]")
    ax[1, 0].legend()

    # # cp — melt (model + optional exp)
    # if cp_melt_df is not None:
    #     ax[1, 1].scatter(cp_melt_df["T"].to_numpy(), cp_melt_df["cp"].to_numpy(),
    #                      s=16, alpha=0.7, label="cp exp (melt)")
    # ax[1, 1].plot(Ts_melt, cp_melt_line, lw=2, label="cp model (melt)")
    # ax[1, 1].axvline(T_t, ls="--", lw=1, color="gray")
    # ax[1, 1].set_title(r"$c_p$ along melting")
    # ax[1, 1].set_xlabel("T [K]")
    # ax[1, 1].legend()

    plt.tight_layout()
    plt.show()

    # --- optional: Argon identity sanity check ---
    # check_cp_identity(PARAMS_ARGON, p_sub,  [2,5,10,20,21,22,80],   label="Argon — sub")
    # check_cp_identity(PARAMS_ARGON, p_melt, [90,120,200,300,350],   label="Argon — melt")

if __name__ == "__main__":
    main()
