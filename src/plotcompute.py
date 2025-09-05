# plot_argon_vm_two_models_debug.py
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# expects (T [K], p [MPa], params) -> props
from computethermoprops import compute_thermo_props
# expects (T [K], p [MPa]) -> props
from testcompute import compute_thermo_props_for_result

# ----------------- constants & p(T)
T_t, P_t = 83.8058, 0.0689
E1_sub, E2_sub, E3_sub = -10.763, -1.526, -0.4245
E4, E5, E6, E7 = 1506.54, 1.731, 4677.16, 0.98493


# triple point (just for completeness here)
T_t = 83.8058   # K
P_t = 0.0689    # MPa

E1_sub, E2_sub, E3_sub = -10.763, -1.526, -0.4245
E4, E5, E6, E7 = 1506.54, 1.731, 4677.16, 0.98493


def p_sub(T):
    """Sublimation pressure (MPa) for T<=T_t. Returns scalar or array to match input."""
    T = np.asarray(T, dtype=float)
    out = np.full_like(T, np.nan, dtype=float)

    valid = np.isfinite(T) & (T > 0.0) & (T <= T_t)
    if np.any(valid):
        theta = 1.0 - (T[valid] / T_t)
        expo = (T_t / T[valid]) * (E1_sub*theta +
                                   E2_sub*theta**1.5 + E3_sub*theta**5)
        out[valid] = P_t * np.exp(expo)

    # return scalar if input was scalar
    return float(out) if out.ndim == 0 else out


def p_melt(T):
    """Melting pressure (MPa) for T>=T_t. Returns scalar or array to match input."""
    T = np.asarray(T, dtype=float)
    out = np.full_like(T, np.nan, dtype=float)

    valid = np.isfinite(T) & (T >= T_t)
    if np.any(valid):
        x = (T[valid] / T_t) - 1.0
        out[valid] = P_t * (1.0 + E4 * x**E5 + E6 * x**E7)

    return float(out) if out.ndim == 0 else out
 


# ----------------- Argon params (paper layout your EOS expects)
PARAMS = np.zeros(31, float)
PARAMS[0] = 22.555
PARAMS[1:4] = [2656.5, 7298.0, 10.0]
PARAMS[9] = 86.44
PARAMS[15] = 2.68
PARAMS[21] = 0.0024
PARAMS[27:30] = [0.0128, 0.388, 7.85]

# set 1e6 if your EOS wants Pa
P_TO_EOS = 1.0

HERE = Path(__file__).resolve().parent
FILE_SUB = HERE / "Vm_sublimation_all_data.txt"
FILE_MELT = HERE / "Vm_melting_all_data.txt"


def load_vm_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    T = pd.to_numeric(df["Temperature"], errors="coerce")
    Vm = pd.to_numeric(df["CellVolume"],  errors="coerce")
    out = pd.DataFrame({"T": T, "Vm": Vm}).dropna(
    ).sort_values("T").reset_index(drop=True)
    # sanity
    assert len(out["T"]) == len(out["Vm"]), "T and Vm column lengths differ"
    return out


def safe_vm_eos(T: float, p: float):
    try:
        props = compute_thermo_props(float(T), float(p)*P_TO_EOS, PARAMS)
        # <- if this isn’t Vm, you’ll see nonsense
        v = float(np.asarray(props)[0])
        return v
    except Exception as e:
        print(f"[EOS] fail at T={T:.3f} K, p={p:.6g} MPa: {e}")
        return np.nan


def safe_vm_res(T: float, p: float):
    try:
        props = compute_thermo_props_for_result(float(T), float(p))
        v = float(np.asarray(props)[0])
        return v
    except Exception as e:
        print(f"[Result] fail at T={T:.3f} K, p={p:.6g} MPa: {e}")
        return np.nan


def vm_line(func, Ts, branch: str):
    ps = p_sub(Ts) if branch == "sub" else p_melt(Ts)
    out = np.full_like(Ts, np.nan, float)
    for i, (T, p) in enumerate(zip(Ts, ps)):
        if np.isfinite(p):
            out[i] = func(T, p)
    return out


def main():
    sub = load_vm_table(FILE_SUB)   # cols: T, Vm
    melt = load_vm_table(FILE_MELT)  # cols: T, Vm

    # experimental arrays (avoid df.T!)
    x_sub = sub["T"].to_numpy()
    y_sub = sub["Vm"].to_numpy()
    x_melt = melt["T"].to_numpy()
    y_melt = melt["Vm"].to_numpy()

    tmin_sub,  tmax_sub = float(np.nanmin(x_sub)),  float(np.nanmax(x_sub))
    tmin_melt, tmax_melt = float(np.nanmin(x_melt)), float(np.nanmax(x_melt))

    Ts_sub = np.linspace(max(1.0, tmin_sub),   min(T_t, tmax_sub),  500)
    Ts_melt = np.linspace(max(T_t, tmin_melt),   tmax_melt,          500)

    # quick probes (safe for scalars now)
    for T in (10.0, 40.0, 80.0, 90.0, 200.0, 350.0):
        p = p_sub(T) if T <= T_t else p_melt(T)
        if np.isfinite(p):
            v0 = safe_vm_eos(T, p)
            v1 = safe_vm_res(T, p)
            print(
                f"T={T:6.1f} K  p={p:9.6g} MPa  Vm(EOS)={v0:8.4f}  Vm(Result)={v1:8.4f}")

    # build model lines
    Vm_sub_eos = vm_line(safe_vm_eos, Ts_sub,  "sub")
    Vm_sub_res = vm_line(safe_vm_res, Ts_sub,  "sub")
    Vm_melt_eos = vm_line(safe_vm_eos, Ts_melt, "melt")
    Vm_melt_res = vm_line(safe_vm_res, Ts_melt, "melt")

    print("max |ΔVm| sub :", np.nanmax(np.abs(Vm_sub_eos - Vm_sub_res)))
    print("max |ΔVm| melt:", np.nanmax(np.abs(Vm_melt_eos - Vm_melt_res)))

    # plot
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    m0 = np.isfinite(Vm_sub_eos)
    m1 = np.isfinite(Vm_sub_res)
    ax[0].scatter(x_sub, y_sub, s=16, alpha=0.7, label="Exp (sub)")
    ax[0].plot(Ts_sub[m0], Vm_sub_eos[m0], lw=2, label="Model (EOS)")
    ax[0].plot(Ts_sub[m1], Vm_sub_res[m1], lw=2,
               ls="--", label="Model (Result)")
    ax[0].axvline(T_t, ls="--", lw=1, color="gray")
    ax[0].set_title("Vm along sublimation")
    ax[0].set_xlabel("T [K]")
    ax[0].set_ylabel(r"$V_m$ [cm$^3$/mol]")
    ax[0].legend()

    m2 = np.isfinite(Vm_melt_eos)
    m3 = np.isfinite(Vm_melt_res)
    ax[1].scatter(x_melt, y_melt, s=16, alpha=0.7, label="Exp (melt)")
    ax[1].plot(Ts_melt[m2], Vm_melt_eos[m2], lw=2, label="Model (EOS)")
    ax[1].plot(Ts_melt[m3], Vm_melt_res[m3], lw=2,
               ls="--", label="Model (Result)")
    ax[1].axvline(T_t, ls="--", lw=1, color="gray")
    ax[1].set_title("Vm along melting")
    ax[1].set_xlabel("T [K]")
    ax[1].legend()

    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    main()
