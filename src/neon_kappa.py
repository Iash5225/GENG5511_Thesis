# python
import numpy as np
from thermopropsv2 import compute_thermo_props
from constants import PARAMS_INIT


def hand_target(T0, p0, v0, BetaT_exp, BetaS_exp=None):
    Kt_exp = 1.0 / BetaT_exp
    dpdv_needed = -1.0 / (v0 * Kt_exp)
    kout = {"KappaT_exp": Kt_exp, "dpdv_needed": dpdv_needed}
    if BetaS_exp is not None:
        Ks_exp = 1.0 / BetaS_exp
        kout["cv_over_cp_needed"] = Ks_exp / Kt_exp
    return kout


def model_point(params, T0, p0):
    props = compute_thermo_props(float(T0), float(p0), params)
    # indices: KappaT=1, KappaS=2, cp=4, cv=5, v=0
    v, Kt, Ks, _, cp, cv = props[0], props[1], props[2], props[3], props[4], props[5]
    dpdv_model = -1.0 / (v * Kt)
    return dict(v=float(v), KappaT=float(Kt), KappaS=float(Ks), cp=float(cp), cv=float(cv), dpdv=dpdv_model)


def sensitivity(params, T0, p0, deltas=None):
    if deltas is None:
        deltas = {"Th0": 1e-2, "g0": 1e-3, "q0": 1e-3,
                  "aa": 1e-4, "bb": 1e-3, "cc": 1e-3}
    # param indices in your layout:
    idx = {"Th0": 9, "g0": 15, "q0": 21, "aa": 27, "bb": 28, "cc": 29}
    base = model_point(params, T0, p0)
    sens = {}
    for name, delta in deltas.items():
        p2 = params.copy()
        p2[idx[name]] += delta
        m2 = model_point(p2, T0, p0)
        sens[name] = {
            k: (m2[k] - base[k]) / delta for k in ("KappaT", "KappaS", "dpdv", "cv", "cp")}
    return base, sens


# Example usage:
T0, p0, v0, BetaT_exp, BetaS_exp = 20.0, 0.002770165, 13.391, 0.00139, 0.00116
print(hand_target(T0, p0, v0, BetaT_exp, BetaS_exp))
base, sens = sensitivity(PARAMS_INIT, T0, p0)
print(base)
for k,v in sens.items(): print(k, v)
