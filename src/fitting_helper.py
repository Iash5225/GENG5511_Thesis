from scipy.optimize import minimize
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from read import load_all_gas_data
from thermal_script import plot_all_gas_properties
from computethermoprops import *
from p_functions import pmelt, psub
from constants import CUSTOMCOLORS, CUSTOMMARKERS, PARAMS_INIT, LOWER_BOUND, UPPER_BOUND, KRYPTON_P_t, KRYPTON_T_t, KRYPTON_REFERENCE_ENTROPY, KRYPTON_REFERENCE_ENTHALPY, PARAM_LABELS
from scipy.integrate import quad

# solid_EOS_fitting.py



def rms(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    return 0.0 if n == 0 else np.sqrt(np.sum(x*x) / n)


def _mean_sq(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 0.0 if x.size == 0 else np.mean(x*x)

