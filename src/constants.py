import numpy as np

# Universal gas constant [J K^-1 mol^-1]
R = 8.31451

# Universal constants
ZERO = 0.0
HALF = 0.5
ONE = 1.0
TWO = 2.0
THREE = 3.0
FOUR = 4.0
EIGHT = 8.0
TWENTY = 20.0
PT375 = 0.375

W_VM_SUB = 0
W_VM_MELT = 0
W_VM_HIGHP = 0
W_CP_SUB = 0
W_ALPHA_SUB =0
W_BETAT_SUB = 0.0
W_BETAS_SUB = 0.0
W_H_SOLID_SUB = 1.0
W_H_SOLID_MELT = 1.0
W_P_SUB = 0.0
W_P_MELT = 0.0
W_GAMMA_T = 0.0


FUNCTION_TOL = 1e-15
GRADIENT_TOL = 1e-12
MAX_ITERATIONS = 1000
EPS = 1e-9
MAXLS = 100
N_OUTER = 1

READ_FROM_EXCEL = False  # Flag to read data from Excel files
# READ_FROM_TXT = True  # Flag to read data from text files
DISPLAY_PLOTS = True  # Flag to display plots

PARAM_LABELS = [
    ("c1", "MPa"),
    ("c2", "MPa"),
    ("c3", "MPa"),
    ("Theta_D,0", "K"),
    ("gamma_D,0", ""),
    ("q_D", ""),
    ("b1", ""),
    ("b2", ""),
    ("b3", ""),
    ("S_m(g, T_t, p_t)", "J mol-1 K-1"),
]

# Machine-dependent constants
SMALLEST_POSITIVE = 2.2251E-308  # Smallest positive number
SMALLEST_SPACING = 1.11E-16  # Smallest relative spacing
MAX_LIMIT = np.log(np.sqrt(np.finfo(float).max)) - 1  # Limit for large X

# REFPROP NIST Constants
KRYPTON_REFERENCE_ENTROPY = 76.99 # kJ/kgK
KRYPTON_REFERENCE_ENTHALPY = 8911/1000  #J/kg

NEON_REFERENCE_ENTROPY = 69.25 # kJ/kgK
NEON_REFERENCE_ENTHALPY = 1695/1000  # J/kg

XENON_REFERENCE_ENTROPY = 77.43 # kJ/kgK
XENON_REFERENCE_ENTHALPY = 12494.09/1000  # J/kg


# Constants for melting pressure equation
KRYPTON_P_t = 0.072954 # MPa
KRYPTON_T_t = 115.8     # K
XENON_P_t = 0.08177
XENON_T_t = 161.405
NEON_P_t = 0.043332
NEON_T_t = 24.5561

# krypton
KRYPTON_E_4 = 1487.56
KRYPTON_E_5 = 1.819
KRYPTON_E_6 = 5158.09
KRYPTON_E_7 = 0.98822

XENON_E_4 = 4136.681829018883
XENON_E_5 = 1.280037747122412
XENON_E_6 = 2209.571430062904
XENON_E_7 = 0.8558320170107692

NEON_E_4 = 4151.27
NEON_E_5 = 1.355
NEON_E_6 = 558.60
NEON_E_7 = 0.16851


# Sublimation pressure constants
KRYPTON_E_1_SUB = -9.975928773420845
KRYPTON_E_2_SUB = -2.8993773847096307
KRYPTON_E_3_SUB = 2.2848747077114253

XENON_E_1_SUB = -18.910506919795417
XENON_E_2_SUB = 11.285160377819482
XENON_E_3_SUB = -12.280765945817006

NEON_E_1_SUB = -14.883845098313778
NEON_E_2_SUB = 5.821759117974304
NEON_E_3_SUB = -0.39834285015300536

gas_params = {
    'krypton': (KRYPTON_T_t, KRYPTON_P_t),
    'xenon': (XENON_T_t, XENON_P_t),
    'neon': (NEON_T_t, NEON_P_t)
}
IMG_OUTPUT_FOLDER = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\output\v3"
#FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\literature_data.xlsx"
OUTPUT_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\output\v2"
XE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Xenon Literature Review.xlsx"
KR_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Krypton Literature Review.xlsx"
NE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Neon Literature Review.xlsx"
TXT_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\txt"

MELTING_SHEETNAME = "Melting"
# Replace with the actual sheet name
THERMAL_COEFFICIENT_SHEETNAME = "Thermal Expansion Coefficient"
HEAT_CAPACITY_SHEETNAME = "Heat Capacity"
CELL_VOLUME_SHEETNAME = "Cell Volume "
BULK_MODULUS_SHEETNAME = "Compressibility Bulk Modulus"
SUBLIMATION_SHEETNAME = "Sublimation"
FUSION_SHEETNAME = "Heat of Fusion"
HEAT_OF_SUBLIMATION_SHEETNAME = "Heat of Sublimation"

EOS_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\eos"

# Define custom colors (normalized, excluding white)
CUSTOMCOLORS = np.array([
    [0, 0, 0], [112, 48, 160], [192, 0, 0], [1, 175, 146], [222, 110, 38],
    [0, 0, 255], [150, 150, 150], [95, 58, 91], [72, 113, 57], [27, 71, 116],
    [222, 110, 38], [139, 44, 42], [0, 200, 0], [255, 0, 240], [92, 103, 177],
    [71, 30, 118], [100, 200, 0], [239, 144, 42], [120, 100, 255], [55, 200, 80],
    [200, 20, 150], [25, 105, 88], [88, 10, 198], [100, 55, 241], [254, 120, 62],
    [165, 158, 171], [224, 21, 138], [
        155, 100, 8], [84, 184, 93], [193, 233, 41],
    [250, 199, 64], [200, 175, 41], [127, 217, 16], [255, 0, 0], [0, 0, 255],
    [95, 58, 91], [72, 113, 57], [27, 71, 116], [222, 110, 38], [139, 44, 42],
    [0, 200, 0], [255, 0, 240], [92, 103, 177], [71, 30, 118], [100, 200, 0],
    [239, 144, 42], [120, 100, 255], [
        55, 200, 80], [200, 20, 150], [25, 105, 88],
    [88, 10, 198], [100, 55, 241], [254, 120, 62], [
        165, 158, 171], [224, 21, 138],
    [155, 100, 8], [84, 184, 93], [193, 233, 41], [250, 199, 64], [200, 175, 41],
    [127, 217, 16]
]) / 256  # Normalize to range [0,1]

# Define only filled markers
CUSTOMMARKERS = ['o', 's', 'D', '^', 'v', '>', '<', 'p', 'h', '*']

MELTING_DEVIATION_Y_MIN=-20
MELTING_DEVIATION_Y_MAX=20

# Axis Limits
MELTING_KRYPTON_X_MIN=100
MELTING_KRYPTON_X_MAX=300
MELTING_KRYPTON_Y_MIN=1
MELTING_KRYPTON_Y_MAX = 10**3

MELTING_NEON_X_MIN = 0
MELTING_NEON_X_MAX = 350
MELTING_NEON_Y_MIN = 10**1
MELTING_NEON_Y_MAX = 10**4


MELTING_XENON_X_MIN = 150
MELTING_XENON_X_MAX = 300
MELTING_XENON_Y_MIN = 1
MELTING_XENON_Y_MAX = 10**3

MARKERSIZE=50
AXIS_FONT_SIZE=14


PARAMS_INIT_NEON = np.array([
    13.38,   # 0: v00
    997.79,   # 1: a1
    469.13,    # 2: a2
    114.76,     # 3: a3
    0, 0, 0, 0, 0,           # 4–8 unused
    65, 0, 0, 0, 0, 0,     # 9–14: Theta_D,0 (smaller -> more T effect)
    2.69,  0, 0, 0, 0, 0,     # 15–20: gamma_D,0 (bigger -> more α)
    0.2,  0, 0, 0, 0, 0,     # 21–26: q_D (moderate)
    0, 0, 0,           # 27–29: aa, bb, cc (keep 0 for Vm-only)
    NEON_REFERENCE_ENTROPY
])


LOWER_BOUND_NEON = np.array([
    13.00,   # 0: v00
    997.00,   # 1: a1
    469.00,    # 2: a2
    114.00,     # 3: a3
    0, 0, 0, 0, 0,           # 4–8 unused
    64, 0, 0, 0, 0, 0,     # 9–14: Theta_D,0 (smaller -> more T effect)
    2.65,  0, 0, 0, 0, 0,     # 15–20: gamma_D,0 (bigger -> more α)
    0.1,  0, 0, 0, 0, 0,     # 21–26: q_D (moderate)
    -1, -1,-1,           # 27–29: aa, bb, cc (keep 0 for Vm-only)
    NEON_REFERENCE_ENTROPY
])

UPPER_BOUND_NEON = np.array([
    15,     # v00
    1100,    # a1
    600,    # a2
    300,     # a3
    0, 0, 0, 0, 0,
    75, 0, 0, 0, 0, 0,   # ΘD0
    3,  0, 0, 0, 0, 0,   # γ0
    2,  0, 0, 0, 0, 0,   # q0
    100, 100, 100,
    NEON_REFERENCE_ENTROPY
])
# ...existing code...
# Xenon
# PARAMS_INIT_XENON = np.array([
#     34.85, 2739.48, 7328.48, 122.62,
#     0, 0, 0, 0, 0,
#     79.69, 0, 0, 0, 0, 0,
#     2.62, 0, 0, 0, 0, 0,
#     -1.54, 0, 0, 0, 0, 0,
#     0.0, 5.89, 6.16,
#     XENON_REFERENCE_ENTROPY  # S* (entropy reference)
# ])

PARAMS_INIT_XENON = np.array([
    34.53, 3386.08, 7310.03, 123.09,
    0, 0, 0, 0, 0,
    58.44, 0, 0, 0, 0, 0,
    2.68, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0.0, 5.89, 6.16,
    XENON_REFERENCE_ENTROPY  # S* (entropy reference)
])

# Lower bounds
LOWER_BOUND_XENON = np.array([
    30, 0, 0, 0, 0, 0, 0, 0, 0,
    20, 0, 0, 0, 0, 0,
    -15, 0, 0, 0, 0, 0,
    -10, 0, 0, 0, 0, 0,
    0, 0, 5, 5
])

# Upper bounds
UPPER_BOUND_XENON = np.array([
    40, 10000, 10000, 200, 0, 0, 0, 0, 0,
    90, 0, 0, 0, 0, 0,
    15, 0, 0, 0, 0, 0,
    10, 0, 0, 0, 0, 0,
    2, 6,7, 200
])

# Upper bounds
PARAMS_INIT = np.array([
    27.13, 3087.52, 7316.55, 123.13,
    0, 0, 0, 0, 0,
    67.74, 0, 0, 0, 0, 0,
    2.64, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0.0, 1.93, 7.93,
    69.13165833  # S* (entropy reference)
])
# UPPER_BOUND = np.array([
#     40, 10000, 10000, 10000, 0, 0, 0, 0, 0,
#     300, 0, 0, 0, 0, 0,
#     5, 0, 0, 0, 0, 0,
#     10, 0, 0, 0, 0, 0,
#     100, 100, 100, 200
# ])

# # Lower bounds
# LOWER_BOUND = np.array([
#     20, 0, 0, 0, 0, 0, 0, 0, 0,
#     20, 0, 0, 0, 0, 0,
#     -5, 0, 0, 0, 0, 0,
#     -10, 0, 0, 0, 0, 0,
#     0, 0, 0, 60
# ])


# Lower bounds
LOWER_BOUND = np.array([
    27.13, 3087.52, 7316.55, 123.13,
    0, 0, 0, 0, 0,
    67.74, 0, 0, 0, 0, 0,
    2.64, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0.0, 1.93, 7.93,
    69  # S* (entropy reference)
])
UPPER_BOUND = np.array([
    27.13, 3087.52, 7316.55, 123.13,
    0, 0, 0, 0, 0,
    67.74, 0, 0, 0, 0, 0,
    2.64, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
    0.0, 1.93, 7.93,
    75  # S* (entropy reference)
])

GAMMA_POS_SLOPE_OFFSET = 45.0
GAMMA_POS_SLOPE_MULT =200.0
GAMMA_NEG_SLOPE_MULT = 50.0

PERCENT_SCALE = 100.0

PERCENT_SCALE = 100.0


# --- Weights for each deviation term ---
PMELT_EXTRA_WEIGHT_T_K = 500.0
PMELT_EXTRA_FACTOR = 5.0
CP_TEMP_THRESHOLD_K = 40
CP_WEIGHT_BELOW = 1
CP_WEIGHT_ABOVE = 1




history = {
    "total": [],
    "Vm_sub": [],
    "Vm_melt": [],
    "Vm_highp": [],
    "cp_sub": [],
    "alpha_sub": [],
    "BetaT_sub": [],
    "BetaS_sub": [],
    "H_solid_sub": [],
    "H_solid_melt": [],
    "p_sub": [],
    "p_melt": [],
    "Gamma_T": []
}


# ---------- Physical-constraint knobs ----------
CP_SPLIT_K = 25.0           # K: split between "low" and "high" heat-capacity weighting
CP_W_BELOW = 3.0             # stronger weight where cp -> 0
CP_W_ABOVE = 1.0

PMELT_T_HI_K = 500.0         # K: emphasize melting pressures at high T
# extra weight multiplier for T >= 500 K (on top of log-RMS)
PMELT_W_HI = 4.0
GAMMA_DEBYE_T_MAX = 30.0     # K: range over which gamma -> gamma_D0
W_GAMMA_LOW_T = 3.0      # weight for (gamma - gamma_D0) penalty in 0–30 K

W_GAMMA_MONO = 3.5      # already in your code; keep or tune

