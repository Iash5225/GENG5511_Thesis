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

# Machine-dependent constants
SMALLEST_POSITIVE = 2.2251E-308  # Smallest positive number
SMALLEST_SPACING = 1.11E-16  # Smallest relative spacing
MAX_LIMIT = np.log(np.sqrt(np.finfo(float).max)) - 1  # Limit for large X

# Specific constants for Debye calculations
DEBINF = 5.13299112734217E-02

# Chebyshev coefficients for Debye3 and derivatives
ADEB3 = [
    2.70773706832744, 0.340068135211092, -1.29451501844409E-02, 7.96375538017382E-04,
    -5.46360009590824E-05, 3.92430195988049E-06, -
    2.8940328235386E-07, 2.173176139625E-08,
    -1.65420999498E-09, 1.2727961892E-10, -9.87963459E-12, 7.725074E-13,
    -6.077972E-14, 4.80759E-15, -3.8204E-16, 3.048E-17, -2.44E-18, 2E-19, -2E-20
]

# Constants for melting pressure equation

# krypton
KRYPTON_E_4 = 1487.56
KRYPTON_E_5 = 1.819
KRYPTON_E_6 = 5158.09
KRYPTON_E_7 = 0.98822
KRYPTON_P_t = 0.072954  
KRYPTON_T_t = 115.8

# XENON_E_4 = 4053.39
# XENON_E_5 = 1.291
# XENON_E_6 = 2231.27
# XENON_E_7 = 0.84562
XENON_P_t = 0.08177
XENON_T_t = 161.405

XENON_E_4 = 4136.681829018883
XENON_E_5 = 1.280037747122412
XENON_E_6 = 2209.571430062904
XENON_E_7 = 0.8558320170107692


NEON_E_4 = 4151.27
NEON_E_5 = 1.355
NEON_E_6 = 558.60
NEON_E_7 = 0.16851
NEON_P_t = 0.043332
NEON_T_t = 24.5561

gas_params = {
    'krypton': (KRYPTON_T_t, KRYPTON_P_t),
    'xenon': (XENON_T_t, XENON_P_t),
    'neon': (NEON_T_t, NEON_P_t)
}

#FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\literature_data.xlsx"
OUTPUT_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\output\v2"
XE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Xenon Literature Review.xlsx"
KR_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Krypton Literature Review.xlsx"
NE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Neon Literature Review.xlsx"

MELTING_SHEETNAME = "Melting"
# Replace with the actual sheet name
THERMAL_COEFFICIENT_SHEETNAME = "Thermal Expansion Coefficient"
HEAT_CAPACITY_SHEETNAME = "Heat Capacity"
CELL_VOLUME_SHEETNAME = "Cell Volume"
BULK_MODULUS_SHEETNAME = "Compressibility Bulk Modulus"
SUBLIMATION_SHEETNAME = "Sublimation"



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

