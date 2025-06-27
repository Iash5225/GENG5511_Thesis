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

READ_FROM_EXCEL = True  # Flag to read data from Excel files
READ_FROM_TXT = False  # Flag to read data from text files
DISPLAY_PLOTS = True  # Flag to display plots


# Machine-dependent constants
SMALLEST_POSITIVE = 2.2251E-308  # Smallest positive number
SMALLEST_SPACING = 1.11E-16  # Smallest relative spacing
MAX_LIMIT = np.log(np.sqrt(np.finfo(float).max)) - 1  # Limit for large X

# Constants for melting pressure equation

# krypton
KRYPTON_E_4 = 1487.56
KRYPTON_E_5 = 1.819
KRYPTON_E_6 = 5158.09
KRYPTON_E_7 = 0.98822
KRYPTON_P_t = 0.072954  
KRYPTON_T_t = 115.8

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

#FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\literature_data.xlsx"
OUTPUT_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\output\v2"
XE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Xenon Literature Review.xlsx"
KR_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Krypton Literature Review.xlsx"
NE_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\MASTER Neon Literature Review.xlsx"
TXT_DATA_FILEPATH = r"C:\Users\iashb\OneDrive - UWA\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\txt"


MELTING_SHEETNAME = "Melting"
# Replace with the actual sheet name
THERMAL_COEFFICIENT_SHEETNAME = "Thermal Expansion Coefficient"
HEAT_CAPACITY_SHEETNAME = "Heat Capacity"
CELL_VOLUME_SHEETNAME = "Cell Volume"
BULK_MODULUS_SHEETNAME = "Compressibility Bulk Modulus"
SUBLIMATION_SHEETNAME = "Sublimation"
FUSION_SHEETNAME = "Heat of Fusion"
HEAT_OF_SUBLIMATION_SHEETNAME = "Heat of Sublimation"



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

