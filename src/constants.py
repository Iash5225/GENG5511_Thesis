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
# E_4 = 1506.54
# E_5 = 1.731
# E_6 = 4677.16
# E_7 = 0.98493

# krypton
E_4 = 1411.48
E_5 = 1.851
E_6 = 5233.98
E_7 = 0.99137
P_t = 0.072954  # Reference pressure, adjust as needed
T_t = 115.8  # Reference temperature, adjust as needed

FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\data\literature_data.xlsx"
OUTPUT_FILEPATH = r"C:\Users\iashb\OneDrive - The University of Western Australia\UWA\05. Year 5\Semester 1\GENG5511 MPE Engineering Research Project\Project\Research_Project\GENG5511_Thesis\img\output"

MELTING_SHEETNAME = "melting"
# Replace with the actual sheet name
THERMAL_COEFFICIENT_SHEETNAME = "thermal expansion coefficient"


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
