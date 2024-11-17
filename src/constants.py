import numpy as np

# Universal gas constant [J K^-1 mol^-1]
R = 8.31451

# Universal constants
ZERO = 0.0
HALF = 0.5
ONE = 1.0
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
