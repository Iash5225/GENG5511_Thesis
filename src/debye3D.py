import numpy as np
# from constants import ZERO, PT375, HALF, ONE, THREE, FOUR, EIGHT, DEBINF, ADEB3, SMALLEST_POSITIVE, SMALLEST_SPACING
from cheval import cheval


import numpy as np

# assumes your cheval(a, T) is defined above / imported


def Debye3D(x):
    """
    First derivative of the Debye-3 function:
        D3'(x) = d/dx [ 3/x^3 ∫_0^x t^3/(exp(t)-1) dt ],  for x >= 0.

    Strategy:
      - small x: series  (-3/8) + x/10
      - mid x in [0,4]: Chebyshev derivative with chain rule: (x/4) * cheval(D, t) - 3/8
                         where t = x^2/8 - 1, and D are Chebyshev coeffs of d/dt of series
      - large x: asymptotic -3*Ainf/x^4 with exp(-x) tail when representable
      - x < 0: NaN (domain)
    """
    x = np.asarray(x, dtype=np.float64)

    # constants
    zeta4 = (np.pi**4) / 90.0
    Ainf = 18.0 * zeta4            # 1/DEBINF in your MATLAB
    epsr = np.finfo(float).eps
    rmin = np.finfo(float).tiny

    # Chebyshev coefficients for Debye3 on t ∈ [-1,1] (your ADEB3)
    A = np.array([
        2.70773706832744,
        0.340068135211092,
        -1.29451501844409e-02,
        7.96375538017382e-04,
        -5.46360009590824e-05,
        3.92430195988049e-06,
        -2.89403282353860e-07,
        2.17317613962500e-08,
        -1.65420999498000e-09,
        1.27279618920000e-10,
        -9.87963459000000e-12,
        7.72507400000000e-13,
        -6.07797200000000e-14,
        4.80759000000000e-15,
        -3.82040000000000e-16,
        3.04800000000000e-17,
        -2.44000000000000e-18,
        2.00000000000000e-19,
        -2.00000000000000e-20
    ], dtype=np.float64)

    # Trim negligible tail as in MATLAB (abs(a_k) > eps/100)
    tol = epsr / 100.0
    last = np.where(np.abs(A) > tol)[0]
    if last.size:
        A = A[:last[-1] + 1]

    # Build Chebyshev coefficients for derivative wrt t (your D array)
    # If A corresponds to degrees k=0..n, then:
    #   D[n] = 0, D[n-1] = 2n * A[n], D[k-1] = D[k+1] + 2k * A[k]
    n = len(A) - 1
    D = np.zeros(n+1, dtype=np.float64)
    if n >= 1:
        D[-2] = 2.0*n * A[-1]   # D[n-1]
        # fill downward: k = n-1 .. 1  (index in A is k)
        for k in range(n-1, 0, -1):
            # D[k-1] = D[k+1] + 2k * A[k]
            D[k-1] = (D[k+1] if (k+1) <= n else 0.0) + 2.0*k * A[k]
    # (D[0]..D[n-1]) valid; D[n] is unused / zero

    y = np.empty_like(x)

    # masks
    xneg = x < 0
    x0 = x == 0
    x_small = (x > 0) & (x < np.sqrt(8.0 * epsr))
    x_mid = (x >= np.sqrt(8.0 * epsr)) & (x <= 4.0)

    # large-x thresholds
    x_upper = -np.log(2.0 * epsr)            # exp(-x) ~ eps
    x_lim1 = -np.log(rmin)                  # exp(-x) underflows
    # where Ainf/x^4 underflows toward zero
    x_big = (Ainf / rmin)**(1.0/4.0)

    x_pure_asymp = (x > 4.0) & (x >= x_lim1)             # exp(-x) underflows
    x_asymp_tail = (x > 4.0) & (x < x_lim1) & (x >= x_upper)
    x_asymp_full = (x > 4.0) & (x < x_upper)
    x_underflowA = x > x_big

    # 1) domain
    y[xneg] = np.nan

    # 2) x=0: limit of derivative is -3/8
    y[x0] = -0.375

    # 3) small x: series  (-3/8) + x/10
    if np.any(x_small):
        xs = x[x_small]
        y[x_small] = -0.375 + 0.1*xs

    # 4) mid x: (x/4) * cheval(D, t) - 3/8, with t = x^2/8 - 1
    if np.any(x_mid):
        xm = x[x_mid]
        t = (xm**2)/8.0 - 1.0
        y[x_mid] = 0.25 * xm * cheval(D, t) - 0.375

    # 5) extremely huge x where Ainf/x^4 underflows -> 0
    y[x_underflowA] = 0.0

    # 6) large x: pure asymptotic -3*Ainf/x^4
    mask = x_pure_asymp & ~x_underflowA
    if np.any(mask):
        xa = x[mask]
        y[mask] = -3.0 * Ainf / (xa**4)

    # 7) large x with exp(-x) tail representable
    #    D3' ≈ -3*Ainf/x^4 + 3*e^{-x} * Q(x)/x^4,
    #    where Q(x) = x^4 + 4x^3 + 9x^2 + 18x + 18  (from differentiating the tail)
    #    Your MATLAB implements this via SUM polynomials; we write it directly.
    for mask in (x_asymp_tail, x_asymp_full):
        mask = mask & ~x_underflowA
        if np.any(mask):
            xl = x[mask]
            with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
                expmx = np.exp(-xl)
                Q = (((xl + 4.0)*xl + 9.0)*xl + 18.0) * \
                    xl + 18.0  # x^4+4x^3+9x^2+18x+18
                y[mask] = (-3.0*Ainf / (xl**4)) + 3.0 * expmx * (Q / (xl**4))

    return y
