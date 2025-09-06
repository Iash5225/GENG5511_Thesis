import numpy as np
from cheval import cheval

# paste your cheval() here or import it
# from my_chebyshev import cheval


def Debye3(x):
    """
    Debye function of order 3:
        D3(x) = 3/x^3 ∫_0^x t^3 / (exp(t) - 1) dt, for x >= 0.
    Uses:
      - small-x series
      - Chebyshev (your cheval) on x ∈ [0,4] with t = x^2/8 - 1
      - large-x asymptotic with exp(-x) tail
    Returns NaN for x < 0 (domain).
    """
    x = np.asarray(x, dtype=np.float64)

    # constants
    zeta4 = (np.pi**4) / 90.0
    Ainf = 18.0 * zeta4               # coefficient in Ainf / x^3
    epsr = np.finfo(float).eps
    rmin = np.finfo(float).tiny

    # Chebyshev coefficients (ADEB3)
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

    # Trim negligible tail once (mirrors abs(a_k) > eps/100)
    tol = epsr / 100.0
    idx = np.where(np.abs(A) > tol)[0]
    if idx.size:
        A = A[:idx[-1] + 1]

    y = np.empty_like(x)

    # masks
    xneg = x < 0
    x0 = x == 0
    x_small = (x > 0) & (x < np.sqrt(8.0 * epsr))
    x_mid = (x >= np.sqrt(8.0 * epsr)) & (x <= 4.0)

    # large-x thresholds (avoid needless exp)
    x_upper = -np.log(2.0 * epsr)              # exp(-x) ~ eps
    x_lim1 = -np.log(rmin)                    # exp(-x) underflows
    x_big = (Ainf / rmin) ** (1.0 / 3.0)     # Ainf/x^3 underflows

    x_pure_asymp = (x > 4.0) & (x >= x_lim1)           # exp(-x) underflows
    x_asymp_tail = (x > 4.0) & (x < x_lim1) & (x >= x_upper)
    x_asymp_full = (x > 4.0) & (x < x_upper)
    x_underflowA = x > x_big

    # 1) domain
    y[xneg] = np.nan

    # 2) x=0
    y[x0] = 1.0

    # 3) small x: 1 - 3x/8 + x^2/20
    if np.any(x_small):
        xs = x[x_small]
        y[x_small] = 1.0 - (3.0/8.0)*xs + (1.0/20.0)*xs**2

    # 4) mid x in [0,4]: Chebyshev with t = x^2/8 - 1, then subtract 3/8 x
    if np.any(x_mid):
        xm = x[x_mid]
        t = (xm**2)/8.0 - 1.0
        y[x_mid] = cheval(A, t) - 0.375 * xm

    # 5) extremely huge x where Ainf/x^3 underflows -> 0
    y[x_underflowA] = 0.0

    # 6) large x: pure Ainf/x^3 (exp(-x) underflow)
    mask = x_pure_asymp & ~x_underflowA
    if np.any(mask):
        xa = x[mask]
        y[mask] = Ainf / (xa**3)

    # 7) large x with exp tail representable
    #    D3 ≈ Ainf/x^3 - 3 e^{-x} (x^3+3x^2+6x+6)/x^3
    for mask in (x_asymp_tail, x_asymp_full):
        mask = mask & ~x_underflowA
        if np.any(mask):
            xl = x[mask]
            with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
                expmx = np.exp(-xl)
                poly = ((xl + 3.0)*xl + 6.0)*xl + 6.0
                y[mask] = (Ainf / (xl**3)) - 3.0 * (poly / (xl**3)) * expmx

    return y
