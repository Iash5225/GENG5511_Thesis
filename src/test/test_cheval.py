import numpy as np
import pytest
# Replace with actual filename, e.g., from cheval import Cheval
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),"..")))
from cheval import Cheval


def test_zero_coefficients():
    """All coefficients zero => result should be zero"""
    a = np.zeros(5)
    T = 0.5
    result = Cheval(len(a), a, T)
    assert result == pytest.approx(0.0)


def test_single_coefficient():
    """Only a[0] non-zero => result should be a[0]"""
    a = np.array([3.0])
    T = 0.3
    result = Cheval(len(a), a, T)
    assert result == pytest.approx(3.0)


def test_linear_case():
    """Test T approximation when a[0] and a[1] are set"""
    a = np.array([1.0, 2.0])  # Should be 1 + 2*T for Chebyshev T0 + T1
    T = 0.4
    expected = 1 + 2*T
    result = Cheval(len(a), a, T)
    assert result == pytest.approx(expected, abs=1e-6)


def test_symmetry_even_terms():
    """Check even function property for even-degree-only coefficients"""
    a = np.array([1.0, 0.0, 1.0])  # T0 + T2: symmetric
    T = 0.3
    result_pos = Cheval(len(a), a, T)
    result_neg = Cheval(len(a), a, -T)
    assert result_pos == pytest.approx(result_neg, abs=1e-6)


def test_large_T_vs_numpy():
    """Compare Cheval with NumPy's Chebyshev evaluator at |T| >= 0.6"""
    a = np.array([1.0, 2.0, 3.0])
    T = 0.9
    result = Cheval(len(a), a, T)

    from numpy.polynomial.chebyshev import Chebyshev
    expected = Chebyshev(a)(T)

    assert result == pytest.approx(expected, abs=1e-6)
