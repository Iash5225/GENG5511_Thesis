from src.debye import Debye
from src.debye3 import Debye3
from src.debye3D import Debye3D


def test_debye():
    assert abs(Debye(3.0) - 8.846913398) < 1e-6
    assert abs(Debye(0.1) - 0.400663) < 1e-6
    assert abs(Debye(1e10) - 1.2e10) < 1e0


def test_debye3():
    assert abs(Debye3(3.0) - 0.531914) < 1e-6
    assert abs(Debye3(0.1) - 0.33367) < 1e-5
    assert abs(Debye3(1e10)) < 1e-6


def test_debye3d():
    assert abs(Debye3D(3.0) - -0.0628) < 1e-4
    assert abs(Debye3D(0.1) - -0.7416) < 1e-4
    assert abs(Debye3D(1e10)) < 1e-6


def test_edge_cases():
    assert abs(Debye(0.0)) < 1e-6
    assert abs(Debye3(0.0) - 1.0) < 1e-6
    assert abs(Debye3D(0.0)) < 1e-6
