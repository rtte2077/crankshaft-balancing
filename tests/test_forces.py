import numpy as np
import pytest

from crankshaft_balancing import (centrifugal_force_magnitude,
                                  centrifugal_force_vectors, radial_distance)


def test_force_magnitude_units():
    # 1 кг на радиусе 1000 мм (=1 м) при omega=2 рад/с -> 4 Н
    assert centrifugal_force_magnitude(np.array([1000.0]), np.array([1.0]), 2.0)[0] == pytest.approx(4.0)


def test_vectors_point_along_center():
    centers = np.array([[0, 0, 3.0], [0, -4.0, 0]])
    v = centrifugal_force_vectors(centers, np.array([10.0, 5.0]))
    assert np.allclose(v, [[0, 10], [-5, 0]])


def test_radial_distance_ignores_x():
    assert radial_distance(np.array([[100.0, 3.0, 4.0]]))[0] == pytest.approx(5.0)
