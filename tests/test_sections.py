import numpy as np

from crankshaft_balancing import section_grid, total_distribution


def test_grid_includes_end():
    z = section_grid(10.0, step=1.0)
    assert z[0] == 0 and z[-1] >= 10.0


def test_single_force_gives_linear_moment():
    z = np.array([0.0, 1.0, 2.0, 5.0])
    m = total_distribution(z, [1.0], [10.0], [], [])
    assert np.allclose(m, [0, 0, 10, 40])
