import random

import numpy as np
import pytest

from crankshaft_balancing import (EngineConfig, ShaftComponent, build_geometry,
                                  derive_crankshaft_dimensions)


def test_order_of_input_is_irrelevant(shaft):
    _, components = shaft
    ref = build_geometry(components)
    shuffled = components.copy()
    random.Random(0).shuffle(shuffled)
    g = build_geometry(shuffled)
    assert np.allclose(ref.cheek_centers, g.cheek_centers)
    assert np.allclose(ref.crankpin_centers, g.crankpin_centers)


def test_unknown_kind_rejected():
    with pytest.raises(ValueError):
        build_geometry([ShaftComponent("bolt", 0, 0, 0, 1)])


def test_example_length_and_pitch(shaft):
    engine, components = shaft
    length, pitch = derive_crankshaft_dimensions(build_geometry(components), engine)
    assert length == pytest.approx(524.0)
    assert pitch == pytest.approx(131.0)


def test_empty_geometry_rejected():
    with pytest.raises(ValueError):
        derive_crankshaft_dimensions(build_geometry([]), EngineConfig())
