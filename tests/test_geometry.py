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


def _one_crank(cheek_width=None, pin_width=None):
    return [ShaftComponent("cheek", 10, 0, 1, 1, width=cheek_width),
            ShaftComponent("crankpin", 30, 0, 1, 1, width=pin_width),
            ShaftComponent("cheek", 50, 0, 1, 1)]


def test_defaults_used_without_individual_widths():
    engine = EngineConfig(main_journal_width=20, crankpin_journal_width=30, cheek_width=10)
    length, _ = derive_crankshaft_dimensions(build_geometry(_one_crank()), engine)
    assert length == pytest.approx(10 + 10 + 30 + 10 + 10)


def test_individual_cheek_and_crankpin_widths_override_defaults():
    engine = EngineConfig(main_journal_width=20, crankpin_journal_width=30, cheek_width=10)
    g = build_geometry(_one_crank(cheek_width=15, pin_width=40))
    length, _ = derive_crankshaft_dimensions(g, engine)
    assert length == pytest.approx(10 + 15 + 40 + 10 + 10)


def test_individual_main_journal_widths(shaft):
    engine, components = shaft
    base = EngineConfig(main_journal_widths={"0": 47, "2": 30})
    length, _ = derive_crankshaft_dimensions(build_geometry(components), base)
    assert length == pytest.approx(524.0 + (47 - 37) / 2 + (30 - 37))


def test_main_journal_index_out_of_range_rejected():
    engine = EngineConfig(main_journal_widths={3: 10})
    with pytest.raises(ValueError):
        derive_crankshaft_dimensions(build_geometry(_one_crank()), engine)


def test_non_positive_width_rejected():
    with pytest.raises(ValueError):
        build_geometry(_one_crank(cheek_width=0))
