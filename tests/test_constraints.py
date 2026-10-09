from pathlib import Path

import numpy as np
import pytest

from crankshaft_balancing import (CounterweightConstraints, analyze,
                                  load_constraints, resolve_constraints,
                                  verify_kkt_conditions,
                                  verify_optimum_via_duality)
from crankshaft_balancing.cli import main

ROOT = Path(__file__).resolve().parent.parent
CONSTRAINED = ROOT / "examples" / "example_shaft_constrained.json"
W = None  # заполняется из engine в тестах


def run(shaft, constraints, step=0.05):
    engine, components = shaft
    return analyze(engine, components, step=step, verbose=False, constraints=constraints)


def balance_and_end_moment(res, axis="y"):
    f, r, w = res.forces, getattr(res, f"result_{axis}"), res.engine.angular_speed
    Fc, Fp = ((f.Fy_cheeks, f.Fy_crankpins) if axis == "y"
              else (f.Fx_cheeks, f.Fx_crankpins))
    return Fc.sum() + Fp.sum() + r["s"].sum() * w ** 2, r["M"][-1] * 1e-3


@pytest.fixture(scope="module")
def free(shaft):
    return run(shaft, None)


def test_zero_removes_counterweights_and_keeps_balance(shaft, free):
    res = run(shaft, {"all": CounterweightConstraints.from_dict({"zero": [2, 5]})})
    s = res.result_y["s"]
    assert s[2] == 0 and s[5] == 0
    force, m_end = balance_and_end_moment(res)
    assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)
    # ограничение не может улучшить оптимум
    assert res.max_moment_y >= free.max_moment_y - 1e-6


def test_manual_value_is_fixed_and_others_compensate(shaft):
    res = run(shaft, {"y": CounterweightConstraints.from_dict({"fixed": {"0": 0.1}})})
    assert res.result_y["s"][0] == pytest.approx(0.1)
    force, m_end = balance_and_end_moment(res)
    assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)


def test_bounds_are_respected(shaft):
    res = run(shaft, {"y": CounterweightConstraints.from_dict(
        {"bounds": {"0": [-0.05, 0.05], "7": [None, 0.0]}})})
    s = res.result_y["s"]
    assert -0.05 - 1e-9 <= s[0] <= 0.05 + 1e-9
    assert s[7] <= 1e-9


def test_axis_specific_overrides_all():
    c = resolve_constraints(
        {"all": CounterweightConstraints.from_dict({"fixed": {"1": 0.2}}),
         "y": CounterweightConstraints.from_dict({"fixed": {"1": 0.3}})}, "y")
    assert c.fixed == {1: 0.3}
    c_x = resolve_constraints(
        {"all": CounterweightConstraints.from_dict({"fixed": {"1": 0.2}}),
         "y": CounterweightConstraints.from_dict({"fixed": {"1": 0.3}})}, "x")
    assert c_x.fixed == {1: 0.2}


def test_certificates_hold_with_constraints(shaft):
    res = run(shaft, {"y": CounterweightConstraints.from_dict(
        {"zero": [2], "bounds": {"0": [-0.15, 0.15]}})})
    assert verify_optimum_via_duality(res.result_y, verbose=False)["confirmed"]
    assert verify_kkt_conditions(res.result_y, verbose=False)["confirmed"]


def test_too_few_free_coefficients_is_infeasible(shaft):
    # одна свободная переменная при двух условиях (сила и момент) — решения нет
    fixed = {str(i): 0.05 * (i + 1) for i in range(7)}
    with pytest.raises(RuntimeError):
        run(shaft, {"y": CounterweightConstraints.from_dict({"fixed": fixed})})


def test_all_fixed_just_evaluates_given_values(shaft):
    fixed = {str(i): 0.1 for i in range(8)}
    res = run(shaft, {"y": CounterweightConstraints.from_dict({"fixed": fixed})})
    assert res.result_y["result"] is None
    assert np.allclose(res.result_y["s"], 0.1)


@pytest.mark.parametrize("data, exc", [
    ({"zero": [8]}, IndexError),
    ({"zero": [-1]}, IndexError),
    ({"fixed": {"0": 0.1}, "bounds": {"0": [0, 1]}}, ValueError),
    ({"bounds": {"0": [1, 0]}}, ValueError),
    ({"zero": [0], "fixed": {"0": 0.5}}, ValueError),
    ({"unknown": 1}, ValueError),
])
def test_invalid_constraints_rejected(shaft, data, exc):
    with pytest.raises(exc):
        run(shaft, {"y": CounterweightConstraints.from_dict(data)})


def test_load_constraints_from_json():
    c = load_constraints(CONSTRAINED)
    assert c["all"].fixed == {2: 0.0, 5: 0.0}
    assert c["y"].bounds[0] == (-0.15, 0.15)


def test_cli_with_constrained_example(tmp_path):
    assert main([str(CONSTRAINED), "--out", str(tmp_path), "--step", "0.1"]) == 0


def test_axis_fixed_overrides_all_bounds_without_error(shaft):
    # "all" ограничивает щеку 0 диапазоном, а ось y фиксирует её — побеждает fixed
    res = run(shaft, {"all": CounterweightConstraints.from_dict({"bounds": {"0": [-1, 1]}}),
                      "y": CounterweightConstraints.from_dict({"fixed": {"0": 0.1}})})
    assert res.result_y["s"][0] == pytest.approx(0.1)


def test_custom_counterweight_positions(shaft):
    engine, components = shaft
    xs = [31.469, 99.531, 262.0, 424.469, 492.531]   # три на щёках + одна между колен
    res = analyze(engine, components, step=0.05, verbose=False, counterweight_positions=xs)
    assert len(res.result_y["s"]) == 5
    assert np.allclose(res.L_counterweights, engine.crankshaft_length - np.array(sorted(xs)))
    force, m_end = balance_and_end_moment(res)
    assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)


def test_positions_out_of_shaft_rejected(shaft):
    engine, components = shaft
    with pytest.raises(ValueError):
        analyze(engine, components, step=0.5, verbose=False,
                counterweight_positions=[10.0, 9999.0])


def test_constraint_index_checked_against_positions(shaft):
    engine, components = shaft
    with pytest.raises(IndexError):   # позиций 3, а номер 5
        analyze(engine, components, step=0.5, verbose=False,
                counterweight_positions=[31.469, 262.0, 492.531],
                constraints={"all": CounterweightConstraints.from_dict({"zero": [5]})})


def test_cli_sensitivity_and_convergence(tmp_path):
    assert main([str(CONSTRAINED), "--out", str(tmp_path), "--step", "0.1",
                 "--sensitivity", "--convergence"]) == 0
    names = {p.name for p in tmp_path.iterdir()}
    assert "sensitivity_y.png" in names
    assert "convergence_y.png" in names or "convergence_x.png" in names
