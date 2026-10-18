import csv
import json
from pathlib import Path

import numpy as np
import pytest

from crankshaft_balancing import (CounterweightConstraints, analyze, load_smooth_percent,
                                  optimize_counterweights, smooth_sweep,
                                  verify_kkt_conditions, verify_optimum_via_duality,
                                  with_group2_percent, with_smooth_percent)
from crankshaft_balancing.cli import main
from crankshaft_balancing.sections import NMM_TO_NM

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "example_shaft.json"
PERCENTS = list(range(0, 101, 10))


@pytest.fixture(scope="module")
def base(shaft):
    engine, components = shaft
    return analyze(engine, components, step=0.05, verbose=False)


@pytest.fixture(scope="module")
def sweep(base):
    out, levels = {}, base
    for q in PERCENTS:
        out[q] = with_smooth_percent(levels, q)
        if q > 0:
            levels = out[q]
    return out


def balance_and_end_moment(res, axis):
    f, w = res.forces, res.engine.angular_speed
    Fc, Fp = ((f.Fy_cheeks, f.Fy_crankpins) if axis == "y"
              else (f.Fx_cheeks, f.Fx_crankpins))
    r = getattr(res, f"result_{axis}")
    return Fc.sum() + Fp.sum() + r["s"].sum() * w ** 2, r["M"][-1] * NMM_TO_NM


def test_default_is_unsmoothed(base):
    assert base.smooth_percent == 0.0
    r = base.result_y
    assert np.isnan(r["moment_range_min"])
    assert r["moment_range"] == pytest.approx(np.ptp(r["M"]))
    assert r["max_moment_unsmoothed"] == pytest.approx(np.max(np.abs(r["M"])))
    assert base.max_moment_y == pytest.approx(135.5, abs=0.5)


def test_q0_equals_unsmoothed(base, sweep):
    assert np.allclose(sweep[0].result_y["s"], base.result_y["s"])


def test_q100_reaches_minimal_range(sweep):
    r = sweep[100].result_y
    assert r["moment_range"] == pytest.approx(r["moment_range_min"], rel=1e-6)
    assert r["moment_range_min"] < r["moment_range_max"]
    # в примере минимальный размах достигается односторонней эпюрой
    assert sweep[100].moment_range_y == pytest.approx(203.2, abs=0.3)


def test_range_is_linear_and_max_moment_monotone(sweep):
    r = sweep[100].result_y
    R_hi, R_lo = r["moment_range_max"], r["moment_range_min"]
    values = []
    for q in PERCENTS:
        expected = R_hi - q / 100 * (R_hi - R_lo)
        assert sweep[q].result_y["moment_range"] == pytest.approx(expected, rel=1e-6)
        values.append(sweep[q].max_moment_y)
    assert all(b >= a * (1 - 1e-9) for a, b in zip(values, values[1:]))


def test_smoothing_costs_least_max_moment(sweep):
    # на q = 100 нет решения с тем же размахом и меньшим max|M|: max|M| >= размах/2,
    # а здесь эпюра односторонняя (M(0) = 0), поэтому max|M| = размах
    assert sweep[100].max_moment_y == pytest.approx(sweep[100].moment_range_y, rel=1e-6)


def test_equalities_hold_for_all_q(sweep):
    for q in PERCENTS:
        for axis in ("x", "y"):
            force, m_end = balance_and_end_moment(sweep[q], axis)
            assert force == pytest.approx(0, abs=1e-6)
            assert m_end == pytest.approx(0, abs=1e-6)


def test_certificates_with_smoothing(sweep):
    for r in (sweep[50].result_y, sweep[50].result_x):
        assert r["lp"].has_range
        assert verify_optimum_via_duality(r, verbose=False)["confirmed"]
        assert verify_kkt_conditions(r, verbose=False)["confirmed"]


def test_smoothing_respects_group2_budget(shaft):
    engine, components = shaft
    res = analyze(engine, components, step=0.1, verbose=False, group2_percent=90,
                  smooth_percent=50)
    r = res.result_y
    B_p = r["static_moment_min"] + 0.9 * (r["static_moment_max"] - r["static_moment_min"])
    assert r["static_moment_sum"] <= B_p * (1 + 1e-6)
    assert r["moment_range"] < r["moment_range_max"]
    # смена p с сохранением сглаживания
    other = with_group2_percent(res, 100)
    assert other.smooth_percent == 50
    for axis in ("x", "y"):
        force, m_end = balance_and_end_moment(other, axis)
        assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)


def test_smoothing_with_constraints(shaft):
    engine, components = shaft
    cons = {"all": CounterweightConstraints.from_dict({"zero": [2, 5]})}
    prev = None
    for q in (0, 50, 100):
        res = analyze(engine, components, step=0.1, verbose=False, constraints=cons,
                      smooth_percent=q)
        assert res.result_y["s"][2] == 0 and res.result_y["s"][5] == 0
        force, m_end = balance_and_end_moment(res, "y")
        assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)
        if prev is not None:
            assert res.moment_range_y <= prev * (1 + 1e-9)
        prev = res.moment_range_y


def test_degenerate_range_does_not_fail(base, capsys):
    # по оси X силы нулевые: размах уже 0
    f, w = base.forces, base.engine.angular_speed
    r = optimize_counterweights(base.Z, f.L_cheeks, f.Fx_cheeks, f.L_crankpins,
                                f.Fx_crankpins, w, [1.0] * 8, axis="x",
                                smooth_percent=60)
    assert r["moment_range"] == pytest.approx(0, abs=1e-6)
    assert "размах уже минимален" in capsys.readouterr().out


def test_all_fixed_reports_range(base):
    f, w = base.forces, base.engine.angular_speed
    r = optimize_counterweights(base.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                                f.Fy_crankpins, w, ["F"] * 8, verbose=False,
                                smooth_percent=40)
    assert r["result"] is None
    assert r["moment_range"] == r["moment_range_min"] == r["moment_range_max"]


@pytest.mark.parametrize("q", [-5, 101])
def test_percent_out_of_range_rejected(base, q):
    f, w = base.forces, base.engine.angular_speed
    with pytest.raises(ValueError):
        optimize_counterweights(base.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                                f.Fy_crankpins, w, [1.0] * 8, verbose=False,
                                smooth_percent=q)


def test_smooth_sweep_rows(base):
    rows = smooth_sweep(base, percents=[0, 50, 100])
    assert [r["smooth_percent"] for r in rows] == [0, 50, 100]
    assert rows[0]["moment_range_y"] > rows[1]["moment_range_y"] > rows[2]["moment_range_y"]
    assert rows[0]["max_moment_y"] <= rows[1]["max_moment_y"] <= rows[2]["max_moment_y"]


def test_json_smooth_percent(tmp_path):
    data = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    assert load_smooth_percent(EXAMPLE) is None
    data["counterweights"] = {"smooth_percent": 25, "all": {"zero": [2]}}
    path = tmp_path / "s.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    assert load_smooth_percent(path) == 25
    data["counterweights"]["smooth_percent"] = -1
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError):
        load_smooth_percent(path)


def test_cli_smooth_validation(capsys):
    with pytest.raises(SystemExit):
        main([str(EXAMPLE), "--no-plots", "--smooth-percent", "120"])
    assert "0..100" in capsys.readouterr().err


def test_cli_smooth_sweep(tmp_path, capsys):
    assert main([str(EXAMPLE), "--out", str(tmp_path), "--step", "0.5",
                 "--smooth-percent", "30", "--smooth-sweep"]) == 0
    for name in ("smooth_tradeoff.png", "smooth_moments.png", "smooth_tradeoff.csv"):
        assert (tmp_path / name).exists()
    data = list(csv.DictReader((tmp_path / "smooth_tradeoff.csv").open(encoding="utf-8")))
    assert [float(r["smooth_percent"]) for r in data] == PERCENTS
    assert "Сглаживание: 30 %" in capsys.readouterr().out
