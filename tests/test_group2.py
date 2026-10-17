import csv
from pathlib import Path

import numpy as np
import pytest

from crankshaft_balancing import (CounterweightConstraints, analyze, counterweight_rows,
                                  format_counterweight_table, group2_sweep,
                                  load_group2_percent, minimal_moment_budget,
                                  minimal_static_moment, optimize_counterweights,
                                  save_counterweight_csv, static_moment,
                                  verify_kkt_conditions, verify_optimum_via_duality,
                                  with_group2_percent)
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
    return {p: with_group2_percent(base, p) for p in PERCENTS}


def axis_args(res, axis):
    f = res.forces
    Fc, Fp = ((f.Fy_cheeks, f.Fy_crankpins) if axis == "y"
              else (f.Fx_cheeks, f.Fx_crankpins))
    return res.Z, f.L_cheeks, Fc, f.L_crankpins, Fp, res.engine.angular_speed


def balance_and_end_moment(res, axis):
    _, _, Fc, _, Fp, w = axis_args(res, axis)
    r = getattr(res, f"result_{axis}")
    return Fc.sum() + Fp.sum() + r["s"].sum() * w ** 2, r["M"][-1] * NMM_TO_NM


def test_default_is_full_group2_and_regression(shaft):
    engine, components = shaft
    res = analyze(engine, components, step=0.01, verbose=False)
    assert res.group2_percent == 100.0
    assert res.max_moment_y == pytest.approx(135.5133, abs=1e-4)
    assert res.result_y["max_moment_opt"] * NMM_TO_NM == pytest.approx(135.5133, abs=1e-4)


def test_new_result_keys(base):
    for r in (base.result_x, base.result_y):
        for key in ("s", "F3", "M", "success", "result", "lp", "s_group1", "s_group2",
                    "static_moment_sum", "static_moment_min", "static_moment_max",
                    "group2_percent", "max_moment_opt"):
            assert key in r
        assert np.allclose(r["s_group1"] + r["s_group2"], r["s"], atol=1e-12)


def test_p0_reaches_minimal_static_moment(sweep):
    for axis in ("x", "y"):
        r0, r100 = getattr(sweep[0], f"result_{axis}"), getattr(sweep[100], f"result_{axis}")
        assert np.abs(r0["s"]).sum() == pytest.approx(r0["static_moment_min"], abs=1e-9)
        assert np.abs(r0["s"]).sum() <= np.abs(r100["s"]).sum() + 1e-9


def test_static_moment_is_linear_in_p(sweep):
    r = sweep[0].result_y
    B_min, B_max = r["static_moment_min"], r["static_moment_max"]
    assert B_max > B_min
    for p in PERCENTS:
        expected = B_min + p / 100 * (B_max - B_min)
        assert sweep[p].result_y["static_moment_sum"] == pytest.approx(expected, rel=1e-6,
                                                                       abs=1e-9)


def test_max_moment_monotone_in_p(sweep):
    values = [sweep[p].max_moment_y for p in PERCENTS]
    assert all(b <= a * (1 + 1e-9) + 1e-9 for a, b in zip(values, values[1:]))
    assert values[-1] == pytest.approx(sweep[100].result_y["max_moment_opt"] * NMM_TO_NM,
                                       rel=1e-8)


def test_equalities_hold_for_all_p(sweep):
    for p in PERCENTS:
        for axis in ("x", "y"):
            force, m_end = balance_and_end_moment(sweep[p], axis)
            assert force == pytest.approx(0, abs=1e-6)
            assert m_end == pytest.approx(0, abs=1e-6)


def test_certificates_with_budget(sweep):
    for r in (sweep[50].result_y, sweep[50].result_x):
        assert r["lp"].n_u > 0                    # проверяется именно задача с бюджетом
        assert verify_optimum_via_duality(r, verbose=False)["confirmed"]
        assert verify_kkt_conditions(r, verbose=False)["confirmed"]


def test_step_functions_agree_with_result(base):
    args = axis_args(base, "y")
    s1, B_min = minimal_static_moment(*args, [1.0] * 8, L3=base.L_counterweights)
    s3, B_max, M_opt = minimal_moment_budget(*args, [1.0] * 8, L3=base.L_counterweights)
    assert B_min == pytest.approx(base.result_y["static_moment_min"], abs=1e-9)
    assert B_max == pytest.approx(base.result_y["static_moment_max"], rel=1e-6)
    assert M_opt == pytest.approx(base.result_y["max_moment_opt"], rel=1e-9)
    assert static_moment(s3) == pytest.approx(B_max)


@pytest.mark.parametrize("constraints", [
    {"all": {"zero": [2, 5]}, "y": {"bounds": {"0": [-0.15, 0.15], "7": [-0.15, 0.15]}}},
    {"y": {"fixed": {"0": 0.05}}},
])
def test_constraints_are_compatible_with_budget(shaft, constraints):
    engine, components = shaft
    cons = {k: CounterweightConstraints.from_dict(v) for k, v in constraints.items()}
    prev = None
    for p in (0, 50, 100):
        res = analyze(engine, components, step=0.1, verbose=False, constraints=cons,
                      group2_percent=p)
        for axis in ("x", "y"):
            force, m_end = balance_and_end_moment(res, axis)
            assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)
        s = res.result_y["s"]
        for i, v in cons.get("all", CounterweightConstraints()).fixed.items():
            assert s[i] == v
        for i, v in cons.get("y", CounterweightConstraints()).fixed.items():
            assert s[i] == v
            assert res.result_y["s_group1"][i] == pytest.approx(v)
        for i, (lo, hi) in cons.get("y", CounterweightConstraints()).bounds.items():
            assert lo - 1e-9 <= s[i] <= hi + 1e-9
        r = res.result_y
        expected = (r["static_moment_min"]
                    + p / 100 * (r["static_moment_max"] - r["static_moment_min"]))
        assert r["static_moment_sum"] == pytest.approx(expected, rel=1e-6, abs=1e-9)
        if prev is not None:
            assert res.max_moment_y <= prev * (1 + 1e-9)
        prev = res.max_moment_y


def test_custom_positions_with_budget(shaft):
    engine, components = shaft
    res = analyze(engine, components, step=0.1, verbose=False, group2_percent=30,
                  counterweight_positions=[20, 120, 260, 400, 500])
    for axis in ("x", "y"):
        force, m_end = balance_and_end_moment(res, axis)
        assert force == pytest.approx(0, abs=1e-6) and m_end == pytest.approx(0, abs=1e-6)
    assert len(res.result_y["s_group1"]) == 5


def test_linf_norm(base):
    args = axis_args(base, "y")
    r0 = optimize_counterweights(*args, [1.0] * 8, verbose=False, group2_percent=0,
                                 norm="linf")
    r100 = optimize_counterweights(*args, [1.0] * 8, verbose=False, norm="linf")
    assert np.abs(r0["s"]).max() == pytest.approx(r0["static_moment_min"], abs=1e-9)
    assert r100["static_moment_sum"] == pytest.approx(r100["static_moment_max"], rel=1e-6)
    assert np.max(np.abs(r100["M"])) == pytest.approx(r100["max_moment_opt"], rel=1e-8)


def test_infeasible_equalities_raise(base):
    # вал в примере уравновешен сам (сумма сил 0), а все s > 0 — баланс сил невыполним
    args = axis_args(base, "y")
    with pytest.raises(RuntimeError, match="два равенства"):
        optimize_counterweights(*args, [1.0] * 8, s_bounds=[(0.1, 0.2)] * 8,
                                verbose=False, group2_percent=50)


def test_all_fixed_ignores_percent(base):
    args = axis_args(base, "y")
    r = optimize_counterweights(*args, ["F"] * 8, verbose=False, group2_percent=30)
    assert r["result"] is None
    assert np.array_equal(r["s_group1"], r["s"]) and not np.any(r["s_group2"])


def test_degenerate_budget_does_not_fail(base, capsys):
    # по оси X силы нулевые: B_min = B_max = 0
    args = axis_args(base, "x")
    r = optimize_counterweights(*args, [1.0] * 8, verbose=True, group2_percent=40, axis="x")
    assert r["static_moment_min"] == r["static_moment_max"] == pytest.approx(0, abs=1e-9)
    assert "B_max = B_min" in capsys.readouterr().out


@pytest.mark.parametrize("p", [-1, 100.5])
def test_percent_out_of_range_rejected(base, p):
    with pytest.raises(ValueError):
        optimize_counterweights(*axis_args(base, "y"), [1.0] * 8, verbose=False,
                                group2_percent=p)


def test_group2_sweep_rows(base):
    rows = group2_sweep(base, percents=[0, 100])
    assert [r["group2_percent"] for r in rows] == [0, 100]
    assert rows[0]["max_moment_y"] >= rows[1]["max_moment_y"]


def test_report_has_group_columns(sweep, tmp_path):
    rows = counterweight_rows(sweep[50])
    for r, s1x, s1y in zip(rows, sweep[50].result_x["s_group1"],
                           sweep[50].result_y["s_group1"]):
        assert r.has_groups
        assert r.s1_x == pytest.approx(s1x) and r.s1_y == pytest.approx(s1y)
        assert r.s1_x + r.s2_x == pytest.approx(r.s_x, abs=1e-9)
        assert r.s1_y + r.s2_y == pytest.approx(r.s_y, abs=1e-9)
    assert "гр.1" in format_counterweight_table(rows)
    assert "gr.2" in format_counterweight_table(rows, lang="en")
    path = tmp_path / "t.csv"
    save_counterweight_csv(rows, path)
    data = list(csv.DictReader(path.open(encoding="utf-8")))
    assert {"s1_x_kg_m", "s1_y_kg_m", "s2_x_kg_m", "s2_y_kg_m"} <= set(data[0])


def test_json_percent(tmp_path):
    import json
    data = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    assert load_group2_percent(EXAMPLE) is None
    data["counterweights"] = {"group2_percent": 40, "all": {"zero": [2]}}
    path = tmp_path / "s.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    assert load_group2_percent(path) == 40
    data["counterweights"]["group2_percent"] = 140
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError):
        load_group2_percent(path)


def test_cli_percent_validation(capsys):
    with pytest.raises(SystemExit):
        main([str(EXAMPLE), "--no-plots", "--group2-percent", "150"])
    assert "0..100" in capsys.readouterr().err


def test_cli_group2_sweep(tmp_path, capsys):
    assert main([str(EXAMPLE), "--out", str(tmp_path), "--step", "0.5",
                 "--group2-percent", "40", "--group2-sweep"]) == 0
    assert (tmp_path / "group2_tradeoff.png").exists()
    data = list(csv.DictReader((tmp_path / "group2_tradeoff.csv").open(encoding="utf-8")))
    assert [float(r["group2_percent"]) for r in data] == PERCENTS
    assert "Группа 2: 40 %" in capsys.readouterr().out
