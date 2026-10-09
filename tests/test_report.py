import csv
import math
from pathlib import Path

import numpy as np
import pytest

from crankshaft_balancing import (analyze, counterweight_rows,
                                  format_counterweight_table,
                                  rows_from_components, save_counterweight_csv)
from crankshaft_balancing.cli import main

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "example_shaft.json"


def test_angle_convention_from_z_axis():
    # (s_x, s_y): (0,+) -> 0°, (+,0) -> 90°, (0,-) -> 180°, (-,0) -> 270°
    rows = rows_from_components([0, 1, 2, 3], [0, 1, 0, -1], [1, 0, -1, 0])
    assert [round(r.angle_deg, 6) for r in rows] == [0.0, 90.0, 180.0, 270.0]
    assert all(r.s_abs == pytest.approx(1.0) for r in rows)


def test_diagonal_component_and_modulus():
    (r,) = rows_from_components([5.0], [3.0], [4.0])
    assert r.s_abs == pytest.approx(5.0)
    assert r.angle_deg == pytest.approx(math.degrees(math.atan2(3, 4)))


def test_zero_counterweight_has_no_angle():
    (r,) = rows_from_components([5.0], [0.0], [0.0])
    assert r.s_abs == 0 and math.isnan(r.angle_deg)


def test_solver_noise_does_not_flip_angle():
    # шум 1e-17 в нулевой составляющей не должен давать ни 359.99999.., ни 360
    rows = rows_from_components([0, 1], [-1e-17, 1e-17], [0.12, 0.12])
    assert [r.angle_deg for r in rows] == [0.0, 0.0]


def test_length_mismatch_rejected():
    with pytest.raises(ValueError):
        rows_from_components([1, 2], [1], [1, 2])


def test_rows_for_real_result(shaft):
    engine, components = shaft
    res = analyze(engine, components, step=0.1, verbose=False)
    rows = counterweight_rows(res)
    assert len(rows) == 8
    assert [r.x for r in rows] == pytest.approx(
        list(res.geometry.cheek_centers[:, 0]))
    for r, sx, sy in zip(rows, res.result_x["s"], res.result_y["s"]):
        assert r.s_abs == pytest.approx(math.hypot(sx, sy), abs=1e-8)
    # в примере все y = 0: вектор лежит на оси z, т.е. угол 0° или 180°
    assert {round(r.angle_deg) for r in rows} <= {0, 180}
    text = format_counterweight_table(rows)
    assert "угол от z" in text and text.count("\n") >= 8


def test_zeroed_counterweight_shows_dash(shaft):
    from crankshaft_balancing import CounterweightConstraints
    engine, components = shaft
    res = analyze(engine, components, step=0.1, verbose=False,
                  constraints={"all": CounterweightConstraints.from_dict({"zero": [2]})})
    rows = counterweight_rows(res)
    assert math.isnan(rows[2].angle_deg)
    assert "—" in format_counterweight_table(rows)


def test_csv_roundtrip(tmp_path):
    rows = rows_from_components([10.0, 20.0], [0.0, 0.0], [0.5, 0.0])
    path = tmp_path / "t.csv"
    save_counterweight_csv(rows, path)
    data = list(csv.DictReader(path.open(encoding="utf-8")))
    assert data[0]["angle_from_z_deg"] == "0.000000"
    assert data[1]["angle_from_z_deg"] == ""        # противовеса нет
    assert float(data[0]["s_abs_kg_m"]) == pytest.approx(0.5)


def test_cli_prints_table_and_saves_csv(tmp_path, capsys):
    assert main([str(EXAMPLE), "--out", str(tmp_path), "--step", "0.1"]) == 0
    out = capsys.readouterr().out
    assert "Статические моменты противовесов" in out
    assert (tmp_path / "counterweights.csv").exists()


def test_cli_explicit_csv_without_plots(tmp_path):
    target = tmp_path / "my.csv"
    assert main([str(EXAMPLE), "--no-plots", "--step", "0.1", "--csv", str(target)]) == 0
    assert target.exists()
