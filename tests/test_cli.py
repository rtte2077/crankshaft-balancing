from pathlib import Path

from crankshaft_balancing.cli import main

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "example_shaft.json"


def test_cli_creates_plots(tmp_path):
    assert main([str(EXAMPLE), "--out", str(tmp_path), "--step", "0.1", "--verify"]) == 0
    names = {p.name for p in tmp_path.iterdir()}
    assert {"moment_y_optimized.png", "moment_3d.png"} <= names
