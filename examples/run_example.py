"""Минимальный пример использования библиотеки как Python-пакета.

    python examples/run_example.py
"""
from pathlib import Path

from crankshaft_balancing import (analyze, load_shaft, verify_kkt_conditions,
                                  verify_optimum_via_duality)
from crankshaft_balancing.plotting import plot_distribution
from crankshaft_balancing.analysis import raw_moments

engine, components = load_shaft(Path(__file__).with_name("example_shaft.json"))
res = analyze(engine, components, step=0.01)

for axis, r in (("y", res.result_y), ("x", res.result_x)):
    verify_optimum_via_duality(r, axis=axis)
    verify_kkt_conditions(r, axis=axis)

f = res.forces
plot_distribution(
    res.Z, res.result_y["M"] * 1e-3, engine.crankshaft_length,
    groups=[{"L": f.L_cheeks, "F": f.Fy_cheeks, "name": "Cheek (F1)"},
            {"L": f.L_crankpins, "F": f.Fy_crankpins, "name": "crankpin (F2)"},
            {"L": f.L_cheeks, "F": res.result_y["F3"], "name": "F3 (оптим.)"}],
    title="Момент M(z) после оптимизации F3 (ось Y)", ylabel="Момент, Н·м",
)
