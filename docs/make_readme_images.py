"""Пересоздаёт картинки для README (русские и английские версии).

    python docs/make_readme_images.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from crankshaft_balancing import (CounterweightConstraints, analyze,  # noqa: E402
                                  counterweight_rows, load_constraints, load_shaft)
from crankshaft_balancing.plotting import (plot_counterweight_table,  # noqa: E402
                                           plot_distribution)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "img"
OUT.mkdir(parents=True, exist_ok=True)

# 1. Эпюра момента (основной пример, без ограничений)
engine, components = load_shaft(ROOT / "examples" / "example_shaft.json")
res = analyze(engine, components, step=0.01, verbose=False)
f, L = res.forces, res.engine.crankshaft_length
M = res.result_y["M"] * 1e-3

plot_distribution(
    res.Z, M, L,
    groups=[{"L": f.L_cheeks, "F": f.Fy_cheeks, "name": "Cheek"},
            {"L": f.L_crankpins, "F": f.Fy_crankpins, "name": "crankpin"},
            {"L": res.L_counterweights, "F": res.result_y["F3"], "name": "F3 (оптим.)"}],
    title="Момент M(z) после оптимизации F3 (ось Y)", ylabel="Момент, Н·м",
    save_path=OUT / "moment_y_optimized.png", show=False)

plot_distribution(
    res.Z, M, L,
    groups=[{"L": f.L_cheeks, "F": f.Fy_cheeks, "name": "Crank web (F1)"},
            {"L": f.L_crankpins, "F": f.Fy_crankpins, "name": "Crankpin (F2)"},
            {"L": res.L_counterweights, "F": res.result_y["F3"],
             "name": "Counterweight (F3)"}],
    title="Bending moment M(z) after counterweight optimization (y axis)",
    ylabel="Bending moment, N·m",
    xlabel="Position along the shaft, mm (0 = left end, L = right end)",
    curve_label="Bending moment",
    save_path=OUT / "moment_y_optimized_en.png", show=False)

# 2. Таблица противовесов (пример с ограничениями: видны убранные противовесы)
cfg = ROOT / "examples" / "example_shaft_constrained.json"
engine, components = load_shaft(cfg)
res_c = analyze(engine, components, step=0.01, verbose=False,
                constraints=load_constraints(cfg))
rows = counterweight_rows(res_c)
plot_counterweight_table(rows, save_path=OUT / "counterweights_table.png",
                         show=False, lang="ru")
plot_counterweight_table(rows, save_path=OUT / "counterweights_table_en.png",
                         show=False, lang="en")
print("Готово:", *sorted(p.name for p in OUT.iterdir()), sep="\n  ")
