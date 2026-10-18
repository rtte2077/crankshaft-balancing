"""Пересоздаёт картинки для README (русские и английские версии).

    python docs/make_readme_images.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from crankshaft_balancing import (CounterweightConstraints, analyze,  # noqa: E402
                                  counterweight_rows, load_constraints, load_shaft,
                                  with_smooth_percent)
from crankshaft_balancing.plotting import (plot_counterweight_table,  # noqa: E402
                                           plot_distribution, plot_smooth_moments)

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

# 2. Таблица противовесов, основной пример без ограничений: 8 щёк — 8 противовесов
rows = counterweight_rows(res)
plot_counterweight_table(rows, save_path=OUT / "counterweights_table.png",
                         show=False, lang="ru")
plot_counterweight_table(rows, save_path=OUT / "counterweights_table_en.png",
                         show=False, lang="en")

# 3. Та же таблица для примера с ограничениями: на щеках 2 и 5 противовесов нет
cfg = ROOT / "examples" / "example_shaft_constrained.json"
engine_c, components_c = load_shaft(cfg)
res_c = analyze(engine_c, components_c, step=0.01, verbose=False,
                constraints=load_constraints(cfg))
rows_c = counterweight_rows(res_c)
plot_counterweight_table(rows_c, save_path=OUT / "counterweights_table_constrained.png",
                         show=False, lang="ru")
plot_counterweight_table(rows_c, save_path=OUT / "counterweights_table_constrained_en.png",
                         show=False, lang="en")
# 4. Сглаживание эпюры: q = 0, 50, 100 (основной пример)
smooth, levels = {}, res
for q in (0, 50, 100):
    smooth[q] = with_smooth_percent(levels, q)
    if q:
        levels = smooth[q]
for q, r in smooth.items():
    Mq = r.result_y["M"] * 1e-3
    plot_distribution(
        r.Z, Mq, L,
        groups=[{"L": f.L_cheeks, "F": f.Fy_cheeks, "name": "Cheek"},
                {"L": f.L_crankpins, "F": f.Fy_crankpins, "name": "crankpin"},
                {"L": r.L_counterweights, "F": r.result_y["F3"], "name": "F3 (оптим.)"}],
        title=(f"Сглаживание q = {q} %: max|M| = {r.max_moment_y:.1f} Н·м, "
               f"размах = {r.moment_range_y:.1f} Н·м (ось Y)"),
        ylabel="Момент, Н·м", save_path=OUT / f"smooth_q{q}.png", show=False)
    plot_distribution(
        r.Z, Mq, L,
        groups=[{"L": f.L_cheeks, "F": f.Fy_cheeks, "name": "Crank web (F1)"},
                {"L": f.L_crankpins, "F": f.Fy_crankpins, "name": "Crankpin (F2)"},
                {"L": r.L_counterweights, "F": r.result_y["F3"],
                 "name": "Counterweight (F3)"}],
        title=(f"Smoothing q = {q} %: max|M| = {r.max_moment_y:.1f} N·m, "
               f"range = {r.moment_range_y:.1f} N·m (y axis)"),
        ylabel="Bending moment, N·m",
        xlabel="Position along the shaft, mm (0 = left end, L = right end)",
        curve_label="Bending moment", save_path=OUT / f"smooth_q{q}_en.png", show=False)
curves = [(q, r.Z, r.result_x["M"], r.result_y["M"]) for q, r in smooth.items()]
plot_smooth_moments(curves, L, axes=("y",), save_path=OUT / "smooth_overlay.png",
                    show=False)
plot_smooth_moments(curves, L, axes=("y",), save_path=OUT / "smooth_overlay_en.png",
                    show=False, title="Bending moment diagrams at different smoothing degrees",
                    xlabel="Position along the shaft, mm (0 = left end)",
                    ylabel="M ({axis} axis), N·m")
for q, r in smooth.items():
    print(f"q = {q}: max|My| = {r.max_moment_y:.2f}, размах = {r.moment_range_y:.2f}, "
          f"sum|s_y| = {r.result_y['static_moment_sum']:.4f}")
    for row in counterweight_rows(r):
        print(f"   {row.number} {row.x:.3f} s_y = {row.s_y:+.4f} |s| = {row.s_abs:.4f} "
              f"угол = {row.angle_deg}")

print("Готово:", *sorted(p.name for p in OUT.iterdir()), sep="\n  ")
