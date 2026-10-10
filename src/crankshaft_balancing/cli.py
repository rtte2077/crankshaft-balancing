"""Командная строка: ``crankbalance shaft.json [--out output] [--verify] ...``."""
import argparse
import sys
from pathlib import Path

import numpy as np

from .analysis import analyze, raw_moments
from .constraints import resolve_constraints
from .loader import load_constraints, load_counterweight_positions, load_shaft
from .optimization import optimize_counterweights_history
from .report import (counterweight_rows, format_counterweight_table,
                     save_counterweight_csv)
from .sections import NMM_TO_NM, section_grid
from .verification import (sweep_sensitivity, verify_kkt_conditions,
                           verify_optimum_via_duality)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="crankbalance",
        description="Расчёт балансировки коленвала: эпюры момента и подбор противовесов.")
    p.add_argument("config", help="JSON-файл с параметрами и элементами вала")
    p.add_argument("--out", default="output", help="папка для графиков (по умолчанию: output)")
    p.add_argument("--step", type=float, default=0.01,
                   help="шаг сетки сечений, мм (по умолчанию 0.01)")
    p.add_argument("--csv", metavar="ФАЙЛ",
                   help="сохранить таблицу статических моментов противовесов в CSV "
                        "(без этого ключа CSV кладётся в папку --out вместе с графиками)")
    p.add_argument("--table-lang", choices=["ru", "en"], default="ru",
                   help="язык подписей таблицы противовесов в консоли и на картинке "
                        "(остальные графики и сообщения — на русском)")
    p.add_argument("--no-plots", action="store_true", help="только расчёт, без графиков")
    p.add_argument("--verify", action="store_true",
                   help="проверить глобальность оптимума (двойственность + ККТ)")
    p.add_argument("--sensitivity", action="store_true",
                   help="сохранить графики чувствительности max|M| к каждому s_i")
    p.add_argument("--convergence", action="store_true",
                   help="сохранить графики сходимости симплекс-метода (грубая сетка, "
                        "нужен SciPy с linprog(method='simplex'))")
    p.add_argument("--interactive", action="store_true",
                   help="дополнительно сохранить интерактивный 3D-график (нужен plotly)")
    return p


def _convergence_plots(res, constraints, out: Path) -> None:
    """Графики хода симплекс-метода на грубой сетке (как раздел 8 исходного ноутбука).

    Для оси X без пользовательских границ ставятся границы (-5, 5): без них
    вырожденная задача (момент по X нулевой по построению) для симплекса «unbounded».
    """
    from .plotting import plot_convergence

    f, w = res.forces, res.engine.angular_speed
    Zc = section_grid(res.engine.crankshaft_length, step=2.0)
    n = len(res.L_counterweights)
    for axis, Fc, Fp in (("y", f.Fy_cheeks, f.Fy_crankpins),
                         ("x", f.Fx_cheeks, f.Fx_crankpins)):
        s_init, s_fixed, s_bounds = resolve_constraints(constraints, axis).to_solver_args(n)
        if s_bounds is None and axis == "x":
            s_bounds = [(-5, 5)] * n
        try:
            h = optimize_counterweights_history(
                Zc, f.L_cheeks, Fc, f.L_crankpins, Fp, w, s_init, s_fixed=s_fixed,
                L3=res.L_counterweights, s_bounds=s_bounds, axis=axis)
        except (RuntimeError, ValueError) as exc:
            print(f"[{axis}] график сходимости пропущен: {exc}")
            continue
        plot_convergence(h["history"], axis=axis,
                         save_path=out / f"convergence_{axis}.png", show=False)


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    engine, components = load_shaft(args.config)
    constraints = load_constraints(args.config)
    res = analyze(engine, components, step=args.step, constraints=constraints,
                  counterweight_positions=load_counterweight_positions(args.config))

    print(f"\nmax|Mx| = {res.max_moment_x:.4f} Н·м, max|My| = {res.max_moment_y:.4f} Н·м")

    rows = counterweight_rows(res)
    print("\n" + format_counterweight_table(rows, lang=args.table_lang))
    if args.csv:
        save_counterweight_csv(rows, args.csv)
        print(f"Таблица сохранена: {args.csv}")

    if args.verify:
        print()
        for axis, r in (("y", res.result_y), ("x", res.result_x)):
            if r["result"] is None:
                print(f"[{axis}] все коэффициенты зафиксированы — проверять нечего\n")
                continue
            verify_optimum_via_duality(r, axis=axis)
            verify_kkt_conditions(r, axis=axis)
            print()

    if args.no_plots:
        return 0

    import matplotlib
    matplotlib.use("Agg")  # без окон — графики только в файлы
    from .plotting import (plot_counterweight_table, plot_distribution,
                           plot_moment_3d, plot_sensitivity)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if not args.csv:
        save_counterweight_csv(rows, out / "counterweights.csv")
    plot_counterweight_table(rows, save_path=out / "counterweights.png", show=False,
                             lang=args.table_lang)
    f, Z, L = res.forces, res.Z, res.engine.crankshaft_length

    for axis, Fc, Fp, r in (("x", f.Fx_cheeks, f.Fx_crankpins, res.result_x),
                            ("y", f.Fy_cheeks, f.Fy_crankpins, res.result_y)):
        raw = raw_moments(Z, f)[0 if axis == "x" else 1]
        groups = [{"L": f.L_cheeks, "F": Fc, "name": "Cheek"},
                  {"L": f.L_crankpins, "F": Fp, "name": "crankpin"}]
        plot_distribution(Z, raw, L, groups=groups,
                          title=f"Момент M(z) до оптимизации (ось {axis.upper()})",
                          ylabel="Момент, Н·м", save_path=out / f"moment_{axis}_raw.png",
                          show=False)
        plot_distribution(Z, r["M"] * NMM_TO_NM, L,
                          groups=groups + [{"L": res.L_counterweights, "F": r["F3"],
                                            "name": "F3 (оптим.)"}],
                          title=f"Момент M(z) после оптимизации F3 (ось {axis.upper()})",
                          ylabel="Момент, Н·м", save_path=out / f"moment_{axis}_optimized.png",
                          show=False)

    plot_moment_3d(Z, res.result_x, res.result_y, shaft_length=L,
                   save_path=out / "moment_3d.png", show=False, verbose=False)
    if args.interactive:
        plot_moment_3d(Z, res.result_x, res.result_y, shaft_length=L, interactive=True,
                       save_path=str(out / "moment_3d.html"), verbose=False)

    if args.sensitivity:
        for axis, Fc, Fp, r in (("y", f.Fy_cheeks, f.Fy_crankpins, res.result_y),
                                ("x", f.Fx_cheeks, f.Fx_crankpins, res.result_x)):
            if r["result"] is None:
                continue
            curves = sweep_sensitivity(r, Z, f.L_cheeks, Fc, f.L_crankpins, Fp,
                                       res.engine.angular_speed)
            plot_sensitivity(curves, r["s"], axis=axis,
                             save_path=out / f"sensitivity_{axis}.png", show=False)

    if args.convergence:
        _convergence_plots(res, constraints, out)

    import matplotlib.pyplot as plt
    plt.close("all")   # фигуры уже сохранены в файлы — освобождаем память
    print(f"Графики сохранены в {out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
