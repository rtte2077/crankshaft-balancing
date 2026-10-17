"""Командная строка: ``crankbalance shaft.json [--out output] [--verify] ...``."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from .analysis import analyze, group2_sweep, raw_moments
from .constraints import resolve_constraints
from .loader import (load_constraints, load_counterweight_positions,
                     load_group2_percent, load_shaft)
from .optimization import optimize_counterweights_history
from .report import (counterweight_rows, format_counterweight_table,
                     save_counterweight_csv)
from .sections import NMM_TO_NM, section_grid
from .verification import (sweep_sensitivity, verify_kkt_conditions,
                           verify_optimum_via_duality)


def _percent(value: str) -> float:
    """Тип аргумента: число 0..100."""
    try:
        p = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"ожидается число 0..100, получено {value!r}")
    if not 0.0 <= p <= 100.0:
        raise argparse.ArgumentTypeError(
            f"доля группы 2 должна быть в диапазоне 0..100, получено {p:g}")
    return p


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
    p.add_argument("--group2-percent", type=_percent, default=None, metavar="P",
                   help="доля группы 2 статических моментов, 0..100: 0 — минимальный "
                        "суммарный статический момент, 100 — минимальный max|M| "
                        "(по умолчанию: counterweights.group2_percent из JSON, иначе 100)")
    p.add_argument("--norm", choices=["l1", "linf"], default="l1",
                   help="бюджет статического момента по оси: l1 — сумма |s_i| "
                        "(по умолчанию), linf — max |s_i|")
    p.add_argument("--group2-sweep", action="store_true",
                   help="посчитать p = 0, 10, ..., 100 и сохранить в --out "
                        "group2_tradeoff.png и group2_tradeoff.csv")
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


def _group2_sweep(res, out: Path, norm: str, plots: bool = True) -> None:
    """p = 0, 10, ..., 100: таблица group2_tradeoff.csv и график group2_tradeoff.png."""
    rows = group2_sweep(res)
    out.mkdir(parents=True, exist_ok=True)
    keys = ["group2_percent", "max_moment_x", "max_moment_y",
            "static_moment_x", "static_moment_y"]
    header = ["group2_percent", "max_abs_Mx_N_m", "max_abs_My_N_m",
              f"static_moment_x_{norm}_kg_m", f"static_moment_y_{norm}_kg_m"]
    with open(out / "group2_tradeoff.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for r in rows:
            writer.writerow([f"{r[k]:.9g}" for k in keys])
    print("\nГруппа 2, %   max|Mx|, Н·м   max|My|, Н·м   B_x, кг·м   B_y, кг·м")
    for r in rows:
        print(f"{r['group2_percent']:>12g} {r['max_moment_x']:>14.4f} "
              f"{r['max_moment_y']:>14.4f} {r['static_moment_x']:>11.6f} "
              f"{r['static_moment_y']:>11.6f}")
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        from .plotting import plot_group2_tradeoff
        plot_group2_tradeoff(rows, norm=norm, save_path=out / "group2_tradeoff.png",
                             show=False)
    print(f"Компромисс группы 2 сохранён: {(out / 'group2_tradeoff.csv').resolve()}")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    engine, components = load_shaft(args.config)
    constraints = load_constraints(args.config)
    group2_percent = args.group2_percent
    if group2_percent is None:
        group2_percent = load_group2_percent(args.config)
    res = analyze(engine, components, step=args.step, constraints=constraints,
                  counterweight_positions=load_counterweight_positions(args.config),
                  group2_percent=100.0 if group2_percent is None else group2_percent,
                  norm=args.norm)

    print(f"\nmax|Mx| = {res.max_moment_x:.4f} Н·м, max|My| = {res.max_moment_y:.4f} Н·м")
    print(f"Группа 2: {res.group2_percent:g} %; суммарный статический момент "
          f"({args.norm}): x = {res.result_x['static_moment_sum']:.6f}, "
          f"y = {res.result_y['static_moment_sum']:.6f} кг·м")

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

    out = Path(args.out)
    if args.group2_sweep:
        _group2_sweep(res, out, args.norm, plots=not args.no_plots)

    if args.no_plots:
        return 0

    import matplotlib
    matplotlib.use("Agg")  # без окон — графики только в файлы
    from .plotting import (plot_counterweight_table, plot_distribution,
                           plot_moment_3d, plot_sensitivity)

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
