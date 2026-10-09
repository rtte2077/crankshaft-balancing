"""Командная строка: ``crankbalance shaft.json [--out output] [--no-plots]``."""
import argparse
import sys
from pathlib import Path

from .analysis import analyze, raw_moments
from .loader import load_shaft
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
    p.add_argument("--no-plots", action="store_true", help="только расчёт, без графиков")
    p.add_argument("--verify", action="store_true",
                   help="проверить глобальность оптимума (двойственность + ККТ)")
    p.add_argument("--interactive", action="store_true",
                   help="дополнительно сохранить интерактивный 3D-график (нужен plotly)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    engine, components = load_shaft(args.config)
    res = analyze(engine, components, step=args.step)

    print(f"\nmax|Mx| = {res.max_moment_x:.4f} Н·м, max|My| = {res.max_moment_y:.4f} Н·м")

    if args.verify:
        print()
        for axis, r in (("y", res.result_y), ("x", res.result_x)):
            verify_optimum_via_duality(r, axis=axis)
            verify_kkt_conditions(r, axis=axis)
            print()

    if args.no_plots:
        return 0

    import matplotlib
    matplotlib.use("Agg")  # без окон — графики только в файлы
    from .plotting import plot_distribution, plot_moment_3d

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    f, Z, L = res.forces, res.Z, res.engine.crankshaft_length
    nm = 1e-3

    for axis, Fc, Fp, r in (("x", f.Fx_cheeks, f.Fx_crankpins, res.result_x),
                            ("y", f.Fy_cheeks, f.Fy_crankpins, res.result_y)):
        raw = raw_moments(Z, f)[0 if axis == "x" else 1]
        groups = [{"L": f.L_cheeks, "F": Fc, "name": "Cheek"},
                  {"L": f.L_crankpins, "F": Fp, "name": "crankpin"}]
        plot_distribution(Z, raw, L, groups=groups,
                          title=f"Момент M(z) до оптимизации (ось {axis.upper()})",
                          ylabel="Момент, Н·м", save_path=out / f"moment_{axis}_raw.png",
                          show=False)
        plot_distribution(Z, r["M"] * nm, L,
                          groups=groups + [{"L": f.L_cheeks, "F": r["F3"], "name": "F3 (оптим.)"}],
                          title=f"Момент M(z) после оптимизации F3 (ось {axis.upper()})",
                          ylabel="Момент, Н·м", save_path=out / f"moment_{axis}_optimized.png",
                          show=False)

    plot_moment_3d(Z, res.result_x, res.result_y, shaft_length=L,
                   save_path=out / "moment_3d.png", show=False, verbose=False)
    if args.interactive:
        plot_moment_3d(Z, res.result_x, res.result_y, shaft_length=L, interactive=True,
                       save_path=str(out / "moment_3d.html"), verbose=False)

    print(f"Графики сохранены в {out.resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
