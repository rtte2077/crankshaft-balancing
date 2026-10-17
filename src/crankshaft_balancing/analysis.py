"""Сквозной расчёт: геометрия -> силы -> оптимизация противовесов по осям X и Y."""
from dataclasses import dataclass, replace
from typing import List, Mapping, Optional, Sequence

import numpy as np

from .config import EngineConfig
from .constraints import CounterweightConstraints, resolve_constraints
from .forces import (centrifugal_force_magnitude, centrifugal_force_vectors,
                     distance_from_right_end, radial_distance)
from .geometry import (CrankshaftGeometry, ShaftComponent, build_geometry,
                       derive_crankshaft_dimensions)
from .optimization import optimize_counterweights
from .sections import NMM_TO_NM, section_grid, total_distribution


@dataclass
class ShaftForces:
    """Силы и точки их приложения (L — от правого конца вала)."""

    L_cheeks: np.ndarray
    L_crankpins: np.ndarray
    Fx_cheeks: np.ndarray
    Fy_cheeks: np.ndarray
    Fx_crankpins: np.ndarray
    Fy_crankpins: np.ndarray
    F_cheeks: np.ndarray       # модули
    F_crankpins: np.ndarray


@dataclass
class BalanceResult:
    """Итог расчёта балансировки."""

    engine: EngineConfig
    geometry: CrankshaftGeometry
    Z: np.ndarray
    forces: ShaftForces
    module_width: float
    result_x: dict
    result_y: dict
    L_counterweights: np.ndarray   # точки приложения F3 (L от правого конца)
    constraints: Optional[Mapping[str, CounterweightConstraints]] = None
    norm: str = "l1"
    budget_eps: float = 1e-9

    @property
    def max_moment_x(self) -> float:
        """max|Mx| после оптимизации, Н·м."""
        return float(np.max(np.abs(self.result_x["M"])) * NMM_TO_NM)

    @property
    def max_moment_y(self) -> float:
        """max|My| после оптимизации, Н·м."""
        return float(np.max(np.abs(self.result_y["M"])) * NMM_TO_NM)

    @property
    def group2_percent(self) -> float:
        return self.result_y["group2_percent"]


def compute_forces(geometry: CrankshaftGeometry, engine: EngineConfig) -> ShaftForces:
    """Центробежные силы щёк и шатунных шеек и их проекции на оси."""
    r_cheeks = radial_distance(geometry.cheek_centers)
    r_crankpins = radial_distance(geometry.crankpin_centers)

    # К массе шатунной шейки добавляется m2 — вращающаяся часть массы шатуна.
    # Это НЕ противовес: противовесы (F3) находятся оптимизацией.
    F_cheeks = centrifugal_force_magnitude(r_cheeks, geometry.cheek_masses,
                                           engine.angular_speed)
    F_crankpins = centrifugal_force_magnitude(
        r_crankpins, geometry.crankpin_masses + engine.rod_rotating_mass,
        engine.angular_speed)

    v_cheeks = centrifugal_force_vectors(geometry.cheek_centers, F_cheeks)
    v_crankpins = centrifugal_force_vectors(geometry.crankpin_centers, F_crankpins)

    return ShaftForces(
        L_cheeks=distance_from_right_end(geometry.cheek_centers, engine.crankshaft_length),
        L_crankpins=distance_from_right_end(geometry.crankpin_centers,
                                            engine.crankshaft_length),
        Fx_cheeks=v_cheeks[:, 0], Fy_cheeks=v_cheeks[:, 1],
        Fx_crankpins=v_crankpins[:, 0], Fy_crankpins=v_crankpins[:, 1],
        F_cheeks=F_cheeks, F_crankpins=F_crankpins,
    )


def raw_moments(Z: np.ndarray, f: ShaftForces) -> tuple:
    """Эпюры Mx, My до подбора противовесов, Н·м."""
    Mx = total_distribution(Z, f.L_cheeks, f.Fx_cheeks, f.L_crankpins, f.Fx_crankpins)
    My = total_distribution(Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins, f.Fy_crankpins)
    return Mx * NMM_TO_NM, My * NMM_TO_NM


def _optimize_axes(Z, f: ShaftForces, L3, omega: float,
                   constraints: Optional[Mapping[str, CounterweightConstraints]],
                   group2_percent: float, norm: str, budget_eps: float,
                   verbose: bool, levels: Optional[Mapping[str, dict]] = None) -> dict:
    """Оптимизация противовесов по осям y и x; ``levels`` — прежние результаты по
    осям (B_min, B_max, M_opt не пересчитываются)."""
    n = len(L3)
    results = {}
    for axis, Fc, Fp in (("y", f.Fy_cheeks, f.Fy_crankpins),
                         ("x", f.Fx_cheeks, f.Fx_crankpins)):
        s_init, s_fixed, s_bounds = resolve_constraints(constraints, axis).to_solver_args(n)
        results[axis] = optimize_counterweights(
            Z, f.L_cheeks, Fc, f.L_crankpins, Fp, omega,
            s_init, s_fixed=s_fixed, L3=L3, s_bounds=s_bounds, axis=axis,
            verbose=verbose, group2_percent=group2_percent, norm=norm,
            budget_eps=budget_eps, levels=None if levels is None else levels[axis])
    return results


def analyze(engine: EngineConfig, components: Sequence[ShaftComponent],
            step: float = 0.01, verbose: bool = True,
            constraints: Optional[Mapping[str, CounterweightConstraints]] = None,
            counterweight_positions: Optional[Sequence[float]] = None,
            group2_percent: float = 100.0, norm: str = "l1",
            budget_eps: float = 1e-9) -> BalanceResult:
    """Полный расчёт балансировки для заданных параметров и элементов вала.

    ``engine.crankshaft_length`` заполняется автоматически из геометрии.
    ``step`` — шаг сетки сечений, мм.
    ``constraints`` — ограничения на противовесы: словарь с ключами ``"x"``, ``"y"``
    и/или ``"all"`` (обе оси), значения — ``CounterweightConstraints``.
    Номер противовеса = номер щеки вдоль вала по возрастанию x, с 0.
    ``counterweight_positions`` — координаты x (мм, как у элементов вала), где стоят
    противовесы; по умолчанию противовесы стоят на щёках. Если заданы, номер
    противовеса — номер позиции в этом списке (после сортировки по возрастанию).
    ``group2_percent`` — доля группы 2 статических моментов, 0..100 (100 — минимальный
    max|M|, 0 — минимальный суммарный статический момент); ``norm`` — ``"l1"``
    (бюджет ``sum|s|``) или ``"linf"`` (``max|s|``); бюджет считается по каждой оси
    отдельно. См. ``optimize_counterweights``.
    """
    geometry = build_geometry(components)
    engine.crankshaft_length, module_width = derive_crankshaft_dimensions(geometry, engine)
    if verbose:
        print(f"Длина коленвала (из геометрии): {engine.crankshaft_length:.3f} мм")
        print(f"Шаг колена (из геометрии):      {module_width:.3f} мм")
        cheeks = ", ".join(f"{i}: x={x:.3f}"
                           for i, x in enumerate(geometry.cheek_centers[:, 0]))
        print(f"Щёки (номер противовеса: координата, мм): {cheeks}")

    Z = section_grid(engine.crankshaft_length, step=step)
    f = compute_forces(geometry, engine)

    if counterweight_positions is None:
        L3 = f.L_cheeks
    else:
        xs = np.sort(np.asarray(counterweight_positions, dtype=float))
        if xs.size == 0:
            raise ValueError("counterweight_positions не должен быть пустым")
        if xs.min() < 0 or xs.max() > engine.crankshaft_length:
            raise ValueError("Позиции противовесов должны лежать в пределах вала "
                             f"(0..{engine.crankshaft_length:.3f} мм)")
        L3 = engine.crankshaft_length - xs
        if verbose:
            print("Позиции противовесов (номер: x, мм): "
                  + ", ".join(f"{i}: x={x:.3f}" for i, x in enumerate(xs)))
    results = _optimize_axes(Z, f, L3, engine.angular_speed, constraints,
                             group2_percent, norm, budget_eps, verbose)

    return BalanceResult(engine=engine, geometry=geometry, Z=Z, forces=f,
                         module_width=module_width,
                         result_x=results["x"], result_y=results["y"],
                         L_counterweights=np.asarray(L3, dtype=float),
                         constraints=constraints, norm=norm, budget_eps=budget_eps)


def with_group2_percent(result: BalanceResult, group2_percent: float,
                        verbose: bool = False) -> BalanceResult:
    """Тот же расчёт с другим ``group2_percent``. B_min, B_max и M_opt берутся из
    ``result`` (они от p не зависят), решается только задача шага 4."""
    results = _optimize_axes(result.Z, result.forces, result.L_counterweights,
                             result.engine.angular_speed, result.constraints,
                             group2_percent, result.norm, result.budget_eps, verbose,
                             levels={"x": result.result_x, "y": result.result_y})
    return replace(result, result_x=results["x"], result_y=results["y"])


def group2_sweep(result: BalanceResult,
                 percents: Sequence[float] = tuple(range(0, 101, 10))) -> List[dict]:
    """max|M| и суммарный статический момент по осям для каждого p из ``percents``.

    Строки: ``group2_percent``, ``max_moment_x``/``max_moment_y`` (Н·м),
    ``static_moment_x``/``static_moment_y`` (кг·м).
    """
    rows = []
    for p in percents:
        r = with_group2_percent(result, p)
        rows.append({"group2_percent": float(p),
                     "max_moment_x": r.max_moment_x, "max_moment_y": r.max_moment_y,
                     "static_moment_x": r.result_x["static_moment_sum"],
                     "static_moment_y": r.result_y["static_moment_sum"]})
    return rows
