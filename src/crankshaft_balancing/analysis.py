"""Сквозной расчёт: геометрия -> силы -> оптимизация противовесов по осям X и Y."""
from dataclasses import dataclass
from typing import Mapping, Optional, Sequence

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

    @property
    def max_moment_x(self) -> float:
        """max|Mx| после оптимизации, Н·м."""
        return float(np.max(np.abs(self.result_x["M"])) * NMM_TO_NM)

    @property
    def max_moment_y(self) -> float:
        """max|My| после оптимизации, Н·м."""
        return float(np.max(np.abs(self.result_y["M"])) * NMM_TO_NM)


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


def analyze(engine: EngineConfig, components: Sequence[ShaftComponent],
            step: float = 0.01, verbose: bool = True,
            constraints: Optional[Mapping[str, CounterweightConstraints]] = None,
            counterweight_positions: Optional[Sequence[float]] = None
            ) -> BalanceResult:
    """Полный расчёт балансировки для заданных параметров и элементов вала.

    ``engine.crankshaft_length`` заполняется автоматически из геометрии.
    ``step`` — шаг сетки сечений, мм.
    ``constraints`` — ограничения на противовесы: словарь с ключами ``"x"``, ``"y"``
    и/или ``"all"`` (обе оси), значения — ``CounterweightConstraints``.
    Номер противовеса = номер щеки вдоль вала по возрастанию x, с 0.
    ``counterweight_positions`` — координаты x (мм, как у элементов вала), где стоят
    противовесы; по умолчанию противовесы стоят на щёках. Если заданы, номер
    противовеса — номер позиции в этом списке (после сортировки по возрастанию).
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
    n = len(L3)

    results = {}
    for axis, Fc, Fp in (("y", f.Fy_cheeks, f.Fy_crankpins),
                         ("x", f.Fx_cheeks, f.Fx_crankpins)):
        s_init, s_fixed, s_bounds = resolve_constraints(constraints, axis).to_solver_args(n)
        results[axis] = optimize_counterweights(
            Z, f.L_cheeks, Fc, f.L_crankpins, Fp, engine.angular_speed,
            s_init, s_fixed=s_fixed, L3=L3, s_bounds=s_bounds, axis=axis,
            verbose=verbose)

    return BalanceResult(engine=engine, geometry=geometry, Z=Z, forces=f,
                         module_width=module_width,
                         result_x=results["x"], result_y=results["y"],
                         L_counterweights=np.asarray(L3, dtype=float))
