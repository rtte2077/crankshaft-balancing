"""Сквозной расчёт: геометрия -> силы -> оптимизация противовесов по осям X и Y."""
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .config import EngineConfig
from .forces import (centrifugal_force_magnitude, centrifugal_force_vectors,
                     distance_from_right_end, radial_distance)
from .geometry import (CrankshaftGeometry, ShaftComponent, build_geometry,
                       derive_crankshaft_dimensions)
from .optimization import make_initial_coefficients, optimize_counterweights
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
            step: float = 0.01, verbose: bool = True) -> BalanceResult:
    """Полный расчёт балансировки для заданных параметров и элементов вала.

    ``engine.crankshaft_length`` заполняется автоматически из геометрии.
    ``step`` — шаг сетки сечений, мм.
    """
    geometry = build_geometry(components)
    engine.crankshaft_length, module_width = derive_crankshaft_dimensions(geometry, engine)
    if verbose:
        print(f"Длина коленвала (из геометрии): {engine.crankshaft_length:.3f} мм")
        print(f"Шаг колена (из геометрии):      {module_width:.3f} мм")

    Z = section_grid(engine.crankshaft_length, step=step)
    f = compute_forces(geometry, engine)

    result_y = optimize_counterweights(
        Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins, f.Fy_crankpins,
        engine.angular_speed, make_initial_coefficients(len(f.Fy_cheeks)),
        axis="y", verbose=verbose)
    result_x = optimize_counterweights(
        Z, f.L_cheeks, f.Fx_cheeks, f.L_crankpins, f.Fx_crankpins,
        engine.angular_speed, make_initial_coefficients(len(f.Fx_cheeks)),
        axis="x", verbose=verbose)

    return BalanceResult(engine=engine, geometry=geometry, Z=Z, forces=f,
                         module_width=module_width,
                         result_x=result_x, result_y=result_y)
