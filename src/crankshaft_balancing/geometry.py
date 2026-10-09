"""Геометрия вала: щёки, шатунные шейки, вывод длины вала."""
from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

from .config import EngineConfig


@dataclass
class CrankshaftGeometry:
    """Координаты центров масс и массы щёк и шатунных шеек.

    Координаты — 3D-векторы [x, y, z] в мм:
        x — положение вдоль оси вала,
        y, z — положение в плоскости поперечного сечения.
    """

    cheek_centers: np.ndarray
    cheek_masses: np.ndarray
    crankpin_centers: np.ndarray
    crankpin_masses: np.ndarray

    def __post_init__(self):
        self.cheek_centers = np.asarray(self.cheek_centers, dtype=float).reshape(-1, 3)
        self.crankpin_centers = np.asarray(self.crankpin_centers, dtype=float).reshape(-1, 3)
        self.cheek_masses = np.asarray(self.cheek_masses, dtype=float)
        self.crankpin_masses = np.asarray(self.crankpin_masses, dtype=float)

        if len(self.cheek_centers) != len(self.cheek_masses):
            raise ValueError("Число центров масс щёк должно совпадать с числом их масс")
        if len(self.crankpin_centers) != len(self.crankpin_masses):
            raise ValueError("Число центров масс шеек должно совпадать с числом их масс")


@dataclass
class ShaftComponent:
    """Один элемент коленвала — щека или шатунная шейка.

    kind: "cheek" (щека) или "crankpin" (шатунная шейка).
    x, y, z: координаты центра масс, мм (x — вдоль оси вала).
    mass: масса элемента, кг.
    """

    kind: str
    x: float
    y: float
    z: float
    mass: float


def build_geometry(components: Sequence[ShaftComponent]) -> CrankshaftGeometry:
    """Собирает CrankshaftGeometry из списка элементов вала в произвольном порядке.

    Элементы сортируются по ``x`` и раскладываются на группы «щёки»/«шейки».
    """
    allowed_kinds = {"cheek", "crankpin"}
    unknown = {c.kind for c in components} - allowed_kinds
    if unknown:
        raise ValueError(f"Неизвестный тип элемента: {unknown}. "
                         f"Допустимые значения kind: {allowed_kinds}")

    components_sorted = sorted(components, key=lambda c: c.x)

    def _to_arrays(kind: str):
        items = [c for c in components_sorted if c.kind == kind]
        centers = np.array([[c.x, c.y, c.z] for c in items], dtype=float).reshape(-1, 3)
        masses = np.array([c.mass for c in items], dtype=float)
        return centers, masses

    cheek_centers, cheek_masses = _to_arrays("cheek")
    crankpin_centers, crankpin_masses = _to_arrays("crankpin")

    return CrankshaftGeometry(
        cheek_centers=cheek_centers, cheek_masses=cheek_masses,
        crankpin_centers=crankpin_centers, crankpin_masses=crankpin_masses,
    )


def derive_crankshaft_dimensions(geometry: CrankshaftGeometry,
                                 engine: EngineConfig) -> Tuple[float, float]:
    """Восстанавливает длину вала и шаг колена из геометрии.

    Длина получается последовательным суммированием ширин элементов в порядке
    их следования вдоль вала (x-координаты используются только для порядка):
      - вал начинается с половины коренной шейки;
      - щека добавляет ``cheek_width``, шатунная шейка — ``crankpin_journal_width``;
      - если подряд идут две щеки без шатунной шейки — между ними опорная
        коренная шейка (``main_journal_width``);
      - вал заканчивается половиной коренной шейки.

    Возвращает (длина вала, шаг колена) в мм.
    """
    tagged = ([(x, "cheek") for x in geometry.cheek_centers[:, 0]]
              + [(x, "crankpin") for x in geometry.crankpin_centers[:, 0]])
    if not tagged:
        raise ValueError("В геометрии нет ни одной щеки или шейки")
    tagged.sort(key=lambda item: item[0])
    kinds = [kind for _, kind in tagged]

    width_of = {"cheek": engine.cheek_width, "crankpin": engine.crankpin_journal_width}

    length = engine.main_journal_width / 2
    for i, kind in enumerate(kinds):
        length += width_of[kind]
        if i + 1 < len(kinds) and kind == "cheek" and kinds[i + 1] == "cheek":
            length += engine.main_journal_width
    length += engine.main_journal_width / 2

    if len(geometry.crankpin_centers) >= 2:
        pitch_source = np.sort(geometry.crankpin_centers[:, 0])
    else:
        pitch_source = np.sort(geometry.cheek_centers[:, 0])
    module_width = float(np.mean(np.diff(pitch_source))) if len(pitch_source) >= 2 else length

    return length, module_width
