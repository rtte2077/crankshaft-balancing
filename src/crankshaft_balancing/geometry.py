"""Геометрия вала: щёки, шатунные шейки, вывод длины вала."""
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

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
    # Индивидуальные ширины элементов, мм (NaN — взять ширину по умолчанию из EngineConfig).
    cheek_widths: Optional[np.ndarray] = None
    crankpin_widths: Optional[np.ndarray] = None

    def __post_init__(self):
        self.cheek_centers = np.asarray(self.cheek_centers, dtype=float).reshape(-1, 3)
        self.crankpin_centers = np.asarray(self.crankpin_centers, dtype=float).reshape(-1, 3)
        self.cheek_masses = np.asarray(self.cheek_masses, dtype=float)
        self.crankpin_masses = np.asarray(self.crankpin_masses, dtype=float)

        if len(self.cheek_centers) != len(self.cheek_masses):
            raise ValueError("Число центров масс щёк должно совпадать с числом их масс")
        if len(self.crankpin_centers) != len(self.crankpin_masses):
            raise ValueError("Число центров масс шеек должно совпадать с числом их масс")

        self.cheek_widths = _widths_array(self.cheek_widths, len(self.cheek_masses), "щёк")
        self.crankpin_widths = _widths_array(self.crankpin_widths,
                                             len(self.crankpin_masses), "шеек")


def _widths_array(widths, n: int, what: str) -> np.ndarray:
    """Массив индивидуальных ширин длины ``n``; None/NaN — ширина по умолчанию."""
    if widths is None:
        return np.full(n, np.nan)
    arr = np.array([np.nan if w is None else w for w in widths], dtype=float)
    if len(arr) != n:
        raise ValueError(f"Число ширин {what} должно совпадать с числом элементов")
    if np.any(arr[~np.isnan(arr)] <= 0):
        raise ValueError(f"Ширины {what} должны быть положительными")
    return arr


@dataclass
class ShaftComponent:
    """Один элемент коленвала — щека или шатунная шейка.

    kind: "cheek" (щека) или "crankpin" (шатунная шейка).
    x, y, z: координаты центра масс, мм (x — вдоль оси вала).
    mass: масса элемента, кг.
    width: толщина (ширина вдоль оси вала) именно этого элемента, мм. Необязательна:
        если не задана, берётся общая ``cheek_width``/``crankpin_journal_width``
        из ``EngineConfig``.
    """

    kind: str
    x: float
    y: float
    z: float
    mass: float
    width: Optional[float] = None


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
        widths = [c.width for c in items]
        return centers, masses, widths

    cheek_centers, cheek_masses, cheek_widths = _to_arrays("cheek")
    crankpin_centers, crankpin_masses, crankpin_widths = _to_arrays("crankpin")

    return CrankshaftGeometry(
        cheek_centers=cheek_centers, cheek_masses=cheek_masses,
        crankpin_centers=crankpin_centers, crankpin_masses=crankpin_masses,
        cheek_widths=cheek_widths, crankpin_widths=crankpin_widths,
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

    Общие ширины из ``engine`` — значения по умолчанию. Их можно переопределить
    для конкретного элемента: ``ShaftComponent.width`` для щеки или шатунной
    шейки и ``engine.main_journal_widths`` (номер коренной шейки -> ширина) для
    коренных шеек. Коренные шейки нумеруются слева направо с 0: 0 — левая
    концевая, последняя — правая концевая (от концевых берётся половина ширины).

    Возвращает (длина вала, шаг колена) в мм.
    """
    tagged = ([(x, "cheek", w) for x, w in zip(geometry.cheek_centers[:, 0],
                                                  geometry.cheek_widths)]
              + [(x, "crankpin", w) for x, w in zip(geometry.crankpin_centers[:, 0],
                                                      geometry.crankpin_widths)])
    if not tagged:
        raise ValueError("В геометрии нет ни одной щеки или шейки")
    tagged.sort(key=lambda item: item[0])
    kinds = [kind for _, kind, _ in tagged]

    default_width = {"cheek": engine.cheek_width, "crankpin": engine.crankpin_journal_width}

    n_main = 2 + sum(1 for a, b in zip(kinds, kinds[1:]) if a == b == "cheek")
    main_overrides = engine.main_journal_widths or {}
    bad = [i for i in main_overrides if not 0 <= i < n_main]
    if bad:
        raise ValueError(f"Нет коренных шеек с номерами {sorted(bad)}: "
                         f"на валу {n_main} коренных шеек (номера 0..{n_main - 1})")
    if any(w <= 0 for w in main_overrides.values()):
        raise ValueError("Ширины коренных шеек должны быть положительными")

    def main_width(i: int) -> float:
        return float(main_overrides.get(i, engine.main_journal_width))

    main_index = 0
    length = main_width(main_index) / 2
    for i, (_, kind, width) in enumerate(tagged):
        length += default_width[kind] if np.isnan(width) else width
        if i + 1 < len(kinds) and kind == "cheek" and kinds[i + 1] == "cheek":
            main_index += 1
            length += main_width(main_index)
    length += main_width(main_index + 1) / 2

    if len(geometry.crankpin_centers) >= 2:
        pitch_source = np.sort(geometry.crankpin_centers[:, 0])
    else:
        pitch_source = np.sort(geometry.cheek_centers[:, 0])
    module_width = float(np.mean(np.diff(pitch_source))) if len(pitch_source) >= 2 else length

    return length, module_width
