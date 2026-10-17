"""Таблица статических моментов противовесов: компоненты, модуль, угол.

Статический момент противовеса ``s = m * r`` (масса на радиус центра масс, кг·м) —
вектор в плоскости поперечного сечения; он направлен так же, как радиус-вектор центра
масс противовеса (сила F3 = s * omega^2 сонаправлена с ним).

Оси в расчёте и в JSON:
    расчётная ось "x"  ==  поле ``y`` в JSON  (первая координата плоскости сечения)
    расчётная ось "y"  ==  поле ``z`` в JSON  (вторая координата плоскости сечения)

Угол отсчитывается от оси z (= расчётная ось "y") в сторону оси y JSON
(= расчётная ось "x"), в диапазоне 0..360°::

    s = (s_x, s_y):  (0, +) -> 0°,  (+, 0) -> 90°,  (0, -) -> 180°,  (-, 0) -> 270°

Группы статических моментов (см. ``optimize_counterweights``): группа 1 ``s1`` — решение
шага 1 (минимальная уравновешивающая конфигурация), группа 2 ``s2 = s - s1``. Это
декомпозиция отчёта, а не два физически отдельных груза: на щеке стоит один противовес
с итоговым ``s``. В таблице — модули векторов ``|s1|`` и ``|s2|``, в CSV — составляющие.
"""
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Union

import numpy as np

# Подписи таблицы. Язык влияет только на тексты таблицы (консоль и картинка).
LABELS = {
    "ru": {
        "title": "Статические моменты противовесов",
        "columns": ["№", "x, мм", "s_x (ось y)\nкг·м", "s_y (ось z)\nкг·м",
                    "|s|\nкг·м", "угол от оси z\nградусы"],
        "group_columns": ["|s₁| гр. 1\nкг·м", "|s₂| гр. 2\nкг·м"],
        "footnote": ("s = m·r — статический момент противовеса; s_x и s_y — составляющие по "
                     "осям y и z (расчётные «x» и «y»).\n"
                     "Угол отсчитывается от оси z в сторону оси y, 0…360°."),
        "group_footnote": ("\nГруппа 1 (s₁) — минимальная уравновешивающая конфигурация, "
                           "группа 2: s₂ = s − s₁ (модули векторов).\nЭто разложение "
                           "отчёта: на щеке стоит один противовес с итоговым s."),
        "text_title": ["Статические моменты противовесов s = m·r "
                       "(s_x — по оси y JSON, s_y — по оси z JSON;",
                       " угол — от оси z в сторону оси y, 0..360°)"],
        "text_columns": ["№", "x, мм", "s_x (y), кг·м", "s_y (z), кг·м", "|s|, кг·м",
                         "угол от z, °"],
        "text_group_columns": ["|s1| гр.1, кг·м", "|s2| гр.2, кг·м"],
        "text_group_note": ("Группа 1 (s1) — минимальная уравновешивающая конфигурация, "
                            "группа 2: s2 = s - s1 (модули векторов; разложение отчёта, "
                            "на щеке один противовес s)"),
    },
    "en": {
        "title": "Counterweight static moments (unbalance)",
        "columns": ["No.", "x, mm", "s_x (y axis)\nkg·m", "s_y (z axis)\nkg·m",
                    "|s|\nkg·m", "Phase angle\nfrom z axis, deg"],
        "group_columns": ["|s₁| group 1\nkg·m", "|s₂| group 2\nkg·m"],
        "footnote": ("s = m·r — static moment (unbalance) of a counterweight; s_x, s_y — "
                     "components along the y and z axes (calculation axes x and y).\n"
                     "The angle is measured from the z axis toward the y axis, 0…360°."),
        "group_footnote": ("\nGroup 1 (s₁) is the minimal balancing configuration, "
                           "group 2: s₂ = s − s₁ (vector moduli).\nThis is a breakdown "
                           "for the report: each web carries one counterweight with the "
                           "total s."),
        "text_title": ["Counterweight static moments s = m·r "
                       "(s_x — along the JSON y axis, s_y — along the JSON z axis;",
                       " angle — from the z axis toward the y axis, 0..360 deg)"],
        "text_columns": ["No.", "x, mm", "s_x (y), kg·m", "s_y (z), kg·m", "|s|, kg·m",
                         "angle from z, deg"],
        "text_group_columns": ["|s1| gr.1, kg·m", "|s2| gr.2, kg·m"],
        "text_group_note": ("Group 1 (s1) is the minimal balancing configuration, "
                            "group 2: s2 = s - s1 (vector moduli; a report breakdown, "
                            "each web carries one counterweight s)"),
    },
}


def get_labels(lang: str) -> dict:
    """Подписи таблицы для языка ``"ru"`` или ``"en"``."""
    if lang not in LABELS:
        raise ValueError(f"Неизвестный язык {lang!r}; доступны: {sorted(LABELS)}")
    return LABELS[lang]


# Компоненты меньше этой доли от максимальной считаются численным шумом решателя
# (например, 1e-17 вместо нуля) и обнуляются, чтобы угол не «прыгал» между 0° и 360°.
_NOISE_REL_TOL = 1e-9


@dataclass
class CounterweightRow:
    """Одна строка таблицы: один противовес."""

    number: int          # номер противовеса (с 0, по возрастанию x)
    x: float             # положение вдоль вала, мм
    s_x: float           # составляющая по расчётной оси x (поле y в JSON), кг·м
    s_y: float           # составляющая по расчётной оси y (поле z в JSON), кг·м
    s_abs: float         # модуль sqrt(s_x^2 + s_y^2), кг·м
    angle_deg: float     # угол от оси z, градусы 0..360; nan, если противовеса нет
    # Группы статических моментов (nan — разложение не задано):
    s1_x: float = math.nan   # группа 1 (уравновешивающая), составляющие
    s1_y: float = math.nan
    s2_x: float = math.nan   # группа 2 = s - s1
    s2_y: float = math.nan

    @property
    def has_groups(self) -> bool:
        return not math.isnan(self.s1_x)

    @property
    def s1_abs(self) -> float:
        return math.hypot(self.s1_x, self.s1_y)

    @property
    def s2_abs(self) -> float:
        return math.hypot(self.s2_x, self.s2_y)


def rows_from_components(x_positions: Sequence[float], s_x: Sequence[float],
                         s_y: Sequence[float],
                         s1_x: Optional[Sequence[float]] = None,
                         s1_y: Optional[Sequence[float]] = None) -> List[CounterweightRow]:
    """Строит строки таблицы по положениям и составляющим s (без привязки к расчёту).

    ``s1_x``, ``s1_y`` — необязательные составляющие группы 1; группа 2 = s - s1.
    """
    x_positions = np.asarray(x_positions, dtype=float)
    sx = np.asarray(s_x, dtype=float).copy()
    sy = np.asarray(s_y, dtype=float).copy()
    if not (len(x_positions) == len(sx) == len(sy)):
        raise ValueError("Положения, s_x и s_y должны быть одной длины")

    max_abs = float(max(np.max(np.abs(sx), initial=0.0), np.max(np.abs(sy), initial=0.0)))
    tol = _NOISE_REL_TOL * max_abs
    sx[np.abs(sx) <= tol] = 0.0
    sy[np.abs(sy) <= tol] = 0.0

    if (s1_x is None) != (s1_y is None):
        raise ValueError("s1_x и s1_y задаются вместе")
    if s1_x is None:
        groups = np.full((4, len(sx)), np.nan)
    else:
        g1x = np.asarray(s1_x, dtype=float).copy()
        g1y = np.asarray(s1_y, dtype=float).copy()
        if not (len(g1x) == len(g1y) == len(sx)):
            raise ValueError("s1_x и s1_y должны быть той же длины, что и s")
        groups = np.vstack([g1x, g1y, sx - g1x, sy - g1y])
        groups[np.abs(groups) <= tol] = 0.0

    rows = []
    for i, (x, a, b) in enumerate(zip(x_positions, sx, sy)):
        modulus = math.hypot(a, b)
        if modulus == 0.0:
            angle = float("nan")
        else:
            angle = round(math.degrees(math.atan2(a, b)), 9) % 360.0
        g = groups[:, i]
        rows.append(CounterweightRow(number=i, x=float(x), s_x=float(a), s_y=float(b),
                                     s_abs=modulus, angle_deg=angle,
                                     s1_x=float(g[0]), s1_y=float(g[1]),
                                     s2_x=float(g[2]), s2_y=float(g[3])))
    return rows


def counterweight_rows(result) -> List[CounterweightRow]:
    """Строки таблицы для результата ``analyze()`` (``BalanceResult``)."""
    x_positions = result.engine.crankshaft_length - np.asarray(result.L_counterweights)
    rx, ry = result.result_x, result.result_y
    return rows_from_components(x_positions, rx["s"], ry["s"],
                                rx.get("s_group1"), ry.get("s_group1"))


def format_counterweight_table(rows: Sequence[CounterweightRow], lang: str = "ru") -> str:
    """Текстовая таблица для печати в консоль (``lang``: ``"ru"`` или ``"en"``)."""
    labels = get_labels(lang)
    c = labels["text_columns"]
    groups = bool(rows) and all(r.has_groups for r in rows)
    header = (f"{c[0]:>3} {c[1]:>9} {c[2]:>14} {c[3]:>14} {c[4]:>10} {c[5]:>17}")
    if groups:
        g = labels["text_group_columns"]
        header += f" {g[0]:>16} {g[1]:>16}"
    lines = [*labels["text_title"], header, "-" * len(header)]
    for r in rows:
        angle = "—" if math.isnan(r.angle_deg) else f"{r.angle_deg:.2f}"
        line = (f"{r.number:>3} {r.x:>9.3f} {r.s_x:>14.4f} {r.s_y:>14.4f} "
                f"{r.s_abs:>10.4f} {angle:>17}")
        if groups:
            line += f" {r.s1_abs:>16.4f} {r.s2_abs:>16.4f}"
        lines.append(line)
    if groups:
        lines.append(labels["text_group_note"])
    return "\n".join(lines)


def save_counterweight_csv(rows: Sequence[CounterweightRow],
                           path: Union[str, Path]) -> None:
    """Сохраняет таблицу в CSV (UTF-8, разделитель — запятая; нет угла — пустая ячейка).

    Столбцы групп (``s1_*`` — группа 1, ``s2_*`` — группа 2 = s - s1) — составляющие;
    если разложение на группы не задано, ячейки пустые.
    """
    def num(v: float) -> str:
        return "" if math.isnan(v) else f"{v:.9f}"

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["number", "x_mm", "s_x_kg_m", "s_y_kg_m", "s_abs_kg_m",
                         "angle_from_z_deg", "s1_x_kg_m", "s1_y_kg_m",
                         "s2_x_kg_m", "s2_y_kg_m"])
        for r in rows:
            writer.writerow([r.number, f"{r.x:.6f}", f"{r.s_x:.9f}", f"{r.s_y:.9f}",
                             f"{r.s_abs:.9f}",
                             "" if math.isnan(r.angle_deg) else f"{r.angle_deg:.6f}",
                             num(r.s1_x), num(r.s1_y), num(r.s2_x), num(r.s2_y)])
