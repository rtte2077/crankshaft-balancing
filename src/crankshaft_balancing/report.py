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
"""
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Union

import numpy as np

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


def rows_from_components(x_positions: Sequence[float], s_x: Sequence[float],
                         s_y: Sequence[float]) -> List[CounterweightRow]:
    """Строит строки таблицы по положениям и составляющим s (без привязки к расчёту)."""
    x_positions = np.asarray(x_positions, dtype=float)
    sx = np.asarray(s_x, dtype=float).copy()
    sy = np.asarray(s_y, dtype=float).copy()
    if not (len(x_positions) == len(sx) == len(sy)):
        raise ValueError("Положения, s_x и s_y должны быть одной длины")

    max_abs = float(max(np.max(np.abs(sx), initial=0.0), np.max(np.abs(sy), initial=0.0)))
    tol = _NOISE_REL_TOL * max_abs
    sx[np.abs(sx) <= tol] = 0.0
    sy[np.abs(sy) <= tol] = 0.0

    rows = []
    for i, (x, a, b) in enumerate(zip(x_positions, sx, sy)):
        modulus = math.hypot(a, b)
        if modulus == 0.0:
            angle = float("nan")
        else:
            angle = round(math.degrees(math.atan2(a, b)), 9) % 360.0
        rows.append(CounterweightRow(number=i, x=float(x), s_x=float(a), s_y=float(b),
                                     s_abs=modulus, angle_deg=angle))
    return rows


def counterweight_rows(result) -> List[CounterweightRow]:
    """Строки таблицы для результата ``analyze()`` (``BalanceResult``)."""
    x_positions = result.engine.crankshaft_length - np.asarray(result.L_counterweights)
    return rows_from_components(x_positions, result.result_x["s"], result.result_y["s"])


def format_counterweight_table(rows: Sequence[CounterweightRow]) -> str:
    """Текстовая таблица для печати в консоль."""
    header = (f"{'№':>3} {'x, мм':>9} {'s_x (y), кг·м':>14} {'s_y (z), кг·м':>14} "
              f"{'|s|, кг·м':>10} {'угол от z, °':>13}")
    lines = ["Статические моменты противовесов s = m·r "
             "(s_x — по оси y JSON, s_y — по оси z JSON;",
             " угол — от оси z в сторону оси y, 0..360°)",
             header, "-" * len(header)]
    for r in rows:
        angle = "—" if math.isnan(r.angle_deg) else f"{r.angle_deg:.2f}"
        lines.append(f"{r.number:>3} {r.x:>9.3f} {r.s_x:>14.4f} {r.s_y:>14.4f} "
                     f"{r.s_abs:>10.4f} {angle:>13}")
    return "\n".join(lines)


def save_counterweight_csv(rows: Sequence[CounterweightRow],
                           path: Union[str, Path]) -> None:
    """Сохраняет таблицу в CSV (UTF-8, разделитель — запятая; нет угла — пустая ячейка)."""
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["number", "x_mm", "s_x_kg_m", "s_y_kg_m", "s_abs_kg_m",
                         "angle_from_z_deg"])
        for r in rows:
            writer.writerow([r.number, f"{r.x:.6f}", f"{r.s_x:.9f}", f"{r.s_y:.9f}",
                             f"{r.s_abs:.9f}",
                             "" if math.isnan(r.angle_deg) else f"{r.angle_deg:.6f}"])
