"""Метод сечений: эпюра момента M(z) вдоль оси вала.

Момент в сечении ``z`` — сумма ``H(z - L_i) * F_i`` по всем силам ``F_i``,
приложенным в точках ``L_i`` (L отсчитывается от ПРАВОГО конца вала).
Единицы: F — Н, z и L — мм, результат — Н·мм (``NMM_TO_NM`` переводит в Н·м).
"""
from typing import Optional

import numpy as np

# Перевод Н·мм -> Н·м.
NMM_TO_NM = 1e-3


def section_grid(length: float, step: float = 0.01) -> np.ndarray:
    """Сетка координат z для метода сечений: от 0 до ``length`` с шагом ``step``."""
    return np.arange(0, length + step, step)


def heaviside_ramp(diff: np.ndarray) -> np.ndarray:
    """Функция Хевисайда: 0 при diff <= 0, иначе diff."""
    return np.where(diff > 0, diff, 0.0)


def section_contribution(z: np.ndarray, application_points: np.ndarray,
                         forces: np.ndarray) -> np.ndarray:
    """Вклад одной группы сил в результат для каждого сечения z."""
    z = np.asarray(z, dtype=float)
    application_points = np.asarray(application_points, dtype=float)
    forces = np.asarray(forces, dtype=float)

    diff = z[:, None] - application_points[None, :]   # (n_z, n_объектов)
    return heaviside_ramp(diff) @ forces


def total_distribution(z, L1, F1, L2, F2,
                       F3: Optional[np.ndarray] = None, L3=None) -> np.ndarray:
    """result(z) = H(z-L1)*F1 + H(z-L2)*F2 [+ H(z-L3)*F3]

    Если ``L3`` не задана, а ``F3`` задана, точки приложения F3 совпадают с L1.
    Результат в Н·мм.
    """
    total = section_contribution(z, L1, F1) + section_contribution(z, L2, F2)
    if F3 is not None:
        total = total + section_contribution(z, L1 if L3 is None else L3, F3)
    return total
