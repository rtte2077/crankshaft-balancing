"""Центробежные силы инерции щёк и шатунных шеек."""
import numpy as np


def radial_distance(centers: np.ndarray) -> np.ndarray:
    """Расстояние от оси вала до центра масс (без координаты x вдоль оси)."""
    centers = np.asarray(centers)
    return np.linalg.norm(centers[:, 1:], axis=1)


def distance_from_right_end(centers: np.ndarray, shaft_length: float) -> np.ndarray:
    """Расстояние L от ПРАВОГО конца вала до каждой точки (нужно для H(z - L))."""
    return shaft_length - np.asarray(centers)[:, 0]


def centrifugal_force_magnitude(radius: np.ndarray, mass: np.ndarray,
                                omega: float) -> np.ndarray:
    """Модуль центробежной силы F = m * r * omega^2, Н.

    ``radius`` — мм (переводится в м), ``mass`` — кг, ``omega`` — рад/с.
    """
    return radius * mass * 1e-3 * omega ** 2


def centrifugal_force_vectors(centers: np.ndarray, magnitudes: np.ndarray) -> np.ndarray:
    """2D-векторы сил (в плоскости YZ), сонаправленные с центром масс каждого объекта."""
    centers = np.asarray(centers)
    yz = centers[:, 1:3]
    unit_vectors = yz / np.linalg.norm(yz, axis=1, keepdims=True)
    return unit_vectors * np.asarray(magnitudes).reshape(-1, 1)
