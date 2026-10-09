"""Подбор коэффициентов противовесов линейным программированием.

Момент M(z) линеен по коэффициентам ``s_i`` (третья сила F3 = s * omega^2), поэтому
``max|M(z)|`` — выпуклая функция от ``s``, и задача сводится к ЛП:

* минимизировать вспомогательную переменную ``t`` при ``-t <= M(z_k) <= t``
  для всех сечений;
* при условии баланса сил ``sum(F1) + sum(F2) + sum(F3) = 0``;
* и нулевом моменте на конце вала ``M(z[-1]) = 0``.

Это гарантирует глобальный минимум (в отличие от градиентных методов на
негладкой функции ``max(|...|)``).
"""
import warnings
from dataclasses import dataclass
from typing import List, Optional, Sequence, Union

import numpy as np
from scipy.optimize import linprog

from .sections import NMM_TO_NM, heaviside_ramp, total_distribution


def make_initial_coefficients(n: int, default: float = 1.0,
                              overrides: Optional[dict] = None) -> list:
    """Начальные коэффициенты s длиной n.

    ``overrides`` — словарь {индекс: значение}; значение ``"F"`` фиксирует
    коэффициент на нуле и исключает его из оптимизации. Пример::

        make_initial_coefficients(4, overrides={2: "F"}) -> [1.0, 1.0, "F", 1.0]
    """
    s_init = [default] * n
    for idx, val in (overrides or {}).items():
        if not 0 <= idx < n:
            raise IndexError(f"Индекс {idx} вне диапазона (0..{n - 1})")
        s_init[idx] = val
    return s_init


@dataclass
class CounterweightLP:
    """Матрицы ЛП для подбора противовесов + служебные данные.

    Переменные ЛП: [s_free (n_free штук), t].
    """

    A_ub: np.ndarray
    b_ub: np.ndarray
    A_eq: np.ndarray
    b_eq: np.ndarray
    c: np.ndarray
    bounds: list
    free_idx: List[int]
    L3: np.ndarray
    s_const_vals: np.ndarray
    M_fixed: np.ndarray
    omega: float

    @property
    def n_free(self) -> int:
        return len(self.free_idx)

    def assemble_solution(self, x: np.ndarray) -> np.ndarray:
        """Полный вектор s (с зафиксированными значениями) из решения ЛП."""
        s_full = self.s_const_vals.copy()
        for k, i in enumerate(self.free_idx):
            s_full[i] = x[k]
        return s_full


def build_counterweight_lp(z, L1, F1, L2, F2, omega,
                           s_init: Sequence[Union[float, str]],
                           s_fixed: Optional[Sequence[bool]] = None,
                           L3=None, s_bounds=None) -> CounterweightLP:
    """Строит ЛП для подбора коэффициентов противовесов (см. docstring модуля)."""
    z = np.asarray(z, dtype=float)
    s_init = list(s_init)

    L3_use = np.asarray(L1 if L3 is None else L3, dtype=float)
    n = len(L3_use)
    n_z = len(z)

    s_fixed = [False] * n if s_fixed is None else list(s_fixed)
    for i in range(n):
        if s_init[i] == "F":
            s_fixed[i] = True
            s_init[i] = 0.0
    s_const_vals = np.array([float(v) for v in s_init])

    if s_bounds is not None and len(s_bounds) != n:
        raise ValueError(f"Длина s_bounds ({len(s_bounds)}) должна совпадать с {n}")

    free_idx = [i for i in range(n) if not s_fixed[i]]
    n_free = len(free_idx)

    # Момент от фиксированной части (F1, F2 и зафиксированные компоненты F3)
    fixed_F3 = np.where(s_fixed, s_const_vals, 0.0) * omega ** 2
    M_fixed = total_distribution(z, L1, F1, L2, F2, F3=fixed_F3, L3=L3_use)

    # H_free[k, j] = H(z_k - L3[free_idx[j]]) * omega^2  -> (n_z, n_free)
    H_free = heaviside_ramp(z[:, None] - L3_use[None, free_idx]) * omega ** 2

    n_vars = n_free + 1
    t_idx = n_free
    c = np.zeros(n_vars)
    c[t_idx] = 1.0

    # |M(z_k)| <= t  <=>  M - t <= 0  и  -M - t <= 0
    A_ub = np.zeros((2 * n_z, n_vars))
    b_ub = np.zeros(2 * n_z)
    A_ub[0:n_z, 0:n_free] = H_free
    A_ub[0:n_z, t_idx] = -1.0
    b_ub[0:n_z] = -M_fixed
    A_ub[n_z:, 0:n_free] = -H_free
    A_ub[n_z:, t_idx] = -1.0
    b_ub[n_z:] = M_fixed

    # Баланс сил и нулевой момент на конце вала
    A_eq = np.zeros((2, n_vars))
    b_eq = np.zeros(2)
    A_eq[0, 0:n_free] = omega ** 2
    b_eq[0] = -(np.sum(F1) + np.sum(F2) + np.sum(fixed_F3))
    A_eq[1, 0:n_free] = H_free[-1, :]
    b_eq[1] = -M_fixed[-1]

    var_bounds = ([(None, None)] * n_free if s_bounds is None
                  else [tuple(s_bounds[i]) for i in free_idx])
    bounds = var_bounds + [(0, None)]  # t >= 0

    return CounterweightLP(A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, c=c,
                           bounds=bounds, free_idx=free_idx, L3=L3_use,
                           s_const_vals=s_const_vals, M_fixed=M_fixed, omega=omega)


def optimize_counterweights(z, L1, F1, L2, F2, omega,
                            s_init: Sequence[Union[float, str]],
                            s_fixed: Optional[Sequence[bool]] = None,
                            L3=None, s_bounds=None,
                            axis: str = "x", verbose: bool = True) -> dict:
    """Подбирает коэффициенты s_i (F3_i = s_i * omega^2), минимизируя max|M(z)|
    при балансе сил и нулевом моменте на конце вала.

    Параметры:
        z, L1, F1, L2, F2 — как в ``total_distribution``.
        omega — угловая скорость, используемая в F3 = s * omega**2.
        s_init — начальные коэффициенты; элемент ``"F"`` фиксирует объект на нуле.
        s_fixed — список bool, какие s_i держать постоянными (значение из s_init).
        L3 — точки приложения F3 (по умолчанию совпадают с L1).
        s_bounds — границы [(min, max), ...] для каждого s_i.
        axis — только для подписей в выводе.

    Возвращает dict: ``s``, ``F3``, ``M`` (Н·мм), ``success``, ``result`` (linprog)
    и ``lp`` (CounterweightLP — нужен для проверок оптимальности).
    """
    z = np.asarray(z, dtype=float)
    lp = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init,
                                s_fixed=s_fixed, L3=L3, s_bounds=s_bounds)

    if lp.n_free == 0:
        if verbose:
            s_all = lp.s_const_vals
            print(f"[{axis}] все коэффициенты зафиксированы — оптимизация не запускается")
            print(f"[{axis}] max|M| = {np.max(np.abs(lp.M_fixed)) * NMM_TO_NM:.4f} Н·м")
            print(f"[{axis}] баланс сил = "
                  f"{round(np.sum(F1) + np.sum(F2) + np.sum(s_all) * omega ** 2, 6)} Н")
            print(f"[{axis}] момент на конце M(z[-1]) = "
                  f"{round(lp.M_fixed[-1] * NMM_TO_NM, 6)} Н·м")
        return {"s": lp.s_const_vals, "F3": lp.s_const_vals * omega ** 2,
                "M": lp.M_fixed, "success": True, "result": None, "lp": lp}

    res = linprog(lp.c, A_ub=lp.A_ub, b_ub=lp.b_ub, A_eq=lp.A_eq, b_eq=lp.b_eq,
                  bounds=lp.bounds, method="highs")
    if res.x is None:
        raise RuntimeError(f"[{axis}] решатель не нашёл решения: {res.message}")

    s_full = lp.assemble_solution(res.x)
    F3_opt = s_full * omega ** 2
    M_opt = total_distribution(z, L1, F1, L2, F2, F3=F3_opt, L3=lp.L3)

    if verbose:
        print(f"[{axis}] статус: {res.message} | успех: {res.success}")
        print(f"[{axis}] max|M| = {np.max(np.abs(M_opt)) * NMM_TO_NM:.4f} Н·м")
        print(f"[{axis}] s = {np.round(s_full, 4)}")
        print(f"[{axis}] баланс сил = "
              f"{round(np.sum(F1) + np.sum(F2) + np.sum(s_full) * omega ** 2, 6)} Н")
        print(f"[{axis}] момент на конце M(z[-1]) = "
              f"{round(M_opt[-1] * NMM_TO_NM, 6)} Н·м")

    return {"s": s_full, "F3": F3_opt, "M": M_opt, "success": res.success,
            "result": res, "lp": lp}


def optimize_counterweights_history(z, L1, F1, L2, F2, omega,
                                    s_init, s_fixed=None, L3=None,
                                    s_bounds=None, axis: str = "",
                                    maxiter: int = 5000,
                                    verbose: bool = True) -> dict:
    """Как ``optimize_counterweights``, но возвращает историю итераций
    (max|M(z)| на каждом шаге) — для графика сходимости.

    Использует устаревший метод ``simplex`` (единственный в ``scipy.linprog`` с
    поддержкой callback), он на порядки медленнее HiGHS — берите грубую сетку
    ``z`` (шаг ~2 мм). Если установленный SciPy уже не содержит этот метод,
    будет выброшено RuntimeError.
    """
    z = np.asarray(z, dtype=float)
    lp = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init,
                                s_fixed=s_fixed, L3=L3, s_bounds=s_bounds)
    if lp.n_free == 0:
        raise ValueError("Нечего оптимизировать — все s зафиксированы")

    history = {"iteration": [], "max_abs_M": []}

    def _callback(res):
        s_cur = lp.assemble_solution(res["x"])
        M_cur = total_distribution(z, L1, F1, L2, F2,
                                   F3=s_cur * omega ** 2, L3=lp.L3)
        history["iteration"].append(int(res["nit"]))
        history["max_abs_M"].append(float(np.max(np.abs(M_cur))))

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            res = linprog(lp.c, A_ub=lp.A_ub, b_ub=lp.b_ub, A_eq=lp.A_eq,
                          b_eq=lp.b_eq, bounds=lp.bounds, method="simplex",
                          callback=_callback,
                          options={"maxiter": maxiter, "autoscale": True})
    except (ValueError, TypeError) as exc:
        raise RuntimeError(
            "Установленная версия SciPy не поддерживает linprog(method='simplex') "
            "с callback. История итераций недоступна; основной расчёт "
            "(optimize_counterweights, HiGHS) от этого не зависит."
        ) from exc

    s_full = lp.assemble_solution(res.x)
    F3_opt = s_full * omega ** 2
    M_opt = total_distribution(z, L1, F1, L2, F2, F3=F3_opt, L3=lp.L3)

    if verbose:
        print(f"[{axis}] статус: {res.message} | успех: {res.success} | итераций: {res.nit}")
        print(f"[{axis}] итоговый max|M| (грубая сетка) = "
              f"{np.max(np.abs(M_opt)) * NMM_TO_NM:.4f} Н·м")

    return {"s": s_full, "F3": F3_opt, "M": M_opt, "success": res.success,
            "history": history, "result": res, "lp": lp}
