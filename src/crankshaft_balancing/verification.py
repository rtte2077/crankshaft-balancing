"""Независимые проверки того, что найденный минимум — глобальный.

Ни одна проверка не меняет ``optimize_counterweights``:

1. Двойственный сертификат: по теореме сильной двойственности оптимумы прямой и
   двойственной задач совпадают тогда и только тогда, когда решение глобально оптимально.
2. Условия ККТ по множителям Лагранжа, которые уже вернул HiGHS.
3. Эмпирическая проверка выпуклости: вариация каждого s_i вокруг оптимума.
"""
from typing import Optional

import numpy as np
from scipy.optimize import linprog

from .optimization import CounterweightLP
from .sections import total_distribution


def verify_optimum_via_duality(primal_result: dict, lp: Optional[CounterweightLP] = None,
                               axis: str = "", tol: float = 1e-6,
                               verbose: bool = True) -> dict:
    """Строит и решает ДВОЙСТВЕННУЮ задачу и сравнивает её оптимум с прямой.

    ``lp`` по умолчанию берётся из ``primal_result["lp"]``.
    """
    lp = lp if lp is not None else primal_result["lp"]
    A_ub, b_ub = lp.A_ub, lp.b_ub
    A_eq, b_eq, c = lp.A_eq, lp.b_eq, lp.c
    n_vars = len(c)

    # Границы переменных -> дополнительные строки A_ub x <= b_ub, чтобы все
    # переменные стали свободными и двойственная задача строилась по одной формуле.
    extra_rows, extra_b = [], []
    for j, (lo, hi) in enumerate(lp.bounds):
        if hi is not None:
            row = np.zeros(n_vars)
            row[j] = 1.0
            extra_rows.append(row)
            extra_b.append(hi)
        if lo is not None:
            row = np.zeros(n_vars)
            row[j] = -1.0
            extra_rows.append(row)
            extra_b.append(-lo)

    if extra_rows:
        A_ub_full = np.vstack([A_ub, np.array(extra_rows)])
        b_ub_full = np.concatenate([b_ub, np.array(extra_b)])
    else:
        A_ub_full, b_ub_full = A_ub, b_ub

    m1, m2 = A_ub_full.shape[0], A_eq.shape[0]

    # Стационарность: A_ub^T y_ub + A_eq^T y_eq = -c
    A_dual_eq = np.hstack([A_ub_full.T, A_eq.T])
    b_dual_eq = -c
    c_dual = np.concatenate([b_ub_full, b_eq])
    dual_bounds = [(0, None)] * m1 + [(None, None)] * m2

    dual_res = linprog(c_dual, A_eq=A_dual_eq, b_eq=b_dual_eq,
                       bounds=dual_bounds, method="highs")

    primal_obj = primal_result["result"].fun
    dual_obj = -dual_res.fun
    gap = primal_obj - dual_obj
    confirmed = dual_res.success and abs(gap) < tol * max(1.0, abs(primal_obj))

    if verbose:
        print(f"[{axis}] прямая задача (primal)     = {primal_obj:.6f}")
        print(f"[{axis}] двойственная задача (dual)  = {dual_obj:.6f}")
        print(f"[{axis}] зазор двойственности        = {gap:.3e}  -> "
              f"{'ГЛОБАЛЬНЫЙ ОПТИМУМ ПОДТВЕРЖДЁН' if confirmed else 'ВНИМАНИЕ: зазор великоват'}")

    return {"primal_objective": primal_obj, "dual_objective": dual_obj,
            "gap": gap, "confirmed": bool(confirmed), "dual_result": dual_res}


def verify_kkt_conditions(primal_result: dict, lp: Optional[CounterweightLP] = None,
                          axis: str = "", tol: float = 1e-6,
                          verbose: bool = True) -> dict:
    """Проверяет условия Каруша-Куна-Таккера по множителям Лагранжа HiGHS.

    Проверяются допустимость прямой и двойственной задач, стационарность и
    дополняющая нежёсткость.
    """
    lp = lp if lp is not None else primal_result["lp"]
    res = primal_result["result"]
    A_ub, A_eq, c = lp.A_ub, lp.A_eq, lp.c

    y_ineq = res.ineqlin.marginals
    y_eq = res.eqlin.marginals
    lower_marg = res.lower.marginals
    upper_marg = res.upper.marginals
    slack = res.ineqlin.residual

    primal_feasible = bool(res.success)
    dual_feasible = bool(np.all(y_ineq <= tol))

    stationarity = c - A_ub.T @ y_ineq - A_eq.T @ y_eq - lower_marg - upper_marg
    stationarity_residual = float(np.max(np.abs(stationarity)))
    complementary_slackness = float(np.max(np.abs(y_ineq * slack)))

    confirmed = (primal_feasible and dual_feasible
                 and stationarity_residual < tol and complementary_slackness < tol)

    if verbose:
        print(f"[{axis}] допустимость прямой задачи:        {primal_feasible}")
        print(f"[{axis}] допустимость двойственной задачи:   {dual_feasible}")
        print(f"[{axis}] невязка стационарности:             {stationarity_residual:.3e}")
        print(f"[{axis}] невязка доп. нежёсткости:           {complementary_slackness:.3e}")
        print(f"[{axis}] условия ККТ выполнены -> "
              f"{'ГЛОБАЛЬНЫЙ ОПТИМУМ ПОДТВЕРЖДЁН' if confirmed else 'ВНИМАНИЕ: нарушены'}")

    return {"primal_feasible": primal_feasible, "dual_feasible": dual_feasible,
            "stationarity_residual": stationarity_residual,
            "complementary_slackness": complementary_slackness,
            "confirmed": bool(confirmed)}


def sweep_sensitivity(primal_result: dict, z, L1, F1, L2, F2, omega,
                      lp: Optional[CounterweightLP] = None,
                      span: Optional[float] = None, n_points: int = 41) -> list:
    """Варьирует каждый свободный s_i вокруг оптимума и пересчитывает max|M(z)|
    напрямую по формуле, без решателя.

    Возвращает список ``(i, s_values, max_abs_M)`` по каждому свободному s_i.
    """
    lp = lp if lp is not None else primal_result["lp"]
    s_opt = np.asarray(primal_result["s"], dtype=float)

    if span is None:
        span = max(0.5 * np.max(np.abs(s_opt)), 0.05)

    curves = []
    for i in lp.free_idx:
        sweep_vals = np.linspace(s_opt[i] - span, s_opt[i] + span, n_points)
        max_abs_M = np.empty(n_points)
        for k, val in enumerate(sweep_vals):
            s_trial = s_opt.copy()
            s_trial[i] = val
            M_trial = total_distribution(z, L1, F1, L2, F2,
                                         F3=s_trial * omega ** 2, L3=lp.L3)
            max_abs_M[k] = np.max(np.abs(M_trial))
        curves.append((i, sweep_vals, max_abs_M))
    return curves
