"""Подбор коэффициентов противовесов линейным программированием.

Момент M(z) линеен по коэффициентам ``s_i`` (третья сила F3 = s * omega^2), поэтому
``max|M(z)|`` — выпуклая функция от ``s``, и задача сводится к ЛП:

* минимизировать вспомогательную переменную ``t`` при ``-t <= M(z_k) <= t``
  для всех сечений;
* при условии баланса сил ``sum(F1) + sum(F2) + sum(F3) = 0``;
* и нулевом суммарном внешнем моменте ``M(z[-1]) = 0``.

Это гарантирует глобальный минимум (в отличие от градиентных методов на
негладкой функции ``max(|...|)``).

Две группы статических моментов (``group2_percent``)
----------------------------------------------------
Бюджет ``B`` — суммарный статический момент по оси: ``sum_i |s_i|`` (``norm="l1"``)
или ``max_i |s_i|`` (``norm="linf"``). Модуль линеаризуется вспомогательными
переменными ``u_i >= s_i``, ``u_i >= -s_i``, поэтому все задачи остаются ЛП:

1. ``B_min`` — минимальный бюджет, при котором выполняются два равенства
   (эпюра M не ограничивается) — ``minimal_static_moment``;
2. ``M_opt`` — минимально возможный ``max|M|`` (исходная задача);
3. ``B_max`` — минимальный бюджет среди решений с ``max|M| <= M_opt`` и наименьшим
   при этом колебанием эпюры (с допуском ``budget_eps``) — ``minimal_moment_budget``;
4. для процента ``p``: ``min t`` при бюджете
   ``B(p) = B_min + p/100 * (B_max - B_min)`` — ``optimize_counterweights``;
5. среди решений шага 4 (тот же max|M|, тот же бюджет) — наименьшее колебание эпюры
   ``max M - min M``.

Оптимум по max|M| обычно неединственный; шаги 3 и 5 выбирают из равноценных по
max|M| решений: сначала с наименьшим колебанием эпюры, затем самое лёгкое.

Бюджет считается по каждой оси отдельно: минимизируются ``sum|s_x|`` и ``sum|s_y|``
по отдельности, а не ``sum sqrt(s_x^2 + s_y^2)`` — последнее связало бы оси и
вывело задачу из ЛП.
"""
import warnings
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple, Union

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

    Переменные ЛП: [s_free (n_free штук), u (n_u штук), t', (hi', lo')]. ``u`` есть
    только в задачах с бюджетом статического момента: ``n_free`` переменных при
    ``norm="l1"`` (u_i >= |s_i|) или одна при ``norm="linf"`` (u >= |s_i| для всех i).
    ``hi'``, ``lo'`` — только в задачах с колебанием эпюры: ``lo <= M(z_k) <= hi``,
    колебание ``hi - lo`` (в тех же единицах ``/ omega**2``, что и ``t'``).

    Для обусловленности строки эпюры поделены на ``omega**2``, а вместо
    ``t`` (Н·мм) переменная ``t' = t / omega**2``; коэффициент цели при ней —
    ``omega**2``, так что значение цели ``min t`` по-прежнему равно max|M| в Н·мм.
    Само ``t`` — ``t_value(x)``.
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
    n_u: int = 0
    norm: str = "l1"
    budget_const: float = 0.0   # вклад зафиксированных s в бюджет (константа)
    has_range: bool = False     # есть ли переменные hi', lo' (колебание эпюры)

    @property
    def n_free(self) -> int:
        return len(self.free_idx)

    @property
    def u_slice(self) -> slice:
        return slice(self.n_free, self.n_free + self.n_u)

    @property
    def t_idx(self) -> int:
        return self.n_free + self.n_u

    @property
    def t_scale(self) -> float:
        """t = t_scale * t' (Н·мм)."""
        return self.omega ** 2

    def t_value(self, x: np.ndarray) -> float:
        """max|M| (Н·мм) из решения ЛП."""
        return float(x[self.t_idx] * self.t_scale)

    def range_value(self, x: np.ndarray) -> float:
        """Колебание эпюры hi - lo (Н·мм) из решения ЛП."""
        if not self.has_range:
            raise ValueError("В этой ЛП нет переменных колебания")
        return float((x[self.t_idx + 1] - x[self.t_idx + 2]) * self.t_scale)

    def assemble_solution(self, x: np.ndarray) -> np.ndarray:
        """Полный вектор s (с зафиксированными значениями) из решения ЛП."""
        s_full = self.s_const_vals.copy()
        for k, i in enumerate(self.free_idx):
            s_full[i] = x[k]
        return s_full


NORMS = ("l1", "linf")


def static_moment(s, norm: str = "l1") -> float:
    """Суммарный статический момент по оси: ``sum|s_i|`` (l1) или ``max|s_i|`` (linf)."""
    _check_norm(norm)
    a = np.abs(np.asarray(s, dtype=float))
    if a.size == 0:
        return 0.0
    return float(a.sum() if norm == "l1" else a.max())


def _check_norm(norm: str) -> None:
    if norm not in NORMS:
        raise ValueError(f"Неизвестная норма {norm!r}; допустимы: {NORMS}")


def _critical_rows(z: np.ndarray, *points) -> np.ndarray:
    """Номера сечений сетки ``z``, где достаточно проверять ограничения на M(z).

    M(z) линейна между соседними точками приложения сил (и концами вала), поэтому
    на каждом таком участке экстремумы M по узлам сетки — в крайних узлах участка.
    Остаются: первый и последний узлы сетки и узлы, ближайшие слева и справа к
    каждой точке приложения силы. Ограничения ``|M| <= t`` и ``lo <= M <= hi`` в
    этих узлах равносильны ограничениям во всех узлах сетки.
    """
    nodes = np.unique(np.concatenate([np.asarray(p, dtype=float).ravel() for p in points]))
    right = np.searchsorted(z, nodes, side="left")          # первый z_k >= узла
    left = np.searchsorted(z, nodes, side="right") - 1      # последний z_k <= узла
    idx = np.concatenate([[0, len(z) - 1], left, right])
    return np.unique(idx[(idx >= 0) & (idx < len(z))])


def build_counterweight_lp(z, L1, F1, L2, F2, omega,
                           s_init: Sequence[Union[float, str]],
                           s_fixed: Optional[Sequence[bool]] = None,
                           L3=None, s_bounds=None, *,
                           budget_vars: bool = False, norm: str = "l1",
                           objective: str = "t", budget: Optional[float] = None,
                           t_max: Optional[float] = None,
                           moment_rows: bool = True, range_vars: bool = False,
                           range_max: Optional[float] = None) -> CounterweightLP:
    """Строит ЛП для подбора коэффициентов противовесов (см. docstring модуля).

    Без дополнительных аргументов — исходная задача ``min t``. Остальное нужно
    для задач с бюджетом статического момента:

    budget_vars — добавить переменные ``u`` (``u >= s``, ``u >= -s``);
    norm        — ``"l1"`` (бюджет ``sum|s|``) или ``"linf"`` (``max|s|``);
    objective   — ``"t"`` (min max|M|) или ``"budget"`` (min бюджета);
    budget      — верхняя граница бюджета ``B`` (с учётом зафиксированных s);
    t_max       — верхняя граница ``t``, т.е. ``max|M| <= t_max`` (Н·мм);
    moment_rows — включать ли строки ``-t <= M(z_k) <= t`` (для шага 1 не нужны).

    Колебание эпюры ``max M - min M``:

    range_vars  — добавить переменные ``hi``, ``lo`` и строки ``lo <= M(z_k) <= hi``;
    objective   — ``"range"``: минимизировать колебание ``hi - lo``;
    range_max   — верхняя граница колебания ``hi - lo`` (Н·мм).
    """
    _check_norm(norm)
    if objective not in ("t", "budget", "range"):
        raise ValueError(f"Неизвестная цель {objective!r}; допустимы: 't', 'budget', "
                         f"'range'")
    if (objective == "budget" or budget is not None) and not budget_vars:
        raise ValueError("Бюджет требует budget_vars=True")
    if (objective == "range" or range_max is not None) and not range_vars:
        raise ValueError("Колебание эпюры требует range_vars=True")

    z = np.asarray(z, dtype=float)
    s_init = list(s_init)

    L3_use = np.asarray(L1 if L3 is None else L3, dtype=float)
    n = len(L3_use)

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

    # Строки эпюры делятся на omega^2 (t' = t / omega^2), иначе
    # коэффициенты при s ~ omega^2 * z ~ 1e7 против 1 при t — плохая обусловленность.
    w2 = omega ** 2
    # G[k, j] = H(z_k - L3[free_idx[j]])  -> (n_z, n_free);  M = w2 * (G s + M_fixed / w2)
    G = heaviside_ramp(z[:, None] - L3_use[None, free_idx])

    n_u = (n_free if norm == "l1" else 1) if budget_vars and n_free else 0
    t_idx = n_free + n_u
    hi_idx, lo_idx = t_idx + 1, t_idx + 2      # только при range_vars
    n_vars = t_idx + 1 + (2 if range_vars else 0)
    c = np.zeros(n_vars)
    if objective == "t":
        c[t_idx] = w2
    elif objective == "budget":
        c[n_free:t_idx] = 1.0
    else:
        c[hi_idx], c[lo_idx] = w2, -w2

    # Строки по эпюре — только в сечениях, где M может достигать экстремума
    # (см. _critical_rows): ЛП та же, что по всей сетке, но в сотни раз меньше.
    k = _critical_rows(z, L1, L2, L3_use)
    Gk, Mk, n_k = G[k], M_fixed[k], len(k)

    rows, rhs = [], []
    if moment_rows:
        # |M(z_k)| <= t  <=>  (M - t) / w2 <= 0  и  (-M - t) / w2 <= 0
        block = np.zeros((2 * n_k, n_vars))
        block[0:n_k, 0:n_free] = Gk
        block[n_k:, 0:n_free] = -Gk
        block[:, t_idx] = -1.0
        rows.append(block)
        rhs.append(np.concatenate([-Mk, Mk]) / w2)

    if range_vars:
        # M(z_k) <= hi  и  lo <= M(z_k)  (поделено на w2, как строки |M| <= t)
        block = np.zeros((2 * n_k, n_vars))
        block[0:n_k, 0:n_free] = Gk
        block[0:n_k, hi_idx] = -1.0
        block[n_k:, 0:n_free] = -Gk
        block[n_k:, lo_idx] = 1.0
        rows.append(block)
        rhs.append(np.concatenate([-Mk, Mk]) / w2)
        if range_max is not None:
            row = np.zeros((1, n_vars))
            row[0, hi_idx], row[0, lo_idx] = 1.0, -1.0
            rows.append(row)
            rhs.append(np.array([range_max / w2]))

    fixed_abs = np.abs(s_const_vals[[i for i in range(n) if s_fixed[i]]])
    budget_const = static_moment(fixed_abs, norm)
    if n_u:
        # u_i >= |s_i|:  s_i - u_i <= 0  и  -s_i - u_i <= 0
        abs_rows = np.zeros((2 * n_free, n_vars))
        for j in range(n_free):
            u_j = n_free + (j if norm == "l1" else 0)
            abs_rows[2 * j, j], abs_rows[2 * j, u_j] = 1.0, -1.0
            abs_rows[2 * j + 1, j], abs_rows[2 * j + 1, u_j] = -1.0, -1.0
        rows.append(abs_rows)
        rhs.append(np.zeros(2 * n_free))
        if budget is not None:
            # l1:   sum u + sum|s_fixed| <= B;   linf: u <= B (max|s_fixed| <= B
            # обеспечивается тем, что B >= B_min >= max|s_fixed|)
            row = np.zeros((1, n_vars))
            row[0, n_free:t_idx] = 1.0
            rows.append(row)
            rhs.append(np.array([budget - budget_const if norm == "l1" else budget]))

    A_ub = np.vstack(rows) if rows else np.zeros((0, n_vars))
    b_ub = np.concatenate(rhs) if rhs else np.zeros(0)

    # Баланс сил и нулевой суммарный внешний момент
    A_eq = np.zeros((2, n_vars))
    b_eq = np.zeros(2)
    # (в исходных единицах Н и Н·мм: допуск решателя на равенства — абсолютный)
    A_eq[0, 0:n_free] = w2
    b_eq[0] = -(np.sum(F1) + np.sum(F2) + np.sum(fixed_F3))
    A_eq[1, 0:n_free] = G[-1, :] * w2
    b_eq[1] = -M_fixed[-1]

    var_bounds = ([(None, None)] * n_free if s_bounds is None
                  else [tuple(s_bounds[i]) for i in free_idx])
    t_max_scaled = None if t_max is None else t_max / w2
    bounds = var_bounds + [(0, None)] * n_u + [(0, t_max_scaled)]  # u >= 0, 0 <= t <= t_max
    if range_vars:
        bounds += [(None, None), (None, None)]                       # hi, lo свободны

    return CounterweightLP(A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, c=c,
                           bounds=bounds, free_idx=free_idx, L3=L3_use,
                           s_const_vals=s_const_vals, M_fixed=M_fixed, omega=omega,
                           n_u=n_u, norm=norm, budget_const=budget_const,
                           has_range=range_vars)


def _solve(lp: CounterweightLP, axis: str, what: str = ""):
    """linprog/HiGHS по матрицам ``lp``; RuntimeError, если решения нет."""
    kwargs = dict(A_ub=lp.A_ub, b_ub=lp.b_ub, A_eq=lp.A_eq, b_eq=lp.b_eq,
                  bounds=lp.bounds, method="highs")
    res = linprog(lp.c, **kwargs)
    if res.x is None and res.status == 4:
        # «численные трудности»: presolve старых HiGHS (SciPy <= 1.13) иногда не
        # справляется с почти вырожденными ограничениями шагов выбора из равноценных
        # решений; без presolve та же задача решается.
        res = linprog(lp.c, **kwargs, options={"presolve": False})
    if res.x is None:
        raise RuntimeError(f"[{axis}] решатель не нашёл решения: {res.message}{what}")
    return res


def _tolerance(value: float, eps: float) -> float:
    """Допуск ``eps``: относительный при |value| >= 1, иначе абсолютный."""
    return eps * max(abs(value), 1.0)


def minimal_static_moment(z, L1, F1, L2, F2, omega, s_init, s_fixed=None,
                          L3=None, s_bounds=None, norm: str = "l1",
                          axis: str = "x") -> Tuple[np.ndarray, float]:
    """Шаг 1: минимальный бюджет статического момента, при котором выполняются
    баланс сил и нулевой суммарный внешний момент (эпюра M не ограничивается).

    Возвращает ``(s, B_min)``; ``s`` — эталонная «уравновешивающая» конфигурация
    (группа 1). RuntimeError, если два равенства невыполнимы при fixed/bounds.
    """
    lp = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init, s_fixed=s_fixed,
                                L3=L3, s_bounds=s_bounds, budget_vars=True, norm=norm,
                                objective="budget", moment_rows=False)
    if lp.n_free == 0:
        return lp.s_const_vals.copy(), static_moment(lp.s_const_vals, norm)
    res = _solve(lp, axis, ". Не выполняются два равенства (баланс сил и нулевой "
                           "суммарный внешний момент) при заданных fixed/bounds")
    s = lp.assemble_solution(res.x)
    return s, static_moment(s, norm)


def minimal_moment_budget(z, L1, F1, L2, F2, omega, s_init, s_fixed=None,
                          L3=None, s_bounds=None, norm: str = "l1",
                          budget_eps: float = 1e-9,
                          axis: str = "x") -> Tuple[np.ndarray, float, float]:
    """Шаги 2-3: ``M_opt`` — минимально возможный max|M| (Н·мм); затем среди решений с
    ``max|M| <= M_opt + допуск`` — наименьшее колебание эпюры ``R_opt``, и ``B_max`` —
    минимальный бюджет среди решений с ``max|M| <= M_opt`` и колебанием ``R_opt``.

    Допуск ``budget_eps`` относительный (абсолютный, в Н·мм, при ``M_opt < 1``).
    Возвращает ``(s, B_max, M_opt)``.
    """
    common = dict(s_fixed=s_fixed, L3=L3, s_bounds=s_bounds)
    lp = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init, **common)
    if lp.n_free == 0:
        s = lp.s_const_vals.copy()
        return s, static_moment(s, norm), float(np.max(np.abs(lp.M_fixed)))
    M_opt = lp.t_value(_solve(lp, axis).x)

    t_cap = M_opt + _tolerance(M_opt, budget_eps)
    lp_r = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init, **common,
                                  t_max=t_cap, range_vars=True, objective="range")
    R_opt = lp_r.range_value(_solve(lp_r, axis).x)
    lp3 = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init, **common,
                                 budget_vars=True, norm=norm, objective="budget",
                                 t_max=t_cap, range_vars=True,
                                 range_max=R_opt + _tolerance(R_opt, budget_eps))
    s = lp3.assemble_solution(_solve(lp3, axis).x)
    return s, static_moment(s, norm), M_opt


def _zero_noise(*arrays: np.ndarray) -> None:
    """Обнуляет (на месте) значения меньше 1e-9 от максимума — как в report.py."""
    max_abs = max((float(np.max(np.abs(a), initial=0.0)) for a in arrays), default=0.0)
    for a in arrays:
        a[np.abs(a) <= 1e-9 * max_abs] = 0.0


def optimize_counterweights(z, L1, F1, L2, F2, omega,
                            s_init: Sequence[Union[float, str]],
                            s_fixed: Optional[Sequence[bool]] = None,
                            L3=None, s_bounds=None,
                            axis: str = "x", verbose: bool = True,
                            group2_percent: float = 100.0, norm: str = "l1",
                            budget_eps: float = 1e-9,
                            levels: Optional[dict] = None) -> dict:
    """Подбирает коэффициенты s_i (F3_i = s_i * omega^2), минимизируя max|M(z)|
    при балансе сил, нулевом суммарном внешнем моменте и бюджете статического
    момента ``B(p) = B_min + p/100 * (B_max - B_min)``.

    Параметры:
        z, L1, F1, L2, F2 — как в ``total_distribution``.
        omega — угловая скорость, используемая в F3 = s * omega**2.
        s_init — начальные коэффициенты; элемент ``"F"`` фиксирует объект на нуле.
        s_fixed — список bool, какие s_i держать постоянными (значение из s_init).
        L3 — точки приложения F3 (по умолчанию совпадают с L1).
        s_bounds — границы [(min, max), ...] для каждого s_i.
        axis — только для подписей в выводе.
        group2_percent — доля группы 2, 0..100: 0 — минимальный суммарный
            статический момент, 100 (по умолчанию) — минимальный max|M|.
        norm — ``"l1"`` (бюджет ``sum|s_i|``) или ``"linf"`` (``max|s_i|``).
        budget_eps — допуск: в шаге 3 на max|M| (относительный; абсолютный в Н·мм,
            если M_opt < 1) и на бюджет B(p) (относительный).
        levels — результат этой же функции для той же задачи: из него берутся
            ``s_group1``, B_min, B_max, M_opt (шаги 1-3 не пересчитываются).

    Бюджет по каждой оси свой: минимизируется ``sum|s|`` по этой оси, а не
    сумма модулей векторов ``sqrt(s_x^2 + s_y^2)``.

    Среди решений с минимальным max|M| берётся решение с наименьшим колебанием эпюры
    ``max M - min M`` (при ``group2_percent = 100`` — ещё и самое лёгкое из них).

    Возвращает dict: ``s``, ``F3``, ``M`` (Н·мм), ``success``, ``result`` (linprog;
    последняя задача — минимум колебания при найденном max|M|),
    ``lp`` (CounterweightLP — нужен для проверок оптимальности), а также
    ``s_group1`` (решение шага 1), ``s_group2 = s - s_group1`` (декомпозиция
    отчёта, а не отдельные грузы), ``static_moment_sum``, ``static_moment_min``
    (B_min), ``static_moment_max`` (B_max), ``group2_percent``, ``max_moment_opt``
    (M_opt, Н·мм), ``moment_range`` (колебание итоговой эпюры, Н·мм).
    """
    if not 0.0 <= float(group2_percent) <= 100.0:
        raise ValueError(f"group2_percent должен быть в диапазоне 0..100, "
                         f"получено {group2_percent}")
    _check_norm(norm)
    z = np.asarray(z, dtype=float)
    common = dict(s_fixed=s_fixed, L3=L3, s_bounds=s_bounds)
    lp = build_counterweight_lp(z, L1, F1, L2, F2, omega, s_init, **common)

    if lp.n_free == 0:
        if verbose:
            s_all = lp.s_const_vals
            print(f"[{axis}] все коэффициенты зафиксированы — оптимизация не запускается")
            print(f"[{axis}] max|M| = {np.max(np.abs(lp.M_fixed)) * NMM_TO_NM:.4f} Н·м")
            print(f"[{axis}] баланс сил = "
                  f"{round(np.sum(F1) + np.sum(F2) + np.sum(s_all) * omega ** 2, 6)} Н")
            print(f"[{axis}] суммарный внешний момент M(z[-1]) = "
                  f"{round(lp.M_fixed[-1] * NMM_TO_NM, 6)} Н·м")
        s = lp.s_const_vals
        B = static_moment(s, norm)
        return {"s": s, "F3": s * omega ** 2, "M": lp.M_fixed, "success": True,
                "result": None, "lp": lp,
                "s_group1": s.copy(), "s_group2": np.zeros_like(s),
                "static_moment_sum": B, "static_moment_min": B, "static_moment_max": B,
                "group2_percent": float(group2_percent),
                "max_moment_opt": float(np.max(np.abs(lp.M_fixed))),
                "moment_range": float(np.ptp(lp.M_fixed))}

    args = (z, L1, F1, L2, F2, omega, s_init)
    if levels is None:
        s_group1, B_min = minimal_static_moment(*args, **common, norm=norm, axis=axis)
        _, B_max, M_opt = minimal_moment_budget(*args, **common, norm=norm,
                                                budget_eps=budget_eps, axis=axis)
    else:
        s_group1 = np.asarray(levels["s_group1"], dtype=float)
        B_min, B_max = levels["static_moment_min"], levels["static_moment_max"]
        M_opt = levels["max_moment_opt"]
    # численно B_max может оказаться чуть меньше B_min — интервал не бывает пустым
    B_max = max(B_max, B_min)
    degenerate = B_max - B_min <= _tolerance(B_max, budget_eps)

    B_p = B_min + float(group2_percent) / 100.0 * (B_max - B_min)
    budgeted = dict(**common, budget_vars=True, norm=norm, budget=B_p * (1.0 + budget_eps))
    lp4 = build_counterweight_lp(*args, **budgeted)
    T_p = lp4.t_value(_solve(lp4, axis).x)
    # Шаг 5: среди решений с max|M| = T_p (и тем же бюджетом) — наименьшее колебание.
    lp = build_counterweight_lp(*args, **budgeted, t_max=T_p + _tolerance(T_p, budget_eps),
                                range_vars=True, objective="range")
    res = _solve(lp, axis)

    s_full = lp.assemble_solution(res.x)
    s_group1 = s_group1.copy()
    s_group2 = s_full - s_group1
    _zero_noise(s_group1, s_group2)
    F3_opt = s_full * omega ** 2
    M_res = total_distribution(z, L1, F1, L2, F2, F3=F3_opt, L3=lp.L3)

    if verbose:
        print(f"[{axis}] статус: {res.message} | успех: {res.success}")
        print(f"[{axis}] группа 2: {group2_percent:g} %, норма бюджета {norm}: "
              f"B_min = {B_min:.6f}, B_max = {B_max:.6f}, B(p) = {B_p:.6f} кг·м")
        if degenerate:
            print(f"[{axis}] B_max = B_min: группа 2 пуста, все p дают одно решение")
        print(f"[{axis}] max|M| = {np.max(np.abs(M_res)) * NMM_TO_NM:.4f} Н·м "
              f"(минимально возможный {M_opt * NMM_TO_NM:.4f} Н·м), колебание эпюры "
              f"max M - min M = {np.ptp(M_res) * NMM_TO_NM:.4f} Н·м")
        print(f"[{axis}] s = {np.round(s_full, 4)}")
        print(f"[{axis}] баланс сил = "
              f"{round(np.sum(F1) + np.sum(F2) + np.sum(s_full) * omega ** 2, 6)} Н")
        print(f"[{axis}] суммарный внешний момент M(z[-1]) = "
              f"{round(M_res[-1] * NMM_TO_NM, 6)} Н·м")

    return {"s": s_full, "F3": F3_opt, "M": M_res, "success": res.success,
            "result": res, "lp": lp,
            "s_group1": s_group1, "s_group2": s_group2,
            "static_moment_sum": static_moment(s_full, norm),
            "static_moment_min": float(B_min), "static_moment_max": float(B_max),
            "group2_percent": float(group2_percent), "max_moment_opt": float(M_opt),
            "moment_range": float(np.ptp(M_res))}


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
