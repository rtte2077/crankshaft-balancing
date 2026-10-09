"""Ограничения на противовесы: убрать противовес на позиции, зафиксировать
заданное вручную значение, ограничить диапазон.

Противовесы стоят на щёках; номер противовеса — это номер щеки при нумерации
вдоль вала по возрастанию ``x``, начиная с 0.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from .optimization import make_initial_coefficients

Bound = Tuple[Optional[float], Optional[float]]
AXES = ("x", "y")
ALLOWED_KEYS = {"x", "y", "all"}


def _parse_value(v) -> float:
    """Число или "F" (= 0, как в make_initial_coefficients)."""
    return 0.0 if v == "F" else float(v)


@dataclass
class CounterweightConstraints:
    """Ограничения на коэффициенты противовесов ``s`` для одной оси.

    fixed:  {номер: значение} — коэффициент фиксируется и не оптимизируется
            (0 означает «противовеса на этой щеке нет»). Фиксированные значения
            участвуют в балансе сил и момента: оптимизатор подстраивает остальные.
    bounds: {номер: (min, max)} — диапазон для оптимизируемого коэффициента;
            ``None`` вместо границы означает «без ограничения с этой стороны».
    """

    fixed: Dict[int, float] = field(default_factory=dict)
    bounds: Dict[int, Bound] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Mapping) -> "CounterweightConstraints":
        """Из словаря вида (так же выглядит блок в JSON)::

            {"zero": [2, 5],                    # убрать противовесы на щеках 2 и 5
             "fixed": {"0": 0.12},              # зафиксировать значение
             "bounds": {"1": [-0.3, 0.3]}}      # ограничить диапазон
        """
        unknown = set(data) - {"zero", "fixed", "bounds"}
        if unknown:
            raise ValueError(f"Неизвестные ключи в ограничениях: {sorted(unknown)}. "
                             f"Допустимы: zero, fixed, bounds")
        fixed: Dict[int, float] = {}
        for i in data.get("zero", []):
            fixed[int(i)] = 0.0
        for i, v in data.get("fixed", {}).items():
            i, v = int(i), _parse_value(v)
            if i in fixed and fixed[i] != v:
                raise ValueError(f"Щека {i}: противоречащие значения в zero и fixed")
            fixed[i] = v
        bounds = {}
        for i, pair in data.get("bounds", {}).items():
            if len(pair) != 2:
                raise ValueError(f"Щека {i}: bounds должен быть парой [min, max]")
            lo = None if pair[0] is None else float(pair[0])
            hi = None if pair[1] is None else float(pair[1])
            bounds[int(i)] = (lo, hi)
        return cls(fixed=fixed, bounds=bounds)

    def merged_with(self, other: "CounterweightConstraints") -> "CounterweightConstraints":
        """Объединение; значения из ``other`` перекрывают текущие."""
        fixed = {**self.fixed, **other.fixed}
        # если щека зафиксирована в other, СТАРЫЕ границы (из self) для неё теряют
        # смысл; границы, заданные в самом other, оставляем — validate() поймает
        # противоречие "и fixed, и bounds" на одной щеке
        old_bounds = {i: b for i, b in self.bounds.items() if i not in other.fixed}
        bounds = {**old_bounds, **other.bounds}
        return CounterweightConstraints(fixed=fixed, bounds=bounds)

    def validate(self, n: int) -> None:
        """Проверка номеров и границ для ``n`` противовесов."""
        for i in list(self.fixed) + list(self.bounds):
            if not 0 <= i < n:
                raise IndexError(
                    f"Номер щеки {i} вне диапазона 0..{n - 1} (всего щёк: {n})")
        for i, (lo, hi) in self.bounds.items():
            if i in self.fixed:
                raise ValueError(f"Щека {i}: нельзя одновременно фиксировать "
                                 f"значение и задавать границы")
            if lo is not None and hi is not None and lo > hi:
                raise ValueError(f"Щека {i}: min ({lo}) больше max ({hi})")

    def to_solver_args(self, n: int) -> Tuple[list, List[bool], Optional[List[Bound]]]:
        """(s_init, s_fixed, s_bounds) для ``optimize_counterweights``."""
        self.validate(n)
        s_init = make_initial_coefficients(n, overrides=dict(self.fixed))
        s_fixed = [i in self.fixed for i in range(n)]
        s_bounds = ([self.bounds.get(i, (None, None)) for i in range(n)]
                    if self.bounds else None)
        return s_init, s_fixed, s_bounds


def resolve_constraints(constraints: Optional[Mapping[str, CounterweightConstraints]],
                        axis: str) -> CounterweightConstraints:
    """Ограничения для оси: сначала ``"all"`` (обе оси), затем перекрытие ``axis``."""
    constraints = constraints or {}
    unknown = set(constraints) - ALLOWED_KEYS
    if unknown:
        raise ValueError(f"Неизвестные оси в ограничениях: {sorted(unknown)}. "
                         f"Допустимы: {sorted(ALLOWED_KEYS)}")
    base = constraints.get("all", CounterweightConstraints())
    return base.merged_with(constraints.get(axis, CounterweightConstraints()))
