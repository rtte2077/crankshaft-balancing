"""Загрузка конфигурации вала из JSON."""
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .config import EngineConfig
from .constraints import ALLOWED_KEYS, CounterweightConstraints
from .geometry import ShaftComponent

# Ключи блока counterweights, которые не являются ограничениями по осям.
NON_CONSTRAINT_KEYS = {"positions", "group2_percent", "smooth_percent"}


def load_shaft(path: Union[str, Path]) -> Tuple[EngineConfig, List[ShaftComponent]]:
    """Читает JSON вида::

        {
          "engine": {"angular_speed": 240.9, "rod_rotating_mass": 1.86, ...},
          "components": [
            {"kind": "cheek", "x": 31.469, "y": 0, "z": 26.5926, "mass": 1.8496},
            ...
          ]
        }

    Ключи в ``engine`` необязательны — пропущенные берутся из ``EngineConfig``.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if "components" not in data:
        raise ValueError(f"В {path} нет ключа 'components'")

    engine = EngineConfig(**data.get("engine", {}))
    components = [ShaftComponent(**c) for c in data["components"]]
    return engine, components


def load_constraints(path: Union[str, Path]) -> Dict[str, CounterweightConstraints]:
    """Читает необязательный блок ``counterweights`` из того же JSON::

        "counterweights": {
          "all": {"zero": [2]},                       # обе оси
          "y":   {"fixed": {"0": 0.12}, "bounds": {"1": [-0.3, 0.3]}}
        }

    Ключи ``x``, ``y``, ``all``; внутри — ``zero``, ``fixed``, ``bounds``
    (см. ``CounterweightConstraints.from_dict``). Нет блока — пустой словарь.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    block = {k: v for k, v in data.get("counterweights", {}).items()
             if k not in NON_CONSTRAINT_KEYS}
    unknown = set(block) - ALLOWED_KEYS
    if unknown:
        raise ValueError(f"В counterweights неизвестные ключи: {sorted(unknown)}. "
                         f"Допустимы: {sorted(ALLOWED_KEYS)}")
    return {k: CounterweightConstraints.from_dict(v) for k, v in block.items()}


def load_counterweight_positions(path: Union[str, Path]) -> Optional[List[float]]:
    """Необязательный список ``counterweights.positions`` — x-координаты (мм), где
    стоят противовесы. Нет ключа — ``None`` (противовесы на щёках)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    positions = data.get("counterweights", {}).get("positions")
    return None if positions is None else [float(x) for x in positions]


def _load_percent(path: Union[str, Path], key: str) -> Optional[float]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    value = data.get("counterweights", {}).get(key)
    if value is None:
        return None
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"В {path}: counterweights.{key} должен быть "
                         f"в диапазоне 0..100, получено {value}")
    return value


def load_group2_percent(path: Union[str, Path]) -> Optional[float]:
    """Необязательный ``counterweights.group2_percent`` — доля группы 2 статических
    моментов, 0..100. Нет ключа — ``None`` (берётся значение по умолчанию, 100).
    Значение, переданное в ``analyze()`` или CLI явно, перекрывает JSON."""
    return _load_percent(path, "group2_percent")


def load_smooth_percent(path: Union[str, Path]) -> Optional[float]:
    """Необязательный ``counterweights.smooth_percent`` — степень сглаживания эпюры,
    0..100. Нет ключа — ``None`` (по умолчанию 0, без сглаживания). Значение,
    переданное в ``analyze()`` или CLI явно, перекрывает JSON."""
    return _load_percent(path, "smooth_percent")
