"""Загрузка конфигурации вала из JSON."""
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .config import EngineConfig
from .constraints import ALLOWED_KEYS, CounterweightConstraints
from .geometry import ShaftComponent


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
    block = {k: v for k, v in data.get("counterweights", {}).items() if k != "positions"}
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
