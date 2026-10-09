"""Загрузка конфигурации вала из JSON."""
import json
from pathlib import Path
from typing import List, Tuple, Union

from .config import EngineConfig
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
