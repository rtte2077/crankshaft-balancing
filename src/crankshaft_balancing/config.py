"""Параметры режима работы и характерных размеров сечений вала."""
from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class EngineConfig:
    """Параметры режима работы и характерных размеров сечений вала.

    Здесь только то, что реально используется в расчётах: угловая скорость,
    вращающаяся часть массы шатуна и ширины сечений по умолчанию (нужны только для
    автоматического вывода длины вала — см. ``geometry.derive_crankshaft_dimensions``).
    Число цилиндров, порядок работы, тип развала и т.п. на сам расчёт не
    влияют: вся конфигурация вала задаётся геометрией — списком компонентов.

    ``crankshaft_length`` не задаётся вручную: он всегда вычисляется из геометрии.
    """

    angular_speed: float = 240.9      # omega, рад/с
    rod_rotating_mass: float = 1.86   # m2, кг — вращающаяся часть массы шатуна,
                                      # отнесённая к шатунной шейке

    # Ширины сечений вала, мм — нужны только для вывода длины вала.
    main_journal_width: float = 37.0
    crankpin_journal_width: float = 38.0
    cheek_width: float = 28.0

    # Необязательные индивидуальные ширины коренных шеек, мм: {номер: ширина}.
    # Номера слева направо с 0 (0 — левая концевая). Не указанные шейки
    # получают ``main_journal_width``. Ширины щёк и шатунных шеек задаются
    # у самих элементов — ``ShaftComponent.width``.
    main_journal_widths: Optional[Dict[int, float]] = None

    # Заполняется автоматически из геометрии.
    crankshaft_length: float = field(init=False, default=float("nan"))

    def __post_init__(self):
        # В JSON ключи — строки ("0"), приводим к int.
        if self.main_journal_widths is not None:
            self.main_journal_widths = {int(k): float(v)
                                        for k, v in self.main_journal_widths.items()}
