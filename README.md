# crankshaft-balancing

Расчёт балансировки коленчатого вала **любой конфигурации** (любое число цилиндров,
рядный или V-образный, любая раскладка щёк и шеек):

1. центробежные силы инерции щёк и шатунных шеек;
2. эпюра изгибающего момента `M(z)` вдоль оси вала методом сечений;
3. подбор коэффициентов противовесов линейным программированием: минимум `max|M(z)|`
   при нулевой суммарной силе и нулевом моменте на конце вала;
4. независимая проверка, что найденный минимум глобальный (двойственность, ККТ, чувствительность).

![Эпюра момента после оптимизации](docs/img/moment_y_optimized.png)

## Установка

```bash
git clone https://github.com/<ваш-аккаунт>/crankshaft-balancing.git
cd crankshaft-balancing
pip install -e .                 # numpy, scipy, matplotlib
pip install -e ".[interactive]"  # + plotly для интерактивного 3D-графика
pip install -e ".[dev]"          # + pytest
```

## Быстрый старт

### Из командной строки

```bash
crankbalance examples/example_shaft.json --verify --out output
# или: python -m crankshaft_balancing examples/example_shaft.json
```

В папке `output/` появятся графики моментов до и после оптимизации и 3D-график.
Ключи: `--step` (шаг сетки сечений, мм), `--no-plots`, `--verify`, `--interactive`.

### Из Python

```python
from crankshaft_balancing import analyze, load_shaft

engine, components = load_shaft("examples/example_shaft.json")
res = analyze(engine, components)

print(res.max_moment_y)        # max|My| после оптимизации, Н·м
print(res.result_y["s"])       # коэффициенты противовесов по оси Y
```

Полный пример — `examples/run_example.py`, интерактивная демонстрация —
`notebooks/demo.ipynb`.

## Описание своего вала

Вал описывается JSON-файлом (см. `examples/example_shaft.json`):

```json
{
  "engine": {
    "angular_speed": 240.9,
    "rod_rotating_mass": 1.86,
    "main_journal_width": 37.0,
    "crankpin_journal_width": 38.0,
    "cheek_width": 28.0
  },
  "components": [
    {"kind": "cheek",    "x": 31.469, "y": 0, "z": 26.5926, "mass": 1.8496},
    {"kind": "crankpin", "x": 65.5,   "y": 0, "z": 64.0,    "mass": 1.3565}
  ]
}
```

| Поле | Смысл | Единицы |
|---|---|---|
| `angular_speed` | угловая скорость вала | рад/с |
| `rod_rotating_mass` | вращающаяся часть массы шатуна, приведённая к шатунной шейке | кг |
| `main_journal_width`, `crankpin_journal_width`, `cheek_width` | ширины сечений, нужны только для вывода длины вала | мм |
| `kind` | `"cheek"` (щека) или `"crankpin"` (шатунная шейка) | — |
| `x, y, z` | координаты центра масс; `x` — вдоль оси вала | мм |
| `mass` | масса элемента | кг |

Элементы можно перечислять в любом порядке — они сортируются по `x`. Длина вала
не задаётся вручную, а выводится из геометрии (подробнее в [docs/theory.md](docs/theory.md)).

## Структура репозитория

```
src/crankshaft_balancing/
├── config.py        # EngineConfig — режим работы и ширины сечений
├── geometry.py      # ShaftComponent, build_geometry, вывод длины вала
├── forces.py        # центробежные силы
├── sections.py      # метод сечений, M(z)
├── optimization.py  # ЛП для подбора противовесов
├── verification.py  # двойственность, ККТ, чувствительность
├── plotting.py      # все графики
├── analysis.py      # analyze() — сквозной расчёт
├── loader.py        # чтение JSON
└── cli.py           # команда crankbalance
examples/            # пример вала и скрипт
notebooks/demo.ipynb # демонстрация с графиками
tests/               # pytest
docs/theory.md       # постановка задачи и формулы
```

## Тесты

```bash
pytest
```

## Ограничения

- Все силы считаются в плоскости YZ поперечного сечения; противовесы привязаны к
  позициям щёк (`L3 = L1` по умолчанию).
- `optimize_counterweights_history` (график сходимости) использует устаревший метод
  `simplex` из `scipy.optimize.linprog` — он удалён из новых версий SciPy. Основной расчёт
  (HiGHS) от него не зависит; при отсутствии метода будет понятная ошибка.
- В примере момент по оси X нулевой по построению геометрии (все центры масс лежат в
  плоскости `y = 0`) — для него оптимизация тривиальна.

## Лицензия

MIT — см. [LICENSE](LICENSE).
