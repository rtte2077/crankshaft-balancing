"""Графики: эпюры момента, 3D-момент, сходимость, чувствительность."""
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np

from .sections import NMM_TO_NM


def _finish(fig, save_path: Optional[str], show: bool):
    """Общий финал: сохранить (если нужно) и/или показать, вернуть fig."""
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_distribution(z: np.ndarray, result: np.ndarray, shaft_length: float,
                      groups: Optional[Sequence[dict]] = None,
                      title: str = "Результат вдоль вала",
                      ylabel: str = "Результат",
                      save_path: Optional[str] = None, show: bool = True,
                      xlabel: str = ("Координата вдоль вала (0 = левый конец, "
                                     "L0 = правый конец)"),
                      curve_label: str = "Результат"):
    """График result(z) в реальных координатах вала (0 слева, shaft_length справа).

    ``z`` отсчитывается от ПРАВОГО конца, поэтому x = shaft_length - z.
    ``groups`` — список словарей {"L":.., "F":.., "name":..} с точками приложения
    сил, которые подписываются на графике. ``xlabel`` и ``curve_label`` — подписи
    оси и кривой (по умолчанию русские; для английского графика передайте свои).
    """
    z = np.asarray(z, dtype=float)
    x = shaft_length - z

    fig = plt.figure(figsize=(11, 6))
    plt.plot(x, result, linewidth=2, color="#2563eb", label=curve_label, zorder=3)

    if groups:
        colors = ["#dc2626", "#16a34a", "#9333ea", "#ea580c"]

        y_min, y_max = min(np.min(result), 0), np.max(result)
        y_range = (y_max - y_min) or 1.0
        plt.ylim(y_min - 0.05 * y_range, y_max + 0.35 * y_range)

        label_top = y_max + 0.32 * y_range
        label_bottom = y_max + 0.03 * y_range
        row_ys = np.linspace(label_top, label_bottom, 3)

        points = []
        for gi, group in enumerate(groups):
            color = colors[gi % len(colors)]
            name = group.get("name", f"group{gi + 1}")
            for li, f_val in zip(group["L"], group.get("F", [None] * len(group["L"]))):
                points.append((shaft_length - li, f_val, color, name))
        points.sort(key=lambda p: p[0])

        for idx, (x_li, f_val, color, name) in enumerate(points):
            plt.axvline(x=x_li, color=color, linestyle="--", alpha=0.5, zorder=1)
            label = f"{name}\nx={round(x_li, 4)}"
            if f_val is not None:
                label += f"\nF={round(float(f_val), 4)}"
            plt.text(x_li, row_ys[idx % 3], label, color=color, fontsize=7.5,
                     ha="center", va="top",
                     bbox=dict(facecolor="white", edgecolor=color,
                               boxstyle="round,pad=0.2", alpha=0.9),
                     zorder=4)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend(loc="lower right")
    plt.tight_layout()
    return _finish(fig, save_path, show)


def combine_moment_results(result_x: dict, result_y: dict) -> dict:
    """Объединяет результаты оптимизации по X и Y; Mx, My переводятся в Н·м."""
    sx = np.asarray(result_x["s"], dtype=float)
    sy = np.asarray(result_y["s"], dtype=float)
    Mx = np.asarray(result_x["M"], dtype=float) * NMM_TO_NM
    My = np.asarray(result_y["M"], dtype=float) * NMM_TO_NM

    if len(sx) != len(sy):
        raise ValueError(f"Длина sx ({len(sx)}) и sy ({len(sy)}) должна совпадать")
    if len(Mx) != len(My):
        raise ValueError("Mx и My должны быть одной длины (тот же массив z)")

    return {
        "s_combined": np.sqrt(sx ** 2 + sy ** 2),
        "s_vectors": np.column_stack((sx, sy)),
        "M_combined": np.sqrt(Mx ** 2 + My ** 2),
        "Mx": Mx, "My": My,
    }


def plot_moment_3d(z, result_x: dict, result_y: dict,
                   shaft_length: Optional[float] = None,
                   title: str = "Суммарный момент вдоль вала",
                   interactive: bool = False,
                   save_path: Optional[str] = None, show: bool = True,
                   verbose: bool = True) -> dict:
    """3D-график вектора момента (Mx, My) вдоль вала.

    interactive=False -> matplotlib (save_path — картинка, по умолчанию не сохраняется).
    interactive=True  -> plotly (нужен пакет ``plotly``), save_path — .html
    (по умолчанию ``moment_3d.html``).

    Возвращает dict с ``s_combined``, ``s_vectors``, ``M_combined``.
    """
    combined = combine_moment_results(result_x, result_y)
    Mx, My, M_combined = combined["Mx"], combined["My"], combined["M_combined"]
    z = np.asarray(z, dtype=float)
    x_coord = shaft_length - z if shaft_length is not None else z

    if not interactive:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        fig = plt.figure(figsize=(11, 8))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(x_coord, Mx, My, color="#2563eb", linewidth=2, label="(Mx, My)")

        y_min, x_min = np.min(My), np.min(Mx)
        ax.plot(x_coord, Mx, np.full_like(x_coord, y_min), color="#dc2626", alpha=0.3)
        ax.plot(x_coord, np.full_like(x_coord, x_min), My, color="#16a34a", alpha=0.3)

        ax.set_xlabel("Координата вдоль вала, мм")
        ax.set_ylabel("Mx, Н·м")
        ax.set_zlabel("My, Н·м")
        ax.set_title(title)
        ax.legend()
        plt.tight_layout()
        _finish(fig, save_path, show)
    else:
        try:
            import plotly.graph_objects as go
        except ImportError as exc:
            raise ImportError("Для интерактивного графика установите plotly: "
                              "pip install 'crankshaft-balancing[interactive]'") from exc

        save_path = save_path or "moment_3d.html"
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(
            x=x_coord, y=Mx, z=My, mode="lines",
            line=dict(color=M_combined, colorscale="Viridis", width=6,
                      colorbar=dict(title="|M|")),
            name="(Mx, My)",
            hovertemplate=("Координата: %{x:.3f} мм<br>"
                           "Mx: %{y:.3f} Н·м<br>My: %{z:.3f} Н·м<extra></extra>"),
        ))
        y_min, x_min = np.min(My), np.min(Mx)
        fig.add_trace(go.Scatter3d(x=x_coord, y=Mx, z=np.full_like(x_coord, y_min),
                                   mode="lines", line=dict(color="red", width=2),
                                   opacity=0.3, name="проекция Mx", hoverinfo="skip"))
        fig.add_trace(go.Scatter3d(x=x_coord, y=np.full_like(x_coord, x_min), z=My,
                                   mode="lines", line=dict(color="green", width=2),
                                   opacity=0.3, name="проекция My", hoverinfo="skip"))
        fig.update_layout(title=title,
                          scene=dict(xaxis_title="Координата вдоль вала, мм",
                                     yaxis_title="Mx, Н·м", zaxis_title="My, Н·м"),
                          width=900, height=700)
        fig.write_html(save_path)
        print(f"Интерактивный 3D график сохранён: {save_path}")

    if verbose:
        print(f"max(M_combined) = {np.max(M_combined):.4f} Н·м")
        print(f"s_combined (по объектам) = {np.round(combined['s_combined'], 4)}")
        print("s_vectors ([sx, sy] по объектам):")
        print(np.round(combined["s_vectors"], 4))

    return {k: combined[k] for k in ("s_combined", "s_vectors", "M_combined")}


def plot_convergence(history: dict, axis: str = "", scale: float = NMM_TO_NM,
                     save_path: Optional[str] = None, show: bool = True):
    """График max|M(z)|, Н·м, в зависимости от номера итерации LP-решателя."""
    iterations = history["iteration"]
    values = np.array(history["max_abs_M"]) * scale

    fig = plt.figure(figsize=(9, 5))
    plt.plot(iterations, values, marker="o", markersize=3, linewidth=1.5, color="#2563eb")
    plt.xlabel("Итерация симплекс-метода")
    plt.ylabel("max|M(z)|, Н·м")
    plt.title("Сходимость оптимизации" + (f" (ось {axis})" if axis else ""))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    return _finish(fig, save_path, show)


def plot_sensitivity(curves: list, s_opt: np.ndarray, axis: str = "",
                     save_path: Optional[str] = None, show: bool = True):
    """Сетка «чаш» max|M| от каждого s_i — результат ``sweep_sensitivity``."""
    n = len(curves)
    n_cols = min(4, n)
    n_rows = int(np.ceil(n / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.4 * n_cols, 3 * n_rows),
                             squeeze=False)
    for idx, (i, sweep_vals, max_abs_M) in enumerate(curves):
        ax = axes[idx // n_cols][idx % n_cols]
        values_nm = max_abs_M * NMM_TO_NM
        ax.plot(sweep_vals, values_nm, color="#2563eb", linewidth=1.8)
        ax.axvline(s_opt[i], color="#dc2626", linestyle="--", linewidth=1)
        opt_value = np.interp(s_opt[i], sweep_vals, values_nm)
        ax.scatter([s_opt[i]], [opt_value], color="#dc2626", zorder=5, s=30)
        ax.set_title(f"s[{i}]", fontsize=10)
        ax.grid(True, alpha=0.3)

    for idx in range(len(curves), n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].axis("off")

    fig.suptitle(f"Чувствительность max|M| к каждому s_i (ось {axis})"
                 if axis else "Чувствительность max|M| к каждому s_i")
    fig.supxlabel("Значение коэффициента s_i (красная линия — найденный оптимум)")
    fig.supylabel("max|M|, Н·м")
    plt.tight_layout()
    return _finish(fig, save_path, show)


def plot_counterweight_table(rows: Sequence, title: Optional[str] = None,
                             save_path: Optional[str] = None, show: bool = True,
                             lang: str = "ru"):
    """Таблица статических моментов противовесов в виде картинки.

    ``rows`` — результат ``counterweight_rows(res)``. Формат файла определяется
    расширением ``save_path`` (png, pdf, svg). Строки с убранным противовесом
    (``s = 0``) выделены серым, в столбце угла у них прочерк.
    ``lang`` — язык подписей: ``"ru"`` или ``"en"``; ``title`` — свой заголовок.
    """
    import math

    from .report import get_labels

    labels = get_labels(lang)
    groups = bool(rows) and all(r.has_groups for r in rows)
    cells = []
    for r in rows:
        angle = "—" if math.isnan(r.angle_deg) else f"{r.angle_deg:.2f}"
        cells.append([str(r.number), f"{r.x:.3f}", f"{r.s_x:+.4f}", f"{r.s_y:+.4f}",
                      f"{r.s_abs:.4f}", angle])
        if groups:
            cells[-1] += [f"{r.s1_abs:.4f}", f"{r.s2_abs:.4f}"]

    columns = labels["columns"] + (labels["group_columns"] if groups else [])
    footnote = labels["footnote"] + (labels["group_footnote"] if groups else "")
    widths = [0.07, 0.15, 0.2, 0.2, 0.17, 0.21] + ([0.17, 0.17] if groups else [])
    total = sum(widths)
    widths = [w / total for w in widths]
    extra_lines = footnote.count("\n") - 1

    n_rows = max(len(cells), 3)               # минимум 3 строки высоты — иначе шапка сжимается
    height = 1.58 + 0.26 * n_rows + 0.15 * extra_lines   # заголовок + таблица + пояснение
    fig, ax = plt.subplots(figsize=(11.5 if groups else 9, height))
    ax.axis("off")

    table = ax.table(cellText=cells, colLabels=columns, cellLoc="center",
                     loc="upper center", colWidths=widths)
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.6)

    header_color, zebra_color, zero_color = "#2563eb", "#eff4ff", "#e5e7eb"
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#cbd5e1")
        if row == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(color="white", fontweight="bold")
            cell.set_height(cell.get_height() * 1.5)
            continue
        item = rows[row - 1]
        if item.s_abs == 0:
            cell.set_facecolor(zero_color)
            cell.set_text_props(color="#6b7280")
        else:
            cell.set_facecolor(zebra_color if row % 2 == 0 else "white")

    fig.suptitle(title or labels["title"], fontsize=14, fontweight="bold",
                 y=1 - 0.12 / height)
    fig.text(0.5, 0.1 / height, footnote, ha="center", va="bottom",
             fontsize=8.5, color="#475569")
    fig.subplots_adjust(top=1 - 0.55 / height, bottom=(0.62 + 0.15 * extra_lines) / height,
                        left=0.03, right=0.97)
    return _finish(fig, save_path, show)


def plot_group2_tradeoff(rows: Sequence[dict], norm: str = "l1",
                         save_path: Optional[str] = None, show: bool = True):
    """max|M| (левая ось) и суммарный статический момент (правая ось) от доли
    группы 2 ``p``. ``rows`` — результат ``analysis.group2_sweep``."""
    p = [r["group2_percent"] for r in rows]
    fig, ax_m = plt.subplots(figsize=(9, 5))
    ax_b = ax_m.twinx()
    m_color, b_color = "#2563eb", "#ea580c"
    ax_m.plot(p, [r["max_moment_y"] for r in rows], "-o", color=m_color,
              label="max|M|, ось Y")
    ax_m.plot(p, [r["max_moment_x"] for r in rows], "--s", color=m_color, alpha=0.6,
              label="max|M|, ось X")
    ax_b.plot(p, [r["static_moment_y"] for r in rows], "-o", color=b_color,
              label="статический момент, ось Y")
    ax_b.plot(p, [r["static_moment_x"] for r in rows], "--s", color=b_color, alpha=0.6,
              label="статический момент, ось X")
    budget = "Σ|s_i|" if norm == "l1" else "max|s_i|"
    ax_m.set_xlabel("Доля группы 2 p, %")
    ax_m.set_ylabel("max|M|, Н·м", color=m_color)
    ax_b.set_ylabel(f"Суммарный статический момент {budget}, кг·м", color=b_color)
    ax_m.tick_params(axis="y", colors=m_color)
    ax_b.tick_params(axis="y", colors=b_color)
    ax_m.grid(True, alpha=0.3)
    lines = ax_m.get_lines() + ax_b.get_lines()
    ax_m.legend(lines, [ln.get_label() for ln in lines], loc="upper center", fontsize=9)
    ax_m.set_title("Компромисс: max|M| и статический момент противовесов от p")
    fig.tight_layout()
    return _finish(fig, save_path, show)


def plot_smooth_tradeoff(rows: Sequence[dict], save_path: Optional[str] = None,
                         show: bool = True):
    """max|M| (левая ось) и размах эпюры max M - min M (правая ось) от степени
    сглаживания ``q``. ``rows`` — результат ``analysis.smooth_sweep``."""
    q = [r["smooth_percent"] for r in rows]
    fig, ax_m = plt.subplots(figsize=(9, 5))
    ax_r = ax_m.twinx()
    m_color, r_color = "#2563eb", "#16a34a"
    ax_m.plot(q, [r["max_moment_y"] for r in rows], "-o", color=m_color,
              label="max|M|, ось Y")
    ax_m.plot(q, [r["max_moment_x"] for r in rows], "--s", color=m_color, alpha=0.6,
              label="max|M|, ось X")
    ax_r.plot(q, [r["moment_range_y"] for r in rows], "-o", color=r_color,
              label="размах max M − min M, ось Y")
    ax_r.plot(q, [r["moment_range_x"] for r in rows], "--s", color=r_color, alpha=0.6,
              label="размах max M − min M, ось X")
    ax_m.set_xlabel("Степень сглаживания q, %")
    ax_m.set_ylabel("max|M|, Н·м", color=m_color)
    ax_r.set_ylabel("Размах эпюры max M − min M, Н·м", color=r_color)
    ax_m.tick_params(axis="y", colors=m_color)
    ax_r.tick_params(axis="y", colors=r_color)
    ax_m.grid(True, alpha=0.3)
    lines = ax_m.get_lines() + ax_r.get_lines()
    ax_m.legend(lines, [ln.get_label() for ln in lines], loc="center left", fontsize=9)
    ax_m.set_title("Сглаживание: размах эпюры и max|M| от q")
    fig.tight_layout()
    return _finish(fig, save_path, show)


def plot_smooth_moments(curves: Sequence[tuple], shaft_length: float,
                        save_path: Optional[str] = None, show: bool = True,
                        axes: Sequence[str] = ("y", "x"), title: Optional[str] = None,
                        xlabel: Optional[str] = None, ylabel: str = "M (ось {axis}), Н·м"):
    """Эпюры момента при разной степени сглаживания на одном графике по осям.

    ``curves`` — список ``(q, z, M_x, M_y)``: ``z`` от правого конца вала (как ``Z``),
    ``M_x``, ``M_y`` в Н·мм. ``axes`` — какие оси рисовать (``"y"``, ``"x"``);
    ``title``, ``xlabel``, ``ylabel`` (с ``{axis}``) — подписи, например на английском.
    """
    fig, ax_list = plt.subplots(len(axes), 1, figsize=(10, 3.6 * len(axes) + 0.6),
                                sharex=True, squeeze=False)
    ax_list = ax_list[:, 0]
    for ax, name in zip(ax_list, axes):
        idx = 3 if name == "y" else 2
        for item in curves:
            ax.plot(shaft_length - np.asarray(item[1]), np.asarray(item[idx]) * NMM_TO_NM,
                    label=f"q = {item[0]:g} %", lw=1.8)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_ylabel(ylabel.format(axis=name.upper()))
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
    ax_list[-1].set_xlabel(xlabel or "Координата вдоль вала, мм (0 — левый конец)")
    fig.suptitle(title or "Эпюры момента при разной степени сглаживания")
    fig.tight_layout()
    return _finish(fig, save_path, show)
