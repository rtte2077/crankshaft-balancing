import numpy as np
import pytest

from crankshaft_balancing import (analyze, make_initial_coefficients, total_distribution,
                                  optimize_counterweights,
                                  optimize_counterweights_history,
                                  sweep_sensitivity, verify_kkt_conditions,
                                  verify_optimum_via_duality)
from crankshaft_balancing.sections import NMM_TO_NM


@pytest.fixture(scope="module")
def result(shaft):
    engine, components = shaft
    return analyze(engine, components, step=0.05, verbose=False)


def test_force_and_end_moment_are_zero(result):
    f, w = result.forces, result.engine.angular_speed
    for r, Fc, Fp in ((result.result_y, f.Fy_cheeks, f.Fy_crankpins),
                      (result.result_x, f.Fx_cheeks, f.Fx_crankpins)):
        assert r["success"]
        assert Fc.sum() + Fp.sum() + r["s"].sum() * w ** 2 == pytest.approx(0, abs=1e-6)
        assert r["M"][-1] * NMM_TO_NM == pytest.approx(0, abs=1e-6)


def test_regression_max_moment(result):
    # Значение из исходного ноутбука (сетка 0.01 мм) — на сетке 0.05 мм почти совпадает.
    assert result.max_moment_y == pytest.approx(135.5, abs=0.5)
    assert result.max_moment_x == pytest.approx(0.0, abs=1e-6)


def test_mirrored_coefficients_are_also_optimal(result):
    # Вал симметричен, поэтому зеркальный набор s тоже оптимален. Сами коэффициенты
    # могут отличаться между версиями SciPy/HiGHS: оптимум неединственен (при том же
    # max|M|), поэтому сравниваем значение цели, а не векторы s.
    f, w = result.forces, result.engine.angular_speed
    r = result.result_y
    s_mirror = r["s"][::-1]
    M = total_distribution(result.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                           f.Fy_crankpins, F3=s_mirror * w ** 2,
                           L3=result.L_counterweights)
    assert np.max(np.abs(M)) == pytest.approx(np.max(np.abs(r["M"])), rel=1e-6)
    assert f.Fy_cheeks.sum() + f.Fy_crankpins.sum() + s_mirror.sum() * w ** 2 \
        == pytest.approx(0, abs=1e-6)
    assert M[-1] * NMM_TO_NM == pytest.approx(0, abs=1e-6)


def test_global_optimum_certificates(result):
    for r in (result.result_y, result.result_x):
        assert verify_optimum_via_duality(r, verbose=False)["confirmed"]
        assert verify_kkt_conditions(r, verbose=False)["confirmed"]


def test_perturbation_never_beats_optimum(result):
    f, w = result.forces, result.engine.angular_speed
    curves = sweep_sensitivity(result.result_y, result.Z, f.L_cheeks, f.Fy_cheeks,
                               f.L_crankpins, f.Fy_crankpins, w)
    best = np.max(np.abs(result.result_y["M"]))
    # при изменении любого одного s_i (без пересчёта баланса) max|M| не меньше оптимума
    for _, _, values in curves:
        assert values.min() >= best * (1 - 1e-9)


def test_fixed_coefficient_is_respected(result):
    f, w = result.forces, result.engine.angular_speed
    s_init = make_initial_coefficients(len(f.Fy_cheeks), overrides={2: "F"})
    r = optimize_counterweights(result.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                                f.Fy_crankpins, w, s_init, verbose=False)
    assert r["s"][2] == 0.0
    assert np.max(np.abs(r["M"])) >= np.max(np.abs(result.result_y["M"])) - 1e-6


def test_all_fixed_skips_solver(result):
    f, w = result.forces, result.engine.angular_speed
    r = optimize_counterweights(result.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                                f.Fy_crankpins, w, ["F"] * len(f.Fy_cheeks), verbose=False)
    assert r["result"] is None and r["success"]


def test_bad_override_index():
    with pytest.raises(IndexError):
        make_initial_coefficients(3, overrides={5: "F"})


def test_history_converges(result):
    f, w = result.forces, result.engine.angular_speed
    from crankshaft_balancing import section_grid
    Zc = section_grid(result.engine.crankshaft_length, step=2.0)
    try:
        h = optimize_counterweights_history(Zc, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                                            f.Fy_crankpins, w,
                                            make_initial_coefficients(len(f.Fy_cheeks)),
                                            verbose=False)
    except RuntimeError:
        pytest.skip("В этой версии SciPy нет linprog(method='simplex')")
    vals = h["history"]["max_abs_M"]
    assert len(vals) > 1 and vals[-1] <= vals[0]


def _plain_min_max(z, L1, F1, L2, F2, w, n, L3, s_bounds):
    """Произвольное решение исходной задачи min max|M| (без выбора по размаху)."""
    from scipy.optimize import linprog
    from crankshaft_balancing import build_counterweight_lp
    lp = build_counterweight_lp(z, L1, F1, L2, F2, w, [1.0] * n, L3=L3, s_bounds=s_bounds)
    res = linprog(lp.c, A_ub=lp.A_ub, b_ub=lp.b_ub, A_eq=lp.A_eq, b_eq=lp.b_eq,
                  bounds=lp.bounds, method="highs")
    s = lp.assemble_solution(res.x)
    return total_distribution(z, L1, F1, L2, F2, F3=s * w ** 2, L3=L3)


def test_smallest_range_among_min_max_solutions():
    # синтетический вал, где оптимум по max|M| неединственный и решения различаются
    # размахом эпюры max M - min M
    from crankshaft_balancing import section_grid
    z, w = section_grid(300.0, 1.0), 100.0
    L1, F1 = [15.6, 32.9, 53.8, 118.0], [679.0, -470.0, -870.0, 77.0]
    L2, F2 = [16.0, 98.7], [-863.0, 620.0]
    L3 = np.array([36.0, 47.0, 104.6, 220.9, 279.7])
    bounds = [(-0.2, 0.2)] * 5
    r = optimize_counterweights(z, L1, F1, L2, F2, w, [1.0] * 5, L3=L3, s_bounds=bounds,
                                verbose=False)
    M_plain = _plain_min_max(z, L1, F1, L2, F2, w, 5, L3, bounds)
    assert np.max(np.abs(r["M"])) == pytest.approx(np.max(np.abs(M_plain)), rel=1e-7)
    assert r["moment_range"] == pytest.approx(np.ptp(r["M"]))
    assert r["moment_range"] < 0.7 * np.ptp(M_plain)      # 9622 -> 5300
    assert verify_optimum_via_duality(r, verbose=False)["confirmed"]
    assert verify_kkt_conditions(r, verbose=False)["confirmed"]


def test_range_never_worse_than_plain_min_max(result):
    f, w = result.forces, result.engine.angular_speed
    M_plain = _plain_min_max(result.Z, f.L_cheeks, f.Fy_cheeks, f.L_crankpins,
                             f.Fy_crankpins, w, 8, result.L_counterweights, None)
    assert result.moment_range_y * 1e3 <= np.ptp(M_plain) * (1 + 1e-8)
    assert result.max_moment_y * 1e3 == pytest.approx(np.max(np.abs(M_plain)), rel=1e-7)


def test_reduced_rows_match_full_grid(result):
    # строки ЛП берутся только в «критических» сечениях; на всей сетке ограничение
    # |M| <= max|M| из ЛП тоже должно выполняться
    r = result.result_y
    t = r["lp"].t_value(r["result"].x)
    assert np.max(np.abs(r["M"])) <= t * (1 + 1e-7)
    assert r["lp"].A_ub.shape[0] < len(result.Z) / 10
