"""Cross-check the hand-rolled algorithms against scipy / known results."""

import numpy as np
import pytest
from scipy import stats as sps

from earthstats import stats as es


@pytest.fixture
def trending():
    rng = np.random.default_rng(0)
    t = np.arange(60)
    return t, 0.3 * t + rng.normal(0, 2, 60)


def test_mann_kendall_matches_kendalltau(trending):
    t, x = trending
    mk = es.mann_kendall(x)
    tau, p = sps.kendalltau(t, x)
    assert mk.tau == pytest.approx(tau, abs=1e-9)
    assert mk.p == pytest.approx(p, rel=0.05)  # scipy uses exact/other variance for ties
    assert mk.trend == "increasing"


def test_mann_kendall_null():
    rng = np.random.default_rng(1)
    ps = [es.mann_kendall(rng.normal(size=40)).p for _ in range(300)]
    # under H0, ~5% of p-values should be below 0.05
    assert 0.01 < np.mean(np.array(ps) < 0.05) < 0.10


def test_sens_slope_matches_theilslopes(trending):
    t, x = trending
    slope, intercept, lo, hi = es.sens_slope(x, t)
    ref = sps.theilslopes(x, t)
    assert slope == pytest.approx(ref.slope)
    assert intercept == pytest.approx(ref.intercept)
    assert lo == pytest.approx(ref.low_slope, rel=0.02)
    assert hi == pytest.approx(ref.high_slope, rel=0.02)


def test_pettitt_finds_shift():
    rng = np.random.default_rng(2)
    x = np.concatenate([rng.normal(0, 1, 50), rng.normal(2, 1, 50)])
    k, _, p = es.pettitt(x)
    assert 44 <= k <= 54
    assert p < 0.001


def test_ar1_effective_n():
    rng = np.random.default_rng(3)
    x = np.zeros(2000)
    for i in range(1, 2000):
        x[i] = 0.6 * x[i - 1] + rng.normal()
    n_eff, r1 = es.ar1_effective_n(x)
    assert r1 == pytest.approx(0.6, abs=0.05)
    assert n_eff == pytest.approx(2000 * 0.4 / 1.6, rel=0.15)


def test_bootstrap_ci_covers_mean():
    rng = np.random.default_rng(4)
    x = rng.normal(10, 2, 200)
    _, lo, hi, _ = es.bootstrap_ci(x, np.mean, n_boot=2000, rng=0)
    assert lo < 10 < hi
    assert hi - lo == pytest.approx(2 * 1.96 * 2 / np.sqrt(200), rel=0.2)


def test_permutation_test_detects_difference():
    rng = np.random.default_rng(5)
    a, b = rng.normal(0, 1, 40), rng.normal(1, 1, 40)
    obs, p, _ = es.permutation_test(a, b, n_perm=2000, rng=0)
    assert obs < 0 and p < 0.01


def test_morans_i_clustered_vs_random():
    # 10x10 grid; clustered field = x coordinate -> strong positive I
    gx, gy = np.meshgrid(np.arange(10), np.arange(10))
    lat, lon = gy.ravel() * 0.1 + 40, gx.ravel() * 0.1 - 100
    D = es.distance_matrix_km(lat, lon)
    W = es.knn_weights(D, k=4)
    clustered = es.morans_i(gx.ravel().astype(float), W, n_perm=200, rng=0)
    random = es.morans_i(np.random.default_rng(6).normal(size=100), W, n_perm=200, rng=0)
    assert clustered["I"] > 0.8 and clustered["p_perm"] < 0.01
    assert abs(random["I"]) < 0.25


def test_kriging_is_exact_at_data_points():
    rng = np.random.default_rng(7)
    lat, lon = rng.uniform(35, 45, 30), rng.uniform(-100, -90, 30)
    z = lat * 0.5 + rng.normal(0, 0.2, 30)
    D = es.distance_matrix_km(lat, lon)
    params = (0.0, z.var(), 500.0)
    pred, var = es.ordinary_kriging(D, D, z, es.spherical, params)
    assert np.allclose(pred, z, atol=1e-6)
    assert np.all(var < 1e-6 + 1e-9)


def test_idw_exact_and_bounded():
    z = np.array([1.0, 3.0, 5.0])
    D = np.array([[0.0, 1.0, 2.0], [1.0, 1.0, 1.0]])
    out = es.idw(D, z)
    assert out[0] == pytest.approx(1.0)
    assert out[1] == pytest.approx(3.0)


def test_return_level_gumbel_limit():
    # GEV with shape 0 (Gumbel): z_T = loc - scale * ln(-ln(1-1/T))
    T = np.array([2, 10, 100])
    z = es.return_level_gev(0.0, 100.0, 10.0, T)
    expected = 100 - 10 * np.log(-np.log(1 - 1 / T))
    assert np.allclose(z, expected)


def _ar1(rng, n, rho):
    x = np.empty(n)
    x[0] = rng.normal()
    for i in range(1, n):
        x[i] = rho * x[i - 1] + rng.normal()
    return x


def test_hamed_rao_controls_size_under_ar1():
    rng = np.random.default_rng(8)
    raw = hr = 0
    for _ in range(300):
        x = _ar1(rng, 60, 0.6)
        raw += es.mann_kendall(x).p < 0.05
        hr += es.mann_kendall_hamed_rao(x).p < 0.05
    assert raw / 300 > 0.15           # raw MK badly over-rejects
    assert hr / 300 < 0.12            # corrected test is near nominal


def test_hamed_rao_equals_mk_for_white_noise():
    x = np.random.default_rng(9).normal(size=50)
    assert es.mann_kendall_hamed_rao(x).p == pytest.approx(es.mann_kendall(x).p, abs=0.05)


def test_fit_gev_recovers_parameters():
    rng = np.random.default_rng(10)
    x = sps.genextreme.rvs(-0.2, loc=100, scale=20, size=2000, random_state=rng)
    c, loc, scale = es.fit_gev(x)
    assert c == pytest.approx(-0.2, abs=0.06)
    assert loc == pytest.approx(100, abs=2)
    assert scale == pytest.approx(20, abs=2)
