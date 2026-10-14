"""Cross-checks of the hand-rolled algorithms against scipy / statsmodels / known answers."""

import numpy as np
import pytest
from scipy import stats as sps

from statsheet import stats as ms


def test_cell_declustering_downweights_clusters():
    rng = np.random.default_rng(0)
    # 50 spread out, 50 clustered
    x = np.concatenate([rng.uniform(0, 100, 50), rng.uniform(40, 45, 50)])
    y = np.concatenate([rng.uniform(0, 100, 50), rng.uniform(40, 45, 50)])
    w = ms.cell_declustering(x, y, cell_size=10)
    assert w.sum() == pytest.approx(100)
    assert w[:50].mean() > 3 * w[50:].mean()


def test_composite_intervals_conserves_metal():
    frm = np.array([0, 0.5, 1.5, 3.5, 4.0])
    to = np.array([0.5, 1.5, 3.5, 4.0, 6.0])
    g = np.array([1.0, 2.0, 4.0, 8.0, 1.0])
    c = ms.composite_intervals(frm, to, g, length=2.0)
    assert np.allclose(c["to"] - c["from"], 2.0)
    assert np.sum(c["grade"] * 2.0) == pytest.approx(np.sum(g * (to - frm)))


def test_experimental_variogram_white_noise_equals_variance():
    rng = np.random.default_rng(1)
    x, y = rng.uniform(0, 100, 600), rng.uniform(0, 100, 600)
    z = rng.normal(0, 2, 600)
    h, g, n = ms.experimental_variogram(x, y, z, lag=5, n_lags=10)
    assert np.all(np.abs(g - 4.0) < 0.6)


def test_variogram_fit_recovers_model():
    h = np.linspace(2, 60, 30)
    g = ms.spherical(h, 1.0, 5.0, 40.0)
    p = ms.fit_variogram(h, g, ms.spherical)
    assert np.allclose(p, [1.0, 5.0, 40.0], rtol=0.02)


def test_directional_variogram_sees_anisotropy():
    rng = np.random.default_rng(2)
    x, y = rng.uniform(0, 100, 800), rng.uniform(0, 100, 800)
    z = np.sin(y / 8) + 0.1 * rng.normal(size=800)  # varies only along y: continuous E-W
    _, g_ns, _ = ms.experimental_variogram(x, y, z, lag=4, n_lags=5, azimuth=0)
    _, g_ew, _ = ms.experimental_variogram(x, y, z, lag=4, n_lags=5, azimuth=90)
    assert g_ns[-1] > 3 * g_ew[-1]


def test_kriging_exact_and_variance_zero_at_data():
    rng = np.random.default_rng(3)
    xy = rng.uniform(0, 100, (40, 2))
    z = rng.normal(5, 1, 40)
    est, var = ms.ordinary_kriging(xy, z, xy, ms.spherical, (0.0, 1.0, 30.0), n_max=12)
    assert np.allclose(est, z, atol=1e-6)
    assert np.all(var < 1e-6)


def test_kriging_matches_global_solution():
    """Moving-neighbourhood OK with all points == closed-form global OK."""
    rng = np.random.default_rng(4)
    xy = rng.uniform(0, 50, (15, 2))
    z = rng.normal(size=15)
    q = np.array([[25.0, 25.0]])
    params = (0.2, 1.0, 20.0)
    est, _ = ms.ordinary_kriging(xy, z, q, ms.spherical, params, n_max=15)
    D = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    C = params[1] - ms.spherical(D, *params)
    c0 = params[1] - ms.spherical(np.hypot(xy[:, 0] - 25, xy[:, 1] - 25), *params)
    A = np.ones((16, 16))
    A[:15, :15] = C
    A[15, 15] = 0
    w = np.linalg.solve(A, np.append(c0, 1))[:15]
    assert est[0] == pytest.approx(w @ z)


def test_nscore_roundtrip_and_normality():
    rng = np.random.default_rng(5)
    z = rng.lognormal(0, 1, 500)
    y, table = ms.nscore(z)
    assert abs(y.mean()) < 0.05 and abs(y.std() - 1) < 0.05
    assert np.allclose(ms.nscore_back(y, table), z)


def test_sgs_honours_data_and_variance():
    rng = np.random.default_rng(6)
    xy = rng.uniform(0, 40, (30, 2))
    y = rng.normal(size=30)
    grid = np.array([[i, j] for i in range(0, 41, 4) for j in range(0, 41, 4)], float)
    grid = np.vstack([grid, xy[:5]])  # include 5 data locations as nodes
    sims = ms.sequential_gaussian_simulation(xy, y, grid, ms.spherical, (0.0, 1.0, 15.0),
                                             n_real=3, seed=1)
    assert np.allclose(sims[:, -5:], y[:5], atol=1e-3)
    assert 0.6 < sims.std() < 1.4


def test_clr_ilr_inverses():
    rng = np.random.default_rng(7)
    X = ms.closure(rng.lognormal(0, 1, (20, 5)))
    assert np.allclose(ms.clr(X).sum(axis=1), 0)
    assert np.allclose(ms.clr_inv(ms.clr(X)), X)
    assert np.allclose(ms.ilr_inv(ms.ilr(X)), X)
    assert np.allclose(ms.alr_inv(ms.alr(X)), X)
    assert ms.ilr(X).shape == (20, 4)


def test_thompson_howarth_recovers_precision():
    rng = np.random.default_rng(8)
    c = rng.lognormal(0, 0.8, 3000)
    sd = 0.05 + 0.1 * c                       # per-measurement sd grows with concentration
    a = c + rng.normal(0, 1, 3000) * sd
    b = c + rng.normal(0, 1, 3000) * sd       # both members carry the error
    bm, bsd, s0, k = ms.thompson_howarth(a, b, n_bins=15)
    assert k == pytest.approx(0.1, abs=0.03)
    assert s0 == pytest.approx(0.05, abs=0.04)


def test_km_left_censored_no_censoring_equals_sample_mean():
    rng = np.random.default_rng(9)
    v = rng.lognormal(1, 0.5, 200)
    km = ms.km_left_censored(v, np.zeros(200, bool))
    assert km["mean"] == pytest.approx(v.mean(), rel=1e-9)
    assert km["median"] == pytest.approx(np.median(v), rel=0.02)


def test_ros_recovers_lognormal_mean():
    rng = np.random.default_rng(10)
    v = rng.lognormal(0, 1, 2000)
    dl = np.quantile(v, 0.4)
    cens = v < dl
    v_rep = np.where(cens, dl, v)
    ros = ms.ros_lognormal(v_rep, cens)
    assert ros["mean"] == pytest.approx(v.mean(), rel=0.08)


def test_kaplan_meier_matches_empirical_without_censoring():
    t = np.array([1, 2, 3, 4, 5.0])
    times, S, _ = ms.kaplan_meier(t, np.ones(5, bool))
    assert np.allclose(S, [1, .8, .6, .4, .2, 0])


def test_weibull_mle_recovers_parameters_with_censoring():
    rng = np.random.default_rng(11)
    t = 500 * rng.weibull(2.0, 3000)
    c = 700.0
    event = t < c
    t_obs = np.minimum(t, c)
    fit = ms.weibull_mle(t_obs, event)
    assert fit.shape == pytest.approx(2.0, abs=0.1)
    assert fit.scale == pytest.approx(500, rel=0.04)
    ref = sps.weibull_min.fit(t, floc=0)  # uncensored reference on the full data
    assert fit.shape == pytest.approx(ref[0], abs=0.15)


def test_shewhart_false_alarm_rate():
    rng = np.random.default_rng(12)
    x = rng.normal(10, 1, 20000)
    _, sigma, lo, hi = ms.shewhart_limits(x)
    assert sigma == pytest.approx(1.0, abs=0.03)
    assert 0.001 < np.mean((x < lo) | (x > hi)) < 0.005  # nominal 0.0027


def test_cusum_detects_one_sigma_shift_fast():
    rng = np.random.default_rng(13)
    x = np.concatenate([rng.normal(0, 1, 60), rng.normal(1, 1, 100)])
    _, _, sig = ms.cusum(x, 0, 1)
    after = sig[sig >= 60]
    assert len(after) and after[0] <= 90   # ARL1 for a 1-sigma shift is ~10 samples


def test_gy_mass_and_variance_are_consistent():
    m = ms.gy_min_mass(1.0, 0.05, 1000.0)
    assert ms.gy_fse_variance(1.0, m, 1000.0) == pytest.approx(0.05)
    m2 = ms.gy_min_mass(0.2, 0.1, 1e7, liberation_size_cm=0.005, alpha=1.5)
    v2 = ms.gy_fse_variance(0.2, m2, 1e7, liberation_size_cm=0.005, alpha=1.5)
    assert v2 == pytest.approx(0.1)
    # liberation factor is capped at 1 once the top size reaches the liberation size
    assert ms.gy_fse_variance(0.001, 10, 1e7, liberation_size_cm=0.005, alpha=1.5) == pytest.approx(
        ms.gy_fse_variance(0.001, 10, 1e7))


def test_fit_ou_recovers_mean_reversion():
    rng = np.random.default_rng(14)
    n, b, mean = 5000, 0.95, 2.0
    x = np.empty(n)
    x[0] = mean
    for i in range(1, n):
        x[i] = mean * (1 - b) + b * x[i - 1] + rng.normal(0, 0.1)
    f = ms.fit_ou(x)
    assert f["b"] == pytest.approx(b, abs=0.01)
    assert f["mean"] == pytest.approx(mean, abs=0.05)


def test_correlated_normals_hits_target_correlation():
    corr = np.array([[1, 0.6], [0.6, 1.0]])
    z = ms.correlated_normals(50000, corr, rng=0)
    assert np.corrcoef(z.T)[0, 1] == pytest.approx(0.6, abs=0.02)


# ------------------------------------------------ censored inference, trend and limits
def test_km_mean_places_mass_below_lowest_limit_at_that_limit():
    # one non-detect '< 1' and detects 2, 3: the only defensible mean is 2 (Helsel 2012)
    km = ms.km_left_censored(np.array([1.0, 2.0, 3.0]), np.array([True, False, False]))
    assert km["mean"] == pytest.approx(2.0)
    mean, se = ms.km_mean_se(np.array([1.0, 2.0, 3.0]), np.array([True, False, False]))
    assert mean == pytest.approx(2.0)


def test_km_mean_tracks_truth_with_several_detection_limits():
    rng = np.random.default_rng(21)
    x = rng.lognormal(1, 0.5, 400)
    dls = rng.choice([1.5, 2.5, 4.0], 400)
    c = x < dls
    xc = np.where(c, dls, x)
    km = ms.km_left_censored(xc, c)
    mean, se = ms.km_mean_se(xc, c)
    assert km["mean"] == pytest.approx(mean)
    assert abs(km["mean"] - x.mean()) < 3 * se
    assert abs(km["mean"] - x.mean()) < abs(xc.mean() - x.mean())  # beats DL substitution
    m0, se0 = ms.km_mean_se(x, np.zeros(400, bool))
    assert m0 == pytest.approx(x.mean())
    assert se0 == pytest.approx(x.std(ddof=1) / np.sqrt(400), rel=0.02)


def test_gehan_equals_mann_whitney_without_censoring():
    rng = np.random.default_rng(22)
    a, b = rng.normal(0, 1, 30), rng.normal(0.6, 1, 25)
    v = np.concatenate([a, b])
    g = np.r_[np.zeros(30), np.ones(25)]
    gh = ms.gehan_test(v, np.zeros(55, bool), g)
    mw = sps.mannwhitneyu(a, b, method="asymptotic", use_continuity=False)
    assert gh["p"] == pytest.approx(mw.pvalue, rel=1e-6)


def test_gehan_ignores_uninformative_censored_pairs():
    # group A: all '< 10'; group B: detects between 1 and 5. No pair is ordered: W must be 0.
    v = np.r_[np.full(10, 10.0), np.linspace(1, 5, 10)]
    c = np.r_[np.ones(10, bool), np.zeros(10, bool)]
    g = np.r_[np.zeros(10), np.ones(10)]
    gh = ms.gehan_test(v, c, g)
    assert gh["W"] == 0


def test_akritas_theil_sen_reduces_to_theil_sen_without_censoring():
    rng = np.random.default_rng(23)
    x = np.arange(40.0)
    y = 0.3 * x + rng.normal(0, 1, 40)
    ats = ms.akritas_theil_sen(x, y, np.zeros(40, bool), n_perm=200, rng=0)
    assert ats["slope"] == pytest.approx(sps.theilslopes(y, x).slope, abs=0.01)
    assert ats["tau"] == pytest.approx(sps.kendalltau(x, y).statistic, abs=1e-9)


def test_censored_regression_matches_ols_and_recovers_truth():
    import statsmodels.api as sm
    rng = np.random.default_rng(24)
    X = rng.normal(0, 1, (300, 2))
    beta = np.array([1.0, 0.8, -0.5])
    y = beta[0] + X @ beta[1:] + rng.normal(0, 0.7, 300)
    full = ms.censored_regression_mle(X, y, np.zeros(300, bool))
    ols = sm.OLS(y, sm.add_constant(X)).fit()
    assert np.allclose(full["beta"], ols.params, atol=1e-4)
    assert np.allclose(full["se"], ols.bse, rtol=0.03)
    dl = np.quantile(y, 0.4)
    c = y < dl
    cens = ms.censored_regression_mle(X, np.where(c, dl, y), c)
    sub = sm.OLS(np.where(c, dl / 2, y), sm.add_constant(X)).fit()
    assert np.all(np.abs(cens["beta"] - beta) < 2.5 * cens["se"])
    assert np.abs(cens["beta"] - beta).sum() < np.abs(sub.params - beta).sum()


def test_seasonal_kendall_single_season_is_mann_kendall_and_blocks_seasonality():
    rng = np.random.default_rng(25)
    years = np.arange(1980, 2020)
    v = 0.05 * (years - 1980) + rng.normal(0, 1, 40)
    sk = ms.seasonal_kendall(v, years, np.ones(40))
    kt = sps.kendalltau(years, v)
    assert sk["tau"] == pytest.approx(kt.statistic)
    assert sk["p"] == pytest.approx(kt.pvalue, rel=0.05)
    assert sk["sen"] == pytest.approx(sps.theilslopes(v, years).slope)
    # strong seasonality, no trend: the seasonal test must not reject at the median seed
    months = np.tile(np.arange(12), 40)
    yrs = np.repeat(years, 12)
    ps = []
    for seed in range(20):
        r = np.random.default_rng(seed)
        vals = 3 * np.sin(months / 12 * 2 * np.pi) + r.normal(0, 1, 480)
        ps.append(ms.seasonal_kendall(vals, yrs, months)["p"])
    assert 0.2 < np.median(ps) < 0.8


def test_dunn_two_groups_matches_mann_whitney():
    rng = np.random.default_rng(26)
    a, b = rng.normal(0, 1, 30), rng.normal(0.6, 1, 25)
    d = ms.dunn_test([a, b])
    mw = sps.mannwhitneyu(a, b, method="asymptotic", use_continuity=False)
    assert d.p.iloc[0] == pytest.approx(mw.pvalue, rel=1e-6)
    d3 = ms.dunn_test([a, b, rng.normal(0, 1, 20)], labels=list("ABC"), adjust="holm")
    assert len(d3) == 3 and (d3.p_adj >= d3.p).all()


def test_ucl_student_t_and_recommendation():
    rng = np.random.default_rng(27)
    v = rng.normal(10, 1, 60)
    tab = ms.ucl_table(v, rng=0)
    assert tab["Student-t UCL95"] == pytest.approx(
        v.mean() + sps.t.ppf(0.95, 59) * v.std(ddof=1) / np.sqrt(60))
    assert tab["Chebyshev UCL95"] > tab["Student-t UCL95"]
    assert ms.ucl_recommend(v)["choice"] == "Student-t UCL95"
    skewed = rng.lognormal(0, 1.2, 60)
    assert ms.ucl_recommend(skewed)["choice"] != "Student-t UCL95"
    c = skewed < np.quantile(skewed, 0.3)
    censored = np.where(c, np.quantile(skewed, 0.3), skewed)
    assert ms.ucl_recommend(censored, c)["choice"].startswith("KM")


def test_prediction_and_tolerance_limits_cover_at_nominal_rate():
    rng = np.random.default_rng(28)
    hits_pl, hits_tl = 0, 0
    for _ in range(400):
        bg = rng.normal(0, 1, 20)
        pl = ms.prediction_limit(bg, k_future=1)
        hits_pl += rng.normal() > pl["normal"]
        tl = ms.tolerance_limit(bg, coverage=0.95, confidence=0.95)
        hits_tl += sps.norm.sf(tl) > 0.05   # the limit fails to cover 95 % of the population
    assert 0.02 < hits_pl / 400 < 0.09      # nominal 5 %
    assert hits_tl / 400 < 0.10             # nominal 5 % confidence failure
