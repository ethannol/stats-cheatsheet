# %% [markdown]
# # 07 Resampling: bootstrap, block bootstrap, jackknife, permutation, Monte Carlo
#
# **Cheat sheet**
# - Bootstrap: resample with replacement, recompute the statistic, read off the spread.
#   Works for any statistic with no formula (median, trimmed mean, ratio, quantile, slope).
# - Intervals: percentile (simple, can under-cover for skewed/biased statistics), basic
#   (reflected), BCa (bias-corrected and accelerated; default choice), studentised.
# - It assumes iid. For time series use block bootstrap (moving / stationary blocks) or
#   resample model residuals (residual bootstrap) or AR-sieve. iid bootstrap on autocorrelated
#   data under-estimates uncertainty, often badly.
# - Regression: case (pairs) bootstrap is robust to heteroscedasticity; residual bootstrap
#   assumes the model is right and errors are exchangeable.
# - Jackknife: leave-one-out; gives bias and SE estimates; fails for non-smooth statistics
#   like the median.
# - Permutation: shuffle labels under H0 of exchangeability to get an exact null.
# - Monte Carlo: simulate from a fitted null model (e.g. AR(1)) to calibrate any test.
# - Set the seed. Report the number of resamples. 2,000 is fine for SEs, 10,000 for tails.

# %%
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import ar1_effective_n, block_bootstrap, bootstrap_ci

setup()
SEC = "07_resampling"
rng = np.random.default_rng(42)

cp = data.ghcn_daily()
wet = cp.PRCP[(cp.PRCP > 0) & (cp.index.year == 2015)].values  # one year of wet days, n ~ 120

# %% ----------------------------------------------------------- bootstrap CIs for a skewed statistic
banner(f"Bootstrap CI for the mean and median of wet-day precipitation (2015, n={len(wet)})")
for stat_name, fn in [
    ("mean", np.mean),
    ("median", np.median),
    ("90th pct", lambda x: np.quantile(x, 0.9)),
]:
    est, lo, hi, boots = bootstrap_ci(wet, fn, n_boot=5000, method="percentile", rng=0)
    _, lo_b, hi_b, _ = bootstrap_ci(wet, fn, n_boot=5000, method="basic", rng=0)
    bca = sps.bootstrap(
        (wet,), fn, n_resamples=5000, method="BCa", random_state=0
    ).confidence_interval
    print(
        f"  {stat_name:9s} est={est:6.2f}  percentile=({lo:5.2f},{hi:5.2f})  "
        f"basic=({lo_b:5.2f},{hi_b:5.2f})  BCa=({bca.low:5.2f},{bca.high:5.2f})"
    )
t_ci = sps.t.interval(0.95, len(wet) - 1, loc=wet.mean(), scale=sps.sem(wet))
print(f"  mean: classical t-interval = ({t_ci[0]:.2f}, {t_ci[1]:.2f})")

est, lo, hi, boots = bootstrap_ci(wet, np.mean, n_boot=5000, rng=0)
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(wet, bins=30, alpha=0.8)
axes[0].set_title("Wet-day precipitation 2015 (skewed)")
axes[0].set_xlabel("mm")
axes[1].hist(boots, bins=60, alpha=0.8, density=True)
for v, lab in [(lo, "2.5%"), (hi, "97.5%")]:
    axes[1].axvline(v, color="C3", ls="--", label=lab)
axes[1].axvline(est, color="k", label="estimate")
axes[1].set_title("Bootstrap distribution of the mean (CLT at work)")
axes[1].legend()
savefig(fig, "bootstrap_mean", SEC)

# %% ----------------------------------------------------------- coverage simulation
banner("Coverage check: which interval actually has 95% coverage? (lognormal, n=30)")
true_mu, true_sig = 1.4, 1.2
true_mean = np.exp(true_mu + true_sig**2 / 2)
true_median = np.exp(true_mu)
n_rep = 400
cover = {k: 0 for k in ["mean:t", "mean:percentile", "mean:BCa", "median:percentile", "median:BCa"]}
for r in range(n_rep):
    x = rng.lognormal(true_mu, true_sig, 30)
    lo, hi = sps.t.interval(0.95, 29, loc=x.mean(), scale=sps.sem(x))
    cover["mean:t"] += lo <= true_mean <= hi
    _, lo, hi, _ = bootstrap_ci(x, np.mean, n_boot=1000, rng=r)
    cover["mean:percentile"] += lo <= true_mean <= hi
    b = sps.bootstrap(
        (x,), np.mean, n_resamples=1000, method="BCa", random_state=r
    ).confidence_interval
    cover["mean:BCa"] += b.low <= true_mean <= b.high
    _, lo, hi, _ = bootstrap_ci(x, np.median, n_boot=1000, rng=r)
    cover["median:percentile"] += lo <= true_median <= hi
    b = sps.bootstrap(
        (x,), np.median, n_resamples=1000, method="BCa", random_state=r
    ).confidence_interval
    cover["median:BCa"] += b.low <= true_median <= b.high
for k, v in cover.items():
    print(f"  {k:18s} coverage={v / n_rep:.3f}")
note("""For the mean of a heavily skewed variable at n=30, every interval under-covers; BCa is
the least bad. The t-interval is not terrible either (the CLT is forgiving). The
bootstrap is not magic at small n with heavy tails: it only sees the tail the sample
happened to catch.""")

# %% ----------------------------------------------------------- block bootstrap for a trend
banner("Trend SE: iid bootstrap vs block bootstrap vs HAC (GISTEMP monthly 1975-2025)")
gt = data.gistemp().anomaly["1975":"2025"]
tt = np.arange(len(gt)) / 12.0
X = sm.add_constant(tt)
ols = sm.OLS(gt.values, X).fit()
hac = sm.OLS(gt.values, X).fit(cov_type="HAC", cov_kwds={"maxlags": 24})
resid = ols.resid
n_eff, r1 = ar1_effective_n(resid)


def slope_of(idx_or_vals):
    """slope of resid-resampled series: fitted + resampled residuals."""
    return sm.OLS(ols.fittedvalues + idx_or_vals, X).fit().params[1]


iid = np.array([slope_of(rng.choice(resid, len(resid))) for _ in range(1000)])
blk = block_bootstrap(resid, slope_of, block_len=24, n_boot=1000, rng=1)
# case (pairs) bootstrap for comparison, iid
case = np.array(
    [
        sm.OLS(gt.values[i], X[i]).fit().params[1]
        for i in (rng.integers(0, len(gt), len(gt)) for _ in range(1000))
    ]
)
print(f"slope = {ols.params[1]:.4f} C/yr   residual r1={r1:.2f}   n_eff={n_eff:.0f} of {len(gt)}")
print(f"  SE OLS (naive)             {ols.bse[1]:.5f}")
print(f"  SE iid residual bootstrap  {iid.std():.5f}")
print(f"  SE iid case bootstrap      {case.std():.5f}")
print(f"  SE block bootstrap (24 mo) {blk.std():.5f}")
print(f"  SE HAC Newey-West (24 lag) {hac.bse[1]:.5f}")
print(f"  SE AR(1) n_eff adjusted    {ols.bse[1] * np.sqrt(len(gt) / n_eff):.5f}")
note("""The iid bootstraps simply reproduce the naive OLS SE because resampling destroys the
autocorrelation that makes the naive SE wrong. Block bootstrap keeps 2-year chunks intact
and lands in the same place as HAC and the AR(1) adjustment. If the resampling scheme does
not preserve the dependence structure, the bootstrap inherits the same mistake as the
formula.""")

fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(iid * 10, bins=40, alpha=0.6, density=True, label=f"iid bootstrap SD={iid.std() * 10:.3f}")
ax.hist(
    blk * 10, bins=40, alpha=0.6, density=True, label=f"block bootstrap SD={blk.std() * 10:.3f}"
)
ax.set_xlabel("trend (C / decade)")
ax.set_title("Bootstrap distributions of the 1975-2025 warming trend")
ax.legend()
savefig(fig, "block_bootstrap", SEC)

# %% ----------------------------------------------------------- jackknife
banner("Jackknife: bias and SE of the skewness of annual precipitation totals")
ann_p = cp.PRCP.resample("YE").sum()
ann_p = ann_p[ann_p.index.year <= 2025].values
n = len(ann_p)
theta = sps.skew(ann_p)
loo = np.array([sps.skew(np.delete(ann_p, i)) for i in range(n)])
jack_bias = (n - 1) * (loo.mean() - theta)
jack_se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2))
_, blo, bhi, bb = bootstrap_ci(ann_p, sps.skew, n_boot=5000, rng=0)
print(
    f"skewness={theta:.3f}  jackknife bias={jack_bias:+.3f}  jackknife SE={jack_se:.3f}  "
    f"bootstrap SE={bb.std():.3f}"
)
print(f"bias-corrected skewness = {theta - jack_bias:.3f}")
note("""The jackknife is a deterministic, n-evaluation approximation to the bootstrap. Sample
skewness is biased toward zero at small n; the jackknife quantifies and corrects it. Do
not jackknife a median or other non-smooth statistic; its SE estimate is inconsistent.""")

# %% ----------------------------------------------------------- Monte Carlo null for a trend test
banner("Monte Carlo: calibrate a trend t-test against an AR(1) null (GISTEMP annual 2000-2025)")
gta = data.gistemp().anomaly.resample("YE").mean()
gta = gta[(gta.index.year >= 2000) & (gta.index.year <= 2025)]
y = gta.values
t_idx = np.arange(len(y), dtype=float)
Xa = sm.add_constant(t_idx)
f = sm.OLS(y, Xa).fit()
t_obs = f.tvalues[1]
res = f.resid
phi = np.corrcoef(res[:-1], res[1:])[0, 1]
sig = res.std(ddof=2) * np.sqrt(1 - phi**2)
n_sim = 5000
t_null = np.empty(n_sim)
for s in range(n_sim):
    e = rng.normal(0, sig, len(y))
    x = np.empty(len(y))
    x[0] = e[0] / np.sqrt(1 - phi**2)
    for k in range(1, len(y)):
        x[k] = phi * x[k - 1] + e[k]
    t_null[s] = sm.OLS(x, Xa).fit().tvalues[1]
p_mc = np.mean(np.abs(t_null) >= abs(t_obs))
print(
    f"slope={f.params[1] * 10:+.3f} C/decade  t={t_obs:.2f}  naive p={f.pvalues[1]:.4f}  "
    f"residual phi={phi:.2f}"
)
print(f"Monte Carlo AR(1) null: p={p_mc:.4f}  (null SD of t = {t_null.std():.2f} vs 1.0 for iid)")
note("""Simulating from a fitted null model is the most general calibration tool: fit AR(1)
to the residuals, generate thousands of trend-free series with the same persistence, and
ask how often the t-statistic is as large as observed. The null SD of t exceeds 1 by
roughly sqrt(n/n_eff), the same factor as the analytic adjustment in section 05.""")
