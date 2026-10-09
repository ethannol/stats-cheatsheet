# %% [markdown]
# # 05 Trend detection and change points
#
# Data: GISTEMP annual global anomaly; Potomac annual mean flow; Central Park annual TMIN.
#
# **Cheat sheet**
# - Parametric trend: OLS slope. SE must account for autocorrelation (AR(1) effective n,
#   HAC). Report slope per decade with CI, not just "significant".
# - Non-parametric trend: Mann-Kendall test (monotonic, rank-based, robust to outliers and
#   non-normality) + Sen's slope (median of pairwise slopes). The hydrology/climatology default.
# - Autocorrelation inflates MK false positives; use trend-free pre-whitening (TFPW), block
#   bootstrap, or the Hamed-Rao variance correction. Or test annual values only.
# - Change points: Pettitt (single shift in median, rank-based); CUSUM; piecewise/hinge
#   regression with the breakpoint estimated by profile likelihood; PELT/binary segmentation
#   for multiple change points (ruptures library). Distinguish shift vs gradual trend.
# - Start/end-point sensitivity: trends over short windows are dominated by noise and
#   endpoint choice (the "hiatus" lesson). Show the trend as a function of window.
# - Simulate the null: fit AR(1) to the data, generate many surrogate series, and see how
#   often the test fires. Cheap, honest, interview-friendly.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from earthstats import data
from earthstats.plotting import banner, note, savefig, setup
from earthstats.stats import (
    ar1_effective_n,
    mann_kendall,
    mann_kendall_hamed_rao,
    pettitt,
    prewhiten_simple,
    prewhiten_tfpw,
    sens_slope,
)

setup()
SEC = "05_trends"

gt = data.gistemp().anomaly.resample("YE").mean()
gt = gt[gt.index.year <= 2025]
gt.index = gt.index.year
q = data.usgs_daily().flow_cfs
wy = q.index.year + (q.index.month >= 10)  # water year Oct-Sep
q_ann = q.groupby(wy).mean().loc[1931:2025]
cp = data.ghcn_daily()
tmin_ann = cp.TMIN.resample("YE").mean()
tmin_ann = tmin_ann[tmin_ann.index.year <= 2025]
tmin_ann.index = tmin_ann.index.year


# %% ----------------------------------------------------------- OLS trend with proper SE
def ols_trend(y: pd.Series, label: str):
    X = sm.add_constant(y.index.values.astype(float))
    fit = sm.OLS(y.values, X).fit()
    hac = sm.OLS(y.values, X).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
    n_eff, r1 = ar1_effective_n(fit.resid)
    se_adj = fit.bse[1] * np.sqrt(len(y) / n_eff)
    print(f"{label:32s} slope={fit.params[1] * 10:+.3f}/decade  "
          f"SE naive={fit.bse[1] * 10:.3f}  SE AR1-adj={se_adj * 10:.3f}  SE HAC={hac.bse[1] * 10:.3f}  "
          f"r1(resid)={r1:.2f}  n_eff={n_eff:.0f}/{len(y)}")
    return fit


banner("OLS trend per decade, three ways to get the standard error")
ols_trend(gt, "GISTEMP 1880-2025 (C)")
ols_trend(gt.loc[1975:], "GISTEMP 1975-2025 (C)")
ols_trend(tmin_ann, "Central Park TMIN 1950-2025 (C)")
ols_trend(q_ann / 1000, "Potomac mean flow (kcfs)")
note("""The full-record GISTEMP trend is not linear (flat to 1910, up, flat, up), so the
residuals are strongly autocorrelated and the naive SE is badly optimistic. Over 1975-2025
the residual autocorrelation is modest. Always print r1 of the residuals next to any trend.""")

# %% ----------------------------------------------------------- Mann-Kendall + Sen
banner("Mann-Kendall + Sen's slope (non-parametric)")
for label, y in [("GISTEMP 1975-2025", gt.loc[1975:]), ("Central Park TMIN", tmin_ann),
                 ("Potomac mean flow (kcfs)", q_ann / 1000)]:
    mk = mann_kendall(y.values)
    slope, _, lo, hi = sens_slope(y.values, y.index.values.astype(float))
    pw, r1 = prewhiten_tfpw(y.values)
    hr = mann_kendall_hamed_rao(y.values)
    print(f"{label:26s} MK z={mk.z:+.2f} p={mk.p:.4f} | Hamed-Rao z={hr.z:+.2f} p={hr.p:.4f} "
          f"(r1={r1:.2f}) | Sen {slope * 10:+.3f}/decade (95% {lo * 10:+.3f}, {hi * 10:+.3f})")
note("""Potomac annual flow shows no monotonic trend: hydrologic series are noisy, and this is
the normal outcome. For GISTEMP and TMIN the trend survives the autocorrelation
correction because the signal is large relative to the persistence. Sen's slope is nearly
identical to OLS here (no outliers); it would differ on a series with a few extreme years.""")

# %% ----------------------------------------------------------- Monte Carlo on the null
banner("Autocorrelation and Mann-Kendall: size (no trend) and power (trend) of four variants")
rng = np.random.default_rng(0)
n, n_sim = 50, 1000


def simulate(rho, slope):
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0] / np.sqrt(1 - rho ** 2)
    for t in range(1, n):
        x[t] = rho * x[t - 1] + e[t]
    return x + slope * np.arange(n)


variants = {"raw MK": lambda x: mann_kendall(x).p,
            "simple PW": lambda x: mann_kendall(prewhiten_simple(x)[0]).p,
            "TFPW (Yue)": lambda x: mann_kendall(prewhiten_tfpw(x)[0]).p,
            "Hamed-Rao": lambda x: mann_kendall_hamed_rao(x).p}
rows = []
for rho in [0.0, 0.3, 0.6]:
    for slope, kind in [(0.0, "size"), (0.03, "power")]:
        rej = {k: 0 for k in variants}
        for _ in range(n_sim):
            x = simulate(rho, slope)
            for k, fn in variants.items():
                rej[k] += fn(x) < .05
        rows.append(dict(rho=rho, test=kind, **{k: v / n_sim for k, v in rej.items()}))
sim = pd.DataFrame(rows).set_index(["rho", "test"])
print(sim.round(3).to_string())
note("""Size rows should read 0.05. Raw MK rejects a true null ~30% of the time at rho=0.6.
Simple pre-whitening fixes the size but, as the power rows show, throws away part of a
real trend. Yue's trend-free pre-whitening (TFPW) was designed to keep the power, but
under the null it adds a spurious Sen slope back onto whitened noise and over-rejects even
more than raw MK (Bayazit & Onoz 2007; Hamed 2009). The Hamed-Rao variance correction
holds the size near nominal with a smaller power cost. Knowing which 'standard fix' is
actually broken is the kind of thing that separates practitioners from tutorial-readers.""")

# %% ----------------------------------------------------------- change points
banner("Change points: Pettitt, CUSUM, hinge regression")
for label, y in [("GISTEMP annual", gt), ("Potomac mean flow", q_ann), ("Central Park TMIN", tmin_ann)]:
    k, K, p = pettitt(y.values)
    print(f"{label:20s} Pettitt change after {y.index[k]}  K={K:.0f}  p={p:.2e}")


def hinge_fit(y: pd.Series):
    """Continuous piecewise-linear fit with one breakpoint chosen by minimum SSE."""
    t = y.index.values.astype(float)
    best = None
    for bp in t[5:-5]:
        X = np.column_stack([np.ones_like(t), t, np.maximum(t - bp, 0)])
        fit = sm.OLS(y.values, X).fit()
        if best is None or fit.ssr < best[1].ssr:
            best = (bp, fit)
    return best


bp, hf = hinge_fit(gt)
lin = sm.OLS(gt.values, sm.add_constant(gt.index.values.astype(float))).fit()
print(f"\nGISTEMP hinge: breakpoint={bp:.0f}, slope before={hf.params[1] * 10:+.3f}/decade, "
      f"after={(hf.params[1] + hf.params[2]) * 10:+.3f}/decade")
print(f"AIC linear={lin.aic:.1f}  AIC hinge={hf.aic:.1f}  (hinge uses 2 extra params incl. breakpoint)")
note("""Pettitt finds the year that best splits the series into two medians; on a series with a
gradual trend it will always 'find' a change near the middle, so it cannot by itself
distinguish shift from trend. The hinge (broken-stick) model estimates where the slope
changed; for global temperature the acceleration lands in the 1970s. The breakpoint is
itself estimated, so count it as a parameter when comparing AIC.""")

# CUSUM
cs = (gt - gt.loc[1951:1980].mean()).cumsum()
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(gt.index, gt.values, lw=1, label="annual anomaly")
tt = gt.index.values.astype(float)
axes[0].plot(tt, hf.fittedvalues, color="C3", lw=2, label=f"hinge fit, break {bp:.0f}")
k, *_ = pettitt(gt.values)
axes[0].axvline(gt.index[k], color="C1", ls="--", label=f"Pettitt {gt.index[k]}")
axes[0].legend()
axes[0].set_title("GISTEMP: hinge regression and Pettitt change point")
axes[1].plot(cs.index, cs.values)
axes[1].axhline(0, color="k", lw=.8)
axes[1].set_title("CUSUM of anomalies vs 1951-80 mean (slope change = regime change)")
savefig(fig, "changepoints", SEC)

# %% ----------------------------------------------------------- endpoint sensitivity
banner("Endpoint sensitivity: trend to 2025 as a function of start year")
rows = []
for start in range(1950, 2016):
    y = gt.loc[start:]
    X = sm.add_constant(y.index.values.astype(float))
    f = sm.OLS(y.values, X).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
    rows.append(dict(start=start, slope=f.params[1] * 10, lo=f.conf_int()[1, 0] * 10,
                     hi=f.conf_int()[1, 1] * 10, n=len(y)))
es = pd.DataFrame(rows).set_index("start")
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(es.index, es.slope, color="C0")
ax.fill_between(es.index, es.lo, es.hi, alpha=.25)
ax.axhline(0, color="k", lw=.8)
ax.set_xlabel("start year (series ends 2025)")
ax.set_ylabel("trend (C / decade)")
ax.set_title("Trend estimate and HAC 95% CI widen explosively for short windows")
savefig(fig, "endpoint_sensitivity", SEC)
print(es.loc[[1950, 1975, 1998, 2005, 2010, 2015]].round(3).to_string())
note(f"""Starting in the 1998 El Nino year was the basis of the 'warming hiatus' claim (which
has since filled in: the 1998-2025 trend is now the steepest in the table). For the
11-year window the HAC interval runs {es.loc[2015, 'lo']:.2f} to {es.loc[2015, 'hi']:.2f} C/decade,
from well below to nearly triple the 1950-2025 rate. Short-window trends are endpoint
noise; the honest answer is the CI, and the honest design fixes the window before looking.""")

# %% ----------------------------------------------------------- hydrology flavour
banner("Potomac: quantile trends (is the distribution changing even if the mean is not?)")
qs = {}
for name, fn in [("annual mean", "mean"), ("annual 10th pct (low flow)", lambda s: s.quantile(.1)),
                 ("annual max", "max")]:
    series = q.groupby(wy).agg(fn).loc[1931:2025]
    mk = mann_kendall(series.values)
    slope, *_ = sens_slope(series.values, series.index.values.astype(float))
    qs[name] = series
    print(f"  {name:28s} MK p={mk.p:.3f}  Sen slope={slope * 10 / series.mean() * 100:+.1f}% of mean per decade")
note("""Trend tests on the mean miss changes in the tails. Testing several annual statistics
(mean, low-flow quantile, maximum) is standard in hydrology; remember it is also multiple
testing (section 02).""")
