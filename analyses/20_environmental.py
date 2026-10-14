# %% [markdown]
# # 20 Environmental monitoring: a mining-impacted river
#
# Data: Clark Fork River at Deer Lodge, Montana (USGS 12324200), downstream of the Butte and
# Anaconda copper district, the largest Superfund complex in the US. Metals and pH from the
# Water Quality Portal (1969-2023), daily discharge from NWIS.
#
# **Cheat sheet**
# - Concentration depends on flow. Plot C vs Q on log-log axes: slope < 0 is dilution (point
#   sources), slope > 0 is mobilisation (tailings flushed by high flow). Flow-adjust before
#   testing for trends or the hydrology masquerades as a remediation effect.
# - Metals are heavy-tailed and seasonal; use logs, medians and quantiles.
# - Trend tests: Mann-Kendall / Kendall's tau on annual medians or on flow-adjusted residuals.
#   Seasonal Kendall if sampling is monthly. Report the Sen slope as %/yr on the log scale.
# - Compliance: fraction of samples exceeding a criterion per period with a binomial CI; a
#   criterion that depends on hardness (copper, zinc) needs the paired hardness value.
# - Change points: CUSUM of annual medians; relate to known interventions (remediation).
# - Non-detects (section 14) and changing detection limits bias simple trend tests toward
#   finding decreases; use censored methods for low-level analytes.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import km_left_censored

setup()
SEC = "20_environmental"
wq = data.clark_fork_wq()
flow = data.clark_fork_flow()
cu = wq[(wq.characteristic == "Copper") & (wq.fraction == "Recoverable") & (wq.unit == "ug/l")].copy()
cu = cu.set_index("date").join(flow, how="left").reset_index()
cu = cu.dropna(subset=["flow_cfs"])
cu["year"] = cu.date.dt.year
cu["logc"], cu["logq"] = np.log10(cu.value), np.log10(cu.flow_cfs)
print(f"{len(cu)} total-recoverable copper samples {cu.year.min()}-{cu.year.max()} with matched daily flow; "
      f"{cu.nondetect.sum()} non-detects")

# %% ----------------------------------------------------------- series
banner("Copper and flow through time")
ann = cu.groupby("year").agg(median_cu=("value", "median"), n=("value", "size"), p90=("value", lambda s: s.quantile(.9)))
print(ann[ann.n >= 4].iloc[::5].round(1).to_string())
fig, ax = plt.subplots(figsize=(12, 4.5))
ax.semilogy(cu.date, cu.value, ".", ms=4, alpha=.6, label="total recoverable Cu (ug/L)")
ax.set_ylabel("Cu (ug/L)")
ax2 = ax.twinx()
ax2.semilogy(flow.index, flow.flow_cfs, color="C7", lw=.4, alpha=.6)
ax2.set_ylabel("daily flow (cfs)", color="C7")
ax2.grid(False)
ax.set_title("Clark Fork at Deer Lodge: copper spikes ride on spring runoff")
ax.legend(loc="upper right")
savefig(fig, "timeseries", SEC)

# %% ----------------------------------------------------------- C-Q relationship
banner("Concentration-discharge relationship")
fit = sm.OLS(cu.logc, sm.add_constant(cu.logq)).fit()
slope = fit.params.iloc[1]
print(f"log10 C = {fit.params.iloc[0]:.2f} + {slope:.2f} log10 Q   (R2={fit.rsquared:.2f}, n={len(cu)})")
print(f"a doubling of flow multiplies copper by {2 ** slope:.2f}")
cu["month"] = cu.date.dt.month
monthly = cu.groupby("month").value.median()
print("monthly median Cu (ug/L): " + ", ".join(f"{m}:{v:.0f}" for m, v in monthly.items()))
note(f"""The C-Q slope is {slope:+.2f}: concentration *rises* with flow. That is mobilisation, not
dilution: high water scours copper-laden tailings and streambank sediments deposited by a
century of mining, the defining signature of this Superfund river. A point-source discharge
would show the opposite sign. Any trend analysis must first remove this flow dependence.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
sc = axes[0].scatter(cu.flow_cfs, cu.value, c=cu.year, cmap="viridis", s=12)
xs = np.logspace(np.log10(cu.flow_cfs.min()), np.log10(cu.flow_cfs.max()), 50)
axes[0].plot(xs, 10 ** (fit.params.iloc[0] + slope * np.log10(xs)), "C3")
axes[0].set_xscale("log")
axes[0].set_yscale("log")
axes[0].set_xlabel("daily flow (cfs)")
axes[0].set_ylabel("Cu (ug/L)")
axes[0].set_title(f"C-Q: slope {slope:+.2f} (mobilisation)")
fig.colorbar(sc, ax=axes[0], label="year")
cu.boxplot(column="value", by="month", ax=axes[1], showfliers=False)
axes[1].set_yscale("log")
axes[1].set_title("Seasonality: spring runoff months")
axes[1].set_ylabel("Cu (ug/L)")
fig.suptitle("")
fig.tight_layout()
savefig(fig, "concentration_discharge", SEC)

# %% ----------------------------------------------------------- flow-adjusted trend
banner("Trend: raw vs flow-adjusted")
cu["resid"] = fit.resid
ann_med = cu.groupby("year").agg(raw=("logc", "median"), adj=("resid", "median"), n=("logc", "size"))
ann_med = ann_med[ann_med.n >= 3]
for col, label in [("raw", "raw log10 Cu"), ("adj", "flow-adjusted residual")]:
    tau, p = sps.kendalltau(ann_med.index, ann_med[col])
    sen = sps.theilslopes(ann_med[col], ann_med.index)
    print(f"{label:26s} Kendall tau={tau:+.2f} p={p:.4f}  Sen slope={sen.slope * 100:+.2f} log10-units/century "
          f"= {(10 ** sen.slope - 1) * 100:+.1f}%/yr")
# break the record into eras
eras = pd.cut(cu.year, [1968, 1990, 2005, 2015, 2024], labels=["1969-90", "1991-05", "2006-15", "2016-23"])
era_tab = cu.groupby(eras, observed=True).agg(n=("value", "size"), median_cu=("value", "median"),
                                              median_flow=("flow_cfs", "median"),
                                              adj_median=("resid", "median"))
era_tab["adj_median"] = 10 ** era_tab.adj_median
print("\nby era (adj_median = flow-adjusted multiplier vs long-term C-Q line):\n" + era_tab.round(2).to_string())
note("""Raw and flow-adjusted trends agree in sign here, so the decline is not a flow artefact;
in drier or wetter decades they can disagree completely. The era table is the plain-language
version: the flow-adjusted copper level in the most recent decade relative to the long-term
relationship. Sampling frequency and detection limits also changed across eras, which a
single trend test silently ignores.""")

# %% ----------------------------------------------------------- exceedances
banner("Compliance statistics against an aquatic-life criterion")
CRIT = 13.0  # ug/L: EPA chronic (CCC) copper criterion at 100 mg/L hardness, illustrative only
cu["exceed"] = cu.value > CRIT
rows = []
for era, g in cu.groupby(eras, observed=True):
    k, n = int(g.exceed.sum()), len(g)
    ci = sps.binomtest(k, n).proportion_ci(method="wilson")
    rows.append(dict(era=era, n=n, exceed=k, frac=k / n, lo=ci.low, hi=ci.high))
ex = pd.DataFrame(rows).set_index("era")
print(ex.round(3).to_string())
hi_q = cu.flow_cfs.quantile(.75)
share_hi = cu[cu.exceed].flow_cfs.ge(hi_q).mean()
print(f"\nfraction of exceedances occurring on top-quartile flow days: {share_hi:.0%} "
      f"(exceedance rate on those days {cu[cu.flow_cfs >= hi_q].exceed.mean():.0%} vs {cu[cu.flow_cfs < hi_q].exceed.mean():.0%} otherwise)")
note(f"""Exceedance fractions with Wilson intervals are the honest compliance summary; 'number of
exceedances this year' is a count without a denominator. The decline from {ex.frac.iloc[0]:.0%} to
{ex.frac.iloc[-1]:.0%} is real but the river still exceeds this benchmark in most samples at any
flow: a chronic sediment-sourced problem, not an event-driven one, so the C-Q slope tells you
when it is worst, not when it is compliant. The criterion ({CRIT} ug/L) is the hardness-100
chronic value for illustration; the real standard moves with each sample's hardness, so
compliance assessment needs paired hardness measurements.""")

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.bar(ex.index.astype(str), ex.frac, yerr=[ex.frac - ex.lo, ex.hi - ex.frac], capsize=4, alpha=.8)
ax.set_ylabel(f"fraction of samples > {CRIT:.0f} ug/L")
ax.set_title("Copper exceedance fraction by era (Wilson 95% CI)")
savefig(fig, "exceedances", SEC)

# %% ----------------------------------------------------------- change point & censored analyte
banner("CUSUM of annual flow-adjusted medians, and a censored low-level analyte")
am = ann_med.adj - ann_med.adj.mean()
cs = am.cumsum()
cp_year = cs.idxmax()
print(f"CUSUM of flow-adjusted annual medians peaks in {cp_year}: regime shift from above-average to below-average residuals")
cd = wq[(wq.characteristic == "Cadmium") & (wq.unit == "ug/l")]
cd_eras = pd.cut(cd.date.dt.year, [1968, 1990, 2005, 2024], labels=["1969-90", "1991-05", "2006-23"])
for era, g in cd.groupby(cd_eras, observed=True):
    km = km_left_censored(g.value.values, g.nondetect.values)
    print(f"  cadmium {era}: n={len(g)}, non-detects={g.nondetect.mean():.0%}, KM median={km['median']:.2f} ug/L, "
          f"KM mean={km['mean']:.2f}, naive mean (DL as value)={g.value.mean():.2f}")
note(f"""The CUSUM turns in {cp_year}: flow-adjusted copper ran above its long-term level before and
below it after. Upstream, the Warm Springs Ponds treatment system (which settles metals out
of Silver Bow Creek before it becomes the Clark Fork) was rebuilt and expanded in the early-
to-mid 1990s; the 2010s floodplain tailings removals are a second candidate for the later
decline. A change point found after the fact is a hypothesis to test against the intervention
record, not a proof. For cadmium the non-detect fraction was 60-70% in the early eras with
changing detection limits, so only censored estimators (Kaplan-Meier) give comparable
statistics across eras; the naive 'DL as value' means are 2-3x too high.""")

fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(cs.index, cs.values, "o-", ms=3)
ax.axvline(cp_year, color="C3", ls="--", label=f"CUSUM max {cp_year}")
ax.set_ylabel("cumulative sum of (annual median residual - mean)")
ax.set_title("CUSUM of flow-adjusted copper: slope change marks a regime shift")
ax.legend()
savefig(fig, "cusum_changepoint", SEC)
