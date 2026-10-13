# %% [markdown]
# # 10 Censored inference, seasonal trends and compliance limits
#
# Data: the same Clark Fork River record as section 09 (USGS 12324200 at Deer Lodge, Montana,
# 1969-2023 metals with non-detects, daily flow). Section 03 estimated a mean from censored
# data; this section does the inference that regulators and consultants actually ask for.
#
# **Cheat sheet**
# - Seasonal Kendall (Hirsch, Slack & Smith 1982): Mann-Kendall within each season, summed.
#   Blocks out seasonality instead of letting it inflate the variance; the standard USGS
#   water-quality trend test. Flow-adjust first (LOWESS of C on Q) or the hydrology leaks in.
# - Comparing groups with non-detects: never substitute and run a t-test. Gehan's generalised
#   Wilcoxon scores only the pairs whose order is certain; it reduces to Mann-Whitney without
#   censoring. Peto-Prentice is the alternative when censoring patterns differ between groups.
# - Regression with a censored response: maximum likelihood (Tobit) uses the normal CDF for
#   each non-detect. Akritas-Theil-Sen is the nonparametric line (censored Kendall's S = 0).
# - Upper confidence limit on the mean (ProUCL): the distribution and the skew pick the UCL.
#   Normal -> t; gamma -> gamma UCL; skewed -> Chebyshev at 95/97.5/99 %; non-detects -> KM.
#   Land's H-UCL for lognormal data is unstable and no longer recommended.
# - Prediction and tolerance limits (EPA Unified Guidance) are the compliance-monitoring
#   question: is a new sample consistent with background? They need independent background
#   data; river samples are seasonal and autocorrelated, so treat the example as mechanics.
# - More than two groups: Kruskal-Wallis then Dunn pairwise with Holm or BH; Welch ANOVA /
#   Tukey on logs if the groups are roughly lognormal. Eras of a river are not a designed
#   experiment: flow, method and detection limits all changed with them.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.nonparametric.smoothers_lowess import lowess

from miningstats import data
from miningstats.plotting import banner, note, savefig, setup
from miningstats.stats import (
    akritas_theil_sen,
    censored_regression_mle,
    dunn_test,
    gehan_test,
    km_left_censored,
    prediction_limit,
    seasonal_kendall,
    tolerance_limit,
    ucl_recommend,
    ucl_table,
)

setup()
SEC = "10_censored_inference"
wq = data.clark_fork_wq()
flow = data.clark_fork_flow()


def analyte(char, frac):
    d = wq[(wq.characteristic == char) & (wq.fraction == frac) & (wq.unit == "ug/l")].copy()
    d = d.set_index("date").join(flow, how="left").reset_index().dropna(subset=["flow_cfs"])
    d["year"] = d.date.dt.year
    d["month"] = d.date.dt.month
    d["t"] = d.year + (d.date.dt.dayofyear - 0.5) / 365.25
    d["logc"], d["logq"] = np.log10(d.value), np.log10(d.flow_cfs)
    return d


ERAS = [1968, 1990, 2005, 2015, 2024]
ERA_LABELS = ["1969-90", "1991-05", "2006-15", "2016-23"]


def era_of(years):
    return pd.cut(years, ERAS, labels=ERA_LABELS)


cu = analyte("Copper", "Recoverable")
pb = analyte("Lead", "Dissolved")
as_ = analyte("Arsenic", "Total")
print(f"recoverable copper {len(cu)} samples ({cu.nondetect.sum()} non-detects), dissolved lead {len(pb)} "
      f"({pb.nondetect.sum()} non-detects), total arsenic {len(as_)} ({as_.nondetect.sum()} non-detects)")

# %% ----------------------------------------------------------- seasonal Kendall
banner("Seasonal Kendall on recoverable copper: raw and flow-adjusted")
months_with_data = cu.groupby("month").year.nunique()
print("years sampled per calendar month: " + ", ".join(f"{m}:{n}" for m, n in months_with_data.items()))
# flow adjustment: LOWESS of log C on log Q (Helsel et al. 2020 section 12.5), residuals carry the trend
fitted = lowess(cu.logc, cu.logq, frac=0.5, return_sorted=False)
cu["adj"] = cu.logc - fitted
results = {}
for label, col in [("raw log10 Cu", "logc"), ("flow-adjusted", "adj")]:
    sk = seasonal_kendall(cu[col], cu.year, cu.month)
    results[label] = sk
    ann = cu.groupby("year")[col].median()
    mk = sps.kendalltau(ann.index, ann.values)
    sen_ann = sps.theilslopes(ann.values, ann.index).slope
    print(f"{label:15s} seasonal Kendall tau={sk['tau']:+.3f} z={sk['z']:+.2f} p={sk['p']:.2e} "
          f"Sen={(10 ** sk['sen'] - 1) * 100:+.2f}%/yr   |   annual-median Mann-Kendall tau={mk.statistic:+.3f} "
          f"p={mk.pvalue:.2e} Sen={(10 ** sen_ann - 1) * 100:+.2f}%/yr")
q_trend = sps.kendalltau(cu.t, cu.logq)
print(f"\ndid the sampled flows drift? Kendall tau of log10 Q on time = {q_trend.statistic:+.3f} (p={q_trend.pvalue:.2g})")
per = results["flow-adjusted"]["per_season"].set_index("season")
per["sen_pct"] = (10 ** per.sen - 1) * 100
print("\nflow-adjusted, by season (month):\n" + per.round(3).to_string())
sk_raw, sk_adj = results["raw log10 Cu"], results["flow-adjusted"]
note(f"""Seasonal Kendall compares March with March and June with June, so the runoff cycle never
enters the pairwise comparisons; the annual-median test throws away within-year
information and the whole-series Kendall would let seasonality inflate the variance.
Here every month with enough years agrees on the sign, which is the pattern of a real trend
rather than a seasonal artefact. Flow adjustment moves the evidence from p = {sk_raw['p']:.1g} to
p = {sk_adj['p']:.1g}: copper rides on discharge (section 09), so the flow the sampler happened to
catch is noise around the trend, and removing it is what sharpens the test. The sampled flows
themselves {'drifted' if q_trend.pvalue < 0.05 else 'did not drift'} over the record (tau {q_trend.statistic:+.2f}), so
{'part of the raw trend is a flow artefact' if q_trend.pvalue < 0.05 else 'the raw slope is not a flow artefact, only a noisier estimate'}.
Report the flow-adjusted slope ({(10 ** sk_adj['sen'] - 1) * 100:+.2f}%/yr) with the raw slope
({(10 ** sk_raw['sen'] - 1) * 100:+.2f}%/yr) beside it. Modern USGS practice (WRTDS, Hirsch et al.
2010) fits the same three effects, time, discharge and season, as a weighted regression and
reports flow-normalised concentration; the seasonal Kendall is the hypothesis test it
replaced for estimation.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
ax = axes[0]
for m, g in cu.groupby("month"):
    if months_with_data[m] >= 10:
        ax.plot(g.year, g.adj, ".", ms=4, alpha=.5)
        ann_m = g.groupby("year").adj.median()
        ax.plot(ann_m.index, ann_m.values, "-", lw=.8, alpha=.6, label=f"month {m}")
ax.set_ylabel("flow-adjusted log10 Cu residual")
ax.set_title("Within-season series that the test compares")
ax.legend(ncol=2, fontsize=8)
ax = axes[1]
ax.bar(per.index.astype(str), per.sen_pct, color=["C0" if n >= 10 else "C7" for n in per.n_years])
ax.axhline((10 ** sk_adj["sen"] - 1) * 100, color="C3", ls="--", label="pooled Sen slope")
ax.set_xlabel("season (month); grey = fewer than 10 years")
ax.set_ylabel("Sen slope (%/yr)")
ax.set_title("Per-season slopes agree in sign")
ax.legend()
fig.tight_layout()
savefig(fig, "seasonal_kendall", SEC)

# %% ----------------------------------------------------------- censored two-group test
banner("Two groups with non-detects: dissolved lead, first era vs most recent")
pb["era"] = era_of(pb.year)
for era, g in pb.groupby("era", observed=True):
    dls = sorted(g.loc[g.nondetect, "value"].unique())
    print(f"  {era}: n={len(g)}, non-detects={g.nondetect.mean():.0%}, detection limits {dls}")
two = pb[pb.era.isin([ERA_LABELS[0], ERA_LABELS[-1]])].copy()
gh = gehan_test(two.value.values, two.nondetect.values, two.era.astype(str).values)
higher = gh["groups"][0] if gh["z"] > 0 else gh["groups"][1]
rows = [dict(test="Gehan (censored)", z=gh["z"], p=gh["p"])]
early_dls = sorted(two.loc[(two.era == ERA_LABELS[0]) & two.nondetect, "value"].unique())
early_nd = two.loc[two.era == ERA_LABELS[0], "nondetect"].mean()
late_med = km_left_censored(two.loc[two.era == ERA_LABELS[-1], "value"].values,
                            two.loc[two.era == ERA_LABELS[-1], "nondetect"].values)["median"]
for label, sub in [("DL", two.value), ("DL/2", np.where(two.nondetect, two.value / 2, two.value)),
                   ("0", np.where(two.nondetect, 0, two.value))]:
    a = np.asarray(sub)[two.era.values == ERA_LABELS[0]]
    b = np.asarray(sub)[two.era.values == ERA_LABELS[-1]]
    mw = sps.mannwhitneyu(a, b, method="asymptotic", use_continuity=False)
    tt = sps.ttest_ind(np.log10(np.maximum(a, 0.01)), np.log10(np.maximum(b, 0.01)), equal_var=False)
    rows.append(dict(test=f"Mann-Whitney, substitute {label}", z=np.nan, p=mw.pvalue))
    rows.append(dict(test=f"Welch t on logs, substitute {label}", z=tt.statistic, p=tt.pvalue))
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3g}"))
print(f"Gehan: {higher} is the higher group (positive scores favour the first label)")
km_by_era = {era: km_left_censored(g.value.values, g.nondetect.values)
             for era, g in pb.groupby("era", observed=True)}
print("\nKM medians by era: " + ", ".join(f"{e}: {k['median']:.2f} ug/L" for e, k in km_by_era.items()))
note(f"""The first era reported {early_nd:.0%} of lead results as below {' or '.join(f'{d:g}' for d in early_dls)} ug/L;
the latest era detects almost everything, with a median of {late_med:.2f} ug/L. Substituting the
detection limit turns the early era into a column of 1s and 5s and makes the decrease look
overwhelming (p ~ 1e-15); substituting zero turns the same results into a column of zeros and
reverses the sign. Gehan's test only uses pairs whose order is certain: an early '< 1' against
a recent 0.3 is not evidence either way, because the early value could have been 0.1, while an
early detected 3 against a recent 0.3 is. What survives is a decrease with p = {gh['p']:.2g}: real,
but with far less evidence than substitution claims. When censoring fractions differ
between groups, as here, the Peto-Prentice score test weights the pairs more evenly than Gehan;
both are in Helsel (2012) chapter 9 and the NADA package.""")

fig, ax = plt.subplots(figsize=(8, 4.5))
for era, km in km_by_era.items():
    ax.step(km["x"], km["F"], where="post", label=f"{era} (KM)")
for dl in pb.loc[pb.nondetect, "value"].unique():
    ax.axvline(dl, color="C7", lw=.5, alpha=.4)
ax.set_xscale("log")
ax.set_xlabel("dissolved lead (ug/L); grey lines = detection limits used")
ax.set_ylabel("cumulative probability")
ax.set_title("Censored CDFs by era: the curves start where the detection limits allow")
ax.legend()
savefig(fig, "censored_groups", SEC)

# %% ----------------------------------------------------------- censored regression
banner("Regression and trend with a censored response: dissolved lead")
X = np.column_stack([pb.logq, pb.t - 2000])
mle = censored_regression_mle(X, pb.logc, pb.nondetect.values)
names = ["intercept", "log10 flow", "year - 2000"]
print("Tobit / MLE regression of log10 Pb on log10 flow and year (non-detects as 'below DL'):")
for nm, b, se in zip(names, mle["beta"], mle["se"]):
    print(f"  {nm:12s} {b:+.4f}  (SE {se:.4f})")
print(f"  residual sigma {mle['sigma']:.3f}; converged {mle['converged']}")
import statsmodels.api as sm  # noqa: E402

det = ~pb.nondetect.values
ols_det = sm.OLS(pb.logc[det], sm.add_constant(X[det])).fit()
y_half = np.where(pb.nondetect, np.log10(pb.value / 2), pb.logc)
ols_half = sm.OLS(y_half, sm.add_constant(X)).fit()
comp = pd.DataFrame({"MLE (censored)": mle["beta"], "OLS detects only": ols_det.params.values,
                     "OLS substitute DL/2": np.asarray(ols_half.params)}, index=names)
print("\ncoefficients by method:\n" + comp.round(4).to_string())
pct = lambda b: (10 ** b - 1) * 100  # noqa: E731
print(f"\ntime trend, holding flow constant: MLE {pct(mle['beta'][2]):+.2f}%/yr, detects-only "
      f"{pct(ols_det.params.iloc[2]):+.2f}%/yr, DL/2 substitution {pct(np.asarray(ols_half.params)[2]):+.2f}%/yr")
ats = akritas_theil_sen(pb.t.values, pb.logc.values, pb.nondetect.values, n_perm=1000, rng=0)
mle_t = censored_regression_mle(pb.t.values - 2000, pb.logc.values, pb.nondetect.values)
print(f"\ntrend without the flow term: Akritas-Theil-Sen slope {pct(ats['slope']):+.2f}%/yr (censored tau "
      f"{ats['tau']:+.3f}, permutation p={ats['p']:.3g}); MLE slope {pct(mle_t['beta'][1]):+.2f}%/yr "
      f"(SE {pct(mle_t['se'][1]):.2f})")
by_era = pb.groupby("era", observed=True).apply(
    lambda g: km_left_censored(g.value.values, g.nondetect.values)["median"], include_groups=False)
print("KM median by era: " + ", ".join(f"{e}: {m:.2f}" for e, m in by_era.items()))
note(f"""Dropping the non-detects keeps only the early samples that happened to be high enough to
measure, so the early years look high and the decline is exaggerated
({pct(ols_det.params.iloc[2]):+.2f}%/yr against the MLE's {pct(mle['beta'][2]):+.2f}%/yr); substituting half
the limit injects the laboratory's history into the response. The MLE treats each non-detect
as the probability that a normal residual falls below the limit, and with log flow in the
model it gives a flow-adjusted decline of {pct(mle['beta'][2]):+.2f}%/yr. Without the flow term both the
MLE ({pct(mle_t['beta'][1]):+.2f}%/yr) and the distribution-free Akritas-Theil-Sen line
({pct(ats['slope']):+.2f}%/yr) see much less: the record is a drop between the first two eras and a flat
or slightly rising tail as detection limits fell, which no single straight line fits well.
The disagreement is the finding: the trend depends on the flow term, and a monotone trend is
the wrong model for a step. Section 09's change-point tools are the follow-up.""")

fig, ax = plt.subplots(figsize=(10, 4.5))
ax.semilogy(pb.t[det], pb.value[det], ".", ms=4, alpha=.6, label="detected")
ax.semilogy(pb.t[~det], pb.value[~det], "v", ms=4, color="C3", alpha=.6, label="non-detect (at DL)")
tt = np.linspace(pb.t.min(), pb.t.max(), 100)
q_med = np.median(pb.logq)
for lab, b, style in [("MLE, censored", mle["beta"], "-"), ("OLS detects only", ols_det.params.values, "--"),
                      ("OLS DL/2", np.asarray(ols_half.params), ":")]:
    ax.semilogy(tt, 10 ** (b[0] + b[1] * q_med + b[2] * (tt - 2000)), style, lw=1.6, label=f"{lab} at median flow")
ax.set_ylabel("dissolved lead (ug/L)")
ax.set_title("Three fits to the same censored record")
ax.legend(fontsize=8)
savefig(fig, "censored_regression", SEC)

# %% ----------------------------------------------------------- UCL on the mean
banner("Upper confidence limit on the mean (ProUCL logic): total arsenic, recent era")
as_["era"] = era_of(as_.year)
recent = as_[as_.era == ERA_LABELS[-1]]
v = recent.value.values
print(f"total arsenic {ERA_LABELS[-1]}: n={len(v)}, mean={v.mean():.2f}, median={np.median(v):.2f}, "
      f"sd={v.std(ddof=1):.2f}, skew={sps.skew(v):.2f}, sd of logs={np.std(np.log(v), ddof=1):.2f}")
tab = ucl_table(v, rng=0)
rec = ucl_recommend(v)
print("goodness of fit: " + ", ".join(f"{k}={p:.3f}" for k, p in rec["fits"].items()))
print(pd.Series(tab).round(2).to_string())
print(f"\nrecommended: {rec['choice']}  ({rec['rationale']})")
cd = analyte("Cadmium", "Dissolved")
cd["era"] = era_of(cd.year)
cd_era = cd[cd.era == ERA_LABELS[1]]
cd_tab = ucl_table(cd_era.value.values, cd_era.nondetect.values)
cd_rec = ucl_recommend(cd_era.value.values, cd_era.nondetect.values)
print(f"\ndissolved cadmium {ERA_LABELS[1]}: n={len(cd_era)}, non-detects={cd_era.nondetect.mean():.0%}")
print(pd.Series(cd_tab).round(3).to_string())
print(f"recommended: {cd_rec['choice']}  ({cd_rec['rationale']})")
note(f"""The UCL95 is the exposure-point concentration in a risk assessment: a one-sided limit that
the true mean exceeds with 5% probability if the sampling is representative. The spread
between the Student-t ({tab['Student-t UCL95']:.2f}) and the Chebyshev 99% ({tab['Chebyshev UCL99']:.2f})
is the cost of not knowing the distribution. ProUCL's decision table picks by fit and skew
because the t-interval under-covers for skewed data and Land's lognormal H-UCL can exceed
the maximum observation. With non-detects the mean and its standard error come from
Kaplan-Meier and the same choice between t and Chebyshev follows. One caveat the software
does not print: a UCL assumes independent samples, and {len(v)} river samples over eight years
are seasonal and serially correlated, so the effective n is smaller than it looks.""")

fig, ax = plt.subplots(figsize=(9, 4.2))
names_u = [k for k in tab if k != "mean"]
vals_u = [tab[k] for k in names_u]
colors = ["C3" if k == rec["choice"] else "C0" for k in names_u]
ax.barh(names_u, vals_u, color=colors, alpha=.85)
ax.axvline(tab["mean"], color="k", lw=1, label=f"sample mean {tab['mean']:.2f}")
ax.axvline(v.max(), color="C7", ls="--", lw=1, label=f"maximum observation {v.max():.1f}")
ax.set_xlabel("total arsenic (ug/L); red = ProUCL-style recommendation")
ax.set_title("The same mean, eight upper confidence limits")
ax.legend(loc="lower right")
ax.invert_yaxis()
savefig(fig, "ucl_comparison", SEC)

# %% ----------------------------------------------------------- prediction and tolerance limits
banner("Prediction and tolerance limits (Unified Guidance mechanics): is a new sample unusual?")
bg = as_[as_.era == ERA_LABELS[2]].value.values
new = as_[as_.era == ERA_LABELS[3]].value.values
lbg = np.log(bg)
pl = prediction_limit(lbg, k_future=1)
pl_k = prediction_limit(lbg, k_future=len(new))
tl = tolerance_limit(lbg, coverage=0.95, confidence=0.95)
print(f"background = {ERA_LABELS[2]} total arsenic, n={len(bg)} (log scale; Shapiro p on logs = {sps.shapiro(lbg).pvalue:.3f})")
print(f"  upper prediction limit, 1 future sample:    {np.exp(pl['normal']):.2f} ug/L   (nonparametric: background max "
      f"{np.exp(pl['nonparametric']):.2f}, confidence {pl['nonparametric_confidence']:.0%})")
print(f"  upper prediction limit, {len(new)} future samples (Bonferroni): {np.exp(pl_k['normal']):.2f} ug/L")
print(f"  upper tolerance limit, 95% coverage / 95% confidence: {np.exp(tl):.2f} ug/L")
exc1 = (new > np.exp(pl["normal"])).sum()
exck = (new > np.exp(pl_k["normal"])).sum()
print(f"  {ERA_LABELS[3]} samples (n={len(new)}): {exc1} exceed the one-sample limit, {exck} exceed the "
      f"{len(new)}-sample limit, {(new > np.exp(tl)).sum()} exceed the tolerance limit")
note(f"""A prediction limit bounds the next sample(s); a tolerance limit bounds a stated fraction of
the population. Compliance monitoring under the Unified Guidance tests each new downgradient
sample against a prediction limit built from independent background, with the number of
future comparisons folded into the limit (Bonferroni here; the guidance uses exact kappa
tables) and usually a retest rule. Testing {len(new)} new samples one at a time against a
one-sample limit would expect about {len(new) * 0.05:.1f} false exceedances; the k-sample limit is wider
for that reason. River samples are not independent background, so {exc1} exceedance{'s' if exc1 != 1 else ''} of the
one-sample limit {'are' if exc1 != 1 else 'is'} a demonstration of the mechanics, not a finding.""")

# %% ----------------------------------------------------------- k groups
banner("More than two groups: recoverable copper across four eras")
cu["era"] = era_of(cu.year)
groups = [g.logc.values for _, g in cu.groupby("era", observed=True)]
labels = [str(e) for e, _ in cu.groupby("era", observed=True)]
kw = sps.kruskal(*groups)
ag = sps.alexandergovern(*groups)
print(f"Kruskal-Wallis H={kw.statistic:.1f} p={kw.pvalue:.2e};  Welch-type (Alexander-Govern) ANOVA on logs p={ag.pvalue:.2e}")
dn = dunn_test(groups, labels, adjust="holm")
tk = sps.tukey_hsd(*groups)
dn["tukey_p (logs)"] = [tk.pvalue[i, j] for i in range(4) for j in range(i + 1, 4)]
dn["median ratio"] = [10 ** (np.median(groups[j]) - np.median(groups[i])) for i in range(4) for j in range(i + 1, 4)]
print(dn.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
note("""Kruskal-Wallis says at least one era differs; Dunn's pairwise z-tests on the pooled ranks say
which, with Holm keeping the family-wise error at 5% across six comparisons. Tukey on the logs
agrees here because log copper is close to normal within eras. The usual mistakes are
running six Mann-Whitney tests with no correction, or reading a significant omnibus test as
'every era differs'. Eras are an observational grouping: flow regime, laboratory method and
detection limits also changed between them, which is why section 09 flow-adjusts first.""")
