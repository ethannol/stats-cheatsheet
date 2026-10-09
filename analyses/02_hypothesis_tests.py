# %% [markdown]
# # 02 Hypothesis testing
#
# Question: are Central Park summer (JJA) nights warmer in 1991-2020 than in 1951-1980?
# Contrast: summer daytime highs (TMAX) at this station show essentially no change, a known
# local effect (tree canopy growth around the Central Park instrument shelter).
#
# **Cheat sheet**
# - p-value = P(data at least this extreme | H0). It is NOT P(H0 | data) and NOT effect size.
#   A null result is "no evidence of a difference", not "evidence of no difference".
# - Pick the unit of analysis so observations are approximately independent. Daily values
#   within a summer are autocorrelated; use seasonal/annual means (n ~ 30 per period).
# - Two-sample means: Welch t-test (unequal variances; default choice), Student t only if
#   variances are equal. Mann-Whitney U for ordinal / heavy-tailed data (tests stochastic
#   dominance, not medians). Permutation test when you want zero distributional assumptions.
# - Paired designs (same station, before/after): paired t / Wilcoxon signed-rank.
# - Always report: effect size (Cohen's d, difference with CI), not just p.
# - Variances: Levene (robust) rather than F-test. Distributions: KS / Anderson-Darling 2-sample.
# - Counts/proportions: chi-square test of independence, Fisher exact for small tables,
#   two-proportion z-test.
# - Multiple comparisons: Bonferroni (controls FWER, conservative), Benjamini-Hochberg
#   (controls FDR, standard for many-tests screening).
# - Power: decide n before collecting. Power = 1 - beta = P(reject | H1 true).

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower
from statsmodels.stats.proportion import proportions_ztest

from earthstats import data
from earthstats.plotting import banner, note, savefig, setup
from earthstats.stats import cohens_d, permutation_test

setup()
SEC = "02_hypothesis_tests"
cp = data.ghcn_daily()
jja = cp[cp.index.month.isin([6, 7, 8])]


def period_means(var):
    s = jja[var].resample("YE").mean()
    s.index = s.index.year
    return s.loc[1951:1980], s.loc[1991:2020]


# %% ----------------------------------------------------------- two-sample tests
banner("Summer mean TMIN (nights): 1951-1980 vs 1991-2020 (n=30 summers each)")
early, late = period_means("TMIN")
diff = late.mean() - early.mean()
print(f"early: mean={early.mean():.2f} sd={early.std():.2f} | late: mean={late.mean():.2f} "
      f"sd={late.std():.2f} | diff={diff:+.2f} C")

t_welch = sps.ttest_ind(late, early, equal_var=False)
t_student = sps.ttest_ind(late, early, equal_var=True)
mw = sps.mannwhitneyu(late, early, alternative="two-sided")
lev = sps.levene(late, early)
ks = sps.ks_2samp(late, early)
obs, p_perm, null = permutation_test(late.values, early.values, n_perm=20000, rng=0)
d = cohens_d(late, early)
ci = t_welch.confidence_interval()

print(f"Welch t:        t={t_welch.statistic:.2f}  p={t_welch.pvalue:.4f}  "
      f"95% CI diff=({ci.low:.2f}, {ci.high:.2f})")
print(f"Student t:      t={t_student.statistic:.2f}  p={t_student.pvalue:.4f}")
print(f"Mann-Whitney U: U={mw.statistic:.0f}  p={mw.pvalue:.4f}")
print(f"Permutation:    diff={obs:+.2f}  p={p_perm:.4f}  ({len(null)} shuffles)")
print(f"Levene (equal var?): W={lev.statistic:.2f} p={lev.pvalue:.3f}")
print(f"KS 2-sample:    D={ks.statistic:.3f} p={ks.pvalue:.4f}")
print(f"Cohen's d = {d:.2f}")
note(f"""All four tests agree, the usual outcome when n is decent and the effect is real. Report
the test whose assumptions you can defend plus the effect size: summer nights are
{diff:.1f} C warmer with d = {d:.2f} (large). Levene finds no variance difference so
Student and Welch coincide; prefer Welch by default anyway, it costs almost nothing.""")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(early, bins=12, alpha=.6, label="1951-1980")
axes[0].hist(late, bins=12, alpha=.6, label="1991-2020")
axes[0].set_xlabel("JJA mean TMIN (C)")
axes[0].legend()
axes[0].set_title("Two samples of 30 summers")
axes[1].hist(null, bins=60, alpha=.8, density=True)
axes[1].axvline(obs, color="C3", lw=2, label=f"observed diff={obs:+.2f}")
axes[1].axvline(-obs, color="C3", lw=1, ls="--")
axes[1].set_title(f"Permutation null distribution (p={p_perm:.4f})")
axes[1].set_xlabel("mean difference under label shuffling")
axes[1].legend()
savefig(fig, "two_sample", SEC)

# %% ----------------------------------------------------------- a null result
banner("Contrast: summer mean TMAX (daytime highs)")
e2, l2 = period_means("TMAX")
t2 = sps.ttest_ind(l2, e2, equal_var=False)
ci2 = t2.confidence_interval()
print(f"diff={l2.mean() - e2.mean():+.2f} C  Welch p={t2.pvalue:.3f}  "
      f"95% CI=({ci2.low:.2f}, {ci2.high:.2f})  d={cohens_d(l2, e2):.2f}  "
      f"Levene p={sps.levene(l2, e2).pvalue:.3f}")
note(f"""A null result. The CI ({ci2.low:.2f}, {ci2.high:.2f}) is the honest summary: the data
are compatible with anything from a small cooling to a ~0.4 C warming. 'Not significant'
does not mean 'no change'. Note also Levene p={sps.levene(l2, e2).pvalue:.2f}: summer
TMAX variance increased even though the mean did not, which a t-test is blind to.""")

# %% ----------------------------------------------------------- the wrong unit of analysis
banner("Pitfall: testing on daily values instead of seasonal means")
early_d = jja.loc["1951":"1980", "TMIN"].dropna()
late_d = jja.loc["1991":"2020", "TMIN"].dropna()
td = sps.ttest_ind(late_d, early_d, equal_var=False)
r1 = jja.TMIN.autocorr(1)
print(f"Daily values: n={len(early_d)} vs {len(late_d)}, t={td.statistic:.1f}, p={td.pvalue:.1e}")
print(f"Seasonal means: n=30 vs 30, t={t_welch.statistic:.1f}, p={t_welch.pvalue:.1e}")
print(f"lag-1 autocorrelation of JJA daily TMIN: {r1:.2f}")
note("""Same data, same mean difference, but the daily-level p-value is ~14 orders of magnitude
smaller because it pretends to have 2,760 independent observations per period. Days
within a summer are autocorrelated (and share the summer's anomaly), so the effective
sample size is far lower. Aggregate to the unit that is plausibly independent, or model
the dependence explicitly (mixed model, GLS).""")

# %% ----------------------------------------------------------- bootstrap CI for the difference
rng = np.random.default_rng(0)
boots = np.array([rng.choice(late, len(late)).mean() - rng.choice(early, len(early)).mean()
                  for _ in range(5000)])
print(f"Bootstrap 95% CI for difference in means: ({np.percentile(boots, 2.5):.2f}, "
      f"{np.percentile(boots, 97.5):.2f})  vs Welch CI ({ci.low:.2f}, {ci.high:.2f})")

# %% ----------------------------------------------------------- proportions / chi-square
banner("Proportions: warm nights (TMIN >= 21.1 C / 70 F) by period")
warm_e = (early_d >= 21.1).sum()
warm_l = (late_d >= 21.1).sum()
z, p = proportions_ztest([warm_l, warm_e], [len(late_d), len(early_d)])
table = np.array([[warm_l, len(late_d) - warm_l], [warm_e, len(early_d) - warm_e]])
chi2, p_chi, dof, _ = sps.chi2_contingency(table)
odds_ratio = (table[0, 0] * table[1, 1]) / (table[0, 1] * table[1, 0])
print(f"warm-night fraction: early={warm_e / len(early_d):.3f} late={warm_l / len(late_d):.3f}")
print(f"two-proportion z={z:.2f} p={p:.2e} | chi2={chi2:.1f} dof={dof} p={p_chi:.2e}")
print(f"odds ratio (late vs early) = {odds_ratio:.2f}")
warm_per_summer = (jja.TMIN >= 21.1).resample("YE").sum()
warm_per_summer.index = warm_per_summer.index.year
a, b = warm_per_summer.loc[1951:1980], warm_per_summer.loc[1991:2020]
mw2 = sps.mannwhitneyu(b, a)
print(f"Warm nights per summer: early median={a.median():.0f} late median={b.median():.0f}  "
      f"Mann-Whitney p={mw2.pvalue:.4f}  (n=30 vs 30 summers)")
note("""Chi-square and the two-proportion z-test are the same test for a 2x2 table (z^2 = chi2).
Both inherit the daily-independence problem. The cleaner design is one count per summer
compared with Mann-Whitney or a Poisson / negative-binomial GLM with a period term.""")

# %% ----------------------------------------------------------- multiple comparisons
banner("Multiple comparisons: monthly TMIN change, 12 tests")
rows = []
monthly = cp.TMIN.resample("ME").mean()
for m in range(1, 13):
    s = monthly[monthly.index.month == m]
    s.index = s.index.year
    a, b = s.loc[1951:1980], s.loc[1991:2020]
    t = sps.ttest_ind(b, a, equal_var=False)
    rows.append(dict(month=m, diff=b.mean() - a.mean(), p_raw=t.pvalue))
mt = pd.DataFrame(rows)
mt["p_bonferroni"] = multipletests(mt.p_raw, method="bonferroni")[1]
mt["p_bh_fdr"] = multipletests(mt.p_raw, method="fdr_bh")[1]
mt["sig_raw"] = mt.p_raw < .05
mt["sig_bonf"] = mt.p_bonferroni < .05
mt["sig_bh"] = mt.p_bh_fdr < .05
print(mt.round(4).to_string(index=False))
print(f"\nsignificant months: raw={mt.sig_raw.sum()} Bonferroni={mt.sig_bonf.sum()} "
      f"BH-FDR={mt.sig_bh.sum()}")
note("""Running 12 tests at alpha=0.05 gives a 46% chance of at least one false positive under
H0. Bonferroni divides alpha by m (controls the family-wise error rate; conservative when
tests are correlated). Benjamini-Hochberg controls the expected fraction of false
discoveries among rejections; it is the default for screening many hypotheses. Months
that survive BH but not Bonferroni are the typical grey zone.""")

# %% ----------------------------------------------------------- power analysis
banner("Power analysis: how many summers to detect a given TMIN shift?")
sd = np.sqrt((early.var() + late.var()) / 2)
pw = TTestIndPower()
for delta in [0.25, 0.5, 1.0]:
    n = pw.solve_power(effect_size=delta / sd, alpha=.05, power=.8, alternative="two-sided")
    print(f"  shift {delta:.2f} C (d={delta / sd:.2f}): n={np.ceil(n):.0f} per group for 80% power")
achieved = pw.power(effect_size=abs(d), nobs1=30, alpha=.05)
print(f"  achieved power for observed d={d:.2f} with n=30: {achieved:.3f}")

fig, ax = plt.subplots(figsize=(7, 4))
ns = np.arange(5, 120)
for delta in [0.25, 0.5, 1.0]:
    ax.plot(ns, pw.power(effect_size=delta / sd, nobs1=ns, alpha=.05), label=f"shift {delta} C")
ax.axhline(.8, color="k", ls="--", lw=.8)
ax.set_xlabel("n per group (years)")
ax.set_ylabel("power")
ax.set_title("Power curves, two-sample t, alpha=0.05")
ax.legend()
savefig(fig, "power_curves", SEC)
note("""Power scales with (effect size / SD) and sqrt(n). A 0.25 C shift needs >100 summers per
group at this station's interannual variability: the signal-to-noise problem is why
climate detection relies on spatial averaging and long records rather than single sites.
Post-hoc 'observed power' is just a transform of the p-value; plan power in advance.""")
