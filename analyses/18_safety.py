# %% [markdown]
# # 18 Mine safety: rates, exposure, overdispersion, small numbers
#
# Data: MSHA open data, every reportable injury at US mines since 2000 with annual employee
# hours per mine. Real data; rates here are computed directly and will differ slightly from
# MSHA's published figures (which apply additional inclusion rules and contractor adjustments).
#
# **Cheat sheet**
# - Counts mean nothing without exposure. Incidence rate = injuries x 200,000 / hours worked
#   (200,000 h = 100 full-time workers for a year). Fatalities use the same denominator.
# - Rates are Poisson-like but overdispersed across mines (heterogeneous risk); the negative
#   binomial, or a quasi-Poisson, gives honest standard errors. Check var/mean of counts after
#   accounting for exposure.
# - Trend and comparisons: Poisson/NB GLM with log(hours) as an offset; coefficients are log
#   rate ratios. Rate ratio CI from exp(beta +- 1.96 SE).
# - Small mines have few injuries: a single incident makes their rate the worst in the country,
#   and zero injuries make it the best. Shrink small-sample rates toward the fleet mean
#   (empirical Bayes, Gamma-Poisson) before ranking anyone.
# - u-charts for rates with variable exposure have limits that widen in years with fewer hours.
# - Lost days are heavy-tailed: report medians and the fraction above thresholds, not means.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps
from scipy.optimize import minimize

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import u_chart_limits

setup()
SEC = "18_safety"
acc = data.msha_accidents()
prod = data.msha_production()
mines = data.msha_mines()
YEARS = (2000, 2024)
INJURY_CODES = [1, 2, 3, 4, 5, 6]  # fatal .. no days lost; excludes 'accident only' and illnesses
acc = acc[acc.year.between(*YEARS) & acc.degree_injury_cd.isin(INJURY_CODES)]
prod = prod[prod.year.between(*YEARS)]
hours = (
    prod.groupby(["mine_id", "year"])
    .agg(
        hours=("annual_hrs", "sum"),
        empl=("avg_annual_empl", "sum"),
        sector=("coal_metal_ind", "first"),
    )
    .reset_index()
)
hours = hours[hours.hours > 0]
inj = (
    acc.groupby(["mine_id", "year"])
    .agg(
        injuries=("degree_injury_cd", "size"), fatal=("degree_injury_cd", lambda s: (s == 1).sum())
    )
    .reset_index()
)
my = hours.merge(inj, on=["mine_id", "year"], how="left").fillna({"injuries": 0, "fatal": 0})
print(
    f"{len(my):,} mine-years with hours, {my.injuries.sum():,.0f} injuries, {my.fatal.sum():,.0f} fatalities, "
    f"{my.hours.sum() / 2000 / 1e6:.2f} million worker-years"
)

# %% ----------------------------------------------------------- national rates
banner("National injury and fatality rates per 200,000 hours")
nat = (
    my.groupby(["year", "sector"])
    .agg(injuries=("injuries", "sum"), fatal=("fatal", "sum"), hours=("hours", "sum"))
    .reset_index()
)
nat["rate"] = nat.injuries * 2e5 / nat.hours
nat["fatal_rate"] = nat.fatal * 2e5 / nat.hours
piv = nat.pivot(index="year", columns="sector", values="rate").rename(
    columns={"C": "coal", "M": "metal/nonmetal"}
)
print(piv.round(2).iloc[[0, 5, 10, 15, 20, -1]].to_string())
# Poisson GLM trend with exposure offset, per sector
for sec, g in nat.groupby("sector"):
    X = sm.add_constant(g.year - 2000)
    fit = sm.GLM(g.injuries, X, family=sm.families.Poisson(), offset=np.log(g.hours)).fit()
    nb = sm.GLM(
        g.injuries, X, family=sm.families.NegativeBinomial(alpha=1e-3), offset=np.log(g.hours)
    ).fit()
    rr = np.exp(fit.params.iloc[1])
    lo, hi = np.exp(fit.conf_int().iloc[1])
    print(
        f"{'coal' if sec == 'C' else 'metal/nonmetal':15s} injury rate ratio per year={rr:.3f} ({lo:.3f}-{hi:.3f}); "
        f"deviance/df={fit.deviance / fit.df_resid:.1f} (Poisson, should be ~1)"
    )
note(f"""Injury rates fell by roughly {(1 - piv.iloc[-1] / piv.iloc[0]).mean() * 100:.0f}% over the period in both
sectors. The national Poisson GLM reports a deviance/df far above 1: even aggregated to the
national level the year-to-year counts are more variable than Poisson, so the printed CI is
too narrow. Use a quasi-Poisson scale or negative binomial for any inference about trend.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for sec, lab in [("C", "coal"), ("M", "metal/nonmetal")]:
    g = nat[nat.sector == sec]
    axes[0].plot(g.year, g.rate, "o-", ms=3, label=lab)
    axes[1].plot(g.year, g.fatal_rate, "o-", ms=3, label=lab)
    # Poisson exact CI on the fatality rate
    lo = sps.chi2.ppf(0.025, 2 * g.fatal) / 2 * 2e5 / g.hours
    hi = sps.chi2.ppf(0.975, 2 * (g.fatal + 1)) / 2 * 2e5 / g.hours
    axes[1].fill_between(g.year, lo, hi, alpha=0.2)
axes[0].set_ylabel("injuries per 200,000 h")
axes[0].set_title("Reportable injury rate, US mines")
axes[0].legend()
axes[1].set_ylabel("fatalities per 200,000 h")
axes[1].set_title("Fatality rate with exact Poisson 95% CI (small counts)")
axes[1].legend()
fig.tight_layout()
savefig(fig, "national_rates", SEC)

# %% ----------------------------------------------------------- overdispersion across mines
banner("Across mines: Poisson vs negative binomial")
big = my[(my.year == 2019) & (my.hours >= 50000)].copy()
X = sm.add_constant((big.sector == "C").astype(float).rename("coal"))
pois = sm.GLM(big.injuries, X, family=sm.families.Poisson(), offset=np.log(big.hours)).fit()
# estimate NB alpha by maximum likelihood via statsmodels' discrete NB
nbm = sm.NegativeBinomial(big.injuries, X, offset=np.log(big.hours)).fit(disp=0)
print(
    f"2019, mines with >= 50,000 h: n={len(big)}; Poisson deviance/df={pois.deviance / pois.df_resid:.2f}; "
    f"NB alpha={nbm.params['alpha']:.2f}"
)
print(
    f"coal vs metal/nonmetal rate ratio: Poisson {np.exp(pois.params['coal']):.2f} "
    f"(SE {pois.bse['coal']:.3f}) | NB {np.exp(nbm.params['coal']):.2f} (SE {nbm.bse['coal']:.3f})"
)
note(f"""The rate ratio is the same under both models but the negative-binomial standard error is
{nbm.bse["coal"] / pois.bse["coal"]:.1f}x the Poisson one. With alpha={nbm.params["alpha"]:.2f}, a mine with an
expected 10 injuries has a variance of 10 + alpha x 100: mines differ in underlying risk
(commodity, method, management), and pretending they are identical Poisson draws produces
spuriously precise comparisons.""")

# %% ----------------------------------------------------------- small numbers & shrinkage
banner("Ranking mines by rate: the small-numbers problem and empirical-Bayes shrinkage")
recent = (
    my[my.year.between(2020, 2024)]
    .groupby("mine_id")
    .agg(injuries=("injuries", "sum"), hours=("hours", "sum"), sector=("sector", "first"))
    .reset_index()
)
recent = recent[recent.hours > 2000]
recent["rate"] = recent.injuries * 2e5 / recent.hours
# Gamma-Poisson empirical Bayes: injuries_i ~ Poisson(lambda_i * E_i), lambda_i ~ Gamma(a, b)
E = recent.hours.values / 2e5
k_obs = recent.injuries.values
# marginal likelihood of a Gamma(a, b) prior on the rate: k_i ~ NegBin(a, b / (b + E_i))


def nll(p):
    a, b = np.exp(p)
    return -np.sum(sps.nbinom.logpmf(k_obs, a, b / (b + E)))


opt = minimize(nll, [0.0, -1.0], method="Nelder-Mead")
a, b = np.exp(opt.x)
m, v = a / b, a / b**2
recent["rate_eb"] = (recent.injuries + a) / (E + b)
top_raw = recent.sort_values("rate", ascending=False).head(5)
top_eb = recent.sort_values("rate_eb", ascending=False).head(5)
print(
    f"Gamma-Poisson prior by maximum likelihood: mean rate {m:.2f}, sd {np.sqrt(v):.2f} per 200k h "
    f"(a={a:.2f}, b={b:.2f}); a mine's rate is shrunk toward {m:.2f} with weight b/(b+E) = {b / (b + np.median(E)):.2f} at median exposure"
)
print(
    "\nworst 5 by raw rate (2020-24):\n"
    + top_raw[["mine_id", "injuries", "hours", "rate", "rate_eb"]].round(1).to_string(index=False)
)
print(
    "\nworst 5 by shrunken rate:\n"
    + top_eb[["mine_id", "injuries", "hours", "rate", "rate_eb"]].round(1).to_string(index=False)
)
frac_zero_small = (recent[recent.hours < 20000].injuries == 0).mean()
print(
    f"\n{frac_zero_small * 100:.0f}% of mines with < 20,000 h over five years recorded zero injuries (rate 0.0)"
)
note("""The raw league table is a list of tiny operations with one or two injuries over a few
thousand hours. Shrinkage pulls them toward the prior in proportion to how little exposure
they have, and the shrunken ranking surfaces large mines with persistently high rates, which
is what an inspector or an insurer actually wants. Zero-injury small mines are not safe;
they are uninformative.""")

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(recent.hours / 2e5, recent.rate, s=5, alpha=0.3, label="raw rate")
ax.scatter(recent.hours / 2e5, recent.rate_eb, s=5, alpha=0.3, label="empirical-Bayes rate")
ax.axhline(m, color="k", lw=0.8)
ax.set_xscale("log")
ax.set_yscale("symlog", linthresh=1)
ax.set_ylim(-0.3, None)
ax.set_xlabel("exposure (units of 200,000 h, 2020-24)")
ax.set_ylabel("injuries per 200,000 h")
ax.set_title("Funnel: raw rates explode at low exposure; shrinkage tames them")
ax.legend()
savefig(fig, "shrinkage_funnel", SEC)

# %% ----------------------------------------------------------- u-chart for one mine
banner("u-chart: one large mine's annual injury rate against Poisson limits")
big_mine = my.groupby("mine_id").hours.sum().idxmax()
g = my[my.mine_id == big_mine].sort_values("year")
name = mines.set_index("mine_id").current_mine_name.get(big_mine, big_mine)
u, lo, hi = u_chart_limits(g.injuries, g.hours, scale=2e5)
g = g.assign(rate=g.injuries * 2e5 / g.hours, lo=lo, hi=hi)
out = g[(g.rate > g.hi) | (g.rate < g.lo)]
print(
    f"mine {big_mine} ({name}): centre line {u:.2f} per 200k h; {len(out)} of {len(g)} years outside Poisson limits: {out.year.tolist()}"
)
mr_sigma = np.abs(np.diff(g.rate)).mean() / 1.128
print(
    f"Poisson sigma at median exposure={np.sqrt(u / (g.hours.median() / 2e5)):.2f}; observed year-to-year sigma (I-chart)={mr_sigma:.2f}"
)
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(g.year, g.rate, "o-", ms=4)
ax.step(g.year, g.hi, where="mid", color="C3", ls="--", lw=0.8)
ax.step(g.year, g.lo, where="mid", color="C3", ls="--", lw=0.8)
ax.axhline(u, color="k", lw=0.8)
ax.set_ylabel("injuries per 200,000 h")
ax.set_title(f"u-chart, {name[:40]} (limits widen when hours fall)")
savefig(fig, "u_chart", SEC)
note(f"""With {len(out)} of {len(g)} years outside the limits, this chart is not flagging {len(out)} special
causes: at tens of millions of hours per year the Poisson limits are so tight (sigma
{np.sqrt(u / (g.hours.median() / 2e5)):.2f}) that ordinary year-to-year variation ({mr_sigma:.2f} by the moving range)
breaches them. Rates at a large mine are overdispersed relative to Poisson (changing workforce,
reporting practice, operations), so use I-chart limits from the observed variability or a
negative-binomial u-chart. The u-chart with Poisson limits is right for small counts, wrong
for huge ones.""")

# %% ----------------------------------------------------------- severity
banner("Severity: lost days are heavy-tailed; experience at injury")
dl = acc.days_lost.dropna()
dl = dl[dl > 0]
print(
    f"days lost (where > 0): n={len(dl):,}  mean={dl.mean():.1f}  median={dl.median():.0f}  P90={dl.quantile(0.9):.0f}  "
    f"P99={dl.quantile(0.99):.0f}  max={dl.max():.0f}"
)
print(
    f"share of all lost days from the worst 10% of cases: {dl[dl >= dl.quantile(0.9)].sum() / dl.sum():.0%}"
)
exp_bins = pd.cut(acc.job_exper, [0, 1, 2, 5, 10, 20, 60], right=False)
print(
    "\ninjuries by years of experience in the job:\n"
    + (exp_bins.value_counts(normalize=True).sort_index() * 100).round(1).to_string()
)
note("""The mean lost-days figure is dominated by a few catastrophic cases; report the median
and the tail. The experience table shows where injuries happen, but without hours worked by
experience band it is not a rate: new workers may be over-represented simply because there
are many of them. Always ask what the denominator is.""")
