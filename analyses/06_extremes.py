# %% [markdown]
# # 06 Extreme value analysis
#
# Data: Potomac annual peak daily flow (floods); Central Park daily precipitation (POT);
# Central Park annual maximum TMIN (nonstationary extremes).
#
# **Cheat sheet**
# - Block maxima -> GEV (Fisher-Tippett). Shape xi: xi>0 Frechet (heavy tail, floods, rain),
#   xi=0 Gumbel, xi<0 Weibull (bounded, temperature). scipy genextreme uses c = -xi.
# - Return level z_T: exceeded on average once per T years, i.e. annual exceedance prob 1/T.
#   100-yr event has a 1% chance each year and ~63% chance over 100 years.
# - Peaks-over-threshold -> Generalised Pareto (GPD). Uses more data than annual maxima. Pick
#   the threshold where the mean-residual-life plot is linear and GPD parameters stabilise.
#   Decluster first so exceedances are independent events, not one storm counted 4 times.
# - Return levels from POT: z_T = u + sigma/xi [ (lambda T)^xi - 1 ] with lambda = events/yr.
# - Uncertainty: the shape parameter is poorly constrained with <50 years; bootstrap or
#   profile likelihood. 100-yr return level CIs are typically +-30% or worse.
# - US flood practice: Bulletin 17C uses log-Pearson III on annual peaks. Shown for comparison.
# - Nonstationarity: let location (and/or scale) depend on time or a covariate; compare with
#   a likelihood-ratio test. "100-yr event" is then a moving target.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps
from scipy.optimize import minimize

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import fit_gev, plotting_positions, return_level_gev

setup()
SEC = "06_extremes"
T = np.array([2, 5, 10, 25, 50, 100, 200, 500])

# %% ----------------------------------------------------------- annual maxima: Potomac
banner("Potomac annual peak daily flow: GEV vs Gumbel vs log-Pearson III")
q = data.usgs_daily().flow_cfs
wy = q.index.year + (q.index.month >= 10)
amax = q.groupby(wy).max().loc[1931:2025].dropna() / 1000  # kcfs
print(
    f"n={len(amax)} water years; median peak {amax.median():.0f} kcfs; record "
    f"{amax.max():.0f} kcfs in WY{amax.idxmax()}"
)

c, loc, scale = fit_gev(amax)
xi = -c
g_loc, g_scale = sps.gumbel_r.fit(amax)
ll_gev = sps.genextreme.logpdf(amax, c, loc, scale).sum()
ll_gum = sps.gumbel_r.logpdf(amax, g_loc, g_scale).sum()
lrt = 2 * (ll_gev - ll_gum)
p_lrt = sps.chi2.sf(lrt, df=1)
print(
    f"GEV: xi={xi:+.3f} loc={loc:.1f} scale={scale:.1f}   Gumbel: loc={g_loc:.1f} scale={g_scale:.1f}"
)
print(f"AIC GEV={-2 * ll_gev + 6:.1f}  AIC Gumbel={-2 * ll_gum + 4:.1f}  LRT p={p_lrt:.3f}")

# log-Pearson III (Bulletin 17C flavour, method of moments on log10 Q)
lq = np.log10(amax)
skew, mu, sd = sps.skew(lq, bias=False), lq.mean(), lq.std(ddof=1)
lp3_levels = 10 ** sps.pearson3.ppf(1 - 1 / T, skew, loc=mu, scale=sd)

levels = pd.DataFrame(
    {
        "T (yr)": T,
        "GEV": return_level_gev(c, loc, scale, T),
        "Gumbel": sps.gumbel_r.ppf(1 - 1 / T, g_loc, g_scale),
        "LP3": lp3_levels,
    }
).set_index("T (yr)")

# parametric bootstrap CI on GEV return levels
rng = np.random.default_rng(0)
boot, boot_xi, failed = [], [], 0
while len(boot) < 1000:
    sample = sps.genextreme.rvs(c, loc, scale, size=len(amax), random_state=rng)
    fit = fit_gev(sample)
    if fit is None:
        failed += 1
        continue
    boot.append(return_level_gev(*fit, T))
    boot_xi.append(-fit[0])
boot = np.array(boot)
levels["GEV 2.5%"] = np.percentile(boot, 2.5, axis=0)
levels["GEV 97.5%"] = np.percentile(boot, 97.5, axis=0)
xi_lo, xi_hi = np.percentile(boot_xi, [2.5, 97.5])
print(
    f"parametric bootstrap: 1000 fits ({failed} discarded as implausible); "
    f"xi 95% CI = ({xi_lo:+.2f}, {xi_hi:+.2f})"
)
print("\nReturn levels (kcfs):\n" + levels.round(0).to_string())
half_width = (
    (levels.loc[100, "GEV 97.5%"] - levels.loc[100, "GEV 2.5%"]) / 2 / levels.loc[100, "GEV"] * 100
)
note(f"""The shape xi={xi:+.2f} is positive (heavy tail), the likelihood-ratio test prefers GEV
over Gumbel (p={p_lrt:.3f}), yet the bootstrap CI on xi runs {xi_lo:+.2f} to {xi_hi:+.2f}:
{len(amax)} years constrain the tail only loosely. The three models agree to ~10 years and
diverge beyond 100, where the GEV 100-year level has a roughly +-{half_width:.0f}% interval.
Extrapolating to 500 years from {len(amax)} years is a modelling assumption, not a measurement.
Say so in the report.""")

x_sorted, T_emp = plotting_positions(amax, a=0.44)  # Gringorten
fig, ax = plt.subplots(figsize=(9, 5))
Tg = np.logspace(0.01, 2.8, 200)
ax.plot(Tg, return_level_gev(c, loc, scale, Tg), label=f"GEV (xi={xi:+.2f})")
ax.plot(Tg, sps.gumbel_r.ppf(1 - 1 / Tg, g_loc, g_scale), ls="--", label="Gumbel")
ax.plot(
    Tg, 10 ** sps.pearson3.ppf(1 - 1 / Tg, skew, loc=mu, scale=sd), ls=":", label="log-Pearson III"
)
ax.fill_between(T, levels["GEV 2.5%"], levels["GEV 97.5%"], alpha=0.2, label="GEV 95% bootstrap")
ax.scatter(T_emp, x_sorted, s=14, color="k", zorder=3, label="observed (Gringorten positions)")
ax.set_xscale("log")
ax.set_xlabel("return period (years)")
ax.set_ylabel("annual peak daily flow (kcfs)")
ax.set_title("Potomac at Little Falls: flood frequency curve")
ax.legend()
savefig(fig, "flood_frequency", SEC)

# %% ----------------------------------------------------------- POT / GPD: precipitation
banner("Central Park daily precipitation: peaks over threshold")
pr = data.ghcn_daily().PRCP.dropna()
pr = pr[pr.index.year <= 2025]
years = pr.index.year.nunique()

# mean residual life plot
us = np.arange(10, 80, 2.0)
mrl = [(pr[pr > u] - u).mean() for u in us]
mrl_se = [(pr[pr > u] - u).std() / np.sqrt((pr > u).sum()) for u in us]
# parameter stability: fit GPD at each threshold
shapes, scales_mod = [], []
for u in us:
    exc = pr[pr > u] - u
    cu, _, su = sps.genpareto.fit(exc, floc=0)
    shapes.append(cu)
    scales_mod.append(su - cu * u)  # modified scale should be constant above a valid u
fig, axes = plt.subplots(1, 3, figsize=(14, 4))
axes[0].errorbar(us, mrl, yerr=np.array(mrl_se) * 1.96, fmt="o-", ms=3)
axes[0].set_xlabel("threshold u (mm)")
axes[0].set_ylabel("mean excess")
axes[0].set_title("Mean residual life (want linear)")
axes[1].plot(us, shapes, "o-", ms=3)
axes[1].set_title("GPD shape vs threshold (want flat)")
axes[1].set_xlabel("u (mm)")
axes[2].plot(us, scales_mod, "o-", ms=3)
axes[2].set_title("modified scale vs threshold (want flat)")
axes[2].set_xlabel("u (mm)")
savefig(fig, "threshold_selection", SEC)

u = 30.0
# decluster: an exceedance starts a new cluster only if >= 2 dry-of-exceedance days separate it
exceed = pr[pr > u]
gaps = exceed.index.to_series().diff().dt.days.fillna(999)
cluster_id = (gaps >= 3).cumsum()
peaks = exceed.groupby(cluster_id).max()
lam = len(peaks) / years
cg, _, sg = sps.genpareto.fit(peaks - u, floc=0)
xi_g = cg
print(
    f"threshold u={u:.0f} mm: {len(exceed)} exceedance days -> {len(peaks)} declustered events "
    f"({lam:.2f}/yr); GPD xi={xi_g:+.3f} sigma={sg:.1f}"
)
rl_pot = u + sg / xi_g * ((lam * T) ** xi_g - 1)
# compare with annual-maxima GEV on the same data
pr_amax = pr.groupby(pr.index.year).max()
ca, la, sa = fit_gev(pr_amax)
rl_am = return_level_gev(ca, la, sa, T)
cmp = pd.DataFrame({"T": T, "POT/GPD": rl_pot, "AnnMax/GEV": rl_am}).set_index("T")
print("\nDaily precipitation return levels (mm):\n" + cmp.round(0).to_string())
print(f"GEV on annual maxima: xi={-ca:+.3f} (n={len(pr_amax)})")
note(f"""POT uses {len(peaks)} events instead of {len(pr_amax)} annual maxima, so the tail is
better constrained. Both approaches put the 100-year daily rainfall near
{rl_pot[T == 100][0]:.0f}-{rl_am[T == 100][0]:.0f} mm. Declustering matters: without it, a 3-day
storm counts as 3 'independent' events and lambda is inflated. The positive xi is the
usual heavy tail of precipitation.""")

# %% ----------------------------------------------------------- nonstationary GEV
banner("Nonstationary GEV: Central Park annual maximum TMIN with a linear trend in location")
tmin = data.ghcn_daily().TMIN.dropna()
tmin = tmin[tmin.index.year <= 2025]
am = tmin.groupby(tmin.index.year).max()
t = (am.index.values - 1950) / 10.0  # decades since 1950
y = am.values


def nll_stationary(p):
    mu, lsig, xi = p
    return -sps.genextreme.logpdf(y, -xi, loc=mu, scale=np.exp(lsig)).sum()


def nll_trend(p):
    mu0, mu1, lsig, xi = p
    return -sps.genextreme.logpdf(y, -xi, loc=mu0 + mu1 * t, scale=np.exp(lsig)).sum()


p0 = [y.mean(), np.log(y.std()), -0.1]
fs = minimize(
    nll_stationary, p0, method="Nelder-Mead", options=dict(maxiter=5000, xatol=1e-6, fatol=1e-8)
)
ft = minimize(
    nll_trend,
    [fs.x[0], 0.0, fs.x[1], fs.x[2]],
    method="Nelder-Mead",
    options=dict(maxiter=10000, xatol=1e-6, fatol=1e-8),
)
lrt = 2 * (fs.fun - ft.fun)
p_lrt = sps.chi2.sf(lrt, 1)
mu0, mu1, lsig, xi_t = ft.x
print(f"stationary: mu={fs.x[0]:.2f} sigma={np.exp(fs.x[1]):.2f} xi={fs.x[2]:+.2f}")
print(f"trend:      mu0={mu0:.2f} mu1={mu1:+.3f} C/decade sigma={np.exp(lsig):.2f} xi={xi_t:+.2f}")
print(f"likelihood-ratio stat={lrt:.2f}  p={p_lrt:.4f}")
rl20_1950 = sps.genextreme.ppf(1 - 1 / 20, -xi_t, loc=mu0, scale=np.exp(lsig))
rl20_2025 = sps.genextreme.ppf(1 - 1 / 20, -xi_t, loc=mu0 + mu1 * t[-1], scale=np.exp(lsig))
p_now = 1 - sps.genextreme.cdf(rl20_1950, -xi_t, loc=mu0 + mu1 * t[-1], scale=np.exp(lsig))
print(
    f"20-yr return level: {rl20_1950:.1f} C in 1950 -> {rl20_2025:.1f} C in 2025; the 1950 "
    f"20-yr night is now a 1-in-{1 / p_now:.0f}-year event"
)
note(f"""Fitting the GEV by writing the negative log-likelihood and calling an optimiser is
worth knowing: it generalises to any covariate. The negative shape (xi={xi_t:+.2f}) is typical
for temperature: a bounded upper tail. Letting location drift with time is strongly
preferred by the likelihood-ratio test; the warmest night of the year has shifted by
{mu1 * t[-1]:.1f} C since 1950, which turns a 1950s 20-year event into a roughly
1-in-{1 / p_now:.0f}-year event today.""")

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.scatter(am.index, y, s=14, color="k", label="annual max TMIN")
yrs = am.index.values
for Tq, ls in [(2, "-"), (20, "--"), (100, ":")]:
    ax.plot(
        yrs,
        sps.genextreme.ppf(1 - 1 / Tq, -xi_t, loc=mu0 + mu1 * t, scale=np.exp(lsig)),
        ls=ls,
        color="C3",
        label=f"{Tq}-yr level (nonstationary)",
    )
ax.set_ylabel("deg C")
ax.set_title("Nonstationary GEV: return levels drift upward with the location parameter")
ax.legend(ncol=2)
savefig(fig, "nonstationary_gev", SEC)
