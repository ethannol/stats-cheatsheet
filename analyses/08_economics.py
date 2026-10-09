# %% [markdown]
# # 08 Commodity prices and project economics under uncertainty
#
# Data: World Bank monthly metal prices 1960-2024 (real); a SYNTHETIC copper project model.
#
# **Cheat sheet**
# - Work with log returns r = ln(P_t/P_{t-1}): additive, roughly symmetric, annualise volatility
#   with sqrt(12) for monthly data. Price levels are non-stationary; returns are not.
# - GBM (random walk with drift) vs mean reversion (Ornstein-Uhlenbeck): test the log price for
#   a unit root (ADF); estimate the half-life from the AR(1) coefficient. Mining analysts often
#   assume long-run reversion to incentive price; GBM gives far wider long-horizon fans.
# - Monte Carlo NPV: specify distributions for grade, recovery, price, costs, capex, schedule;
#   impose correlations (price-cost via energy; grade-recovery) with a Cholesky factor; report
#   the NPV distribution, P(NPV < 0), and expected value vs the deterministic base case.
# - Tornado (one-at-a-time) sensitivity ranks drivers; it ignores interactions and correlations,
#   which the Monte Carlo includes.
# - Break-even cut-off grade g_c = cost per tonne ore / (price x recovery x payability): the link
#   between economics and the grade-tonnage curve of section 04.
# - Real options: flexibility (defer, expand, abandon) has value that static NPV ignores; a
#   simple simulation of 'wait and decide' makes the point without Black-Scholes.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from miningstats import data
from miningstats.plotting import banner, note, savefig, setup
from miningstats.stats import correlated_normals, fit_ou, gbm_paths

setup()
SEC = "08_economics"
px = data.commodity_prices()["1990":]
METALS = ["copper", "gold", "iron_ore", "nickel", "zinc"]

# %% ----------------------------------------------------------- returns & volatility
banner("Log returns, volatility, correlation (monthly, 1990-2024)")
lr = np.log(px[METALS]).diff().dropna()
vol = lr.std() * np.sqrt(12)
print("annualised volatility:\n" + (vol * 100).round(1).astype(str).add("%").to_string())
print("\nreturn correlation:\n" + lr.corr().round(2).to_string())
print(f"\ncopper monthly return skew={lr.copper.skew():.2f}, excess kurtosis={lr.copper.kurt():.2f} "
      f"(normal = 0; fat tails)")
note(f"""Base metals run at {vol[['copper', 'nickel', 'zinc']].min() * 100:.0f}-{vol[['copper', 'nickel', 'zinc']].max() * 100:.0f}%
annualised volatility, gold at {vol['gold'] * 100:.0f}%. Returns are fat-tailed, so a normal
model understates the frequency of crashes; copper and zinc move together (macro demand)
while gold is nearly uncorrelated with the industrial metals, which is why a diversified
miner's cash-flow volatility is lower than any single commodity's.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
(px[METALS] / px[METALS].iloc[0]).plot(ax=axes[0], lw=1, logy=True)
axes[0].set_title("Nominal prices, indexed to Jan 1990 (log scale)")
roll = lr.copper.rolling(36).std() * np.sqrt(12)
axes[1].plot(roll.index, roll * 100)
axes[1].set_ylabel("36-month rolling vol (%)")
axes[1].set_title("Copper volatility is itself volatile (clustering)")
fig.tight_layout()
savefig(fig, "prices_volatility", SEC)

# %% ----------------------------------------------------------- GBM vs mean reversion
banner("Random walk or mean reversion? Copper and gold")
for m in ["copper", "gold"]:
    ou = fit_ou(np.log(px[m].dropna()))
    if np.isfinite(ou["half_life_years"]):
        print(f"{m:8s} ADF p (unit root)={ou['p_unit_root']:.3f}  AR(1) b={ou['b']:.4f}  "
              f"half-life={ou['half_life_years']:.1f} yr  long-run mean={np.exp(ou['mean']):,.0f}")
    else:
        print(f"{m:8s} ADF p (unit root)={ou['p_unit_root']:.3f}  AR(1) b={ou['b']:.4f} >= 1: no mean reversion "
              f"estimable (pure random walk or trending)")
note("""Neither series rejects a unit root at conventional levels over 35 years of monthly data: the
power to distinguish slow mean reversion (half-life of several years) from a random walk is
low. Yet the two assumptions give very different 10-year price fans and NPV distributions, so
state the assumption, show both, and justify reversion with the cost curve (incentive price)
rather than with a statistical test that cannot settle it.""")

cu = px.copper.dropna()
sigma_cu = lr.copper.std()
rng = np.random.default_rng(0)
H = 120
gbm = gbm_paths(cu.iloc[-1], 0.0, sigma_cu * np.sqrt(12), H, 2000, dt=1 / 12, rng=rng)
ou = fit_ou(np.log(cu))
b_m = ou["b"] if np.isfinite(ou["b"]) else 0.99
mean_m = ou["mean"]
sig_e = np.log(cu).diff().dropna().std()
ou_paths = np.empty((2000, H))
x = np.full(2000, np.log(cu.iloc[-1]))
for t in range(H):
    x = mean_m * (1 - b_m) + b_m * x + rng.normal(0, sig_e, 2000)
    ou_paths[:, t] = np.exp(x)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
fut = pd.date_range(cu.index[-1], periods=H + 1, freq="MS")[1:]
for ax, paths, title in [(axes[0], gbm, "GBM (random walk, zero drift)"), (axes[1], ou_paths, "Mean reversion fitted by AR(1)")]:
    for q, a in [(0.05, .15), (0.25, .3)]:
        ax.fill_between(fut, np.quantile(paths, q, axis=0), np.quantile(paths, 1 - q, axis=0), alpha=a, color="C0")
    ax.plot(fut, np.median(paths, axis=0), color="C3", label="median")
    ax.plot(cu["2015":].index, cu["2015":], color="k", lw=1, label="history")
    ax.set_title(f"{title}: 10-yr fan (50% / 90%)")
    ax.legend()
axes[0].set_ylabel("copper $/t")
fig.tight_layout()
savefig(fig, "price_fans", SEC)
print(f"10-year 90% band width: GBM {np.quantile(gbm[:, -1], .95) / np.quantile(gbm[:, -1], .05):.1f}x, "
      f"mean reversion {np.quantile(ou_paths[:, -1], .95) / np.quantile(ou_paths[:, -1], .05):.1f}x")

# %% ----------------------------------------------------------- Monte Carlo NPV
banner("SYNTHETIC copper project: deterministic vs Monte Carlo NPV")
base = dict(tonnes_per_year=20e6, life=12, grade_pct=0.55, recovery=0.88, price=cu.iloc[-1],
            opex_per_t=22.0, capex=1.8e9, discount=0.08, payability=0.96, capex_years=2)


def npv(p):
    t = np.arange(1, p["life"] + 1)
    metal_t = p["tonnes_per_year"] * p["grade_pct"] / 100 * p["recovery"]
    rev = metal_t * p["price"] * p["payability"]
    cf = (rev - p["tonnes_per_year"] * p["opex_per_t"]) / (1 + p["discount"]) ** (t + p["capex_years"])
    capex = p["capex"] / p["capex_years"] * np.sum(1 / (1 + p["discount"]) ** np.arange(1, p["capex_years"] + 1))
    return cf.sum() - capex


npv0 = npv(base)
print(f"base case NPV = ${npv0 / 1e9:.2f} billion at copper {base['price']:,.0f} $/t")
# tornado
rows = []
for k in ["grade_pct", "recovery", "price", "opex_per_t", "capex", "discount", "tonnes_per_year"]:
    lo, hi = dict(base), dict(base)
    lo[k], hi[k] = base[k] * 0.9, base[k] * 1.1
    rows.append(dict(input=k, low=npv(lo) / 1e9, high=npv(hi) / 1e9))
tor = pd.DataFrame(rows).set_index("input")
tor["swing"] = (tor.high - tor.low).abs()
tor = tor.sort_values("swing")
print("\ntornado (+-10% one at a time, $bn):\n" + tor.round(2).to_string())
# Monte Carlo with correlated inputs
N = 20000
names = ["grade", "recovery", "price", "opex", "capex"]
corr = np.array([[1, 0.4, 0, 0, 0],
                 [0.4, 1, 0, 0, 0],
                 [0, 0, 1, 0.5, 0.3],   # price-opex via energy/inputs; price-capex via cycle
                 [0, 0, 0.5, 1, 0.4],
                 [0, 0, 0.3, 0.4, 1]], float)
z = correlated_normals(N, corr, rng=1)
sims = []
for i in range(N):
    p = dict(base)
    p["grade_pct"] = base["grade_pct"] * np.exp(0.10 * z[i, 0] - 0.005)            # ~10% CV, lognormal
    p["recovery"] = np.clip(base["recovery"] + 0.03 * z[i, 1], 0.6, 0.97)
    p["price"] = base["price"] * np.exp(sigma_cu * np.sqrt(12) * np.sqrt(6) * z[i, 2] * 0.5)  # life-average price
    p["opex_per_t"] = base["opex_per_t"] * np.exp(0.15 * z[i, 3])
    p["capex"] = base["capex"] * np.exp(0.20 * z[i, 4] + 0.05)                     # capex overruns skew up
    sims.append(npv(p))
sims = np.array(sims) / 1e9
print(f"\nMonte Carlo ({N:,} runs): mean NPV=${sims.mean():.2f}bn, median=${np.median(sims):.2f}bn, "
      f"P10=${np.quantile(sims, .1):.2f}bn, P90=${np.quantile(sims, .9):.2f}bn, P(NPV<0)={np.mean(sims < 0):.1%}")
note(f"""The tornado says price and grade dominate; costs and capex matter less because they are
spread over the life or discounted. The Monte Carlo adds what the tornado cannot: joint
uncertainty and correlation (high prices come with high costs, so the upside is clipped),
and asymmetric inputs (capex overruns). The result is a mean NPV of ${sims.mean():.2f}bn against a
deterministic ${npv0 / 1e9:.2f}bn and a {np.mean(sims < 0):.0%} chance of destroying value, which
is the number the board needs.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
axes[0].barh(tor.index, tor.high - npv0 / 1e9, left=npv0 / 1e9, color="C2", label="+10%")
axes[0].barh(tor.index, tor.low - npv0 / 1e9, left=npv0 / 1e9, color="C3", label="-10%")
axes[0].axvline(npv0 / 1e9, color="k", lw=.8)
axes[0].set_xlabel("NPV ($bn)")
axes[0].set_title("Tornado: one-at-a-time +-10%")
axes[0].legend()
axes[1].hist(sims, bins=80, alpha=.8)
axes[1].axvline(0, color="k", lw=.8)
axes[1].axvline(npv0 / 1e9, color="C3", ls="--", label="deterministic base case")
axes[1].axvline(sims.mean(), color="C1", ls="-", label="Monte Carlo mean")
axes[1].set_xlabel("NPV ($bn)")
axes[1].set_title(f"NPV distribution, P(NPV<0) = {np.mean(sims < 0):.0%}")
axes[1].legend()
fig.tight_layout()
savefig(fig, "npv", SEC)

# %% ----------------------------------------------------------- cut-off grade
banner("Break-even cut-off grade")
gc = base["opex_per_t"] / (base["price"] * base["recovery"] * base["payability"]) * 100
print(f"g_c = opex / (price x recovery x payability) = {gc:.3f}% Cu at {base['price']:,.0f} $/t; "
      f"at -30% price: {gc / 0.7:.3f}%")
note("""The cut-off grade converts the grade-tonnage curve (section 04) into ore tonnes and
metal, and it moves with price: a price fall raises the cut-off, shrinks the ore inventory
and lowers the mine life, which is why resource statements quote the price assumption.""")

# %% ----------------------------------------------------------- option to defer
banner("Real-option flavour: value of waiting two years for price information")
# decide today vs observe the (simulated) price in 2 years and only build if NPV > 0 then
price_2y = gbm_paths(cu.iloc[-1], 0.0, sigma_cu * np.sqrt(12), 24, 20000, dt=1 / 12, rng=rng)[:, -1]
npv_now_dist = np.array([npv({**base, "price": p}) for p in price_2y]) / 1e9   # NPV as it would turn out
npv_wait = np.maximum(npv_now_dist, 0) / (1 + base["discount"]) ** 2
print(f"commit now: expected NPV ${npv_now_dist.mean():.2f}bn (P(<0)={np.mean(npv_now_dist < 0):.0%}); "
      f"wait 2 years and build only if positive: ${npv_wait.mean():.2f}bn -> option value ${npv_wait.mean() - npv_now_dist.mean():+.2f}bn")
note("""Waiting trades two years of discounting against the right not to build into a price
collapse. When volatility is high and the project is marginal the option is worth a lot;
for a deeply profitable project it is worth little and delay is pure cost. This is the
intuition behind why miners sequence projects and hold undeveloped resources.""")
