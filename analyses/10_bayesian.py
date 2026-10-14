# %% [markdown]
# # 10 Bayesian inference: conjugate updates, grid approximation, Metropolis-Hastings
#
# **Cheat sheet**
# - posterior ∝ likelihood x prior. The prior is a modelling choice to defend, like a link
#   function; weakly informative priors are the default, flat priors are not 'objective'.
# - Credible interval: P(theta in interval | data) = 95%. Confidence interval: 95% of
#   intervals built this way would contain theta. People want the first and compute the second.
# - Conjugate pairs give closed forms: Beta-Binomial (proportions), Gamma-Poisson (rates),
#   Normal-Normal (means with known variance), Normal-Inverse-Gamma (regression).
# - MCMC when there is no closed form. Metropolis-Hastings: propose, accept with probability
#   min(1, p(new)/p(old)). Tune the step for ~25-45% acceptance. Run >= 4 chains from
#   dispersed starts. Diagnostics: trace plots, R-hat (< 1.01), effective sample size, burn-in.
# - Posterior predictive checks: simulate data from the posterior, compare to the observed.
# - With lots of data and weak priors, Bayesian and frequentist answers coincide numerically;
#   the Bayesian advantage is small n, hierarchical structure, and honest propagation of
#   uncertainty into derived quantities (e.g. return levels, ratios).

# %%
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup

setup()
SEC = "10_bayesian"
rng = np.random.default_rng(0)
cp = data.ghcn_daily()
cp = cp[cp.index.year <= 2025]

# %% ----------------------------------------------------------- Beta-Binomial
banner("Beta-Binomial: probability that a July day in NYC is wet (>= 1 mm)")
july = cp[cp.index.month == 7].PRCP.dropna()
wet = (july >= 1.0)
theta = np.linspace(0, 1, 1001)
priors = {"flat Beta(1,1)": (1, 1), "weak Beta(2,2)": (2, 2), "wrong-ish Beta(20,5) (expects 80% wet)": (20, 5)}
# sequential update: first 10 days, first season, all data
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
for ax, n_obs in zip(axes, [10, 93, len(wet)]):
    k = wet.iloc[:n_obs].sum()
    for (name, (a, b)), col in zip(priors.items(), ["C0", "C1", "C3"]):
        post = sps.beta(a + k, b + n_obs - k)
        ax.plot(theta, post.pdf(theta), color=col, label=f"{name}: mean {post.mean():.2f}")
        ax.plot(theta, sps.beta(a, b).pdf(theta), color=col, ls=":", lw=.8)
    ax.set_title(f"after n={n_obs} days ({k} wet)")
    ax.set_xlabel("P(wet day)")
    ax.legend(fontsize=7)
axes[0].set_ylabel("density (dotted = prior)")
savefig(fig, "beta_binomial", SEC)
k, n_ = wet.sum(), len(wet)
post = sps.beta(1 + k, 1 + n_ - k)
lo, hi = post.ppf([.025, .975])
freq_ci = sps.binomtest(int(k), n_).proportion_ci()
print(f"all data: {k}/{n_} wet; posterior mean={post.mean():.4f}  95% credible=({lo:.4f}, {hi:.4f})")
print(f"frequentist Clopper-Pearson 95% CI=({freq_ci.low:.4f}, {freq_ci.high:.4f})")
note("""With 10 observations the prior dominates and the three posteriors disagree; by a full
season they have nearly converged; with 2,300 days the prior is irrelevant and the credible
interval equals the frequentist CI to three decimals. The prior matters exactly when you
have little data, which is also when you most need to think about it.""")

# %% ----------------------------------------------------------- Gamma-Poisson
banner("Gamma-Poisson: hot days (TMAX >= 35 C) and warm nights (TMIN >= 23.9 C) per summer")
jja = cp[cp.index.month.isin([6, 7, 8])]
a0, b0 = 1.0, 0.5  # Gamma(shape, rate) prior, mean 2 events/summer, weak
for label, counts in [("hot days TMAX>=35C", (jja.TMAX >= 35).resample("YE").sum()),
                      ("warm nights TMIN>=23.9C", (jja.TMIN >= 23.9).resample("YE").sum())]:
    counts.index = counts.index.year
    e_counts, l_counts = counts.loc[1951:1980], counts.loc[1991:2020]
    post_e = sps.gamma(a0 + e_counts.sum(), scale=1 / (b0 + len(e_counts)))
    post_l = sps.gamma(a0 + l_counts.sum(), scale=1 / (b0 + len(l_counts)))
    samp_e, samp_l = post_e.rvs(50000, random_state=rng), post_l.rvs(50000, random_state=rng)
    ratio = samp_l / samp_e
    print(f"{label}")
    print(f"  early: {e_counts.sum():4d} events / 30 summers -> rate {post_e.mean():.2f} "
          f"({post_e.ppf(.025):.2f}, {post_e.ppf(.975):.2f})   var/mean={e_counts.var() / e_counts.mean():.1f}")
    print(f"  late:  {l_counts.sum():4d} events / 30 summers -> rate {post_l.mean():.2f} "
          f"({post_l.ppf(.025):.2f}, {post_l.ppf(.975):.2f})   var/mean={l_counts.var() / l_counts.mean():.1f}")
    print(f"  rate ratio late/early: median {np.median(ratio):.2f}, 95% credible "
          f"({np.percentile(ratio, 2.5):.2f}, {np.percentile(ratio, 97.5):.2f});  "
          f"P(ratio > 1) = {(ratio > 1).mean():.3f}")
note("""Posterior sampling makes derived quantities trivial: draw rates, divide, and the ratio's
credible interval and P(ratio > 1) fall out with no delta method. Hot days show no change
(consistent with the TMAX null in section 02); warm nights have become far more frequent.
The var/mean column is an overdispersion check: a Poisson has var/mean = 1, but heat waves
cluster events, so the variance is several times the mean and the Gamma-Poisson intervals
are too narrow. A negative-binomial likelihood (Gamma-mixed Poisson) is the honest fix.""")

# %% ----------------------------------------------------------- Metropolis-Hastings regression
banner("Metropolis-Hastings: Bayesian linear trend for GISTEMP 1975-2025 (annual)")
gt = data.gistemp().anomaly.resample("YE").mean()
gt = gt[(gt.index.year >= 1975) & (gt.index.year <= 2025)]
t = (gt.index.year - 2000) / 10.0  # decades, centred
y = gt.values
ols = sm.OLS(y, sm.add_constant(t)).fit()


def log_post(p, prior_slope_sd=1.0):
    a, b, log_s = p
    s = np.exp(log_s)
    ll = sps.norm.logpdf(y, a + b * t, s).sum()
    lp = (sps.norm.logpdf(a, 0, 5) + sps.norm.logpdf(b, 0, prior_slope_sd)  # weakly informative
          + sps.norm.logpdf(log_s, 0, 2))
    return ll + lp


def metropolis(log_post_fn, start, step, n_iter, rng):
    chain = np.empty((n_iter, len(start)))
    cur, cur_lp = np.array(start, float), log_post_fn(start)
    accepted = 0
    for i in range(n_iter):
        prop = cur + rng.normal(0, step)
        lp = log_post_fn(prop)
        if np.log(rng.uniform()) < lp - cur_lp:
            cur, cur_lp = prop, lp
            accepted += 1
        chain[i] = cur
    return chain, accepted / n_iter


def rhat(chains):
    """Gelman-Rubin potential scale reduction (split chains not used; chains already >= 4)."""
    n_ = chains.shape[1]
    means = chains.mean(axis=1)
    B = n_ * means.var(ddof=1)
    W = chains.var(axis=1, ddof=1).mean()
    var_hat = (n_ - 1) / n_ * W + B / n_
    return np.sqrt(var_hat / W)


def ess(x):
    """Effective sample size from the autocorrelation function (initial positive sequence)."""
    x = x - x.mean()
    acf = np.correlate(x, x, "full")[len(x) - 1:] / (x.var() * len(x))
    tau = 1.0
    for k in range(1, len(acf)):
        if acf[k] < 0.05:
            break
        tau += 2 * acf[k]
    return len(x) / tau


n_iter, burn = 20000, 4000
step = np.array([0.02, 0.02, 0.1])
chains, acc = [], []
for _ in range(4):
    start = [rng.normal(0, 1), rng.normal(0, 0.5), rng.normal(0, 1)]  # dispersed starts
    ch, a = metropolis(log_post, start, step, n_iter, rng)
    chains.append(ch[burn:])
    acc.append(a)
chains = np.array(chains)  # (4, n, 3)
post = chains.reshape(-1, 3)
names = ["intercept", "slope (C/decade)", "log sigma"]
print(f"acceptance rates: {np.round(acc, 2)}  (target ~0.25-0.45)")
for j, nm in enumerate(names):
    print(f"  {nm:18s} R-hat={rhat(chains[:, :, j]):.3f}  ESS={sum(ess(chains[c, :, j]) for c in range(4)):.0f} "
          f"of {post.shape[0]}")
b = post[:, 1]
print(f"\nslope posterior: mean={b.mean():.4f}  95% credible=({np.percentile(b, 2.5):.4f}, "
      f"{np.percentile(b, 97.5):.4f})")
print(f"OLS:             est ={ols.params[1]:.4f}  95% CI      =({ols.conf_int()[1, 0]:.4f}, "
      f"{ols.conf_int()[1, 1]:.4f})")
print(f"P(slope > 0.15 C/decade) = {(b > 0.15).mean():.3f};  P(slope > 0.2) = {(b > 0.2).mean():.3f}")

# prior sensitivity: a sceptical prior centred on zero with SD 0.05
ch_s, _ = metropolis(lambda p: log_post(p, prior_slope_sd=0.05), [0, 0, 0], step, n_iter, rng)
b_s = ch_s[burn:, 1]
print(f"with a sceptical prior N(0, 0.05) on the slope: posterior mean={b_s.mean():.4f} "
      f"(shrunk from {b.mean():.4f}); the data still overwhelm it")

fig, axes = plt.subplots(2, 2, figsize=(12, 7))
for c in range(4):
    axes[0, 0].plot(chains[c, :2000, 1], lw=.5, alpha=.8)
axes[0, 0].set_title("Trace plot, slope, 4 chains (first 2000 post-burn-in draws)")
axes[0, 1].hist(b, bins=60, density=True, alpha=.8, label="posterior (weak prior)")
axes[0, 1].hist(b_s, bins=60, density=True, alpha=.5, label="posterior (sceptical prior)")
xs = np.linspace(b.min(), b.max(), 200)
axes[0, 1].plot(xs, sps.norm.pdf(xs, ols.params[1], ols.bse[1]), "k--", label="OLS sampling dist.")
axes[0, 1].set_title("Slope: posterior vs OLS")
axes[0, 1].legend(fontsize=8)
# posterior predictive check
idx = rng.integers(0, len(post), 200)
for i in idx[:100]:
    a_, b_, ls_ = post[i]
    axes[1, 0].plot(gt.index.year, a_ + b_ * t, color="C0", alpha=.05)
axes[1, 0].scatter(gt.index.year, y, s=12, color="k", zorder=3)
axes[1, 0].set_title("100 posterior regression lines over the data")
sim_sd = [np.std(a_ + b_ * t + rng.normal(0, np.exp(ls_), len(t)) - (a_ + b_ * t)) for a_, b_, ls_ in post[idx]]
axes[1, 1].hist(sim_sd, bins=30, alpha=.8, label="replicated residual SD")
axes[1, 1].axvline(ols.resid.std(), color="C3", label="observed residual SD")
axes[1, 1].set_title("Posterior predictive check on residual spread")
axes[1, 1].legend()
fig.tight_layout()
savefig(fig, "metropolis", SEC)
note("""Sixty lines of numpy is a complete MCMC workflow: dispersed starts, burn-in, acceptance
tuning, R-hat near 1 and a healthy ESS, then a posterior that matches OLS because the
prior is weak and n is moderate. The payoff is the last line: P(slope > 0.15 C/decade) is a
direct probability statement, which no p-value provides. In practice use PyMC / Stan /
NumPyro for the sampler (NUTS) and spend your effort on the model.""")
