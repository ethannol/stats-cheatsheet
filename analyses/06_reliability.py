# %% [markdown]
# # 06 Equipment reliability and production SPC
#
# Data: SYNTHETIC haul-truck fleet (failures with right-censoring; daily production with a
# planted step change). Maintenance records are never public; the generating parameters are
# stated in data.equipment_failures so every estimate can be checked.
#
# **Cheat sheet**
# - Time-between-failure data are right-censored: a component still running at study end has a
#   lifetime *at least* that long. Dropping censored units or treating them as failures biases
#   MTBF low. Kaplan-Meier handles censoring non-parametrically; Weibull MLE parametrically.
# - Weibull shape beta: < 1 infant mortality (hazard falls: burn-in, quality problems), ~1 random
#   failures (exponential, constant hazard, maintenance interval does not help), > 1 wear-out
#   (hazard rises: preventive replacement pays). Scale eta = characteristic life (63.2% failed).
# - MTBF = mean of a Weibull = eta Gamma(1 + 1/beta), not the average of observed TBFs when
#   censoring exists. MTTR from repair times (usually lognormal). Inherent availability =
#   MTBF / (MTBF + MTTR).
# - Weibull probability plot: ln(-ln S(t)) vs ln t is a straight line with slope beta.
# - SPC: Shewhart I-MR chart catches big shifts; CUSUM and EWMA catch sustained small shifts
#   (0.5-1.5 sigma) several times faster. Set limits from an in-control phase I, and remove
#   known structure (seasonality, planned downtime) before charting or every winter alarms.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from miningstats import data
from miningstats.plotting import banner, note, savefig, setup
from miningstats.stats import cusum, ewma, kaplan_meier, shewhart_limits, weibull_mle

setup()
SEC = "06_reliability"
eq = data.equipment_failures()
f, prod = eq["failures"], eq["production"]
TRUE = {"engine": (2.4, 2400.0), "hydraulics": (0.85, 900.0), "tyres": (3.0, 1800.0)}

# %% ----------------------------------------------------------- Weibull by component
banner("Component lifetimes: Kaplan-Meier and Weibull MLE with right censoring")
rows = []
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
for comp, g in f.groupby("component"):
    fit = weibull_mle(g.tbf, g.event == 1)
    naive_mtbf = g.loc[g.event == 1, "tbf"].mean()
    as_failures = weibull_mle(g.tbf, np.ones(len(g), bool))
    mttr = g.ttr.dropna()
    avail = fit.mttf / (fit.mttf + mttr.mean())
    rows.append(dict(component=comp, n_fail=int(g.event.sum()), n_censored=int((g.event == 0).sum()),
                     beta=fit.shape, eta=fit.scale, true_beta=TRUE[comp][0], true_eta=TRUE[comp][1],
                     MTBF_weibull=fit.mttf, MTBF_naive=naive_mtbf,
                     MTBF_if_censored_as_failed=as_failures.mttf, MTTR_h=mttr.mean(), availability=avail))
    t, S, _ = kaplan_meier(g.tbf, g.event == 1)
    axes[0].step(t, S, where="post", label=comp)
    tt = np.linspace(1, g.tbf.max(), 200)
    axes[0].plot(tt, fit.reliability(tt), "--", color=axes[0].lines[-1].get_color(), lw=.8)
    # Weibull probability plot on KM survival
    m = (S > 0) & (S < 1) & (t > 0)
    axes[1].plot(np.log(t[m]), np.log(-np.log(S[m])), ".", label=f"{comp}: beta={fit.shape:.2f}")
    axes[2].plot(tt, fit.hazard(tt) * 1000, label=comp)
tab = pd.DataFrame(rows).set_index("component")
print(tab.round(3).to_string())
note(f"""Shape parameters recover the generating values (hydraulics beta<1: early-life/random
failures where scheduled replacement would waste parts; tyres beta=3: wear-out where
preventive replacement at a fraction of eta pays). The naive MTBF (mean of observed TBFs)
and 'censored as failed' both run below the Weibull MTBF because long lives are the ones that
get censored. Inherent availability ranges from {tab.availability.min():.3f} to {tab.availability.max():.3f};
at fleet level the components combine in series so unavailabilities roughly add.""")
axes[0].set_xlabel("operating hours")
axes[0].set_ylabel("reliability R(t)")
axes[0].set_title("Kaplan-Meier (steps) and Weibull fit (dashed)")
axes[0].legend()
axes[1].set_xlabel("ln(t)")
axes[1].set_ylabel("ln(-ln S(t))")
axes[1].set_title("Weibull probability plot: slope = beta")
axes[1].legend(fontsize=8)
axes[2].set_xlabel("operating hours")
axes[2].set_ylabel("hazard (failures per 1000 h)")
axes[2].set_title("Hazard functions: falling, flat-ish, rising")
axes[2].legend()
fig.tight_layout()
savefig(fig, "weibull", SEC)

# %% ----------------------------------------------------------- preventive replacement
banner("Preventive replacement interval for the wear-out component (tyres)")
g = f[f.component == "tyres"]
fit = weibull_mle(g.tbf, g.event == 1)
cost_pm, cost_fail = 1.0, 6.0  # planned replacement vs in-service failure (downtime, damage)
T = np.linspace(200, 3000, 300)
R = fit.reliability(T)
# cost rate of an age-replacement policy: [c_pm R(T) + c_f (1-R(T))] / expected cycle length
cycle = np.array([np.trapezoid(fit.reliability(np.linspace(0, t, 400)), np.linspace(0, t, 400)) for t in T])
rate = (cost_pm * R + cost_fail * (1 - R)) / cycle
Tstar = T[np.argmin(rate)]
print(f"optimal age-replacement interval ~{Tstar:.0f} h ({Tstar / fit.scale:.2f} eta) at failure/PM cost ratio "
      f"{cost_fail / cost_pm:.0f}; run-to-failure cost rate={cost_fail / fit.mttf:.5f} vs optimum {rate.min():.5f} per h")
note("""Age replacement only makes sense when beta > 1 and failures cost more than planned swaps;
the optimum sits well before the characteristic life. For the hydraulics (beta < 1) the same
curve is monotone decreasing: replace only on failure.""")

# %% ----------------------------------------------------------- fleet homogeneity
banner("Is any truck a 'bad actor'? Poisson homogeneity of failure counts")
counts = f[f.event == 1].groupby("truck").size()
chi2, p = sps.chisquare(counts)
print(f"failures per truck: {counts.to_dict()}")
print(f"chi-square homogeneity: chi2={chi2:.1f} df={len(counts) - 1} p={p:.3f}  (dispersion var/mean={counts.var() / counts.mean():.2f})")
note("""Before blaming a truck or an operator, test whether the spread in counts exceeds what a
common Poisson rate would produce. Here the fleet is homogeneous by construction; a real
fleet with p < 0.05 and var/mean >> 1 has a bad actor or a maintenance-history effect worth
a negative-binomial model with truck-level covariates.""")

# %% ----------------------------------------------------------- production SPC
banner("Production SPC: detecting a sustained 7% drop in daily tonnes")
y = prod.tonnes.values
day = prod.day.values
phase1 = 365
# remove the seasonal cycle estimated from phase I only
Xs = np.column_stack([np.ones(phase1), np.sin(2 * np.pi * day[:phase1] / 365), np.cos(2 * np.pi * day[:phase1] / 365)])
beta = np.linalg.lstsq(Xs, y[:phase1], rcond=None)[0]
Xall = np.column_stack([np.ones(len(day)), np.sin(2 * np.pi * day / 365), np.cos(2 * np.pi * day / 365)])
resid = y - Xall @ beta
centre, sigma, lo, hi = shewhart_limits(resid, phase1=phase1)
shift_day = int(np.argmax(prod.planted_shift.values)) + 1
oob = np.where((resid < lo) | (resid > hi))[0] + 1
first_shew = oob[oob >= shift_day].min() if (oob >= shift_day).any() else None
cp, cm, sig = cusum(resid, centre, sigma, k=0.5, h=5)
cu_first = sig[sig + 1 >= shift_day].min() + 1 if (sig + 1 >= shift_day).any() else None
z, elo, ehi = ewma(resid, lam=0.2, target=centre, sigma=sigma, L=2.7)
ew_sig = np.where((z < elo) | (z > ehi))[0] + 1
first_ewma = ew_sig[ew_sig >= shift_day].min() if (ew_sig >= shift_day).any() else None
shift_sigma = 0.07 * y[:phase1].mean() / sigma
print(f"planted shift at day {shift_day}: -7% = {shift_sigma:.2f} sigma of daily noise")
print(f"Shewhart I-chart: false alarms before shift={np.sum(oob < shift_day)}, first signal after shift: day {first_shew}")
print(f"CUSUM (k=0.5, h=5): false alarms before shift={np.sum(sig + 1 < shift_day)}, first signal after shift: day {cu_first}")
print(f"EWMA (lambda=0.2): false alarms before shift={np.sum(ew_sig < shift_day)}, first signal after shift: day {first_ewma}")
delay = {"Shewhart": first_shew, "CUSUM": cu_first, "EWMA": first_ewma}
note(f"""A 7% production drop is only ~{shift_sigma:.1f} sigma of day-to-day noise. Detection delays
after day {shift_day}: Shewhart {delay['Shewhart'] - shift_day if first_shew else 'none'} days, CUSUM
{delay['CUSUM'] - shift_day if cu_first else 'none'} days, EWMA {delay['EWMA'] - shift_day if first_ewma else 'none'} days.
The accumulating charts win on delay but pay in false alarms over a 600-day in-control run
(ARL0 ~ 465 for these CUSUM settings): the design trade-off is choosing h or lambda for the
false-alarm rate the operation can live with. Charting raw tonnes without removing the
seasonal cycle would alarm every winter and go blind every summer.""")

fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
axes[0].plot(day, resid, lw=.6)
axes[0].axhline(centre, color="k", lw=.8)
axes[0].axhline(lo, color="C3", ls="--")
axes[0].axhline(hi, color="C3", ls="--")
axes[0].axvline(shift_day, color="C1", ls=":", label="planted shift")
axes[0].set_title("Individuals chart of deseasonalised daily tonnes (limits from phase I)")
axes[0].legend()
axes[1].plot(day, cp, label="CUSUM+")
axes[1].plot(day, cm, label="CUSUM-")
axes[1].axhline(5, color="C3", ls="--")
axes[1].axvline(shift_day, color="C1", ls=":")
axes[1].set_title("Tabular CUSUM (restarts after each signal)")
axes[1].legend()
axes[2].plot(day, z, lw=1)
axes[2].plot(day, elo, "C3--", lw=.8)
axes[2].plot(day, ehi, "C3--", lw=.8)
axes[2].axvline(shift_day, color="C1", ls=":")
axes[2].set_title("EWMA (lambda = 0.2)")
axes[2].set_xlabel("day")
savefig(fig, "production_spc", SEC)
