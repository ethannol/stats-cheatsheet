# %% [markdown]
# # 12 Assay EDA: skewed grades, declustering, top-cuts, compositing
#
# Data: Walker Lake (470 samples of V; the 78,000-node exhaustive grid gives us the truth, which
# no real deposit ever does). GSLIB cluster.dat for a reference check of declustering weights.
#
# **Cheat sheet**
# - Grades are positive, right-skewed and often near-lognormal. CV = sd/mean is the first number
#   to quote: < 1 benign, 1-2 needs care, > 2 (nuggety gold) means a few samples carry the metal.
# - Probability plots (log scale) reveal mixed populations (kinks) and the upper tail behaviour.
# - Drilling is preferential: geologists drill where it is rich. The naive mean is biased high.
#   Cell declustering weights samples by 1/(samples per cell); choose the cell size at the
#   minimum of the declustered-mean curve (preferential high-grade sampling) or the sample spacing.
# - Top-cutting (capping) limits the influence of outliers on estimates. Common practice: cap at
#   a probability-plot break, or P97.5-P99, and report metal removed. Capping is a judgment call;
#   document it. Alternatives: restricted search ellipses, indicator methods.
# - Compositing to equal lengths before variography/estimation gives every sample equal support;
#   it conserves metal, lowers variance (volume-variance effect), and removes sample-length bias.
# - Always estimate within geological domains (here the T indicator): stationarity is a modelling
#   decision, not a property of the data.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import (
    cell_declustering,
    composite_intervals,
    declustering_scan,
    top_cut_summary,
)

setup()
SEC = "12_assay_eda"
w = data.walker()
ex = data.walker_exhaustive()
V = w.V.values
truth_mean = ex.V.mean()

# %% ----------------------------------------------------------- distribution
banner("Walker Lake V: distribution and skewness")
desc = pd.Series({"n": len(V), "mean": V.mean(), "median": np.median(V), "sd": V.std(ddof=1),
                  "CV": V.std(ddof=1) / V.mean(), "skew": sps.skew(V), "min": V.min(),
                  "P90": np.quantile(V, .9), "P99": np.quantile(V, .99), "max": V.max(),
                  "frac zero": (V == 0).mean()})
print(desc.round(3).to_string())
print(f"\nexhaustive truth: mean={truth_mean:.1f} median={ex.V.median():.1f} CV={ex.V.std() / truth_mean:.2f}")
pos = V[V > 0]
print(f"lognormal fit to V>0: sigma_ln={np.log(pos).std(ddof=1):.2f}  Shapiro on ln(V>0) p="
      f"{sps.shapiro(np.log(pos)).pvalue:.3f}")
by_T = w.groupby("T").V.agg(["count", "mean", "median", "std"])
by_T["CV"] = by_T["std"] / by_T["mean"]
print("\nby geological type T:\n" + by_T.round(2).to_string())
note(f"""The sample mean ({V.mean():.0f}) sits well above the exhaustive truth ({truth_mean:.0f}):
Walker Lake was sampled preferentially in the high-grade zones, exactly as a drilling campaign
would be. CV of {desc['CV']:.2f} is moderate for a base-metal style variable. The two T types
have very different grades and variability: estimate them separately.""")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
axes[0].hist(V, bins=40, alpha=.8)
axes[0].set_title("Histogram of V (right-skewed)")
axes[0].set_xlabel("V")
axes[1].hist(np.log10(pos), bins=40, alpha=.8)
axes[1].set_title("log10(V > 0): closer to symmetric")
axes[1].set_xlabel("log10 V")
# lognormal probability plot: sorted values vs normal quantiles, log y
srt = np.sort(pos)
pp = (np.arange(1, len(srt) + 1) - 0.5) / len(srt)
axes[2].scatter(sps.norm.ppf(pp), srt, s=8)
axes[2].set_yscale("log")
axes[2].set_xlabel("standard normal quantile")
axes[2].set_ylabel("V (log scale)")
axes[2].set_title("Lognormal probability plot: a straight line = lognormal")
for q in [0.975, 0.99]:
    axes[2].axhline(np.quantile(V, q), color="C3", ls="--", lw=.8)
    axes[2].text(-2.5, np.quantile(V, q) * 1.05, f"P{q * 100:.1f}", color="C3", fontsize=8)
fig.tight_layout()
savefig(fig, "distribution", SEC)

# %% ----------------------------------------------------------- declustering
banner("Cell declustering")
cells = np.arange(2, 60, 2.0)
dm = declustering_scan(w.x.values, w.y.values, V, cells)
best = cells[int(np.argmin(dm))]
wts = cell_declustering(w.x.values, w.y.values, best)
declus_mean = np.average(V, weights=wts)
print(f"naive mean={V.mean():.1f}  declustered mean (cell={best:.0f})={declus_mean:.1f}  "
      f"truth={truth_mean:.1f}")
print(f"weights: min={wts.min():.2f} max={wts.max():.2f} (1 = average sample)")
print(f"bias removed: {(V.mean() - declus_mean) / (V.mean() - truth_mean) * 100:.0f}% of the naive bias")
# reference check on GSLIB cluster.dat, which ships its own declustering weights
gc = data.gslib_cluster()
w_gslib = cell_declustering(gc.x.values, gc.y.values, 5.0)
r = np.corrcoef(w_gslib, gc.declus_weight)[0, 1]
print(f"\nGSLIB cluster.dat: naive mean={gc.primary.mean():.2f}, declustered (ours)="
      f"{np.average(gc.primary, weights=w_gslib):.2f}, declustered (shipped weights)="
      f"{np.average(gc.primary, weights=gc.declus_weight):.2f}; weight correlation r={r:.2f}")
note(f"""The declustered-mean curve bottoms out near cell size {best:.0f} (the spacing of the
regular sampling pass); that minimum is the conventional pick when high grades were
oversampled. Declustering removes most of the bias but not all: the clusters sit on the
richest ground, where even the declustered samples are rich. Declustered weights also feed
the histogram used for simulation and the normal-score transform (section 15).""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
axes[0].plot(cells, dm, "o-", ms=3)
axes[0].axhline(V.mean(), color="k", ls=":", label="naive mean")
axes[0].axhline(truth_mean, color="C3", ls="--", label="exhaustive truth")
axes[0].axvline(best, color="C1", lw=.8)
axes[0].set_xlabel("cell size")
axes[0].set_ylabel("declustered mean of V")
axes[0].set_title("Declustered mean vs cell size")
axes[0].legend()
sc = axes[1].scatter(w.x, w.y, c=V, s=18 + 60 * (wts - wts.min()) / (wts.max() - wts.min()),
                     cmap="viridis", edgecolor="k", lw=.3)
axes[1].set_aspect("equal")
axes[1].set_title("Samples: colour = grade, size = declustering weight")
fig.colorbar(sc, ax=axes[1], label="V")
fig.tight_layout()
savefig(fig, "declustering", SEC)

# %% ----------------------------------------------------------- top cuts
banner("Top-cut (capping) analysis")
rows = [top_cut_summary(V, np.quantile(V, q)) for q in [0.95, 0.975, 0.99, 0.995]]
tc = pd.DataFrame(rows, index=["P95", "P97.5", "P99", "P99.5"])
print(tc.round(3).to_string())
hi_T = w[w["T"] == 2].V.values
print(f"\nWithin T=2 only: P99={np.quantile(hi_T, .99):.0f}, max={hi_T.max():.0f}, "
      f"top 1% of samples hold {np.sort(hi_T)[-max(1, len(hi_T) // 100):].sum() / hi_T.sum() * 100:.1f}% of the metal")
note(f"""Capping at P99 touches {tc.loc['P99', 'n_cut']} samples and removes {tc.loc['P99', 'metal_removed'] * 100:.1f}%
of the metal while lowering the CV from {tc.loc['P99', 'cv_before']:.2f} to {tc.loc['P99', 'cv_after']:.2f}.
For this moderately skewed variable capping barely matters; for nuggety gold the same table
can remove 10-20% of the metal, and the choice becomes the single most consequential
assumption in the resource. Decide caps per domain, from the probability plot, and report
the metal removed.""")

# %% ----------------------------------------------------------- compositing
banner("Compositing irregular samples to 2-unit lengths (drillhole built from Walker Lake column)")
dh = data.drillhole_from_walker()
comp = composite_intervals(dh["from"], dh["to"], dh.grade, length=2.0)
comp4 = composite_intervals(dh["from"], dh["to"], dh.grade, length=4.0)
metal_raw = np.sum(dh.grade * dh.length)
metal_comp = np.sum(comp["grade"] * (comp["to"] - comp["from"]))
sd_raw, sd2, sd4 = dh.grade.std(ddof=1), comp["grade"].std(ddof=1), comp4["grade"].std(ddof=1)
print(f"raw samples: n={len(dh)}, lengths {dh.length.min()}-{dh.length.max()}, mean length={dh.length.mean():.2f}, "
      f"grade sd={sd_raw:.1f}")
print(f"2-unit composites: n={len(comp['grade'])}, grade sd={sd2:.1f};  4-unit composites: n={len(comp4['grade'])}, grade sd={sd4:.1f}")
print(f"metal (grade x length) raw={metal_raw:.0f} composited={metal_comp:.0f} -> conserved")
lw_mean = metal_raw / dh.length.sum()
print(f"naive mean of raw sample grades={dh.grade.mean():.1f} vs length-weighted mean={lw_mean:.1f}")
note(f"""Compositing conserves metal and length-weights the samples. Here the raw mean and the
length-weighted mean agree because these synthetic interval lengths were drawn independently
of grade; in real holes short samples target the richest intervals and the arithmetic mean
of raw assays is biased high. The 2-unit composite sd ({sd2:.0f}) barely differs from the raw sd
({sd_raw:.0f}) because the raw mean length is already ~2, but 4-unit composites drop to {sd4:.0f}:
the volume-variance effect. Variograms and cut-off statistics must be computed on the
support that will be estimated, not on raw assays.""")

fig, ax = plt.subplots(figsize=(12, 3.8))
ax.step(np.append(dh["from"], dh["to"].iloc[-1]), np.append(dh.grade, dh.grade.iloc[-1]), where="post",
        lw=.8, label="raw samples (irregular lengths)")
ax.step(np.append(comp["from"], comp["to"][-1]), np.append(comp["grade"], comp["grade"][-1]), where="post",
        lw=1.8, color="C3", label="2-unit composites")
ax.set_xlabel("downhole depth")
ax.set_ylabel("grade")
ax.set_title("Compositing smooths and length-weights")
ax.legend()
savefig(fig, "compositing", SEC)
