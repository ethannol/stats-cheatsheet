# %% [markdown]
# # 03 Assay QA/QC and censored data
#
# Data: SYNTHETIC gold QA/QC programme (duplicates, CRMs, blanks) with planted problems, so we
# know the answers; then REAL censored data: dissolved copper in the Clark Fork River with
# non-detects.
#
# **Cheat sheet**
# - QA/QC is inserted blind at ~5% each: field duplicates (precision of sampling + prep + assay),
#   certified reference materials (accuracy / lab bias), blanks (contamination / carry-over).
# - Precision: scatter plot, HARD = |a-b|/(a+b) with the rule of thumb '90% of pairs under 20%'
#   for coarse gold duplicates (10% for pulp duplicates), relative difference vs mean plot,
#   Thompson-Howarth sd = s0 + k*c, i.e. precision is concentration-dependent.
# - Accuracy: CRM control chart with certified value +-2SD (warning) and +-3SD (action). A run
#   of points on one side is bias even when none breach the limits: use CUSUM. Mix-ups look like
#   single wild points; remove them from the bias estimate but investigate them.
# - Blanks: plot vs detection limit *and* vs the grade of the preceding sample: carry-over
#   shows as blanks correlated with the previous sample.
# - Non-detects are left-censored, not zeros and not 'DL/2'. Substitution biases means and
#   variances; use Kaplan-Meier (distribution-free), ROS (parametric imputation) or MLE.
#   With > 50-80% censored, report percentiles above the DL only.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from miningstats import data
from miningstats.plotting import banner, note, savefig, setup
from miningstats.stats import cusum, hard, km_left_censored, ros_lognormal, thompson_howarth

setup()
SEC = "03_qaqc"
qa = data.assay_qaqc()
dup, crm, bl, DL = qa["duplicates"], qa["crm"], qa["blanks"], qa["dl"]

# %% ----------------------------------------------------------- duplicates
banner("Field duplicates: precision")
h = hard(dup.original, dup.duplicate)
p90 = np.quantile(h, .9)
print(f"{len(dup)} pairs; HARD median={np.median(h):.3f}, 90th percentile={p90:.3f} "
      f"({'passes' if p90 < .2 else 'fails'} the 90% < 20% criterion)")
n_flag = (h > .5).sum()
n_caught = (dup.planted_swap & (h > .5)).sum()
print(f"pairs with HARD > 0.5: {n_flag}  (planted swaps: {dup.planted_swap.sum()}, of which flagged: {n_caught})")
bm, bsd, s0, k = thompson_howarth(dup.original, dup.duplicate)
print(f"Thompson-Howarth: sd = {s0:.3f} + {k:.3f} * c   -> precision at 1 g/t = {(s0 + k) / 1 * 100:.0f}% "
      f"(generated with sd = sqrt((0.12 c)^2 + 0.02^2), i.e. 12% at high grade, 0.02 g/t floor)")
ratio = dup.duplicate.mean() / dup.original.mean()
print(f"mean duplicate / mean original = {ratio:.3f} (should be ~1; a systematic departure means "
      f"the duplicate stream is treated differently)")
note(f"""HARD summarises precision in one number per pair; its 90th percentile ({p90:.2f}) is
the figure that goes in the QA/QC report. The Thompson-Howarth line recovers the generating
model (a concentration-proportional term plus an absolute floor), which is how one should
state precision: 'X% at grade Y', never a single percentage. Only {n_caught} of the
{dup.planted_swap.sum()} planted swaps show up as HARD > 0.5: a swap between two samples of
similar grade is invisible to any duplicate statistic, which is why sample tracking (barcodes,
chain of custody) matters more than statistics for gross errors.""")

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
axes[0].loglog(dup.original, dup.duplicate, ".", alpha=.6)
lim = [DL, dup[["original", "duplicate"]].max().max() * 1.2]
axes[0].plot(lim, lim, "k-", lw=.8)
axes[0].plot(lim, [v * 1.2 for v in lim], "k:", lw=.8)
axes[0].plot(lim, [v / 1.2 for v in lim], "k:", lw=.8)
axes[0].scatter(dup.original[dup.planted_swap], dup.duplicate[dup.planted_swap], s=60,
                facecolor="none", edgecolor="C3", label="planted swaps")
axes[0].set_xlabel("original (g/t)")
axes[0].set_ylabel("duplicate (g/t)")
axes[0].set_title("Duplicate scatter (dotted = +-20%)")
axes[0].legend()
mean_pair = 0.5 * (dup.original + dup.duplicate)
axes[1].semilogx(mean_pair, h, ".", alpha=.6)
axes[1].axhline(.2, color="C3", ls="--", label="20%")
axes[1].axhline(p90, color="C1", ls=":", label=f"observed P90={p90:.2f}")
axes[1].set_xlabel("pair mean (g/t)")
axes[1].set_ylabel("HARD")
axes[1].set_title("Precision worsens at low grade")
axes[1].legend()
axes[2].plot(bm, bsd, "o")
cc = np.linspace(0, bm.max(), 50)
axes[2].plot(cc, s0 + k * cc, "C3-", label=f"sd = {s0:.3f} + {k:.3f} c")
axes[2].set_xlabel("concentration (g/t)")
axes[2].set_ylabel("sd of a single assay (g/t)")
axes[2].set_title("Thompson-Howarth precision model")
axes[2].legend()
fig.tight_layout()
savefig(fig, "duplicates", SEC)

# %% ----------------------------------------------------------- CRMs
banner("Certified reference materials: accuracy and bias")
fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
for ax, (name, g) in zip(axes, crm.groupby("crm")):
    g = g.sort_values("batch")
    z = (g.assay - g.certified) / g.cert_sd
    out3 = np.abs(z) > 3
    warn = (np.abs(z) > 2) & ~out3
    clean = g[~out3]
    bias = (clean.assay.mean() / clean.certified.iloc[0] - 1) * 100
    cp, cm, sig = cusum(clean.assay.values, clean.certified.iloc[0], clean.cert_sd.iloc[0], k=0.5, h=4)
    first_sig_batch = int(clean.batch.iloc[sig[0]]) if len(sig) else None
    print(f"{name}: certified {g.certified.iloc[0]:.2f}; {out3.sum()} points beyond 3SD (planted mix-ups in "
          f"this CRM: {g.planted_mixup.sum()}); {warn.sum()} warnings; overall bias (ex-outliers)={bias:+.1f}%; "
          f"CUSUM signals at batch {first_sig_batch}")
    ax.plot(g.batch, g.assay, "o-", ms=3, lw=.6)
    for s_, ls in [(2, "--"), (3, "-")]:
        ax.axhline(g.certified.iloc[0] + s_ * g.cert_sd.iloc[0], color="C3", ls=ls, lw=.8)
        ax.axhline(g.certified.iloc[0] - s_ * g.cert_sd.iloc[0], color="C3", ls=ls, lw=.8)
    ax.axhline(g.certified.iloc[0], color="k", lw=.8)
    ax.axvspan(70, 90, color="C1", alpha=.15, label="planted +8% bias period")
    ax.set_ylabel(f"{name} (g/t)")
    if ax is axes[0]:
        ax.legend(loc="upper left")
axes[0].set_title("CRM control charts (dashed +-2SD warning, solid +-3SD action)")
axes[-1].set_xlabel("lab batch")
savefig(fig, "crm_control_charts", SEC)
note("""An 8% bias is below the 3SD action line for most individual CRM inserts, so a limits-only
chart passes a lab that is over-reporting grade by 8%. The CUSUM, or simply 'how many of the
last 10 points are above the line', catches it within a few batches. Bias is estimated per
CRM and per period after removing the obvious mix-ups; a bias that scales with grade
(all three CRMs off by the same percentage) points at calibration, a constant offset at
a blank/contamination problem.""")

# %% ----------------------------------------------------------- blanks
banner("Blanks: contamination and carry-over")
fail = bl.assay > 5 * DL
print(f"{len(bl)} blanks; {bl.below_dl.sum()} below DL; {fail.sum()} above 5xDL "
      f"(planted carry-over events: {bl.planted_contamination.sum()}, flagged: {(fail & bl.planted_contamination).sum()})")
hi_prev = bl.prev_grade > 3
print(f"mean blank after a >3 g/t sample: {bl.assay[hi_prev].mean() * 1000:.1f} ppb vs "
      f"{bl.assay[~hi_prev].mean() * 1000:.1f} ppb otherwise")
fig, axes = plt.subplots(1, 2, figsize=(12, 4.3))
axes[0].semilogy(bl.batch, np.maximum(bl.assay, DL / 10), ".", alpha=.7)
axes[0].axhline(DL, color="k", ls="--", label="DL")
axes[0].axhline(5 * DL, color="C3", ls="-", label="5 x DL action")
axes[0].set_xlabel("batch")
axes[0].set_ylabel("blank assay (g/t)")
axes[0].set_title("Blank chart")
axes[0].legend()
axes[1].loglog(bl.prev_grade, np.maximum(bl.assay, DL / 10), ".", alpha=.7)
axes[1].axhline(5 * DL, color="C3")
axes[1].set_xlabel("grade of preceding sample (g/t)")
axes[1].set_ylabel("blank assay (g/t)")
axes[1].set_title("Carry-over: failing blanks follow high-grade samples")
fig.tight_layout()
savefig(fig, "blanks", SEC)
note("""Every failing blank follows a high-grade sample: that is carry-over in sample preparation
(insufficient cleaning of the pulveriser), not laboratory contamination in general. The
remedy is a prep protocol change and re-assay of the samples that followed high grades,
which is a different (and cheaper) action than rejecting the lab.""")

# %% ----------------------------------------------------------- censored data (real)
banner("Non-detects: trace metals in the Clark Fork River at Deer Lodge (USGS)")
wq = data.clark_fork_wq()
metals = wq[(wq.characteristic != "pH") & (wq.unit == "ug/l")]
nd_counts = metals.groupby(["characteristic", "fraction"]).nondetect.agg(["sum", "size"])
nd_counts["pct"] = nd_counts["sum"] / nd_counts["size"] * 100
print("non-detects by analyte:\n" + nd_counts.sort_values("sum", ascending=False).head(6).round(0).to_string())
char, frac = nd_counts["sum"].idxmax()
cu = metals[(metals.characteristic == char) & (metals.fraction == frac)].copy()
n_cens = cu.nondetect.sum()
print(f"\nusing {frac.lower()} {char.lower()}: {len(cu)} results, {n_cens} non-detects ({n_cens / len(cu) * 100:.0f}%), "
      f"detection limits used: {sorted(cu.loc[cu.nondetect, 'value'].unique())}")
v, c = cu.value.values, cu.nondetect.values
km = km_left_censored(v, c)
ros = ros_lognormal(v, c)
subs = {"substitute 0": np.where(c, 0, v), "substitute DL/2": np.where(c, v / 2, v), "substitute DL": v}
rows = [{"method": k, "mean": s.mean(), "median": np.median(s), "sd": s.std(ddof=1)} for k, s in subs.items()]
rows.append({"method": "Kaplan-Meier", "mean": km["mean"], "median": km["median"], "sd": np.nan})
rows.append({"method": "ROS (lognormal)", "mean": ros["mean"], "median": ros["median"],
             "sd": ros["imputed_values"].std(ddof=1)})
print(pd.DataFrame(rows).round(2).to_string(index=False))
res_tab = pd.DataFrame(rows).set_index("method")
spread = res_tab["mean"].max() - res_tab["mean"].min()
note(f"""With {n_cens / len(cu) * 100:.0f}% non-detects and detection limits spanning three orders of
magnitude (one old method reported '< 100'), substitution alone moves the mean from
{res_tab['mean'].min():.2f} to {res_tab['mean'].max():.2f} ug/L against a Kaplan-Meier estimate of
{res_tab.loc['Kaplan-Meier', 'mean']:.2f}. Kaplan-Meier needs no distributional assumption and is the regulatory default
(EPA ProUCL, Helsel); ROS assumes lognormality and extrapolates below the DL; substitution is
never right and gets worse when detection limits changed over time, as they did here. Above
~50% censoring, stop estimating means and report the upper percentiles only.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
axes[0].semilogy(cu.date[~c], v[~c], ".", ms=4, label="detected")
axes[0].semilogy(cu.date[c], v[c], "v", ms=5, color="C3", label="non-detect (plotted at DL)")
axes[0].set_ylabel(f"{frac.lower()} {char} (ug/L)")
axes[0].set_title(f"Clark Fork {frac.lower()} {char.lower()}: non-detects and changing detection limits")
axes[0].legend()
axes[1].step(km["x"], km["F"], where="post", label="Kaplan-Meier CDF")
srt = np.sort(ros["imputed_values"])
axes[1].plot(srt, np.arange(1, len(srt) + 1) / len(srt), "--", label="ROS completed sample")
axes[1].set_xscale("log")
axes[1].set_xlabel(f"{frac.lower()} {char} (ug/L)")
axes[1].set_ylabel("cumulative probability")
axes[1].set_title("Censored CDF estimates")
axes[1].legend()
fig.tight_layout()
savefig(fig, "nondetects", SEC)
