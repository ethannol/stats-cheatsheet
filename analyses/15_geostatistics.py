# %% [markdown]
# # 15 Geostatistics for resource estimation
#
# Data: Walker Lake V. 470 samples to estimate; 78,000 exhaustive values to check against. This
# is the one dataset where we can score estimators against the truth, so we do.
#
# **Cheat sheet**
# - Variogram gamma(h) = half the mean squared difference of pairs h apart. Nugget (sampling error
#   + micro-scale), sill (~variance), range (correlation length). Compute omnidirectional first,
#   then directional (azimuth cones) to find anisotropy; model it as a stretched ellipse.
# - Fit the short lags carefully: they drive kriging weights. Never fit beyond half the domain.
# - Ordinary kriging = best linear unbiased estimator given the variogram: weights account for
#   distance *and* redundancy (clustered samples share weight, unlike IDW). It returns a kriging
#   variance that depends on data geometry only, not on the values.
# - Cross-validation (leave-one-out) checks the model: mean error ~0, RMSE, and the slope of
#   truth-on-estimate ~1 (conditional unbiasedness). Standardised errors should have variance ~1.
# - Smoothing: kriged/IDW estimates have less variance than the truth; estimates at point
#   support exaggerate the variance of block grades (volume-variance relation). Grade-tonnage
#   curves from point estimates overstate high-grade tonnes: use block support.
# - Simulation (SGS) reproduces the histogram and variogram, giving many equally probable
#   models: use it for uncertainty, probability above cut-off, and grade-tonnage ranges; use
#   kriging for the single best estimate.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.stats import (
    cell_declustering,
    experimental_variogram,
    fit_variogram,
    grade_tonnage,
    idw,
    nearest_neighbour,
    nscore,
    nscore_back,
    ordinary_kriging,
    sequential_gaussian_simulation,
    spherical,
)

setup()
SEC = "15_geostatistics"
w = data.walker()
ex = data.walker_exhaustive()
xy = w[["x", "y"]].to_numpy()
V = w.V.values
BLOCK = 5
# true block grades: average the exhaustive grid in BLOCK x BLOCK cells
ex["bx"] = (ex.x - 1) // BLOCK
ex["by"] = (ex.y - 1) // BLOCK
blocks = ex.groupby(["bx", "by"]).agg(V=("V", "mean"), x=("x", "mean"), y=("y", "mean")).reset_index()
bxy = blocks[["x", "y"]].to_numpy()
truth = blocks.V.values
print(f"{len(V)} samples -> {len(blocks)} blocks of {BLOCK}x{BLOCK}; true block mean={truth.mean():.1f}, "
      f"true block sd={truth.std():.1f}, point sd={ex.V.std():.1f}, sample sd={V.std():.1f}")

# %% ----------------------------------------------------------- variography
banner("Variography: omnidirectional and directional")
lag, nl = 5.0, 16
h_o, g_o, n_o = experimental_variogram(w.x, w.y, V, lag, nl)
p_omni = fit_variogram(h_o, g_o, spherical, counts=n_o, nugget_max=0.6 * V.var())
print(f"omnidirectional spherical: nugget={p_omni[0]:.0f} sill={p_omni[1]:.0f} range={p_omni[2]:.1f}  "
      f"(sample variance {V.var():.0f}); relative nugget={p_omni[0] / p_omni[1]:.2f}")
dirs = {}
for az in [0, 45, 90, 135]:
    h_d, g_d, n_d = experimental_variogram(w.x, w.y, V, lag, nl, azimuth=az, tol_deg=22.5, bandwidth=12)
    # fit range only, nugget & sill fixed at the omnidirectional values
    def wss(r, h_d=h_d, g_d=g_d, n_d=n_d):
        return np.sum(n_d / h_d * (spherical(h_d, p_omni[0], p_omni[1], r) - g_d) ** 2)

    res = minimize_scalar(wss, bounds=(5, 200), method="bounded")
    dirs[az] = dict(h=h_d, g=g_d, n=n_d, range=res.x)
    print(f"  azimuth {az:3d}: range={res.x:.1f}")
major = max(dirs, key=lambda a: dirs[a]["range"])
minor = min(dirs, key=lambda a: dirs[a]["range"])
ratio = dirs[minor]["range"] / dirs[major]["range"]
print(f"anisotropy: major axis azimuth {major} (range {dirs[major]['range']:.0f}), minor {minor} "
      f"(range {dirs[minor]['range']:.0f}), ratio={ratio:.2f}")
aniso_params = (p_omni[0], p_omni[1], dirs[major]["range"])
note(f"""The omnidirectional variogram has a relative nugget of {p_omni[0] / p_omni[1]:.0%} and a range of
~{p_omni[2]:.0f} units. Directionally, continuity is longest along azimuth {major} and shortest along
{minor} (ratio {ratio:.2f}): the mineralisation is elongated, as the sample map suggests, and the
kriging search and model should follow that ellipse. Walker Lake is the textbook case of
geometric anisotropy.""")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
axes[0].scatter(h_o, g_o, s=np.sqrt(n_o), label="omnidirectional (size ~ pairs)")
hh = np.linspace(0, lag * nl, 200)
axes[0].plot(hh, spherical(hh, *p_omni), "C3", label="spherical fit")
axes[0].axhline(V.var(), color="k", ls=":", lw=.8, label="sample variance")
axes[0].set_xlabel("lag distance")
axes[0].set_ylabel("semivariance")
axes[0].set_title("Omnidirectional variogram of V")
axes[0].legend()
for (az, d), c in zip(dirs.items(), ["C0", "C1", "C2", "C4"]):
    axes[1].plot(d["h"], d["g"], "o-", ms=3, lw=.8, color=c, label=f"az {az} (range {d['range']:.0f})")
axes[1].axhline(p_omni[1], color="k", ls=":", lw=.8)
axes[1].set_xlabel("lag distance")
axes[1].set_title("Directional variograms: geometric anisotropy")
axes[1].legend()
fig.tight_layout()
savefig(fig, "variograms", SEC)

# %% ----------------------------------------------------------- estimation vs truth
banner(f"Estimating {BLOCK}x{BLOCK} blocks: OK vs IDW vs nearest neighbour, scored against truth")
# block OK: average of point estimates over a 2x2 discretisation of each block
offs = np.array([[-1.25, -1.25], [-1.25, 1.25], [1.25, -1.25], [1.25, 1.25]])
disc = (bxy[:, None, :] + offs[None, :, :]).reshape(-1, 2)
ok_pts, ok_var_pts = ordinary_kriging(xy, V, disc, spherical, aniso_params, n_max=24,
                                      anisotropy=(major, ratio))
ok_block = ok_pts.reshape(-1, 4).mean(axis=1)
mean_point_var = ok_var_pts.reshape(-1, 4).mean(axis=1)
# This is descriptive point uncertainty, NOT variance of the block mean.
# That requires cross-point error covariances or a direct block-kriging solve.
ok_point, _ = ordinary_kriging(xy, V, bxy, spherical, aniso_params, n_max=24, anisotropy=(major, ratio))
idw2 = idw(xy, V, bxy, power=2)
nn = nearest_neighbour(xy, V, bxy)
wts = cell_declustering(w.x.values, w.y.values, 20)
est = {"OK (block, anisotropic)": ok_block, "OK (point)": ok_point, "IDW p=2": idw2, "nearest neighbour": nn}
rows = []
for name, e in est.items():
    err = e - truth
    slope = np.polyfit(e, truth, 1)[0]
    rows.append(dict(estimator=name, mean_est=e.mean(), bias=err.mean(), RMSE=np.sqrt(np.mean(err ** 2)),
                     r=np.corrcoef(e, truth)[0, 1], sd_est=e.std(), slope_truth_on_est=slope))
res = pd.DataFrame(rows).set_index("estimator")
print(res.round(2).to_string())
print(f"\ntruth: mean={truth.mean():.1f} sd={truth.std():.1f};  declustered sample mean={np.average(V, weights=wts):.1f}")
note(f"""The samples average {V.mean():.0f} against a truth of {truth.mean():.0f}, yet OK and nearest
neighbour are essentially unbiased ({res.loc['OK (block, anisotropic)', 'bias']:+.0f} and
{res.loc['nearest neighbour', 'bias']:+.0f}): both decluster implicitly, kriging through the
redundancy term in its weights, NN through its polygons of influence. IDW does not, and
carries a {res.loc['IDW p=2', 'bias']:+.0f} bias straight from the preferential sampling. OK has the
lowest RMSE and the truth-on-estimate slope closest to 1 (least conditional bias). Every
estimator is smoother than the truth (sd {truth.std():.0f}); NN is the least smooth and the
least accurate. Smoothness is the price of minimising squared error.""")

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
grid_shape = (int(blocks.by.max()) + 1, int(blocks.bx.max()) + 1)


def as_grid(v):
    g = np.full(grid_shape, np.nan)
    g[blocks.by.astype(int), blocks.bx.astype(int)] = v
    return g


vmax = np.quantile(truth, .99)
for ax, (name, v) in zip(axes.ravel(), [("truth (block means)", truth), ("OK block estimate", ok_block),
                                         ("RMS point kriging SD (not block SE)", np.sqrt(np.maximum(mean_point_var, 0))),
                                         ("IDW p=2", idw2), ("nearest neighbour", nn),
                                         ("OK error (estimate - truth)", ok_block - truth)]):
    if "sd" in name:
        im = ax.imshow(as_grid(v), origin="lower", cmap="magma_r")
    elif "error" in name:
        im = ax.imshow(as_grid(v), origin="lower", cmap="RdBu_r", vmin=-400, vmax=400)
    else:
        im = ax.imshow(as_grid(v), origin="lower", cmap="viridis", vmin=0, vmax=vmax)
    ax.scatter((w.x - 1) / BLOCK, (w.y - 1) / BLOCK, s=2, color="w", alpha=.5)
    ax.set_title(name)
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=.8)
fig.tight_layout()
savefig(fig, "estimation_maps", SEC)

# %% ----------------------------------------------------------- cross-validation
banner("Leave-one-out cross-validation of the kriging model (what you can do without a truth)")
loo = np.empty(len(V))
loo_var = np.empty(len(V))
for i in range(len(V)):
    keep = np.arange(len(V)) != i
    e, vv = ordinary_kriging(xy[keep], V[keep], xy[i:i + 1], spherical, aniso_params, n_max=24,
                             anisotropy=(major, ratio))
    loo[i], loo_var[i] = e[0], vv[0]
err = loo - V
z = err / np.sqrt(np.maximum(loo_var, 1e-9))
print(f"LOO: mean error={err.mean():+.1f}  RMSE={np.sqrt(np.mean(err ** 2)):.1f}  r={np.corrcoef(loo, V)[0, 1]:.2f}  "
      f"slope(truth on est)={np.polyfit(loo, V, 1)[0]:.2f}")
print(f"standardised errors: mean={z.mean():+.2f} variance={z.var():.2f} (want ~0 and ~1)")
note(f"""Cross-validation is the only model check available on a real deposit. Standardised
error variance of {z.var():.2f} says the kriging variance is {'about right' if 0.8 < z.var() < 1.25 else 'mis-scaled'}
for the variogram sill used; a value far from 1 means the nugget/sill is wrong. The slope
below 1 is conditional bias: high estimates are too high, low ones too low, which is the
smoothing/selection problem that becomes money at the cut-off grade.""")

# %% ----------------------------------------------------------- volume-variance & grade tonnage
banner("Change of support: variance vs block size, and grade-tonnage curves")
rows = []
for b in [1, 2, 5, 10, 20, 50]:
    g = ex.groupby([(ex.x - 1) // b, (ex.y - 1) // b]).V.mean()
    rows.append(dict(block=b, n_blocks=len(g), variance=g.var(), sd=g.std()))
vv = pd.DataFrame(rows).set_index("block")
print(vv.round(1).to_string())
f_red = vv.loc[BLOCK, "variance"] / vv.loc[1, "variance"]
print(f"\nvariance reduction factor point -> {BLOCK}x{BLOCK}: {f_red:.2f}  (Krige's relation: D2(point,domain) = D2(point,block) + D2(block,domain))")
cutoffs = np.arange(0, 1001, 50)
gt_truth = pd.DataFrame(grade_tonnage(truth, cutoffs))
gt_ok = pd.DataFrame(grade_tonnage(ok_block, cutoffs))
gt_idw = pd.DataFrame(grade_tonnage(idw2, cutoffs))
gt_pts = pd.DataFrame(grade_tonnage(ex.V.values, cutoffs, block_tonnes=1 / BLOCK ** 2))  # point support, same total
c = 500
print(f"\nat cut-off {c}: tonnes above (blocks) truth={gt_truth.set_index('cutoff').loc[c, 'tonnes']:.0f}  "
      f"OK={gt_ok.set_index('cutoff').loc[c, 'tonnes']:.0f}  IDW={gt_idw.set_index('cutoff').loc[c, 'tonnes']:.0f}  "
      f"point-support truth={gt_pts.set_index('cutoff').loc[c, 'tonnes']:.0f}")
note(f"""Variance falls as block size grows (point sd {vv.loc[1, 'sd']:.0f} -> {BLOCK}x{BLOCK} block sd
{vv.loc[BLOCK, 'sd']:.0f}): the volume-variance relation. The grade-tonnage curve therefore depends
on support: at a cut-off of {c}, point-support truth reports {gt_pts.set_index('cutoff').loc[c, 'tonnes'] / gt_truth.set_index('cutoff').loc[c, 'tonnes'] - 1:+.0%}
tonnes relative to the block truth. Smoothed estimates (OK/IDW) err the other way: too few
high-grade tonnes at high cut-offs, too many just above a low cut-off. A resource must state
its selective mining unit, and change-of-support corrections exist for exactly this.""")

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
axes[0].loglog(vv.index, vv.variance, "o-")
axes[0].set_xlabel("block side")
axes[0].set_ylabel("variance of block grades")
axes[0].set_title("Volume-variance relation (exhaustive data)")
for gt, lab, ls in [(gt_truth, "truth (5x5 blocks)", "-"), (gt_ok, "OK block estimates", "--"),
                    (gt_idw, "IDW estimates", ":"), (gt_pts, "truth at point support", "-.")]:
    axes[1].plot(gt.cutoff, gt.tonnes / gt_truth.tonnes.iloc[0], ls, label=lab)
    axes[2].plot(gt.cutoff, gt.grade, ls, label=lab)
axes[1].set_xlabel("cut-off")
axes[1].set_ylabel("fraction of tonnes above cut-off")
axes[1].set_title("Grade-tonnage: tonnes")
axes[1].legend(fontsize=8)
axes[2].set_xlabel("cut-off")
axes[2].set_ylabel("mean grade above cut-off")
axes[2].set_title("Grade-tonnage: grade")
fig.tight_layout()
savefig(fig, "volume_variance_grade_tonnage", SEC)

# %% ----------------------------------------------------------- SGS
banner("Sequential Gaussian simulation: uncertainty instead of a single smoothed map")
y_ns, table = nscore(V, weights=wts)
h_n, g_n, n_n = experimental_variogram(w.x, w.y, y_ns, lag, nl)
p_ns = fit_variogram(h_n, g_n, spherical, counts=n_n, nugget_max=0.5)
p_ns = (p_ns[0], 1.0, p_ns[2])  # standard normal: sill fixed at 1
print(f"normal-score variogram: nugget={p_ns[0]:.2f} range={p_ns[2]:.1f}")
N_REAL = 8
sims_ns = sequential_gaussian_simulation(xy, y_ns, bxy, spherical, p_ns, n_real=N_REAL, n_max=16, seed=0)
sims = np.array([nscore_back(s, table, tail_max=V.max() * 1.1) for s in sims_ns])
# realisations are point-support; average 4 realisations' nodes? No: compare realisations to truth
# honestly at point support by also providing E-type and spread
etype = sims.mean(axis=0)
print(f"realisation sd (mean over {N_REAL}): {sims.std(axis=1).mean():.0f}  vs point truth sd {ex.V.std():.0f}, "
      f"E-type sd {etype.std():.0f}, OK sd {ok_point.std():.0f}")
print(f"realisation mean: {sims.mean():.0f}  (declustered sample mean {np.average(V, weights=wts):.0f}; truth {truth.mean():.0f})")
h_s, g_s, _ = experimental_variogram(bxy[:, 0], bxy[:, 1], sims_ns[0], lag, nl)
print(f"variogram reproduction (realisation 1, normal scores) at lag {h_s[1]:.0f}: {g_s[1]:.2f} vs model {spherical(h_s[1], *p_ns):.2f}")
p_above = (sims >= c).mean(axis=0)
gt_sims = np.array([[r["tonnes"] for r in grade_tonnage(s, cutoffs)] for s in sims])
print(f"tonnes above {c}: realisations {gt_sims[:, cutoffs == c].min():.0f}-{gt_sims[:, cutoffs == c].max():.0f}, "
      f"point truth {gt_pts.set_index('cutoff').loc[c, 'tonnes']:.0f}, block truth {gt_truth.set_index('cutoff').loc[c, 'tonnes']:.0f}")
note(f"""Each realisation has the full variance of its input histogram (sd ~{sims.std(axis=1).mean():.0f}
vs OK's {ok_point.std():.0f}), so it looks like a deposit rather than a smoothed map, and the
E-type average recovers a kriging-like map. The spread across realisations ({gt_sims[:, cutoffs == c].min():.0f}
to {gt_sims[:, cutoffs == c].max():.0f} tonnes above {c}) is the model's uncertainty, but note that
the whole band sits above the point truth ({gt_pts.set_index('cutoff').loc[c, 'tonnes']:.0f}): simulation
reproduces the declustered sample histogram faithfully, and that histogram is still richer
and wider than the deposit (mean {np.average(V, weights=wts):.0f} vs {truth.mean():.0f}). Simulation
quantifies uncertainty *given* the inputs; it cannot repair a biased histogram. Re-block the
realisations to SMU support before reading grade-tonnage for mine planning.""")

fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
for ax, (name, v, cm) in zip(axes, [("realisation 1", sims[0], "viridis"), ("realisation 2", sims[1], "viridis"),
                                     ("E-type (mean of 8)", etype, "viridis"),
                                     (f"P(V > {c})", p_above, "magma")]):
    im = ax.imshow(as_grid(v), origin="lower", cmap=cm, vmin=0, vmax=(vmax if cm == "viridis" else 1))
    ax.set_title(name)
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=.8)
fig.tight_layout()
savefig(fig, "sgs", SEC)
