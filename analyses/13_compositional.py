# %% [markdown]
# # 13 Compositional data analysis (CoDA) for geochemistry
#
# Data: Jura soil metals (Cd, Co, Cr, Cu, Ni, Pb, Zn in ppm) with rock type.
#
# **Cheat sheet**
# - Assays are parts of a whole (ppm, %, wt%). The parts are constrained to sum to a constant
#   (closure), so they cannot vary independently: at least one covariance must be negative.
#   Correlations between raw percentages are partly artefacts (Pearson's 'spurious correlation').
# - Closure also means sub-compositions have different correlations than the full composition
#   (subcompositional incoherence), and Euclidean distance between compositions is meaningless.
# - Log-ratio transforms map the simplex to real space where ordinary statistics apply:
#   alr = log(x_i/x_D) (D-1 coords, asymmetric, depends on the reference part);
#   clr = log(x_i / geometric mean) (D coords summing to zero, singular covariance, good for biplots);
#   ilr = orthonormal rotation of clr (D-1 coords, use for regression, clustering, kriging).
# - The variation matrix var(log(x_i/x_j)) replaces the correlation matrix: near-zero entries mean
#   proportional parts.
# - Zeros/non-detects must be replaced before taking logs (multiplicative replacement at 65% DL).
# - Even when the total is far from 100% (trace elements), the ratios are what carry the geology.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import COLORS, banner, note, savefig, setup
from statsheet.stats import closure, clr, ilr, variation_matrix

setup()
SEC = "13_compositional"
j = data.jura()
METALS = ["cd", "co", "cr", "cu", "ni", "pb", "zn"]
X = j[METALS].to_numpy()
print(f"{len(j)} samples, 7 metals; rock types: {j.rock.value_counts().to_dict()}")

# %% ----------------------------------------------------------- spurious correlation demo
banner("Closure manufactures correlation: a simulation")
rng = np.random.default_rng(0)
A = rng.lognormal([3, 3, 3], 0.6, (2000, 3))  # three *independent* amounts
C = closure(A, 100)  # expressed as percentages
r_raw = np.corrcoef(A.T)[np.triu_indices(3, 1)]
r_closed = np.corrcoef(C.T)[np.triu_indices(3, 1)]
print(f"independent amounts: pairwise r = {np.round(r_raw, 2)}")
print(f"after closing to 100%: pairwise r = {np.round(r_closed, 2)}")
print(
    f"clr of closed data:     r = {np.round(np.corrcoef(clr(C).T)[np.triu_indices(3, 1)], 2)} "
    f"(clr correlations are not interpretable as independence either: rows sum to zero)"
)
note("""Three unrelated quantities become strongly negatively correlated the moment they are
expressed as percentages of a total. Any correlation matrix of raw assay percentages contains
this artefact. The fix is not a different correlation coefficient: it is to analyse ratios.""")

# %% ----------------------------------------------------------- Jura: raw vs log-ratio
banner("Jura metals: raw correlation vs variation matrix")
raw_corr = pd.DataFrame(np.corrcoef(X.T), index=METALS, columns=METALS)
T = pd.DataFrame(variation_matrix(X), index=METALS, columns=METALS)
print("Pearson r on raw ppm:\n" + raw_corr.round(2).to_string())
print("\nVariation matrix var(ln(xi/xj)) (small = proportional):\n" + T.round(2).to_string())
pairs = T.where(np.triu(np.ones_like(T, dtype=bool), 1)).stack().dropna().sort_values()
print(
    "\nmost proportional pairs:",
    ", ".join(f"{a}/{b} ({v:.2f})" for (a, b), v in pairs.head(3).items()),
)
print(
    "least proportional pairs:",
    ", ".join(f"{a}/{b} ({v:.2f})" for (a, b), v in pairs.tail(3).items()),
)
note(f"""In the Jura data the metals are trace elements (sum ~0.03% of the soil), so closure to
the whole is weak and raw correlations are not dominated by the artefact. The variation matrix
still adds information: {pairs.index[0][0]}/{pairs.index[0][1]} is the most proportional pair
(co-located in the same minerals), while {pairs.index[-1][0]}/{pairs.index[-1][1]} vary most
independently. Log-ratio thinking is about geology (which elements travel together), not
just about avoiding an artefact.""")

# %% ----------------------------------------------------------- clr biplot
banner("Compositional biplot (PCA of clr coordinates)")
Z = clr(X)
Zc = Z - Z.mean(axis=0)
U, S, Vt = np.linalg.svd(Zc, full_matrices=False)
expl = S**2 / np.sum(S**2)
scores = U[:, :2] * S[:2]
loads = Vt[:2].T
print(
    f"explained variance PC1={expl[0]:.2f} PC2={expl[1]:.2f} (clr covariance has rank D-1={len(METALS) - 1})"
)
rock_codes = {r: i for i, r in enumerate(sorted(j.rock.unique()))}
fig, ax = plt.subplots(figsize=(8, 7))
for r, i in rock_codes.items():
    m = (j.rock == r).values
    ax.scatter(scores[m, 0], scores[m, 1], s=14, color=COLORS[i], label=r, alpha=0.7)
for k, m in enumerate(METALS):
    ax.arrow(0, 0, loads[k, 0] * 3, loads[k, 1] * 3, color="k", head_width=0.08, alpha=0.7)
    ax.text(loads[k, 0] * 3.3, loads[k, 1] * 3.3, m.capitalize(), fontsize=10, weight="bold")
ax.set_xlabel(f"PC1 ({expl[0]:.0%})")
ax.set_ylabel(f"PC2 ({expl[1]:.0%})")
ax.set_title("clr biplot: ray links ~ log-ratio variance; short link = proportional pair")
ax.legend(fontsize=8)
savefig(fig, "clr_biplot", SEC)
note("""In a compositional biplot the *links* between ray tips matter, not the rays: a short
link between two elements means their log-ratio has low variance (they are proportional).
Rock types separate along the log-ratio axes, which is the geochemical fingerprint a mining
geologist uses for lithology classification from assays.""")

# %% ----------------------------------------------------------- ilr for standard statistics
banner("ilr coordinates: do rock types differ? (MANOVA-style test via per-coordinate ANOVA)")
Yilr = ilr(X)
groups = [Yilr[(j.rock == r).values] for r in rock_codes]
rows = []
for k in range(Yilr.shape[1]):
    f, p = sps.f_oneway(*[g[:, k] for g in groups])
    rows.append(dict(ilr_coord=k + 1, F=f, p=p))
print(pd.DataFrame(rows).round(4).to_string(index=False))
# Aitchison vs Euclidean distance between rock-type centres
cent = {r: X[(j.rock == r).values].mean(axis=0) for r in rock_codes}
geo_cent = {r: np.exp(np.log(X[(j.rock == r).values]).mean(axis=0)) for r in rock_codes}
r1, r2 = list(rock_codes)[:2]
d_euc = np.linalg.norm(cent[r1] - cent[r2])
d_ait = np.linalg.norm(clr(geo_cent[r1][None]) - clr(geo_cent[r2][None]))
print(
    f"\n{r1} vs {r2}: Euclidean distance between mean compositions={d_euc:.1f} ppm (dominated by Zn), "
    f"Aitchison distance={d_ait:.2f} (scale-free)"
)
note("""ilr coordinates are ordinary real numbers: feed them to ANOVA, regression, k-means,
random forests or kriging, then back-transform results to the simplex. The compositional
mean is the geometric mean (closed), and the right distance is Aitchison's, which does not let
the most abundant element (Zn here) dominate the way Euclidean distance on ppm does.""")

# %% ----------------------------------------------------------- ratio maps
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
for ax, (num, den) in zip(axes, [("cu", "zn"), ("pb", "zn")]):
    lr = np.log(j[num] / j[den])
    sc = ax.scatter(
        j.x, j.y, c=lr, cmap="RdBu_r", s=22, vmin=-np.abs(lr).max(), vmax=np.abs(lr).max()
    )
    ax.set_aspect("equal")
    ax.set_title(f"ln({num.capitalize()}/{den.capitalize()}): spatial pattern of a log-ratio")
    fig.colorbar(sc, ax=ax)
fig.tight_layout()
savefig(fig, "logratio_maps", SEC)
