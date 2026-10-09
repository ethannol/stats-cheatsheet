# %% [markdown]
# # 08 Multivariate analysis: correlation structure, PCA, clustering
#
# Data: 1991-2020 climatologies for ~370 US airport stations, 15 variables each.
#
# **Cheat sheet**
# - Standardise before PCA/k-means unless the variables share units and you want the
#   variance weighting (degree-days would otherwise dominate everything).
# - PCA: rotate to orthogonal axes of maximum variance. Loadings tell you what each component
#   means; scores place observations. Scree plot / cumulative variance / Kaiser (eigenvalue
#   > 1) for how many to keep. Signs of components are arbitrary. PCA is linear and sensitive
#   to outliers; it is EOF analysis when applied to a space x time field.
# - k-means: minimises within-cluster variance, assumes roughly spherical equal-size clusters,
#   needs k. Choose k with the elbow (inertia) and silhouette; run several inits (n_init).
# - Hierarchical (Ward) gives a dendrogram, no k needed up front; cut it at the height you can
#   defend. Compare partitions with the adjusted Rand index.
# - Gaussian mixtures when clusters overlap and you want soft membership; DBSCAN for
#   arbitrary shapes / noise points.
# - Clusters are descriptive. Validate against something external (here: geography, Koppen).

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.preprocessing import StandardScaler

from earthstats import data
from earthstats.plotting import COLORS, banner, note, savefig, setup

setup()
SEC = "08_multivariate"

st = data.station_normals()
VARS = ["tavg", "tmax", "tmin", "prcp", "snow", "dp01", "dp10", "dx90", "dx32", "dt32", "dt00",
        "emxt", "emnt", "htdd", "cldd"]
st = st.dropna(subset=VARS).reset_index(drop=True)
X = st[VARS].to_numpy()
Z = StandardScaler().fit_transform(X)
print(f"{len(st)} stations x {len(VARS)} climate variables")

# %% ----------------------------------------------------------- correlation structure
banner("Correlation matrix")
C = pd.DataFrame(np.corrcoef(Z.T), index=VARS, columns=VARS)
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(len(VARS)))
ax.set_xticklabels(VARS, rotation=90)
ax.set_yticks(range(len(VARS)))
ax.set_yticklabels(VARS)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=.8, label="Pearson r")
ax.set_title("Station climatology variables are highly redundant")
savefig(fig, "correlation_matrix", SEC)
pairs = C.where(np.triu(np.ones_like(C, dtype=bool), 1)).stack()
print("strongest |r| pairs:\n" + pairs.reindex(pairs.abs().sort_values(ascending=False).index)
      .head(6).round(3).to_string())
note("""Fifteen variables, but a handful of underlying factors: 'how warm' (tavg, htdd, cldd,
dt32...), 'how wet' (prcp, dp01, dp10), 'how snowy / continental'. PCA will find these.""")

# %% ----------------------------------------------------------- PCA
banner("PCA")
pca = PCA().fit(Z)
ev = pca.explained_variance_ratio_
print("explained variance by component:", np.round(ev[:6], 3), "... cumulative(3) =", round(ev[:3].sum(), 3))
print("Kaiser rule (eigenvalue > 1): keep", int((pca.explained_variance_ > 1).sum()), "components")
load = pd.DataFrame(pca.components_[:3].T, index=VARS, columns=["PC1", "PC2", "PC3"])
print("\nloadings:\n" + load.round(2).to_string())
scores = pca.transform(Z)
st["PC1"], st["PC2"], st["PC3"] = scores[:, 0], scores[:, 1], scores[:, 2]
note("""PC1 is a temperature axis (all warm variables load with one sign, all cold variables
with the other), PC2 is a wetness axis (prcp, dp01, dp10), PC3 is continentality: many hot
days AND many very cold days AND high cooling degree-days load together, which is the
large annual range of the interior. Three components carry >90% of the variance in
fifteen variables. Interpretation is by loadings, and the signs flip freely between
software: always check which direction is 'warm'.""")

fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
axes[0].bar(range(1, 9), ev[:8])
axes[0].plot(range(1, 9), np.cumsum(ev[:8]), "o-", color="C3")
axes[0].set_xlabel("component")
axes[0].set_title("Scree plot and cumulative variance")
sc = axes[1].scatter(st.PC1, st.PC2, c=st.lat, cmap="viridis", s=14)
for v in VARS:
    axes[1].annotate(v, (load.loc[v, "PC1"] * 6, load.loc[v, "PC2"] * 6), fontsize=7, color="C3")
    axes[1].arrow(0, 0, load.loc[v, "PC1"] * 5.5, load.loc[v, "PC2"] * 5.5, color="C3", alpha=.4,
                  head_width=.1)
axes[1].set_xlabel(f"PC1 ({ev[0]:.0%})")
axes[1].set_ylabel(f"PC2 ({ev[1]:.0%})")
axes[1].set_title("Biplot: scores (coloured by latitude) + loadings")
fig.colorbar(sc, ax=axes[1], label="latitude")
sc2 = axes[2].scatter(st.lon, st.lat, c=st.PC1, cmap="RdBu_r", s=16)
axes[2].set_title("PC1 score mapped: the temperature gradient")
axes[2].set_xlabel("lon")
axes[2].set_ylabel("lat")
fig.colorbar(sc2, ax=axes[2], label="PC1")
fig.tight_layout()
savefig(fig, "pca", SEC)

# %% ----------------------------------------------------------- choosing k
banner("k-means: choosing k")
ks = range(2, 11)
inertia, sil = [], []
for k in ks:
    km = KMeans(k, n_init=10, random_state=0).fit(Z)
    inertia.append(km.inertia_)
    sil.append(silhouette_score(Z, km.labels_))
for k, i, s_ in zip(ks, inertia, sil):
    print(f"  k={k:2d} inertia={i:8.1f} silhouette={s_:.3f}")
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
axes[0].plot(list(ks), inertia, "o-")
axes[0].set_xlabel("k")
axes[0].set_title("Elbow: within-cluster sum of squares")
axes[1].plot(list(ks), sil, "o-")
axes[1].set_xlabel("k")
axes[1].set_title("Silhouette (higher = tighter, better separated)")
savefig(fig, "choose_k", SEC)
note("""Neither criterion gives a sharp answer, which is normal for a continuous climate
gradient: there are no 'true' clusters, only a useful discretisation. The silhouette is
highest at small k (the warm/cold split), the elbow suggests 4-6. Pick k by what the
clusters will be used for and report that the choice is a judgement.""")

# %% ----------------------------------------------------------- k=5 regimes, hierarchical, ARI
K = 5
km = KMeans(K, n_init=20, random_state=0).fit(Z)
st["kmeans"] = km.labels_
Lk = linkage(Z, method="ward")
st["ward"] = fcluster(Lk, K, criterion="maxclust") - 1
ari = adjusted_rand_score(st.kmeans, st.ward)
print(f"\nk={K}: adjusted Rand index between k-means and Ward partitions = {ari:.2f}")
prof = st.groupby("kmeans")[["tavg", "prcp", "snow", "dx90", "dt32", "lat", "elev"]].mean().round(1)
prof["n"] = st.kmeans.value_counts().sort_index()
# order clusters warm -> cold for readability
order = prof.sort_values("tavg", ascending=False).index
relabel = {old: new for new, old in enumerate(order)}
st["kmeans"] = st.kmeans.map(relabel)
prof = prof.loc[order].reset_index(drop=True)
prof.index.name = "cluster"
print("\ncluster profiles (means):\n" + prof.to_string())

fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1.3, 1]))
for c in range(K):
    m = st.kmeans == c
    axes[0].scatter(st.lon[m], st.lat[m], s=18, color=COLORS[c],
                    label=f"{c}: tavg {prof.loc[c, 'tavg']:.0f}C, prcp {prof.loc[c, 'prcp']:.0f}mm")
axes[0].legend(fontsize=8, loc="lower left")
axes[0].set_title(f"k-means (k={K}) on 15 standardised climate variables")
axes[0].set_xlabel("lon")
axes[0].set_ylabel("lat")
dendrogram(Lk, ax=axes[1], no_labels=True, color_threshold=Lk[-K + 1, 2])
axes[1].set_title("Ward dendrogram (cut at 5 clusters)")
axes[1].set_ylabel("Ward distance")
fig.tight_layout()
savefig(fig, "clusters_map", SEC)
note(f"""Unsupervised clusters recover a recognisable climate map: warm-dry (Southwest), warm-wet
(Southeast), cool-wet with moderate snow (Northeast, Midwest, Pacific Northwest), cool-dry
at high elevation (interior West) and cold-snowy (northern tier). Two different algorithms
agree substantially (ARI={ari:.2f}) because the structure is real; where they disagree is
along the gradients between regimes.""")
