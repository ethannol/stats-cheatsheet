# %% [markdown]
# # 09 Spatial statistics: autocorrelation, variograms, interpolation
#
# Data: 1991-2020 mean temperature at ~370 US stations.
#
# **Cheat sheet**
# - Tobler's law: near things are more related. Spatial autocorrelation violates the
#   independence assumption of every classical test, just like temporal autocorrelation.
# - Moran's I: global index of spatial autocorrelation on a weights matrix W (kNN, distance
#   band, contiguity). E[I] = -1/(n-1); significance by permutation. Local Moran / Getis-Ord
#   G* for hotspots (not shown).
# - Test the *residuals* of a model for spatial autocorrelation, not just the raw variable. If
#   residuals cluster, SEs are wrong and a spatial model (spatial lag/error, GWR, kriging of
#   residuals) is needed.
# - Variogram gamma(h): half the mean squared difference as a function of separation. Nugget
#   (measurement error + microscale), sill (total variance), range (decorrelation distance).
#   Fit spherical / exponential / Gaussian models to it.
# - Interpolation: IDW (deterministic, no uncertainty), ordinary kriging (BLUE given the
#   variogram, gives a variance), regression kriging / universal kriging (trend from
#   covariates + kriged residuals). Validate with leave-one-out, not R^2 at data points.
# - Use great-circle or projected distances; degrees are not distances. Watch anisotropy.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup
from statsheet.spatial import fit_interpolation_fold
from statsheet.stats import (
    distance_matrix_km,
    exponential,
    fit_variogram,
    haversine_km,
    idw_from_distances,
    knn_weights,
    morans_i,
    ordinary_kriging_from_distances,
    spherical,
    variogram_from_distances,
)

setup()
SEC = "09_spatial"
st = data.station_normals().dropna(subset=["tavg"]).reset_index(drop=True)
st["elev_km"] = st.elev / 1000
n = len(st)
D = distance_matrix_km(st.lat.values, st.lon.values)
print(f"{n} stations; median nearest-neighbour distance {np.sort(D, axis=1)[:, 1].mean():.0f} km")

# %% ----------------------------------------------------------- Moran's I
banner("Moran's I on the raw variable and on regression residuals")
trend = smf.ols("tavg ~ lat + elev_km + lon", data=st).fit()
st["resid"] = trend.resid
W8 = knn_weights(D, k=8)
for label, z in [
    ("tavg (raw)", st.tavg.values),
    ("OLS residual", st.resid.values),
    ("shuffled residual (null check)", np.random.default_rng(0).permutation(st.resid.values)),
]:
    mi = morans_i(z, W8, n_perm=999, rng=0)
    print(
        f"  {label:32s} I={mi['I']:+.3f}  E[I]={mi['EI']:+.3f}  z={mi['z']:+.1f}  p_perm={mi['p_perm']:.3f}"
    )
note("""Raw temperature is massively autocorrelated (of course: latitude). The interesting test
is on the residuals of the lat+elev+lon model: I is still clearly positive, meaning the
model's errors cluster (coastal strips, basins, the Gulf). The regression's standard errors
in section 03 were therefore too small, exactly as with time-series autocorrelation.""")

# %% ----------------------------------------------------------- variogram of residuals
banner("Empirical variogram of the residuals and model fits")
# fit only to ~half the domain width: pair counts thin out and large-scale trend leaks in
bins = np.arange(0, 1300, 75.0)
h, gam, cnt = variogram_from_distances(D, st.resid.values, bins)
sph = fit_variogram(h, gam, spherical, counts=cnt)
expo = fit_variogram(h, gam, exponential, counts=cnt)
print(f"spherical:   nugget={sph[0]:.3f} sill={sph[1]:.3f} range={sph[2]:.0f} km")
print(f"exponential: nugget={expo[0]:.3f} sill={expo[1]:.3f} range={expo[2]:.0f} km (practical)")
print(f"residual variance = {st.resid.var():.3f}")
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.scatter(h, gam, s=np.sqrt(cnt), alpha=0.7, label="empirical (size ~ pair count)")
hh = np.linspace(1, h.max(), 200)
ax.plot(hh, spherical(hh, *sph), label="spherical fit")
ax.plot(hh, exponential(hh, *expo), ls="--", label="exponential fit")
ax.axhline(st.resid.var(), color="k", lw=0.8, ls=":", label="sample variance")
ax.set_xlabel("separation (km)")
ax.set_ylabel("semivariance (C^2)")
ax.set_title("Residual variogram: nugget, rising limb, sill")
ax.legend()
savefig(fig, "variogram", SEC)
note(f"""The nugget ({sph[0]:.2f} C^2) is variance at zero distance: station siting, instrument
differences, microclimate. Semivariance rises toward the sill over ~{sph[2]:.0f} km, the
distance beyond which residuals are uncorrelated. The empirical curve keeps climbing past
the sample variance at large lags: a sign that the lat+elev+lon trend has not removed all
the large-scale structure (non-stationarity). Fit variograms only to roughly half the
domain, and treat an unbounded variogram as a prompt to improve the trend model.""")

# %% ----------------------------------------------------------- interpolation to a grid
banner("Interpolation: IDW vs regression + kriged residuals")
glon, glat = np.meshgrid(np.arange(-125, -66, 0.5), np.arange(24.5, 50, 0.5))
gl = pd.DataFrame({"lon": glon.ravel(), "lat": glat.ravel()})
Dq = haversine_km(
    gl.lat.values[:, None], gl.lon.values[:, None], st.lat.values[None, :], st.lon.values[None, :]
)
mask = Dq.min(axis=1) < 250  # only show grid cells near a station (drops oceans, Canada, Mexico)

idw_t = idw_from_distances(Dq, st.tavg.values, power=2)
# regression kriging: trend from covariates (need elevation on the grid; approximate with IDW of
# station elevation - a real analysis would use a DEM) + ordinary kriging of the residuals
gl["elev_km"] = idw_from_distances(Dq, st.elev_km.values, power=2)
trend_grid = trend.predict(gl)
krig_res, krig_var = ordinary_kriging_from_distances(D, Dq, st.resid.values, spherical, sph)
rk_t = trend_grid.values + krig_res


def grid(v):
    out = np.full(glon.shape, np.nan).ravel()
    out[mask] = v[mask]
    return out.reshape(glon.shape)


fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
for ax, field, title in [
    (axes[0], idw_t, "IDW (power 2) of station tavg"),
    (axes[1], rk_t, "Regression kriging: lat+elev+lon trend + OK residuals"),
]:
    im = ax.pcolormesh(glon, glat, grid(field), cmap="RdYlBu_r", vmin=0, vmax=25, shading="auto")
    ax.scatter(st.lon, st.lat, s=3, color="k", alpha=0.4)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="C")
im = axes[2].pcolormesh(
    glon, glat, grid(np.sqrt(np.maximum(krig_var, 0))), cmap="magma_r", shading="auto"
)
axes[2].scatter(st.lon, st.lat, s=3, color="w", alpha=0.6)
axes[2].set_title("Kriging standard error of the residual (low near stations)")
fig.colorbar(im, ax=axes[2], label="C")
for ax in axes:
    ax.set_aspect(1.3)
    ax.grid(False)
fig.tight_layout()
savefig(fig, "interpolation_maps", SEC)
note("""IDW produces bulls-eyes around stations and cannot extrapolate a trend (mountain cold
spots vanish between stations). Regression kriging uses the physics (lapse rate, latitude)
for the large-scale pattern and kriging for the spatially correlated leftovers, and it
hands you a map of its own uncertainty. The elevation grid here is IDW of station
elevations, which is crude; use a DEM in real work.""")

# %% ----------------------------------------------------------- leave-one-out validation
banner("Leave-one-out cross-validation RMSE (C)")
errs = {
    "IDW p=2": [],
    "IDW p=1": [],
    "OLS trend only": [],
    "OK on tavg": [],
    "regression kriging": [],
}
for i in range(n):
    keep = np.arange(n) != i
    Di = D[np.ix_(keep, keep)]
    Dqi = D[i, keep][None, :]
    errs["IDW p=2"].append(idw_from_distances(Dqi, st.tavg.values[keep], 2)[0] - st.tavg[i])
    errs["IDW p=1"].append(idw_from_distances(Dqi, st.tavg.values[keep], 1)[0] - st.tavg[i])
    f, tavg_sph_fold, resid_sph_fold = fit_interpolation_fold(st, D, keep, bins)
    pred_trend = f.predict(st.iloc[[i]]).iloc[0]
    errs["OLS trend only"].append(pred_trend - st.tavg[i])
    ok, _ = ordinary_kriging_from_distances(Di, Dqi, st.tavg.values[keep], spherical, tavg_sph_fold)
    errs["OK on tavg"].append(ok[0] - st.tavg[i])
    okr, _ = ordinary_kriging_from_distances(Di, Dqi, f.resid.values, spherical, resid_sph_fold)
    errs["regression kriging"].append(pred_trend + okr[0] - st.tavg[i])
loo = pd.DataFrame(
    {k: [np.sqrt(np.mean(np.square(v))), np.mean(np.abs(v))] for k, v in errs.items()},
    index=["RMSE", "MAE"],
).T.sort_values("RMSE")
print(loo.round(3).to_string())
note("""Each leave-one-out fold refits the trend and both variograms on training stations only.
Leave-one-out measures local interpolation; spatial block holdouts are needed to assess
transfer across regions. It is a useful yardstick for interpolators (fit at data points is trivially
perfect for IDW and kriging). Combining a covariate trend with kriged residuals beats
either ingredient alone; ordinary kriging on raw temperature struggles because the field
is not stationary (the mean trends with latitude), which the variogram assumes away.""")
