# %% [markdown]
# # 01 Exploratory data analysis & descriptive statistics
#
# Data: Central Park (NYC) daily weather 1950-present; Mauna Loa CO2.
#
# **Cheat sheet**
# - Always start with: shape, dtypes, missingness pattern, time coverage, units. Plot the raw
#   series before any statistic.
# - Location: mean (sensitive to outliers) vs median (robust). Spread: SD vs IQR vs MAD
#   (MAD * 1.4826 ~ SD for normal data). Shape: skewness, kurtosis.
# - Precipitation is zero-inflated and right-skewed: summarise wet days separately, log or
#   gamma transforms, never report mean +- SD alone.
# - Climatology + anomaly: subtract the seasonal cycle (day-of-year mean, smoothed) before
#   looking for trends or correlations, otherwise seasonality dominates everything.
# - Correlation: Pearson (linear, sensitive to outliers) vs Spearman (rank, monotonic) vs
#   Kendall (pairwise concordance, small n). Seasonality can manufacture correlation.
# - QQ plot beats a normality test for n > ~100: tests reject trivially at large n.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup

setup()
SEC = "01_eda"

cp = data.ghcn_daily()
co2 = data.co2()

# %% ----------------------------------------------------------- data audit
banner("Data audit: Central Park GHCN-Daily")
print(cp.describe().T.round(2))
miss = cp.isna().mean().round(4)
print("\nFraction missing per column:\n", miss.to_string())
gaps = cp.TMAX.isna()
runs = (gaps != gaps.shift()).cumsum()
longest_gap = gaps.groupby(runs).sum().max()
print(f"\nLongest run of missing TMAX days: {longest_gap}")
note("""Missingness is the first thing to characterise: is it random (MCAR), related to
observed covariates (MAR), or related to the missing value itself (MNAR, e.g. a gauge
that fails during floods)? Reindexing to a complete calendar made the gaps visible.""")

# %% ----------------------------------------------------------- central tendency & spread
banner("Location, spread and shape: TMAX vs PRCP")
tmax = cp.TMAX.dropna()
prcp = cp.PRCP.dropna()
wet = prcp[prcp > 0]


def summary(x: pd.Series) -> pd.Series:
    mad = sps.median_abs_deviation(x, scale="normal")  # scaled so MAD ~ SD under normality
    return pd.Series({
        "n": len(x), "mean": x.mean(), "median": x.median(), "sd": x.std(),
        "iqr": x.quantile(.75) - x.quantile(.25), "mad_scaled": mad,
        "skew": sps.skew(x), "excess_kurt": sps.kurtosis(x),
        "p05": x.quantile(.05), "p95": x.quantile(.95), "max": x.max()})


tab = pd.DataFrame({"TMAX (C)": summary(tmax), "PRCP all days (mm)": summary(prcp),
                    "PRCP wet days (mm)": summary(wet), "log PRCP wet (ln mm)": summary(np.log(wet))})
print(tab.round(2).to_string())
print(f"\nFraction of dry days: {(prcp == 0).mean():.3f}")
note("""TMAX: mean ~ median, mild negative skew, platykurtic (it is a mixture of seasons, so
flatter than normal). PRCP: mean >> median, huge positive skew, kurtosis in the hundreds.
On wet days the log transform brings skew near zero: lognormal / gamma are the usual
precipitation models. Report medians and quantiles for skewed variables.""")

# %%
fig, axes = plt.subplots(2, 2, figsize=(11, 7))
axes[0, 0].hist(tmax, bins=60, density=True, alpha=.8)
axes[0, 0].set_title("TMAX daily: bimodal-ish mixture of seasons")
axes[0, 0].set_xlabel("deg C")
axes[0, 1].hist(wet, bins=80, density=True, alpha=.8)
axes[0, 1].set_yscale("log")
axes[0, 1].set_title("Wet-day precipitation (log y): heavy right tail")
axes[0, 1].set_xlabel("mm")
sps.probplot(tmax.sample(3000, random_state=0), dist="norm", plot=axes[1, 0])
axes[1, 0].set_title("QQ plot TMAX vs normal (light tails)")
sps.probplot(np.log(wet.sample(3000, random_state=0)), dist="norm", plot=axes[1, 1])
axes[1, 1].set_title("QQ plot ln(wet-day PRCP) vs normal")
fig.tight_layout()
savefig(fig, "distributions", SEC)

# %% ----------------------------------------------------------- normality tests at large n
banner("Normality tests at large n reject everything")
for name, x in [("TMAX (n=all)", tmax), ("TMAX (n=50 sample)", tmax.sample(50, random_state=1)),
                ("ln wet PRCP (n=all)", np.log(wet))]:
    sw = sps.shapiro(x.sample(min(len(x), 5000), random_state=0))
    da = sps.normaltest(x)  # D'Agostino-Pearson: skew + kurtosis
    print(f"{name:22s} Shapiro W={sw.statistic:.4f} p={sw.pvalue:.2e} | "
          f"D'Agostino K2={da.statistic:.1f} p={da.pvalue:.2e}")
note("""With tens of thousands of points even tiny departures from normality are 'significant'.
Use QQ plots and ask whether the departure matters for the method you plan to use
(the CLT protects means; it does not protect extreme quantiles).""")

# %% ----------------------------------------------------------- climatology & anomalies
banner("Climatology and anomalies")
doy = cp.index.dayofyear
clim = cp.TMAX.groupby(doy).mean()
# smooth the day-of-year climatology with a circular 31-day window (Fourier is the alternative)
clim_s = pd.concat([clim.iloc[-15:], clim, clim.iloc[:15]]).rolling(31, center=True).mean().iloc[15:-15]
anom = cp.TMAX - clim_s.reindex(doy).to_numpy()
cp["TMAX_anom"] = anom

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(clim.index, clim.values, lw=.6, alpha=.6, label="raw day-of-year mean")
axes[0].plot(clim_s.index, clim_s.values, lw=2, label="31-day smoothed climatology")
axes[0].set_xlabel("day of year")
axes[0].set_ylabel("TMAX (C)")
axes[0].legend()
axes[0].set_title("Seasonal climatology")
ann = anom.resample("YE").mean()
axes[1].plot(ann.index, ann.values, marker="o", ms=3)
axes[1].axhline(0, color="k", lw=.8)
axes[1].set_title("Annual mean TMAX anomaly (vs 1950-present climatology)")
axes[1].set_ylabel("deg C")
savefig(fig, "climatology_anomaly", SEC)
print(f"Autocorrelation of daily TMAX, lag 1: {cp.TMAX.autocorr(1):.3f}  "
      f"| of daily anomaly, lag 1: {anom.autocorr(1):.3f}")
note(f"""Raw daily TMAX has a lag-1 autocorrelation of {cp.TMAX.autocorr(1):.2f}, much of it from
the seasonal cycle. Anomalies still have {anom.autocorr(1):.2f}, which is real weather
persistence. Either way, daily values are not independent: n is far smaller than the row
count for inference.""")

# %% ----------------------------------------------------------- box plots by month
fig, ax = plt.subplots(figsize=(10, 4))
cp["month"] = cp.index.month
cp.boxplot(column="PRCP", by="month", ax=ax, showfliers=False)
ax.set_title("Daily precipitation by month (fliers hidden)")
ax.set_ylabel("mm")
fig.suptitle("")
savefig(fig, "prcp_by_month", SEC)

# %% ----------------------------------------------------------- correlation
banner("Correlation: raw vs anomaly, Pearson vs Spearman")
raw = cp[["TMAX", "TMIN", "PRCP"]].dropna()
print("Pearson on raw daily values:\n", raw.corr().round(3))
print("\nSpearman on raw daily values:\n", raw.corr(method="spearman").round(3))
cp["TMIN_anom"] = cp.TMIN - cp.TMIN.groupby(doy).transform("mean")
an = cp[["TMAX_anom", "TMIN_anom", "PRCP"]].dropna()
print("\nPearson on anomalies:\n", an.corr().round(3))
r, p = sps.pearsonr(an.TMAX_anom, an.PRCP)
rho, p2 = sps.spearmanr(an.TMAX_anom, an.PRCP)
print(f"\nTMAX anomaly vs PRCP: Pearson r={r:.3f}, Spearman rho={rho:.3f}")
note("""TMAX-TMIN correlation is 0.95 raw but lower on anomalies: the shared seasonal cycle
inflates it. Rain days are cooler than average (negative TMAX-anomaly correlation),
which the raw Pearson correlation misses because rain is uncorrelated with season here.
Spearman is more negative than Pearson because precipitation is heavy-tailed.""")

# %% ----------------------------------------------------------- CO2 first look
banner("Mauna Loa CO2: first look")
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].plot(co2.index, co2.co2, lw=.8, label="monthly mean")
axes[0].plot(co2.index, co2.deseasonalized, lw=1.2, label="deseasonalised (NOAA)")
axes[0].set_ylabel("ppm")
axes[0].legend()
axes[0].set_title("Trend + seasonal cycle + accelerating growth")
growth = co2.deseasonalized.diff(12).dropna()
axes[1].plot(growth.index, growth.values, lw=.8)
axes[1].plot(growth.rolling(60, center=True).mean(), lw=2)
axes[1].set_ylabel("ppm / yr")
axes[1].set_title("Year-over-year growth rate (and 5-yr mean)")
savefig(fig, "co2_overview", SEC)
print(co2.co2.describe().round(2).to_string())
print(f"First year mean {co2.co2[:'1959'].mean():.1f} ppm, last 12 months mean "
      f"{co2.co2.iloc[-12:].mean():.1f} ppm")
note("""A smooth monotonic trend with a stable annual cycle (~6 ppm amplitude from Northern
Hemisphere photosynthesis) and a growth rate that itself rises from ~0.8 to ~2.5 ppm/yr.
The growth-rate series is where the interesting variability (ENSO) lives; see 04.""")
