# %% [markdown]
# # 04 Time series: decomposition, stationarity, ARIMA, spectral analysis
#
# Data: Mauna Loa monthly CO2 (trend + strong seasonality); GISTEMP for spectral analysis.
#
# **Cheat sheet**
# - Decompose first: trend + seasonal + remainder. STL (loess-based) handles changing
#   seasonality; classical decomposition assumes a fixed seasonal pattern.
# - Stationarity: constant mean/variance/autocovariance. ADF tests H0 = unit root (non-
#   stationary); KPSS tests H0 = stationary. Use both: ADF reject + KPSS not reject = stationary.
# - Differencing removes trend (d=1) or quadratic trend (d=2); seasonal differencing (lag s)
#   removes seasonality. Over-differencing shows up as a strong negative lag-1 ACF.
# - ACF tails off slowly -> AR; cuts off after q lags -> MA. PACF cuts off after p lags -> AR(p).
# - SARIMA(p,d,q)(P,D,Q)s. Check residuals: white noise (Ljung-Box), normal QQ.
# - Forecast evaluation: hold out the end (never shuffle), compare to naive and seasonal-naive
#   baselines, report MAE/RMSE/MAPE and interval coverage.
# - Spectral: periodogram / Welch PSD finds periodicities; detrend first or the trend
#   dominates low frequencies. Red noise (AR(1)) also has power at low frequencies.
# - Cross-correlation at lags shows lead/lag relationships; pre-whiten or you will find
#   spurious correlation between two trending series.

# %%
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.tsa.holtwinters import ExponentialSmoothing
from statsmodels.tsa.seasonal import STL, seasonal_decompose
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.stattools import adfuller, kpss

from statsheet import data
from statsheet.plotting import banner, note, savefig, setup

warnings.filterwarnings("ignore")
setup()
SEC = "04_time_series"

co2 = data.co2().co2.asfreq("MS").interpolate()  # 7 missing months in 1958-64
co2 = co2[:"2025-12"]

# %% ----------------------------------------------------------- decomposition
banner("Decomposition: classical vs STL")
classical = seasonal_decompose(co2, model="additive", period=12)
stl = STL(co2, period=12, robust=True).fit()
fig = stl.plot()
fig.set_size_inches(10, 8)
fig.suptitle("STL decomposition of monthly CO2", y=1.0)
savefig(fig, "stl", SEC)

seas_amp = stl.seasonal.groupby(stl.seasonal.index.year).agg(lambda s: s.max() - s.min())
print(
    f"Seasonal amplitude (peak-to-trough): 1960s {seas_amp.loc[1960:1969].mean():.2f} ppm, "
    f"2010s {seas_amp.loc[2010:2019].mean():.2f} ppm"
)
print(f"Remainder SD: classical {classical.resid.std():.3f}, STL {stl.resid.std():.3f}")
note("""STL lets the seasonal component evolve: the CO2 seasonal amplitude has grown ~15%
since the 1960s (longer growing seasons, more vigorous Northern Hemisphere vegetation).
Classical decomposition would average that away into the remainder.""")

# %% ----------------------------------------------------------- stationarity
banner("Stationarity tests on CO2 and its differences")


def stationarity(x, label):
    adf = adfuller(x.dropna(), autolag="AIC")
    kp = kpss(x.dropna(), regression="c", nlags="auto")
    verdict = (
        "stationary"
        if (adf[1] < 0.05 and kp[1] > 0.05)
        else "non-stationary"
        if (adf[1] >= 0.05 and kp[1] <= 0.05)
        else "ambiguous"
    )
    print(f"{label:34s} ADF p={adf[1]:.3f}  KPSS p={kp[1]:.3f}  -> {verdict}")


stationarity(co2, "CO2 level")
stationarity(co2.diff(), "first difference")
stationarity(co2.diff(12), "seasonal difference (lag 12)")
stationarity(co2.diff().diff(12), "first + seasonal difference")
stationarity(stl.resid, "STL remainder")
note("""The level is obviously non-stationary. A first difference still has the seasonal
cycle (ADF sees stationarity in the mean; KPSS is suspicious); the seasonal difference
still carries the growth trend; both together give white-ish noise. That is the (d=1, D=1)
of the SARIMA below. KPSS has an awkward p-value range because it is table-based.""")

# %% ----------------------------------------------------------- ACF / PACF
fig, axes = plt.subplots(2, 2, figsize=(12, 6))
plot_acf(co2, lags=48, ax=axes[0, 0], title="ACF: CO2 level (slow decay = trend)")
plot_acf(
    co2.diff().dropna(),
    lags=48,
    ax=axes[0, 1],
    title="ACF: first difference (seasonal spikes at 12, 24, 36)",
)
d12 = co2.diff().diff(12).dropna()
plot_acf(d12, lags=48, ax=axes[1, 0], title="ACF: (1-B)(1-B^12) CO2")
plot_pacf(d12, lags=48, ax=axes[1, 1], title="PACF: (1-B)(1-B^12) CO2", method="ywm")
fig.tight_layout()
savefig(fig, "acf_pacf", SEC)
note("""After double differencing, the ACF has a negative spike at lag 1 and at lag 12 and
little else: the signature of an MA(1) x seasonal MA(1), i.e. the classic 'airline model'
SARIMA(0,1,1)(0,1,1)12. Let the data suggest the order, then confirm with AIC.""")

# %% ----------------------------------------------------------- SARIMA with holdout
banner("SARIMA fit, holdout forecast and baselines")
train, test = co2[:"2019-12"], co2["2020-01":]
candidates = {
    "(0,1,1)(0,1,1)12 airline": ((0, 1, 1), (0, 1, 1, 12)),
    "(1,1,1)(1,1,1)12": ((1, 1, 1), (1, 1, 1, 12)),
    "(2,1,0)(0,1,1)12": ((2, 1, 0), (0, 1, 1, 12)),
}
fits = {}
for name, (order, sorder) in candidates.items():
    fits[name] = SARIMAX(train, order=order, seasonal_order=sorder).fit(disp=False)
    print(f"  {name:28s} AIC={fits[name].aic:8.1f}  BIC={fits[name].bic:8.1f}")
best_name = min(fits, key=lambda k: fits[k].aic)
best = fits[best_name]
print(f"best by AIC: {best_name}")
print(best.summary().tables[1])
lb = acorr_ljungbox(best.resid[13:], lags=[12, 24], return_df=True)
print(f"Ljung-Box on residuals p(12)={lb.lb_pvalue.iloc[0]:.3f} p(24)={lb.lb_pvalue.iloc[1]:.3f}")

h = len(test)
fc = best.get_forecast(h)
mean, ci = fc.predicted_mean, fc.conf_int(alpha=0.05)
# baselines
seasonal_naive = pd.Series(
    [
        train.iloc[-12 + (i % 12)]
        + (i // 12 + 1) * (train.iloc[-12:].mean() - train.iloc[-24:-12].mean())
        for i in range(h)
    ],
    index=test.index,
)  # seasonal naive with drift
hw = ExponentialSmoothing(
    train, trend="add", seasonal="add", seasonal_periods=12, damped_trend=False
).fit()
hw_fc = hw.forecast(h)


def metrics(y, yhat):
    e = y - yhat
    return dict(
        MAE=np.abs(e).mean(),
        RMSE=np.sqrt((e**2).mean()),
        MAPE=100 * np.abs(e / y).mean(),
        bias=e.mean(),
    )


res = pd.DataFrame(
    {
        "SARIMA": metrics(test, mean),
        "Holt-Winters": metrics(test, hw_fc),
        "Seasonal naive + drift": metrics(test, seasonal_naive),
    }
).T
print("\nHoldout 2020-2025 (72 months):\n" + res.round(3).to_string())
cover = ((test >= ci.iloc[:, 0]) & (test <= ci.iloc[:, 1])).mean()
print(f"SARIMA 95% interval coverage on holdout: {cover:.2f}")
note(f"""A 6-year horizon is long for a monthly model. The positive bias ({res.loc["SARIMA", "bias"]:.2f}
ppm) says the truth ran above the forecast: CO2 growth accelerated after 2020 (El Nino
2023-24), which a model fit through 2019 cannot know. Interval coverage is {cover:.2f}: the
intervals of a twice-integrated model widen quickly with horizon, so they are wide enough
here; at shorter horizons ARIMA intervals are usually too narrow because they ignore
parameter and model uncertainty. Always report a trivial baseline: SARIMA roughly halves
the seasonal-naive error, which is the honest statement of what the model buys you.""")

fig, ax = plt.subplots(figsize=(11, 4.5))
ax.plot(co2["2012":], lw=1, label="observed")
ax.plot(mean, lw=1.5, label=f"SARIMA {best_name}")
ax.fill_between(ci.index, ci.iloc[:, 0], ci.iloc[:, 1], alpha=0.2)
ax.plot(hw_fc, lw=1, ls="--", label="Holt-Winters")
ax.plot(seasonal_naive, lw=1, ls=":", label="seasonal naive + drift")
ax.axvline(test.index[0], color="k", lw=0.8)
ax.set_ylabel("ppm")
ax.set_title("Out-of-sample forecast, trained through 2019")
ax.legend()
savefig(fig, "sarima_forecast", SEC)

# %% ----------------------------------------------------------- rolling-origin backtest
banner("Rolling-origin (expanding window) backtest, 12-month horizon")
origins = pd.date_range("2005-12-01", "2024-12-01", freq="24MS")
errs = []
for o in origins:
    tr = co2[:o]
    te = co2[o + pd.offsets.MonthBegin(1) :][:12]
    f = SARIMAX(tr, order=(0, 1, 1), seasonal_order=(0, 1, 1, 12)).fit(disp=False)
    errs.append(np.sqrt(((f.forecast(12) - te) ** 2).mean()))
print(
    "12-month RMSE by origin year: "
    + ", ".join(f"{y}: {e:.2f}" for y, e in zip(origins.year, errs))
)
print(f"mean RMSE={np.mean(errs):.2f} ppm")
note("""Time-series cross-validation = refit at successive origins and score the next h steps.
Never random K-fold on time series: the model would see the future.""")

# %% ----------------------------------------------------------- spectral analysis
banner("Spectral analysis: where is the variance in frequency?")
gt = data.gistemp().anomaly.asfreq("MS").interpolate()[:"2025-12"]
gt_detr = signal.detrend(gt.values)
f, P = signal.welch(gt_detr, fs=12, nperseg=360)  # cycles per year
# a linear detrend leaves the quadratic growth in; remove a 2nd-order polynomial instead
tt = np.arange(len(co2))
co2_detr = co2.values - np.polyval(np.polyfit(tt, co2.values, 2), tt)
f2, P2 = signal.periodogram(co2_detr, fs=12)
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].semilogy(f2[1:], P2[1:], lw=0.8)
axes[0].axvline(1, color="C3", ls="--", lw=0.8, label="1 cycle/yr")
axes[0].axvline(2, color="C3", ls=":", lw=0.8, label="2 cycles/yr (harmonic)")
axes[0].set_xlabel("frequency (cycles / year)")
axes[0].set_title("CO2 periodogram (quadratic trend removed)")
axes[0].legend()
axes[1].loglog(f[1:], P[1:], lw=1)
for per in [3.5, 5]:
    axes[1].axvline(1 / per, color="C3", ls="--", lw=0.8)
axes[1].text(1 / 4.5, P[1:].max() * 0.5, "ENSO band\n(3-7 yr)", ha="center")
axes[1].set_xlabel("frequency (cycles / year)")
axes[1].set_title("GISTEMP Welch PSD (detrended): red-noise slope + ENSO bump")
savefig(fig, "spectra", SEC)
peak = f2[1:][np.argmax(P2[1:])]
sec = f2[(f2 > 1.5)][np.argmax(P2[(f2 > 1.5)])]
print(
    f"CO2 dominant frequency: {peak:.3f} cycles/yr (period {1 / peak:.2f} yr); "
    f"strongest above 1.5 cycles/yr: {sec:.3f}"
)
print(
    f"Variance explained by the annual line (f in 0.9-1.1): "
    f"{P2[(f2 > 0.9) & (f2 < 1.1)].sum() / P2[1:].sum():.2f}"
)
note("""CO2's spectrum is a line at 1 cycle/yr plus its harmonic at 2 (the seasonal cycle is
not a pure sinusoid). Global temperature, once detrended, has a 'red' spectrum: power
rising toward low frequencies, as any AR(1)/integrating system does. A bump in the 3-7
year band is ENSO. Red noise is why eyeballing 'cycles' in climate series is dangerous:
a significance test against an AR(1) background spectrum is required.""")

# %% ----------------------------------------------------------- cross-correlation
banner("Cross-correlation: CO2 growth rate vs global temperature anomaly")
growth = co2.diff(12).dropna()
gt_a = gt.reindex(growth.index)
both = pd.DataFrame({"growth": growth, "T": gt_a}).dropna()
# detrend both so the shared trend does not swamp the lagged relationship
g_d = signal.detrend(both.growth.values)
t_d = signal.detrend(both["T"].values)
lags = np.arange(-24, 25)


def lagged_corr(x, y, L):
    """corr(x_t, y_{t+L}); positive L means x leads y."""
    if L >= 0:
        return np.corrcoef(x[: len(x) - L], y[L:])[0, 1]
    return np.corrcoef(x[-L:], y[: len(y) + L])[0, 1]


cc = [lagged_corr(t_d, g_d, L) for L in lags]  # positive lag: temperature leads growth
raw_cc = both.growth.corr(both["T"])
best_lag = lags[int(np.argmax(cc))]
print(f"raw correlation (both trending): r={raw_cc:.2f}")
print(
    f"detrended max cross-correlation r={max(cc):.2f} at lag {best_lag} months "
    f"(positive = temperature leads CO2 growth)"
)
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(lags, cc, width=0.8)
ax.axhline(2 / np.sqrt(len(g_d)), color="C3", ls="--", lw=0.8)
ax.axhline(-2 / np.sqrt(len(g_d)), color="C3", ls="--", lw=0.8)
ax.set_xlabel("lag (months); positive = temperature leads")
ax.set_ylabel("cross-correlation (detrended)")
ax.set_title("Warm years (El Nino) are followed by faster CO2 growth")
savefig(fig, "cross_correlation", SEC)
note(f"""Two trending series always correlate (r={raw_cc:.2f} raw) and it means nothing. Detrended,
the cross-correlation peaks at lag {best_lag:+d} months with temperature leading CO2 growth:
El Nino warmth and drought suppress tropical carbon uptake and raise fire emissions.
Growth is a trailing 12-month difference centred 6 months back, so relative to the window
centre global temperature and CO2 growth are nearly coincident ({best_lag - 6:+d} months);
tropical Pacific SST (Nino 3.4) would lead both by a few more months. The 2/sqrt(n) lines
are the 95% band for white noise; both series are autocorrelated, so the true band is wider.""")
