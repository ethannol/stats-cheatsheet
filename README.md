# Earth-science statistics cheat sheet

A runnable, annotated tour of the statistical toolkit a data scientist uses on environmental
data. Every method is applied to real, freely available observations, and every script prints
an interpretation of what it computed so you can read the output like a worked example.

Two audiences, both me: interview preparation (what is the method, when do you use it, what
are its assumptions, how do you read the output, what goes wrong) and a demonstration of
earth-science data skills (hydrology, climatology, spatial analysis).

Start with [CHEATSHEET.md](CHEATSHEET.md) for the condensed reference, then run the scripts.

## Data

| Dataset | Source | What it is | Used for |
|---|---|---|---|
| Mauna Loa CO2 | [NOAA GML](https://gml.noaa.gov/ccgg/trends/data.html) | Monthly mean CO2, 1958-present | time series, decomposition, SARIMA, spectra |
| GISTEMP v4 | [NASA GISS](https://data.giss.nasa.gov/gistemp/) | Global land-ocean temperature anomaly, monthly 1880-present | trends, change points, regression, Bayesian |
| GHCN-Daily, Central Park | [NOAA NCEI](https://www.ncei.noaa.gov/access/services/data/v1) | Daily TMAX/TMIN/PRCP/SNOW, 1950-present | EDA, hypothesis tests, extremes, bootstrap |
| GHCN-Daily, Washington National | NOAA NCEI | Daily weather for the Potomac basin | rainfall-runoff ML features |
| USGS 01646500 | [USGS NWIS](https://waterservices.usgs.gov/) | Potomac River daily discharge, 1930-present | floods, return periods, hydrology ML |
| GSOY station climatologies | NOAA NCEI | 1991-2020 annual summaries for ~370 US first-order stations | multiple regression, PCA, clustering, spatial |

Everything is downloaded once by `make data` and cached under `data/` (about 20 MB).
Please credit the providers if you reuse the data; NOAA GML in particular asks for it.

## Setup and run

```bash
uv sync --extra dev      # Python >= 3.12; installs numpy/pandas/scipy/statsmodels/sklearn/matplotlib
make data                # download and cache all datasets
make all                 # run every analysis, writing figures/ and printing interpretations
make 06                  # run one topic by its number prefix
make test                # cross-check the hand-rolled algorithms against scipy/statsmodels
make notebooks           # export the scripts to Jupyter notebooks (jupytext)
```

Each script in `analyses/` is a plain Python file with `# %%` cell markers, so it runs
top to bottom, opens as an interactive notebook in VS Code / PyCharm, or converts to `.ipynb`.

## Layout

```
analyses/                        one script per topic; the docstring at the top is the cheat sheet
  01_eda.py                      distributions, robust statistics, climatology/anomalies, correlation
  02_hypothesis_tests.py         t/Welch/Mann-Whitney/permutation, effect size, proportions, FDR, power
  03_regression.py               OLS + diagnostics, HAC/GLS, multicollinearity, multiple regression, GLM
  04_time_series.py              STL, ADF/KPSS, ACF/PACF, SARIMA with holdout + baselines, spectra, CCF
  05_trends_changepoints.py      OLS/MK/Sen trends, autocorrelation corrections, Pettitt, hinge, CUSUM
  06_extremes.py                 GEV/Gumbel/LP3 flood frequency, POT/GPD, nonstationary GEV
  07_resampling.py               bootstrap variants + coverage, block bootstrap, jackknife, Monte Carlo
  08_multivariate.py             correlation structure, PCA, k-means, Ward clustering
  09_spatial.py                  Moran's I, variograms, IDW vs regression kriging, LOO validation
  10_bayesian.py                 Beta-Binomial, Gamma-Poisson, hand-rolled Metropolis-Hastings
  11_ml.py                       rainfall-runoff regression, leakage demo, importance, classification
src/earthstats/
  data.py                        loaders with caching for the six datasets
  stats.py                       hand-rolled: Mann-Kendall (+Hamed-Rao), Sen's slope, pre-whitening,
                                 Pettitt, Moran's I, variogram, ordinary kriging, IDW, bootstraps,
                                 effective sample size, GEV helpers
  plotting.py                    shared style, savefig, banner/note helpers
tests/test_stats.py              the hand-rolled algorithms vs scipy / known answers
figures/<topic>/*.png            written by the scripts
CHEATSHEET.md                    the condensed reference
```

## Headline results (as talking points)

These come straight from the script output; rerun to refresh.

- Central Park summer **nights** warmed 0.7 C between 1951-80 and 1991-2020 (Welch p = 0.0005,
  Cohen's d = 0.95) while summer **daytime highs** did not change (p = 0.74, CI -0.6 to +0.4 C).
  Testing on daily values instead of seasonal means gives p = 1e-17 for the same effect:
  the unit-of-analysis lesson.
- Regressing global temperature on ln(CO2) implies 2.7 C per doubling: a clean fit (R2 = 0.93)
  that overshoots the IPCC transient response because every other forcing is collinear with CO2.
  Residual autocorrelation (DW = 1.4) inflates the naive t-statistic by ~40%.
- The airline SARIMA(0,1,1)(0,1,1)12 fits Mauna Loa CO2; a 6-year holdout forecast has RMSE
  1.0 ppm vs 1.7 ppm for seasonal naive with drift. Global temperature's detrended spectrum is
  red with an ENSO bump; CO2 growth lags warm years.
- Raw Mann-Kendall on AR(1) noise with rho = 0.6 has a 33% false-positive rate; trend-free
  pre-whitening makes it *worse* (42%); the Hamed-Rao variance correction holds it at 9%.
- Potomac 100-year daily flow is ~390 kcfs by GEV with a bootstrap interval of roughly +-37%;
  Gumbel vs GEV cannot be settled by 95 years of data. NYC's warmest night of the year has
  risen 1.8 C since 1950, turning a 1950s 20-year event into a roughly 1-in-2-year event.
- A block bootstrap gives a trend SE 2.6x the iid bootstrap on monthly GISTEMP; the iid
  bootstrap simply reproduces the (wrong) OLS formula.
- Regression kriging (lat + elevation + longitude trend, kriged residuals) interpolates station
  temperature with LOO RMSE 0.8 C vs 2.0 C for IDW. Residuals of the regression have Moran's
  I = 0.42: spatially clustered errors, so the regression SEs are too small.
- Rainfall-runoff: weather-only models all plateau at R2 ~0.66 (information limit of one rain
  gauge), persistence dominates one day ahead, and shuffled K-fold overstates R2 by 0.18 relative
  to a time-aware split.

## Extending

- Swap stations or gauges by ID in `earthstats.data` (`ghcn_daily("USW00023174")` for LAX,
  `usgs_daily("09380000")` for the Colorado at Lees Ferry). Caches are keyed by ID.
- `make refresh` re-downloads after bumping `END_DATE` in `data.py`.
- Candidates not covered here: EOF/PCA on gridded fields (xarray), hierarchical Bayesian models
  (PyMC), spatio-temporal kriging, quantile regression, GAMs, causal inference designs.
