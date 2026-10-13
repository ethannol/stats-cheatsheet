# Earth-science statistics cheat sheet

A runnable, annotated tour of the statistical toolkit a data scientist uses on environmental
data. Every method is applied to real, freely available observations, and every script prints
an interpretation of what it computed so you can read the output like a worked example.

Two audiences, both me: interview preparation (what is the method, when do you use it, what
are its assumptions, how do you read the output, what goes wrong) and a demonstration of
earth-science data skills (hydrology, climatology, spatial analysis).

Start with [CHEATSHEET.md](CHEATSHEET.md) for the condensed reference, then run the scripts.
A question-first decision chart that spans this sheet and its mining companion lives in
[`../mining-stats-cheatsheet/FLOWCHART.md`](../mining-stats-cheatsheet/FLOWCHART.md).

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

## Worked results and reproducibility

[Generated worked results](reports/worked-results.md) contain actual script output rather
than manually refreshed headline numbers. The [run manifest](reports/results.json) records
the source commit, source patch hash, Python/dependency versions, date cutoff, script and
output hashes, and figure hashes. [Input provenance](reports/data-manifest.json) records
request URLs/parameters, SHA-256, sizes and retrieval times. Existing unmanifested caches
are explicitly labelled `legacy-cache`, with unknown retrieval time; metadata is never
backdated.

```bash
make results             # full analysis; run on a workstation/CI, not the laptop
make reproduce           # requires the exact raw-byte hashes in the committed manifest
```

`make results` writes per-topic transcripts and selected portfolio preview figures under
`reports/`. Dates are clipped to `END_DATE` even when an upstream file has newer rows.
A fixed cutoff does not prevent providers revising historical values: retain the raw
snapshot alongside its manifest for exact reproduction. Changed bytes are rejected until
an explicit refresh. Derived caches are rebuilt after loader/cutoff changes.

## Extending

- Swap stations or gauges by ID in `earthstats.data` (`ghcn_daily("USW00023174")` for LAX,
  `usgs_daily("09380000")` for the Colorado at Lees Ferry). Caches are keyed by ID.
- `make refresh` re-downloads after bumping `END_DATE` in `data.py`.
- Candidates not covered here: EOF/PCA on gridded fields (xarray), hierarchical Bayesian models
  (PyMC), spatio-temporal kriging, quantile regression, GAMs, causal inference designs.
