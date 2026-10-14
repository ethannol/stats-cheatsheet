# Statistics cheat sheet

A runnable, annotated tour of the statistical toolkit a data scientist uses on environmental
and mining data. Sections 01-11 cover general methods (EDA, hypothesis tests, regression, time
series, trends and change points, extremes, resampling, multivariate, spatial, Bayesian, machine
learning) worked on climate, river and station data; sections 12-21 cover the methods specific
to mining practice (assay EDA and capping, compositional data, QA/QC, geostatistics, sampling
theory, equipment reliability, safety rates, project economics, environmental monitoring,
censored inference and compliance limits). Every method is applied to real, freely available
observations or to clearly labelled synthetic data, and every script prints an interpretation of
what it computed, so the output reads like a worked example and the narrative follows the numbers.

Two audiences, both me: interview preparation (what is the method, when do you use it, what
are its assumptions, how do you read the output, what goes wrong) and a demonstration of
earth-science and mining data skills (hydrology, climatology, spatial analysis, resource
estimation, QA/QC, compliance).

Start with [CHEATSHEET.md](CHEATSHEET.md) for the condensed reference, or with
[FLOWCHART.md](FLOWCHART.md) to go from a question to the script that answers it, then run the
scripts.

This repository merges two earlier ones, `earth-stats-cheatsheet` (now sections 01-11) and
`mining-stats-cheatsheet` (now sections 12-21). They had been split by repository rather than
by subject, and a reader uses both; their git histories are preserved here.

## Data

| Dataset | Source | Real / synthetic | Used for |
|---|---|---|---|
| Mauna Loa CO2 | [NOAA GML](https://gml.noaa.gov/ccgg/trends/data.html) | real (monthly mean CO2, 1958-present) | time series, decomposition, SARIMA, spectra |
| GISTEMP v4 | [NASA GISS](https://data.giss.nasa.gov/gistemp/) | real (global land-ocean temperature anomaly, monthly 1880-present) | trends, change points, regression, Bayesian |
| GHCN-Daily, Central Park | [NOAA NCEI](https://www.ncei.noaa.gov/access/services/data/v1) | real (daily TMAX/TMIN/PRCP/SNOW, 1950-present) | EDA, hypothesis tests, extremes, bootstrap |
| GHCN-Daily, Washington National | NOAA NCEI | real (daily weather for the Potomac basin) | rainfall-runoff ML features |
| USGS 01646500 | [USGS NWIS](https://waterservices.usgs.gov/) | real (Potomac River daily discharge, 1930-present) | floods, return periods, hydrology ML |
| GSOY station climatologies | NOAA NCEI | real (1991-2020 annual summaries for ~370 US first-order stations) | multiple regression, PCA, clustering, spatial |
| Walker Lake (470 samples + 78,000-node exhaustive grid) | Isaaks & Srivastava 1989, via R `gstat` | real (derived from a Nevada DEM; the standard mining-geostatistics teaching deposit) | assay EDA, declustering, capping, variography, kriging vs truth, change of support, simulation |
| Jura soil metals (359 samples, 7 metals, rock type) | Goovaerts 1997, via R `gstat` | real | compositional data analysis |
| GSLIB `cluster.dat` (140 samples with declustering weights) | Deutsch & Journel, via R `gstat` | real teaching data | declustering cross-check |
| MSHA Mines, MinesProdYearly, Accidents | [MSHA open data](https://arlweb.msha.gov/OpenGovernmentData/OGIMSHA.asp) | real (every US mine, hours and injuries 2000-2025) | safety rates, overdispersion, shrinkage, control charts |
| Clark Fork River at Deer Lodge MT (USGS 12324200) | [Water Quality Portal](https://www.waterqualitydata.us), USGS NWIS | real (1969-2023 metals and pH with non-detects, daily flow) | censored data, C-Q relationships, trends, compliance, seasonal Kendall, censored inference, UCLs and prediction limits |
| World Bank Pink Sheet monthly commodity prices | [World Bank CMO](https://www.worldbank.org/en/research/commodity-markets) | real (1960-2024) | returns, volatility, random walk vs mean reversion |
| Gold assay QA/QC (duplicates, CRMs, blanks) | `statsheet.data.assay_qaqc` | **synthetic**, seeded, with planted problems | precision, bias, carry-over detection |
| Haul-truck component failures and daily production | `statsheet.data.equipment_failures` | **synthetic**, seeded, generating parameters stated | Weibull, Kaplan-Meier, availability, SPC |
| Copper project cash-flow model | `analyses/19_economics.py` | **synthetic** | Monte Carlo NPV, tornado, cut-off grade, option to defer |

No public assay QA/QC database or mine maintenance log exists, so those two are simulated from
stated parameters and every recovered estimate is compared with the value used to generate it.
Everything else is downloaded once by `make data` and cached under `data/` (about 100 MB, most
of it the MSHA accident file). The data are not committed: `reports/data-manifest.json` records
what was fetched, from where, and its hashes. Please credit the providers if you reuse the data;
NOAA GML in particular asks for it.

## Setup and run

```bash
uv sync --extra dev      # Python >= 3.12; installs numpy/pandas/scipy/statsmodels/sklearn/matplotlib/rdata/openpyxl
make data                # download and cache all datasets (MSHA is ~60 MB zipped; be patient)
make all                 # run every analysis, writing figures/ and printing interpretations
make 06                  # run one topic by its number prefix (01-21)
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
  12_assay_eda.py                skewed grades, probability plots, cell declustering vs truth, top-cuts, compositing
  13_compositional.py            closure and spurious correlation, alr/clr/ilr, variation matrix, clr biplot
  14_qaqc.py                     duplicates (HARD, Thompson-Howarth), CRM control charts + CUSUM bias, blanks,
                                 non-detects (substitution vs Kaplan-Meier vs ROS) on real river data
  15_geostatistics.py            omni/directional variograms, anisotropy, OK vs IDW vs NN scored against the truth,
                                 cross-validation, volume-variance, grade-tonnage, sequential Gaussian simulation
  16_sampling_theory.py          Gy's mineralogical constant, FSE nomogram, protocol design, nugget decomposition
  17_reliability.py              Kaplan-Meier, Weibull MLE with censoring, MTBF/MTTR/availability, age replacement, SPC
  18_safety.py                   injuries per 200,000 h, Poisson vs negative binomial, empirical-Bayes shrinkage, u-chart
  19_economics.py                log returns and volatility, GBM vs mean reversion, Monte Carlo NPV, tornado, cut-off
  20_environmental.py            concentration-discharge, flow-adjusted trends, exceedance fractions, CUSUM change point
  21_censored_inference.py       seasonal Kendall (flow-adjusted), Gehan test, censored MLE and Akritas-Theil-Sen
                                 regression, ProUCL-style UCL95 selection, prediction and tolerance limits, Dunn
src/statsheet/
  data.py                        loaders with caching for every dataset; synthetic generators with stated parameters
  stats.py                       hand-rolled: Mann-Kendall (+Hamed-Rao), Sen's slope, pre-whitening, Pettitt,
                                 effective sample size, bootstraps, Moran's I, variograms (from a distance matrix
                                 or directional on planar coordinates), ordinary kriging (global or moving
                                 neighbourhood with anisotropy), IDW, GEV helpers, cell declustering, compositing,
                                 normal scores, sequential Gaussian simulation, grade-tonnage, alr/clr/ilr,
                                 HARD, Thompson-Howarth, Kaplan-Meier (left and right censoring) with SE, ROS,
                                 Gehan, Akritas-Theil-Sen, censored MLE regression, seasonal Kendall, Dunn, UCLs,
                                 prediction and tolerance limits, Weibull MLE, Shewhart/CUSUM/EWMA/u-chart, Gy's FSE,
                                 GBM, OU fit
  spatial.py                     fold-local trend and variogram fitting for the leave-one-out example
  plotting.py                    shared style, savefig, banner/note helpers
tests/                           the hand-rolled algorithms vs scipy / statsmodels / closed forms / known answers,
                                 plus offline fixtures for cache provenance and provider cutoffs
figures/<topic>/*.png            written by the scripts
CHEATSHEET.md                    the condensed reference
FLOWCHART.md                     question-first decision charts (Mermaid), with sources
```

## Worked results and reproducibility

[Generated worked results](reports/worked-results.md) contain actual script output rather
than manually refreshed headline numbers. The [run manifest](reports/results.json) records
the source commit, source patch hash, Python/dependency versions, date cutoff, script and
output hashes, and figure hashes. [Input provenance](reports/data-manifest.json) records
request URLs/parameters, SHA-256, sizes and retrieval times. Existing unmanifested caches
are explicitly labelled `legacy-cache`, with unknown retrieval time; metadata is never
backdated. The small teaching inputs fetched from GitHub use a pinned upstream commit.

```bash
make results             # full analysis; run on a workstation/CI, not the laptop
make reproduce           # requires the exact raw-byte hashes in the committed manifest
```

`make results` writes per-topic transcripts and selected portfolio preview figures under
`reports/`. Dates are clipped to `END_DATE` even when an upstream file has newer rows.
A fixed cutoff does not prevent providers revising historical values: retain the raw
snapshot alongside its manifest for exact reproduction. Changed bytes are rejected until
an explicit refresh. Derived caches are rebuilt after loader/cutoff changes.

## Extending, and what was left out

- Swap stations or gauges by ID in `statsheet.data` (`ghcn_daily("USW00023174")` for LAX,
  `usgs_daily("09380000")` for the Colorado at Lees Ferry). Caches are keyed by ID.
- `make refresh` re-downloads after bumping `END_DATE` in `data.py`.
- Candidates not covered: EOF/PCA on gridded fields (xarray), hierarchical Bayesian models
  (PyMC), spatio-temporal kriging, quantile regression, GAMs, causal inference designs.
- Real drillhole intervals: open state-survey drillhole databases (SARIG, GSQ, GSWA) are multi-GB
  bulk exports; compositing is demonstrated on a 1-D hole constructed from the Walker Lake grid.
- Indicator kriging and change-of-support corrections (affine / discrete Gaussian) are described in
  the cheat sheet but not implemented; SGS plus re-blocking covers the same ground.
- Groundwater / piezometer monitoring: no suitable public mine-site dataset was found; the
  environmental sections use surface-water chemistry instead.

## How AI was used

Both earlier cheat sheets, and this merge of them, were written in Claude Code sessions with
Claude Fable 5.1 from the author's briefs: the agent wrote the loaders, the hand-rolled
implementations, the analysis scripts and their interpretation notes, the tests and these
documents, and the author reviewed the choice of methods, the interpretations and the tests.
No language model is called by any script; every number in the worked results is emitted by a
script and recorded with its source, input and output hashes.
