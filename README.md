# Mining data statistics cheat sheet

A runnable, annotated tour of the statistical methods mining companies actually use: resource
estimation, assay QA/QC, geochemistry, sampling theory, equipment reliability, safety analytics,
project economics and environmental compliance. Every script prints an interpretation of what it
computed, so the output reads like a worked example, and the narrative follows the numbers.

Companion to `../earth-stats-cheatsheet` (same layout and conventions). Start with
[CHEATSHEET.md](CHEATSHEET.md) for the condensed reference, then run the scripts.

## Data

| Dataset | Source | Real / synthetic | Used for |
|---|---|---|---|
| Walker Lake (470 samples + 78,000-node exhaustive grid) | Isaaks & Srivastava 1989, via R `gstat` | real (derived from a Nevada DEM; the standard mining-geostatistics teaching deposit) | EDA, declustering, capping, variography, kriging vs truth, change of support, simulation |
| Jura soil metals (359 samples, 7 metals, rock type) | Goovaerts 1997, via R `gstat` | real | compositional data analysis |
| GSLIB `cluster.dat` (140 samples with declustering weights) | Deutsch & Journel, via R `gstat` | real teaching data | declustering cross-check |
| MSHA Mines, MinesProdYearly, Accidents | [MSHA open data](https://arlweb.msha.gov/OpenGovernmentData/OGIMSHA.asp) | real (every US mine, hours and injuries 2000-2025) | safety rates, overdispersion, shrinkage, control charts |
| Clark Fork River at Deer Lodge MT (USGS 12324200) | [Water Quality Portal](https://www.waterqualitydata.us), [USGS NWIS](https://waterservices.usgs.gov) | real (1969-2023 metals and pH with non-detects, daily flow) | censored data, C-Q relationships, trends, compliance |
| World Bank Pink Sheet monthly commodity prices | [World Bank CMO](https://www.worldbank.org/en/research/commodity-markets) | real (1960-2024) | returns, volatility, random walk vs mean reversion |
| Gold assay QA/QC (duplicates, CRMs, blanks) | `miningstats.data.assay_qaqc` | **synthetic**, seeded, with planted problems | precision, bias, carry-over detection |
| Haul-truck component failures and daily production | `miningstats.data.equipment_failures` | **synthetic**, seeded, generating parameters stated | Weibull, Kaplan-Meier, availability, SPC |
| Copper project cash-flow model | `analyses/08_economics.py` | **synthetic** | Monte Carlo NPV, tornado, cut-off grade, option to defer |

No public assay QA/QC database or mine maintenance log exists, so those two are simulated from
stated parameters and every recovered estimate is compared with the value used to generate it.
Everything else is downloaded once by `make data` and cached under `data/` (about 80 MB, most of
it the MSHA accident file). Please credit the providers if you reuse the data.

## Setup and run

```bash
uv sync --extra dev      # Python >= 3.12
make data                # download and cache all datasets (MSHA is ~60 MB zipped; be patient)
make all                 # run every analysis, writing figures/ and printing interpretations
make 04                  # run one topic by number prefix
make test                # cross-check the hand-rolled algorithms
make notebooks           # export the scripts to Jupyter notebooks (jupytext)
```

Each script in `analyses/` is a plain Python file with `# %%` cell markers: run it top to bottom,
open it as an interactive notebook in VS Code / PyCharm, or convert to `.ipynb`.

## Layout

```
analyses/
  01_assay_eda.py            skewed grades, probability plots, cell declustering vs truth, top-cuts, compositing
  02_compositional.py        closure and spurious correlation, alr/clr/ilr, variation matrix, clr biplot, Aitchison distance
  03_qaqc.py                 duplicates (HARD, Thompson-Howarth), CRM control charts + CUSUM bias, blanks/carry-over,
                             non-detects (substitution vs Kaplan-Meier vs ROS) on real river data
  04_geostatistics.py        omni/directional variograms, anisotropy, OK vs IDW vs NN scored against the exhaustive truth,
                             cross-validation, volume-variance, grade-tonnage, sequential Gaussian simulation
  05_sampling_theory.py      Gy's mineralogical constant, FSE nomogram, protocol design, nugget decomposition, sample counts
  06_reliability.py          Kaplan-Meier, Weibull MLE with censoring, MTBF/MTTR/availability, age replacement, SPC detection
  07_safety.py               injuries per 200,000 h, Poisson vs negative binomial, empirical-Bayes shrinkage, u-chart, severity
  08_economics.py            log returns and volatility, GBM vs mean reversion fans, Monte Carlo NPV, tornado, cut-off, deferral option
  09_environmental.py        concentration-discharge, flow-adjusted trends, exceedance fractions with CIs, CUSUM change point
src/miningstats/
  data.py                    loaders with caching; synthetic generators with stated parameters
  stats.py                   hand-rolled: cell declustering, compositing, directional variograms + models, ordinary kriging
                             (moving neighbourhood, anisotropy), IDW, normal scores, sequential Gaussian simulation,
                             grade-tonnage, alr/clr/ilr, variation matrix, HARD, Thompson-Howarth, Kaplan-Meier (left and
                             right censoring), ROS, Weibull MLE, Shewhart/CUSUM/EWMA/u-chart, Gy's FSE, GBM, OU fit
  plotting.py                shared style and helpers
tests/test_stats.py          20 cross-checks against scipy / closed forms / known answers
figures/<topic>/*.png        written by the scripts
CHEATSHEET.md                the condensed reference
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

## Left out and why

- Real drillhole intervals: open state-survey drillhole databases (SARIG, GSQ, GSWA) are multi-GB
  bulk exports; compositing is demonstrated on a 1-D hole constructed from the Walker Lake grid.
- Indicator kriging and change-of-support corrections (affine / discrete Gaussian) are described in
  the cheat sheet but not implemented; SGS plus re-blocking covers the same ground.
- Groundwater / piezometer monitoring: no suitable public mine-site dataset was found; the
  environmental section uses surface-water chemistry instead.
