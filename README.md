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

## Headline results (as talking points)

From the script output; rerun to refresh.

- Walker Lake's 470 samples average 435 against an exhaustive truth of 278: preferential sampling.
  Cell declustering at the 20-unit cell recovers 287 (94% of the bias removed); the GSLIB weights
  cross-check correlates at r = 0.96 with ours.
- Scored against the true 5x5 block grades, anisotropic ordinary kriging is unbiased with RMSE 115
  and a conditional-bias slope of 0.92; IDW carries a +57 bias straight from the sampling; nearest
  neighbour is unbiased but 30% less accurate. Every estimator is smoother than the truth.
- Variance of block grades falls from 62,400 at point support to 52,300 at 5x5 and 20,100 at 50x50;
  at a cut-off of 500 the point-support truth reports 13% more tonnes than the block truth and the
  kriged model 16% fewer. Simulation reproduces the declustered histogram faithfully, which is still
  richer than the deposit, so its uncertainty band sits above the truth: simulation cannot repair a
  biased histogram.
- Three independent quantities closed to 100% acquire pairwise correlations of about -0.5; in the
  Jura data Co/Ni is the most proportional pair (log-ratio variance 0.08).
- Synthetic QA/QC: Thompson-Howarth recovers the generating precision (sd = 0.01 + 0.10 c vs
  0.02 floor and 12% slope); CUSUM flags the planted +8% lab bias within ~8 batches when 3SD limits
  do not; all 7 carry-over blanks follow samples above 3 g/t; only half the planted sample swaps are
  detectable from duplicate statistics.
- Clark Fork dissolved lead has 36% non-detects with detection limits from 0.04 to 100 ug/L:
  substituting 0, DL/2 or DL gives means of 0.23, 1.01 and 1.79 ug/L against a Kaplan-Meier 0.28.
- Gy's theory: a 1.5 g/t gold sample needs ~57 kg at 1 cm top size for 10% fundamental error (11 t
  under Gy's classic liberation factor) but 37 g at 75 um; splitting coarse gives an 89% protocol
  error, crushing to 2 mm before the first split gives 24%.
- Weibull MLE with right censoring recovers beta 2.27 / 0.91 / 3.08 for generating values 2.4 /
  0.85 / 3.0; the naive MTBF under-states the Weibull MTBF because long lives get censored. A
  1.15-sigma production drop is caught 8 days (EWMA) and 15 days (CUSUM) after the change vs 46
  days for a Shewhart chart.
- US mine injury rates fell from 8.2 (coal) and 5.1 (metal/nonmetal) per 200,000 h in 2000 to 3.4
  and 2.0 in 2024. Across mines the negative-binomial standard error is 1.9x the Poisson one. The
  raw worst-five list is all mines with one injury in ~2,000 hours; Gamma-Poisson shrinkage replaces
  it with large mines at 11-16 injuries per 200,000 h.
- Copper prices: 21% annualised volatility, fat tails (excess kurtosis 4.8), no rejection of a unit
  root over 35 years (estimated half-life 17 years); a 10-year GBM fan is 9x wide vs 6x under the
  fitted mean reversion. The synthetic copper project has a deterministic NPV of $0.9 bn and a 27%
  Monte Carlo probability of negative NPV.
- Clark Fork copper rises with flow (C-Q slope +0.90: tailings mobilisation, not dilution); the
  flow-adjusted decline is 2.3%/yr (Kendall tau -0.51), with a CUSUM regime change in 1996 and
  exceedance of a 13 ug/L benchmark falling from 97% to 85% of samples.

## Left out and why

- Real drillhole intervals: open state-survey drillhole databases (SARIG, GSQ, GSWA) are multi-GB
  bulk exports; compositing is demonstrated on a 1-D hole constructed from the Walker Lake grid.
- Indicator kriging and change-of-support corrections (affine / discrete Gaussian) are described in
  the cheat sheet but not implemented; SGS plus re-blocking covers the same ground.
- Groundwater / piezometer monitoring: no suitable public mine-site dataset was found; the
  environmental section uses surface-water chemistry instead.
