"""statsheet: statistical-methods cheat sheet built on real environmental and mining data.

Sections 01-11 work general methods on climate, river and station data; sections 12-21 work
the methods specific to mining practice on deposit, assay, safety, price and river-chemistry data.

Modules
-------
data        download + cache the datasets (NOAA GML CO2, NASA GISTEMP, NOAA GHCN-Daily, USGS NWIS
            streamflow, NOAA GSOY station climatologies, Walker Lake, Jura, GSLIB cluster, MSHA,
            Water Quality Portal, World Bank commodity prices) and generate the clearly labelled
            synthetic ones (assay QA/QC, equipment failures, a drillhole from the Walker Lake grid)
stats       hand-rolled implementations: Mann-Kendall (+Hamed-Rao), Sen's slope, pre-whitening,
            Pettitt, AR(1) effective sample size, bootstrap variants, permutation test, Moran's I,
            variograms (from a distance matrix, or directional on planar coordinates) and models,
            ordinary kriging (global from distances, or moving neighbourhood with anisotropy), IDW,
            GEV helpers, cell declustering, compositing, normal scores, sequential Gaussian
            simulation, grade-tonnage, log-ratio transforms, Thompson-Howarth, Kaplan-Meier (left
            and right censoring) with SE, ROS, Gehan, censored Kendall / Akritas-Theil-Sen, censored
            MLE regression, seasonal Kendall, Dunn, UCLs, prediction and tolerance limits, Weibull
            MLE, Shewhart / CUSUM / EWMA / u-chart, Gy's fundamental sampling error, GBM, OU fit
spatial     fold-local trend and variogram fitting for the leave-one-out interpolation example
plotting    consistent matplotlib styling and a save helper
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
