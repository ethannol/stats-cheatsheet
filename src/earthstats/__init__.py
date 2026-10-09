"""earthstats: statistical-methods cheat sheet built on real earth-science data.

Modules
-------
data      download + cache the datasets (NOAA GML CO2, NASA GISTEMP, NOAA GHCN-Daily,
          USGS NWIS streamflow, NOAA GSOY station climatologies)
stats     hand-rolled implementations of the methods that matter in environmental
          statistics (Mann-Kendall, Sen's slope, Pettitt, Moran's I, variogram,
          ordinary kriging, IDW, bootstrap variants, AR(1) effective sample size)
plotting  consistent matplotlib styling and a save helper
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
