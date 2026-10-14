"""statsheet: statistical-methods cheat sheet built on mining-industry data problems.

Modules
-------
data      download + cache the datasets (Walker Lake, Jura, GSLIB cluster, MSHA, Water Quality
          Portal, World Bank commodity prices) and generate the clearly labelled synthetic ones
          (assay QA/QC, equipment failures)
stats     hand-rolled implementations: cell declustering, variograms (directional), ordinary
          kriging, IDW, sequential Gaussian simulation, normal-score transform, log-ratio
          transforms, Thompson-Howarth, Kaplan-Meier (left and right censoring), Weibull MLE,
          control charts, Gy's fundamental sampling error
plotting  consistent matplotlib styling and a save helper
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figures"
