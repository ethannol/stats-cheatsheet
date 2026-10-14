"""Small non-rendering leakage regression fixture."""

import numpy as np
import pandas as pd

from statsheet.spatial import fit_interpolation_fold


def test_held_out_temperature_cannot_change_fitted_models():
    rng = np.random.default_rng(29)
    coords = rng.uniform(0, 500, (40, 2))
    distance = np.sqrt(((coords[:, None] - coords[None, :]) ** 2).sum(axis=2))
    st = pd.DataFrame(
        {"lat": coords[:, 0] / 20, "lon": coords[:, 1] / 20, "elev_km": rng.uniform(0, 2, 40)}
    )
    st["tavg"] = 12 - 0.2 * st.lat - st.elev_km + rng.normal(size=40)
    keep = np.arange(40) != 0
    before = fit_interpolation_fold(st, distance, keep, np.arange(0, 801, 100))
    st.loc[0, "tavg"] = 1e9
    after = fit_interpolation_fold(st, distance, keep, np.arange(0, 801, 100))
    np.testing.assert_allclose(before[0].params, after[0].params)
    np.testing.assert_allclose(before[1], after[1])
    np.testing.assert_allclose(before[2], after[2])
