"""Fold-local interpolation fitting, shared with the worked spatial example."""
import numpy as np
import statsmodels.formula.api as smf

from .stats import empirical_variogram, fit_variogram, spherical


def fit_interpolation_fold(stations, distances, keep, bins):
    """Fit trend and raw/residual variograms using only training stations.

    Held-out temperatures cannot influence either fitted variogram. Covariates at query
    stations are known inputs; spatial block holdouts remain a separate transfer test.
    """
    training = stations.loc[keep]
    known_distances = distances[np.ix_(keep, keep)]
    trend = smf.ols("tavg ~ lat + elev_km + lon", data=training).fit()
    raw_params = fit_variogram(
        *empirical_variogram(known_distances, training.tavg.values, bins)[:2], spherical
    )
    residual_params = fit_variogram(
        *empirical_variogram(known_distances, trend.resid.values, bins)[:2], spherical
    )
    return trend, raw_params, residual_params
