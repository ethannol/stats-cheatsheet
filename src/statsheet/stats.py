"""Hand-rolled implementations of statistical methods that come up constantly in
environmental data science. Each function is short enough to read in an interview and is
cross-checked against scipy/statsmodels in tests/test_stats.py.

Trend / change
    mann_kendall, mann_kendall_hamed_rao, sens_slope, prewhiten_tfpw, prewhiten_simple,
    pettitt, ar1_effective_n
Resampling
    bootstrap_ci, block_bootstrap, permutation_test
Spatial
    haversine_km, knn_weights, morans_i, idw, empirical_variogram, fit_variogram,
    ordinary_kriging
Extremes
    return_level_gev, plotting_positions
Effect sizes
    cohens_d
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats as sps
from scipy.optimize import curve_fit


# ============================================================================ trend
@dataclass
class MKResult:
    S: float
    var_S: float
    z: float
    p: float
    tau: float
    n: int

    @property
    def trend(self) -> str:
        if self.p >= 0.05:
            return "no trend (p >= 0.05)"
        return "increasing" if self.z > 0 else "decreasing"

    def __str__(self):
        return (f"Mann-Kendall: S={self.S:.0f} z={self.z:.2f} p={self.p:.4f} "
                f"tau={self.tau:.3f} n={self.n} -> {self.trend}")


def mann_kendall(x) -> MKResult:
    """Non-parametric monotonic trend test (Mann 1945, Kendall 1975).

    S = sum over all pairs i<j of sign(x_j - x_i). Under H0 (iid, no trend) S is
    approximately normal with mean 0 and variance n(n-1)(2n+5)/18, reduced for ties.
    A continuity correction of +-1 is applied to S before standardising.

    Assumes independent observations. Autocorrelation inflates the variance and gives
    false positives; see prewhiten_tfpw() or use annual values.
    """
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    n = len(x)
    # vectorised sign sum over all pairs (upper triangle)
    diff = np.sign(x[None, :] - x[:, None])
    S = np.triu(diff, k=1).sum()
    # tie correction: sum over tie groups of t(t-1)(2t+5)
    _, counts = np.unique(x, return_counts=True)
    ties = counts[counts > 1]
    var_S = (n * (n - 1) * (2 * n + 5) - np.sum(ties * (ties - 1) * (2 * ties + 5))) / 18.0
    if S > 0:
        z = (S - 1) / np.sqrt(var_S)
    elif S < 0:
        z = (S + 1) / np.sqrt(var_S)
    else:
        z = 0.0
    p = 2 * (1 - sps.norm.cdf(abs(z)))
    tau = S / (0.5 * n * (n - 1))
    return MKResult(S=float(S), var_S=float(var_S), z=float(z), p=float(p), tau=float(tau), n=n)


def sens_slope(x, t=None, alpha: float = 0.05):
    """Sen's (Theil-Sen) slope: median of all pairwise slopes. Robust to outliers and
    non-normality, breakdown point ~29%. CI via the rank method (Gilbert 1987) using the
    Mann-Kendall variance. Returns (slope, intercept, lo, hi) with
    intercept = median(x) - slope*median(t) (Conover's method, as in scipy.stats.theilslopes).
    """
    x = np.asarray(x, dtype=float)
    t = np.arange(len(x), dtype=float) if t is None else np.asarray(t, dtype=float)
    ok = ~np.isnan(x)
    x, t = x[ok], t[ok]
    i, j = np.triu_indices(len(x), k=1)
    dt = t[j] - t[i]
    good = dt != 0
    slopes = (x[j] - x[i])[good] / dt[good]
    slopes.sort()
    slope = np.median(slopes)
    intercept = np.median(x) - slope * np.median(t)
    N = len(slopes)
    var_S = mann_kendall(x).var_S
    c = sps.norm.ppf(1 - alpha / 2) * np.sqrt(var_S)
    m1 = int(np.floor((N - c) / 2))
    m2 = int(np.ceil((N + c) / 2))
    lo = slopes[max(m1, 0)]
    hi = slopes[min(m2, N - 1)]
    return float(slope), float(intercept), float(lo), float(hi)


def prewhiten_tfpw(x):
    """Trend-free pre-whitening (Yue et al. 2002) for Mann-Kendall on autocorrelated data.

    1. remove Sen's slope trend, 2. estimate lag-1 autocorrelation r1 of the residual,
    3. remove AR(1) component: y_t = r_t - r1 * r_{t-1}, 4. add the trend back.
    Returns the pre-whitened series (length n-1) and r1. Run mann_kendall() on the result.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    t = np.arange(n, dtype=float)
    slope, _, _, _ = sens_slope(x, t)
    resid = x - slope * t
    r1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
    pw = resid[1:] - r1 * resid[:-1]
    return pw + slope * t[1:], float(r1)


def prewhiten_simple(x):
    """Plain pre-whitening (von Storch 1995): y_t = x_t - r1 x_{t-1} using r1 of the raw
    series. Controls the false-positive rate but removes part of any real trend too, so it
    loses power. Returns (whitened series, r1)."""
    x = np.asarray(x, dtype=float)
    r1 = np.corrcoef(x[:-1], x[1:])[0, 1]
    return x[1:] - r1 * x[:-1], float(r1)


def mann_kendall_hamed_rao(x, alpha: float = 0.05) -> MKResult:
    """Mann-Kendall with the Hamed & Rao (1998) variance correction for autocorrelation.

    Var*(S) = Var(S) * n/n_s*, with
    n/n_s* = 1 + 2/(n(n-1)(n-2)) * sum_k (n-k)(n-k-1)(n-k-2) rho_k
    where rho_k is the lag-k autocorrelation of the *ranks* of x, summed over lags whose
    autocorrelation is significant at level alpha (|rho_k| > z_{1-alpha/2}/sqrt(n)).
    This keeps the test's nominal size under serial dependence without altering the data.
    """
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    base = mann_kendall(x)
    n = base.n
    r = sps.rankdata(x) - (n + 1) / 2
    denom = np.sum(r * r)
    crit = sps.norm.ppf(1 - alpha / 2) / np.sqrt(n)
    corr = 0.0
    for k in range(1, n - 2):
        rho_k = np.sum(r[:-k] * r[k:]) / denom
        if abs(rho_k) > crit:
            corr += (n - k) * (n - k - 1) * (n - k - 2) * rho_k
    factor = max(1.0 + 2.0 / (n * (n - 1) * (n - 2)) * corr, 1e-6)
    var_S = base.var_S * factor
    S = base.S
    z = (S - np.sign(S)) / np.sqrt(var_S) if S != 0 else 0.0
    p = 2 * (1 - sps.norm.cdf(abs(z)))
    return MKResult(S=S, var_S=float(var_S), z=float(z), p=float(p), tau=base.tau, n=n)


def ar1_effective_n(x) -> tuple[float, float]:
    """Effective sample size for an AR(1) series: n_eff = n (1 - r1) / (1 + r1).
    Use it to inflate the standard error of a regression slope: SE_adj = SE * sqrt(n/n_eff)
    (Santer et al. 2000). Returns (n_eff, r1)."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    r1 = np.corrcoef(x[:-1], x[1:])[0, 1]
    r1 = float(np.clip(r1, 0, 0.999))  # negative r1 would *increase* n_eff; cap at n
    n_eff = len(x) * (1 - r1) / (1 + r1)
    return float(n_eff), r1


def pettitt(x):
    """Pettitt (1979) non-parametric single change-point test for a shift in the median.

    U_t = sum_{i<=t} sum_{j>t} sign(x_i - x_j); K = max |U_t|; the change is at argmax.
    Approximate p = 2 exp(-6 K^2 / (n^3 + n^2)). Returns (index, K, p).
    Assumes a single abrupt shift and independent data; a gradual trend will also trigger it.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    # rank-based O(n log n) form: U_t = 2*sum_{i<=t} rank_i - t(n+1)
    r = sps.rankdata(x)
    U = 2 * np.cumsum(r)[:-1] - np.arange(1, n) * (n + 1)
    k = int(np.argmax(np.abs(U)))
    K = float(abs(U[k]))
    p = float(min(1.0, 2 * np.exp(-6 * K**2 / (n**3 + n**2))))
    return k, K, p


# ====================================================================== resampling
def bootstrap_ci(x, statistic=np.mean, n_boot: int = 5000, alpha: float = 0.05,
                 method: str = "percentile", rng=None):
    """Non-parametric bootstrap CI for an iid sample.

    percentile : quantiles of the bootstrap distribution. Simple; biased for skewed stats.
    basic      : 2*theta_hat - percentile bounds (reflected). Better coverage when biased.
    For BCa (bias-corrected, accelerated) use scipy.stats.bootstrap(method='BCa').
    Returns (estimate, lo, hi, boot_samples).
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x)
    x = x[~np.isnan(x)]
    est = statistic(x)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    boots = np.array([statistic(x[i]) for i in idx])
    lo_q, hi_q = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    if method == "percentile":
        lo, hi = lo_q, hi_q
    elif method == "basic":
        lo, hi = 2 * est - hi_q, 2 * est - lo_q
    else:
        raise ValueError(method)
    return float(est), float(lo), float(hi), boots


def block_bootstrap(x, statistic, block_len: int, n_boot: int = 2000, rng=None):
    """Moving-block bootstrap for dependent (autocorrelated) data.

    Resample contiguous blocks of length block_len with replacement and concatenate to the
    original length, preserving short-range dependence. Rule of thumb block_len ~ n^(1/3)
    or a few times the decorrelation time. Returns bootstrap statistics array.
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(x, dtype=float)
    n = len(x)
    n_blocks = int(np.ceil(n / block_len))
    starts = rng.integers(0, n - block_len + 1, size=(n_boot, n_blocks))
    out = np.empty(n_boot)
    offs = np.arange(block_len)
    for b in range(n_boot):
        idx = (starts[b][:, None] + offs).ravel()[:n]
        out[b] = statistic(x[idx])
    return out


def permutation_test(a, b, statistic=lambda a, b: np.mean(a) - np.mean(b),
                     n_perm: int = 10000, rng=None):
    """Two-sample permutation test. Under H0 labels are exchangeable, so shuffle the pooled
    sample and recompute the statistic. Exact in spirit, no distributional assumptions.
    Returns (observed, p_two_sided, null_distribution)."""
    rng = np.random.default_rng(rng)
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    obs = statistic(a, b)
    pooled = np.concatenate([a, b])
    null = np.empty(n_perm)
    for k in range(n_perm):
        rng.shuffle(pooled)
        null[k] = statistic(pooled[:len(a)], pooled[len(a):])
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (n_perm + 1)
    return float(obs), float(p), null


def cohens_d(a, b) -> float:
    """Standardised mean difference with pooled SD. ~0.2 small, 0.5 medium, 0.8 large."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp)


# ========================================================================= spatial
def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance (km) between arrays of points; broadcasts."""
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dlmb = np.radians(lon2) - np.radians(lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dlmb / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def distance_matrix_km(lat, lon):
    lat, lon = np.asarray(lat), np.asarray(lon)
    return haversine_km(lat[:, None], lon[:, None], lat[None, :], lon[None, :])


def knn_weights(D, k: int = 8, row_standardize: bool = True):
    """Binary k-nearest-neighbour spatial weights from a distance matrix."""
    n = D.shape[0]
    W = np.zeros_like(D)
    order = np.argsort(D, axis=1)[:, 1:k + 1]  # skip self
    W[np.arange(n)[:, None], order] = 1.0
    if row_standardize:
        W /= W.sum(axis=1, keepdims=True)
    return W


def morans_i(z, W, n_perm: int = 999, rng=None):
    """Global Moran's I spatial autocorrelation.

    I = (n / S0) * (z' W z) / (z' z) with z mean-centred, S0 = sum(W).
    E[I] under H0 = -1/(n-1). Values > E[I] mean similar values cluster in space.
    Significance by random permutation of values over locations (distribution-free).
    Returns dict(I, EI, z_norm, p_perm).
    """
    rng = np.random.default_rng(rng)
    z = np.asarray(z, float)
    z = z - z.mean()
    n = len(z)
    S0 = W.sum()
    denom = z @ z

    def I_of(v):
        return (n / S0) * (v @ W @ v) / denom

    I = I_of(z)
    EI = -1.0 / (n - 1)
    perms = np.array([I_of(rng.permutation(z)) for _ in range(n_perm)])
    extreme = perms >= I if I > EI else perms <= I
    p = (np.sum(extreme) + 1) / (n_perm + 1)
    z_norm = (I - EI) / perms.std(ddof=1)
    return dict(I=float(I), EI=float(EI), z=float(z_norm), p_perm=float(min(1.0, 2 * p)))


def idw(D_query_known, z_known, power: float = 2.0):
    """Inverse-distance-weighted interpolation. D_query_known: (m, n) distances from each
    query point to each known point. Exact at data points (weight -> inf as d -> 0)."""
    D = np.where(D_query_known == 0, 1e-12, D_query_known)
    w = 1.0 / D ** power
    return (w @ z_known) / w.sum(axis=1)


def empirical_variogram(D, z, bins):
    """Semivariance gamma(h) = 1/(2 N(h)) sum (z_i - z_j)^2 over pairs with distance in bin h.
    Returns (bin_centres, gamma, pair_counts). Use only pairs i<j."""
    z = np.asarray(z, float)
    iu = np.triu_indices(len(z), k=1)
    d = D[iu]
    sq = 0.5 * (z[iu[0]] - z[iu[1]]) ** 2
    which = np.digitize(d, bins) - 1
    centres, gamma, counts = [], [], []
    for b in range(len(bins) - 1):
        m = which == b
        if m.sum() >= 30:
            centres.append(0.5 * (bins[b] + bins[b + 1]))
            gamma.append(sq[m].mean())
            counts.append(int(m.sum()))
    return np.array(centres), np.array(gamma), np.array(counts)


def spherical(h, nugget, sill, rng_):
    h = np.asarray(h, float)
    r = np.clip(h / rng_, 0, 1)
    return nugget + (sill - nugget) * (1.5 * r - 0.5 * r ** 3)


def exponential(h, nugget, sill, rng_):
    return nugget + (sill - nugget) * (1 - np.exp(-3 * np.asarray(h, float) / rng_))


def fit_variogram(h, gamma, model=spherical):
    """Weighted least squares fit of a nugget/sill/range model to an empirical variogram."""
    p0 = [gamma.min() * 0.5, gamma.max(), h.max() * 0.6]
    bounds = ([0, 0, 1e-6], [gamma.max(), gamma.max() * 3, h.max() * 5])
    params, _ = curve_fit(model, h, gamma, p0=p0, bounds=bounds, maxfev=20000)
    return params


def ordinary_kriging(D_known, D_query_known, z_known, model, params):
    """Ordinary kriging: BLUE interpolator honouring the variogram.

    Solve [[Gamma, 1],[1', 0]] [w, mu]' = [gamma_0, 1] for each query point, where Gamma is the
    semivariance between known points and gamma_0 between query and known. The prediction is
    w'z; kriging variance = w'gamma_0 + mu. Returns (pred, var)."""
    n = len(z_known)
    G = model(D_known, *params)
    A = np.ones((n + 1, n + 1))
    A[:n, :n] = G
    A[n, n] = 0.0
    g0 = model(D_query_known, *params)  # (m, n)
    rhs = np.ones((n + 1, g0.shape[0]))
    rhs[:n, :] = g0.T
    sol = np.linalg.solve(A, rhs)  # (n+1, m)
    w, mu = sol[:n, :], sol[n, :]
    pred = w.T @ z_known
    var = np.einsum("nm,mn->m", w, g0) + mu
    return pred, var


# ========================================================================= extremes
def fit_gev(x, max_abs_shape: float = 1.0):
    """GEV maximum likelihood with sensible starting values (L-moment-ish) and a sanity
    bound on the shape. scipy's default start can diverge on bootstrap resamples.
    Returns (c, loc, scale) in scipy's convention (c = -xi) or None if the fit is implausible."""
    x = np.asarray(x, float)
    scale0 = np.sqrt(6) * x.std(ddof=1) / np.pi
    loc0 = x.mean() - 0.5772 * scale0
    c, loc, scale = sps.genextreme.fit(x, 0.0, loc=loc0, scale=scale0)
    if not np.isfinite([c, loc, scale]).all() or abs(c) > max_abs_shape or scale <= 0:
        return None
    return float(c), float(loc), float(scale)


def return_level_gev(shape, loc, scale, T):
    """Return level for return period T (years) from a GEV fitted to annual maxima.
    scipy parameterisation: shape = -xi (scipy's c).
    z_T = loc + scale/xi * [(-ln(1-1/T))^(-xi) - 1]."""
    T = np.asarray(T, float)
    p = 1.0 - 1.0 / T
    return sps.genextreme.ppf(p, shape, loc=loc, scale=scale)


def plotting_positions(x, a: float = 0.0):
    """Empirical return periods for sorted annual maxima. a=0 -> Weibull (i/(n+1)),
    a=0.44 -> Gringorten (unbiased for Gumbel). Returns (sorted_x, return_period)."""
    x = np.sort(np.asarray(x, float))[::-1]
    n = len(x)
    rank = np.arange(1, n + 1)
    P_exceed = (rank - a) / (n + 1 - 2 * a)
    return x, 1.0 / P_exceed
