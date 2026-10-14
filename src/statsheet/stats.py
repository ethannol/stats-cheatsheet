"""Hand-rolled implementations of the methods mining data scientists lean on. Short enough to
read in an interview; cross-checked in tests/test_stats.py.

Resource estimation
    cell_declustering, declustering_scan, composite_intervals, top_cut_summary,
    experimental_variogram, variogram_models (spherical/exponential/gaussian), fit_variogram,
    anisotropic_distance, ordinary_kriging, idw, nearest_neighbour, nscore, nscore_back,
    sequential_gaussian_simulation, grade_tonnage
Compositional data
    closure, alr, alr_inv, clr, clr_inv, ilr, ilr_inv, variation_matrix
QA/QC and censored data
    hard, thompson_howarth, km_left_censored, ros_lognormal
Censored inference, trend and compliance limits
    km_mean_se, gehan_test, censored_kendall_S, akritas_theil_sen, censored_regression_mle,
    seasonal_kendall, dunn_test, ucl_table, ucl_recommend, prediction_limit, tolerance_limit
Reliability / SPC
    kaplan_meier, weibull_mle, shewhart_limits, cusum, ewma, u_chart_limits
Sampling theory
    gy_fse_variance, gy_min_mass
Economics
    gbm_paths, fit_ou, correlated_normals
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats as sps
from scipy.optimize import curve_fit, minimize
from scipy.spatial import cKDTree


# ================================================================ declustering
def cell_declustering(x, y, cell_size: float, n_offsets: int = 5) -> np.ndarray:
    """Cell declustering (Journel 1983, Deutsch 1989): overlay a grid of cells of side cell_size;
    each sample gets weight 1/(number of samples in its cell), normalised to sum n. Averaging over
    several grid origin offsets removes the dependence on where the grid starts.
    Clustered samples share weight; isolated samples count fully."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    w = np.zeros(n)
    for ox in np.linspace(0, cell_size, n_offsets, endpoint=False):
        for oy in np.linspace(0, cell_size, n_offsets, endpoint=False):
            ix = np.floor((x - x.min() + ox) / cell_size).astype(int)
            iy = np.floor((y - y.min() + oy) / cell_size).astype(int)
            _, inv, counts = np.unique(ix * 100003 + iy, return_inverse=True, return_counts=True)
            w += 1.0 / counts[inv]
    w /= n_offsets ** 2
    return w * n / w.sum()


def declustering_scan(x, y, z, cell_sizes):
    """Declustered mean of z as a function of cell size. With preferential sampling of high grades
    the curve dips to a minimum: that cell size is the conventional choice."""
    return np.array([np.average(z, weights=cell_declustering(x, y, c)) for c in cell_sizes])


# ================================================================ compositing
def composite_intervals(frm, to, grade, length: float, min_fraction: float = 0.5):
    """Downhole compositing to fixed length: length-weighted average of the samples overlapping
    each composite. A trailing composite shorter than min_fraction*length is merged into the
    previous one. Returns DataFrame-like dict with from, to, grade, n_samples."""
    frm, to, grade = (np.asarray(a, float) for a in (frm, to, grade))
    start, end = frm.min(), to.max()
    edges = np.arange(start, end, length)
    if end - edges[-1] < min_fraction * length and len(edges) > 1:
        edges = edges[:-1]
    edges = np.append(edges, end)
    out = {"from": [], "to": [], "grade": [], "n_samples": []}
    for a, b in zip(edges[:-1], edges[1:]):
        overlap = np.clip(np.minimum(to, b) - np.maximum(frm, a), 0, None)
        m = overlap > 0
        out["from"].append(a)
        out["to"].append(b)
        out["grade"].append(np.sum(overlap[m] * grade[m]) / overlap[m].sum())
        out["n_samples"].append(int(m.sum()))
    return {k: np.array(v) for k, v in out.items()}


def top_cut_summary(z, cap):
    """Effect of capping at `cap`: fraction of samples cut, metal removed, CV before and after."""
    z = np.asarray(z, float)
    zc = np.minimum(z, cap)
    return dict(cap=cap, n_cut=int((z > cap).sum()), frac_cut=float((z > cap).mean()),
                mean_before=z.mean(), mean_after=zc.mean(),
                metal_removed=float(1 - zc.sum() / z.sum()),
                cv_before=z.std(ddof=1) / z.mean(), cv_after=zc.std(ddof=1) / zc.mean())


# ================================================================== variograms
def experimental_variogram(x, y, z, lag: float, n_lags: int, azimuth: float | None = None,
                           tol_deg: float = 22.5, bandwidth: float | None = None):
    """Semivariogram gamma(h) = 1/(2N) sum (z_i - z_j)^2 for pairs binned by separation.
    azimuth (degrees clockwise from +y, mining convention) restricts pairs to a direction cone;
    None gives the omnidirectional variogram. Returns (lag_centres, gamma, pair_counts)."""
    x, y, z = (np.asarray(a, float) for a in (x, y, z))
    i, j = np.triu_indices(len(z), k=1)
    dx, dy = x[j] - x[i], y[j] - y[i]
    d = np.hypot(dx, dy)
    keep = (d > 0) & (d < lag * (n_lags + 0.5))
    if azimuth is not None:
        az = np.radians(azimuth)
        ux, uy = np.sin(az), np.cos(az)  # unit vector of the direction
        along = dx * ux + dy * uy
        perp = np.abs(dx * uy - dy * ux)
        ang_ok = np.abs(along) / np.maximum(d, 1e-12) >= np.cos(np.radians(tol_deg))
        keep &= ang_ok
        if bandwidth is not None:
            keep &= perp <= bandwidth
    d, sq = d[keep], 0.5 * (z[i[keep]] - z[j[keep]]) ** 2
    bins = np.round(d / lag).astype(int)  # lag k collects pairs within +-lag/2 of k*lag
    centres, gamma, counts = [], [], []
    for k in range(1, n_lags + 1):
        m = bins == k
        if m.sum() >= 20:
            centres.append(d[m].mean())
            gamma.append(sq[m].mean())
            counts.append(int(m.sum()))
    return np.array(centres), np.array(gamma), np.array(counts)


def spherical(h, nugget, sill, rng_):
    h = np.asarray(h, float)
    r = np.clip(h / rng_, 0, 1)
    return nugget + (sill - nugget) * (1.5 * r - 0.5 * r ** 3)


def exponential(h, nugget, sill, rng_):
    return nugget + (sill - nugget) * (1 - np.exp(-3 * np.asarray(h, float) / rng_))


def gaussian(h, nugget, sill, rng_):
    return nugget + (sill - nugget) * (1 - np.exp(-3 * (np.asarray(h, float) / rng_) ** 2))


def fit_variogram(h, gamma, model=spherical, counts=None, nugget_max=None):
    """Weighted least squares fit of (nugget, sill, range). Weights ~ pair count / h so short
    lags (which drive kriging weights) matter most."""
    h, gamma = np.asarray(h, float), np.asarray(gamma, float)
    w = (np.asarray(counts, float) if counts is not None else np.ones_like(h)) / h
    sigma = 1 / np.sqrt(w)
    p0 = [max(gamma.min() * 0.5, 1e-9), gamma.max(), h.max() * 0.5]
    hi_n = nugget_max if nugget_max is not None else gamma.max()
    bounds = ([0, 1e-9, h.min() * 0.5], [hi_n, gamma.max() * 3, h.max() * 4])
    params, _ = curve_fit(model, h, gamma, p0=p0, sigma=sigma, bounds=bounds, maxfev=20000)
    return params


def anisotropic_distance(dx, dy, azimuth_major: float, ratio: float):
    """Rotate separation vectors into the anisotropy frame and stretch the minor axis so that a
    single isotropic model can be used. ratio = minor range / major range (<= 1)."""
    az = np.radians(azimuth_major)
    along = dx * np.sin(az) + dy * np.cos(az)
    perp = -dx * np.cos(az) + dy * np.sin(az)
    return np.hypot(along, perp / ratio)


# ===================================================================== kriging
def _cov_from_vario(h, model, params):
    nugget, sill, rng_ = params
    return sill - model(h, nugget, sill, rng_)


def ordinary_kriging(xy_known, z, xy_query, model, params, n_max: int = 24, radius=None,
                     anisotropy=None):
    """Ordinary kriging with a moving neighbourhood of the n_max nearest samples.
    Solves [C 1; 1' 0][w; mu] = [c0; 1] per query point. Returns (estimate, kriging variance).
    anisotropy = (azimuth_major, ratio) applies anisotropic_distance to all separations."""
    xy_known, xy_query = np.asarray(xy_known, float), np.asarray(xy_query, float)
    z = np.asarray(z, float)
    tree = cKDTree(xy_known)
    k = min(n_max, len(z))
    dist, idx = tree.query(xy_query, k=k)
    if k == 1:
        dist, idx = dist[:, None], idx[:, None]
    est, var = np.full(len(xy_query), np.nan), np.full(len(xy_query), np.nan)
    sill = params[1]

    def D(a, b):
        dx, dy = a[:, None, 0] - b[None, :, 0], a[:, None, 1] - b[None, :, 1]
        if anisotropy is not None:
            return anisotropic_distance(dx, dy, *anisotropy)
        return np.hypot(dx, dy)

    for q in range(len(xy_query)):
        sel = idx[q]
        if radius is not None:
            sel = sel[dist[q] <= radius]
        if len(sel) < 2:
            continue
        pts = xy_known[sel]
        C = _cov_from_vario(D(pts, pts), model, params)
        c0 = _cov_from_vario(D(xy_query[q:q + 1], pts), model, params)[0]
        n = len(sel)
        A = np.ones((n + 1, n + 1))
        A[:n, :n] = C
        A[n, n] = 0
        rhs = np.append(c0, 1.0)
        try:
            sol = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            sol = np.linalg.lstsq(A, rhs, rcond=None)[0]
        w, mu = sol[:n], sol[n]
        est[q] = w @ z[sel]
        var[q] = sill - w @ c0 - mu
    return est, var


def idw(xy_known, z, xy_query, power: float = 2.0, n_max: int = 24):
    xy_known, xy_query = np.asarray(xy_known, float), np.asarray(xy_query, float)
    z = np.asarray(z, float)
    dist, idx = cKDTree(xy_known).query(xy_query, k=min(n_max, len(z)))
    dist = np.maximum(dist, 1e-9)
    w = 1 / dist ** power
    return (w * z[idx]).sum(axis=1) / w.sum(axis=1)


def nearest_neighbour(xy_known, z, xy_query):
    _, idx = cKDTree(np.asarray(xy_known, float)).query(np.asarray(xy_query, float), k=1)
    return np.asarray(z, float)[idx]


# ============================================================== normal scores
def nscore(z, weights=None):
    """Normal-score transform: map each value to the standard-normal quantile of its (weighted)
    cumulative probability. Returns (y, table) where table = (sorted z, sorted y) for back-
    transformation."""
    z = np.asarray(z, float)
    order = np.argsort(z)
    w = np.ones_like(z) if weights is None else np.asarray(weights, float)
    cw = np.cumsum(w[order])
    p = (cw - 0.5 * w[order]) / cw[-1]
    y_sorted = sps.norm.ppf(p)
    y = np.empty_like(z)
    y[order] = y_sorted
    return y, (z[order], y_sorted)


def nscore_back(y, table, tail_max=None):
    """Inverse normal-score transform by interpolation in the table; linear extrapolation in the
    upper tail capped at tail_max (default: table maximum)."""
    zs, ys = table
    zmax = zs[-1] if tail_max is None else tail_max
    return np.interp(y, ys, zs, left=zs[0], right=zmax)


# ======================================================== conditional simulation
def sequential_gaussian_simulation(xy_known, y_known, xy_grid, model, params, n_real: int = 1,
                                   n_max: int = 16, seed: int = 0):
    """Sequential Gaussian simulation in normal-score space with simple kriging (mean 0, sill 1).

    For each realisation: visit grid nodes in a random path; at each node krige from the n_max
    nearest *data + previously simulated nodes*, draw from N(estimate, kriging variance), and add
    the node to the conditioning set. Realisations honour the data, the histogram and the
    variogram, unlike the kriged map which is smoothed. Returns array (n_real, n_grid)."""
    rng = np.random.default_rng(seed)
    xy_known, y_known, xy_grid = (np.asarray(a, float) for a in (xy_known, y_known, xy_grid))
    n_grid = len(xy_grid)
    out = np.empty((n_real, n_grid))
    for r in range(n_real):
        path = rng.permutation(n_grid)
        pts = np.vstack([xy_known, xy_grid[path]])
        vals = np.concatenate([y_known, np.full(n_grid, np.nan)])
        n_known = len(y_known)
        tree = cKDTree(xy_known)
        sim_pts, sim_idx = [], []
        for step, node in enumerate(path):
            q = xy_grid[node]
            _, idx_d = tree.query(q, k=min(n_max, n_known))
            idx_d = np.atleast_1d(idx_d)
            cand = list(idx_d)
            if sim_pts:
                sp = np.array(sim_pts)
                d_s = np.hypot(sp[:, 0] - q[0], sp[:, 1] - q[1])
                near = np.argsort(d_s)[:n_max]
                cand += [n_known + sim_idx[k] for k in near]
            cand = np.array(cand)
            d_all = np.hypot(pts[cand, 0] - q[0], pts[cand, 1] - q[1])
            cand = cand[np.argsort(d_all)[:n_max]]
            P = pts[cand]
            C = _cov_from_vario(np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :,
                1]),
                                model, params)
            c0 = _cov_from_vario(np.hypot(P[:, 0] - q[0], P[:, 1] - q[1]), model, params)
            w = np.linalg.solve(C + 1e-9 * np.eye(len(cand)), c0)
            mean = w @ vals[cand]
            var = max(params[1] - w @ c0, 1e-9)
            v = rng.normal(mean, np.sqrt(var))
            vals[n_known + step] = v
            sim_pts.append(q)
            sim_idx.append(step)
        out[r, path] = vals[n_known:]
    return out


# ================================================================ grade-tonnage
def grade_tonnage(block_grades, cutoffs, block_tonnes: float = 1.0):
    """Tonnes, mean grade and metal above each cut-off for a set of block grades."""
    g = np.asarray(block_grades, float)
    g = g[~np.isnan(g)]
    rows = []
    for c in cutoffs:
        above = g[g >= c]
        t = len(above) * block_tonnes
        rows.append(dict(cutoff=c, tonnes=t, grade=above.mean() if len(above) else np.nan,
                         metal=above.sum() * block_tonnes))
    return rows


# ============================================================ compositional data
def closure(X, total: float = 1.0):
    X = np.asarray(X, float)
    return total * X / X.sum(axis=1, keepdims=True)


def alr(X, denom: int = -1):
    X = np.asarray(X, float)
    d = X[:, denom:denom + 1] if denom != -1 else X[:, -1:]
    keep = np.delete(np.arange(X.shape[1]), denom if denom != -1 else X.shape[1] - 1)
    return np.log(X[:, keep] / d)


def alr_inv(Y, denom_last: bool = True):
    Y = np.asarray(Y, float)
    X = np.hstack([np.exp(Y), np.ones((len(Y), 1))])
    return closure(X)


def clr(X):
    """Centred log-ratio: log(x_i / geometric mean). Rows sum to zero; symmetric but singular."""
    L = np.log(np.asarray(X, float))
    return L - L.mean(axis=1, keepdims=True)


def clr_inv(Y):
    return closure(np.exp(np.asarray(Y, float)))


def _helmert(D):
    """Orthonormal basis (D-1 x D) for the ilr transform."""
    V = np.zeros((D - 1, D))
    for i in range(1, D):
        V[i - 1, :i] = 1.0 / i
        V[i - 1, i] = -1.0
        V[i - 1] *= np.sqrt(i / (i + 1))
    return V


def ilr(X):
    """Isometric log-ratio: clr projected on an orthonormal basis; D-1 unconstrained coordinates
    suitable for any multivariate method."""
    X = np.asarray(X, float)
    return clr(X) @ _helmert(X.shape[1]).T


def ilr_inv(Y):
    Y = np.asarray(Y, float)
    return clr_inv(Y @ _helmert(Y.shape[1] + 1))


def variation_matrix(X):
    """Aitchison variation matrix T_ij = var(log(x_i/x_j)); small = proportional parts."""
    L = np.log(np.asarray(X, float))
    D = L.shape[1]
    T = np.zeros((D, D))
    for i in range(D):
        for j in range(D):
            T[i, j] = np.var(L[:, i] - L[:, j], ddof=1)
    return T


# ======================================================================= QA/QC
def hard(a, b):
    """Half absolute relative difference |a-b| / (a+b). 0 = identical, 1 = one of the pair is 0.
    Industry rule of thumb: 90% of field duplicates with HARD < 0.2 (20%) for gold."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.abs(a - b) / (a + b)


def thompson_howarth(a, b, n_bins: int = 11):
    """Thompson-Howarth precision estimate: group duplicate pairs by mean concentration, take the
    median absolute difference per group (x1.048 for sd), and regress sd on concentration:
    sd = s0 + k*c. Returns (bin_means, bin_sd, s0, k)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    mean = 0.5 * (a + b)
    diff = np.abs(a - b)
    order = np.argsort(mean)
    groups = np.array_split(order, n_bins)
    bm = np.array([mean[g].mean() for g in groups])
    bsd = np.array([np.median(diff[g]) * 1.048 for g in groups])
    k, s0 = np.polyfit(bm, bsd, 1)
    return bm, bsd, s0, k


# ============================================================== censored data
def km_left_censored(values, censored):
    """Kaplan-Meier estimate of the mean and quantiles for left-censored (non-detect) data.
    Flip the data (M - x) so non-detects become right-censored, run the right-censored KM
    estimator, flip back. Distribution-free; the standard for environmental non-detects
    (Helsel 2012). Returns dict(mean, median, p25, p75, surv_x, surv_p)."""
    v, c = np.asarray(values, float), np.asarray(censored, bool)
    M = v.max() + 1.0
    t = M - v
    event = ~c
    order = np.argsort(t)
    t, event = t[order], event[order]
    S = 1.0
    xs, Ss = [], []
    for ti in np.unique(t):
        at_risk = (t >= ti).sum()
        d = (event & (t == ti)).sum()
        if d:
            S *= 1 - d / at_risk
        xs.append(ti)
        Ss.append(S)
    xs, Ss = np.array(xs), np.array(Ss)
    # flip back. S_T(t_i) = P(T > t_i) = P(X < x_i) with x_i = M - t_i, i.e. F(x_i-).
    x_back = M - xs[::-1]
    F_minus = Ss[::-1]                                   # P(X < x_k)
    F = np.concatenate([F_minus[1:], [1.0]])             # P(X <= x_k)
    pmf = F - F_minus                                    # P(X = x_k)
    # the mass below the lowest censoring limit has no estimable location; place it at that
    # limit (Helsel 2012 section 6.4), which makes the KM mean a slight upper bound there
    pmf[0] += F_minus[0]
    mean = float(np.sum(pmf * x_back))

    def quant(p):
        k = np.searchsorted(F, p)
        return float(x_back[min(k, len(x_back) - 1)])

    return dict(mean=mean, median=quant(0.5), p25=quant(0.25), p75=quant(0.75), x=x_back, F=F)


def ros_lognormal(values, censored):
    """Regression on order statistics (Helsel & Cohn 1988): fit a lognormal to the detected
    values' probability plot (plotting positions account for the censored ones), impute each
    censored value from the fitted line, then summarise the completed data.
    Returns dict(mean, median, imputed_values)."""
    v, c = np.asarray(values, float), np.asarray(censored, bool)
    n = len(v)
    order = np.argsort(v)
    v, c = v[order], c[order]
    pp = (np.arange(1, n + 1) - 0.375) / (n + 0.25)  # Blom plotting positions over the whole sample
    det = ~c & (v > 0)  # detected zeros carry no log information; they stay as reported
    slope, intercept = np.polyfit(sps.norm.ppf(pp[det]), np.log(v[det]), 1)
    imputed = v.copy()
    imputed[c] = np.exp(intercept + slope * sps.norm.ppf(pp[c]))
    return dict(mean=float(imputed.mean()), median=float(np.median(imputed)),
                imputed_values=imputed, sigma=slope, mu=intercept)


# ==================================================== censored inference and trend
def km_mean_se(values, censored):
    """Standard error of the Kaplan-Meier mean for left-censored data (Greenwood-type
    formula applied to the flipped, right-censored sample; Helsel 2012 section 6.4):
    Var(mean) = sum over event times of A_i^2 d_i / (n_i (n_i - d_i)), where A_i is the area
    under S(t) to the right of t_i. Returns (mean, se)."""
    v, c = np.asarray(values, float), np.asarray(censored, bool)
    M = v.max() + 1.0
    t, event = M - v, ~c
    order = np.argsort(t)
    t, event = t[order], event[order]
    times, S, n_risk, d_at = [0.0], [1.0], [], []
    s = 1.0
    for ti in np.unique(t[event]):
        n_i = (t >= ti).sum()
        d_i = (event & (t == ti)).sum()
        s *= 1 - d_i / n_i
        times.append(ti)
        S.append(s)
        n_risk.append(n_i)
        d_at.append(d_i)
    times, S = np.array(times), np.array(S)
    # S is defined up to the largest observation; a censored tail is set to the last S
    t_end = t.max()
    knots = np.append(times, t_end)
    areas = np.diff(knots) * S              # area of each step
    mean_t = areas.sum()
    var = 0.0
    for i, (n_i, d_i) in enumerate(zip(n_risk, d_at)):
        A_i = areas[i + 1:].sum()
        if n_i > d_i:
            var += A_i ** 2 * d_i / (n_i * (n_i - d_i))
    return float(M - mean_t), float(np.sqrt(var))


def _censored_pair_sign(x_i, c_i, x_j, c_j):
    """Vectorised 'definitely less / greater' sign for left-censored pairs. An observation
    (x, c=True) means 'below x'. Returns +1 if i > j for sure, -1 if i < j for sure, else 0."""
    both_det = ~c_i & ~c_j
    i_cens_only = c_i & ~c_j
    j_cens_only = ~c_i & c_j
    less = (both_det & (x_i < x_j)) | (i_cens_only & (x_i <= x_j))
    greater = (both_det & (x_i > x_j)) | (j_cens_only & (x_i >= x_j))
    return greater.astype(int) - less.astype(int)


def gehan_test(values, censored, group):
    """Gehan's generalised Wilcoxon test for two groups of left-censored data (Gehan 1965;
    Helsel 2012 chapter 9). Each observation scores +1 against every observation it is
    definitely larger than and -1 against every one it is definitely smaller than; pairs
    involving a non-detect that cannot be ordered score 0. W = sum of scores in group 1;
    under H0 the permutation variance is n1 n2 / (n (n-1)) * sum(score^2).
    Without censoring this is the Mann-Whitney test. Returns dict(W, z, p, scores)."""
    v, c, g = np.asarray(values, float), np.asarray(censored, bool), np.asarray(group)
    n = len(v)
    s = _censored_pair_sign(v[:, None], c[:, None], v[None, :], c[None, :])
    scores = s.sum(axis=1)
    labels = np.unique(g)
    if len(labels) != 2:
        raise ValueError("gehan_test needs exactly two groups")
    g1 = g == labels[0]
    n1, n2 = g1.sum(), n - g1.sum()
    W = scores[g1].sum()
    var = n1 * n2 / (n * (n - 1)) * np.sum(scores ** 2)
    z = W / np.sqrt(var)
    return dict(W=float(W), z=float(z), p=float(2 * sps.norm.sf(abs(z))), scores=scores,
                groups=labels)


def censored_kendall_S(x, y, y_censored):
    """Kendall's S between an uncensored x and a left-censored y, counting only pairs whose
    order in y is certain (Helsel 2012 section 10.2)."""
    x, y, c = np.asarray(x, float), np.asarray(y, float), np.asarray(y_censored, bool)
    sy = _censored_pair_sign(y[:, None], c[:, None], y[None, :], c[None, :])
    sx = np.sign(x[:, None] - x[None, :])
    iu = np.triu_indices(len(x), 1)
    return int(np.sum(sx[iu] * sy[iu]))


def akritas_theil_sen(x, y, y_censored, n_perm: int = 2000, rng=None):
    """Akritas-Theil-Sen line for a left-censored response (Akritas, Murphy & LaValley 1995;
    Helsel 2012 section 10.3): the slope b at which censored Kendall's S between x and the
    residuals y - b x is zero, found by root search over the range of detected pairwise slopes.
    The intercept is the median of y - b x treating non-detects at their limit (an upper bound).
    The p-value for the trend is a permutation p-value of S (x shuffled), which sidesteps the
    tie-and-censoring variance formula. Returns dict(slope, intercept, S, tau, p)."""
    x, y, c = np.asarray(x, float), np.asarray(y, float), np.asarray(y_censored, bool)
    rng = np.random.default_rng(rng)
    S0 = censored_kendall_S(x, y, c)
    n_pairs = len(x) * (len(x) - 1) / 2
    perm = np.array([censored_kendall_S(rng.permutation(x), y, c) for _ in range(n_perm)])
    p = (np.sum(np.abs(perm) >= abs(S0)) + 1) / (n_perm + 1)
    det = ~c
    xd, yd = x[det], y[det]
    dx = xd[:, None] - xd[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        slopes = (yd[:, None] - yd[None, :]) / dx
    slopes = slopes[np.isfinite(slopes)]
    lo, hi = slopes.min() - 1e-9, slopes.max() + 1e-9

    def S_of(b):
        return censored_kendall_S(x, y - b * x, c)

    # S(b) is a non-increasing step function of b; bisection to the sign change
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if S_of(mid) > 0:
            lo = mid
        else:
            hi = mid
    slope = 0.5 * (lo + hi)
    intercept = float(np.median(y - slope * x))
    return dict(slope=float(slope), intercept=intercept, S=S0, tau=S0 / n_pairs, p=float(p))


def censored_regression_mle(X, y, censored, add_constant: bool = True):
    """Maximum-likelihood linear regression with a left-censored response (Tobit-type; Helsel
    2012 chapter 12; ProUCL calls it 'MLE regression'). Each detected y contributes the normal
    density of its residual; each non-detect 'below y' contributes the normal CDF of (y - Xb)/s.
    Standard errors come from the inverse of a finite-difference Hessian.
    Returns dict(beta, se, sigma, loglik, aic, names)."""
    X = np.asarray(X, float)
    X = X[:, None] if X.ndim == 1 else X
    if add_constant:
        X = np.column_stack([np.ones(len(X)), X])
    y, c = np.asarray(y, float), np.asarray(censored, bool)
    k = X.shape[1]
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
    s0 = np.std(y - X @ beta0) + 1e-6
    theta0 = np.append(beta0, np.log(s0))

    def nll(theta):
        b, s = theta[:k], np.exp(theta[k])
        z = (y - X @ b) / s
        ll_det = sps.norm.logpdf(z[~c]) - np.log(s)
        ll_cen = sps.norm.logcdf(z[c])
        return -(ll_det.sum() + ll_cen.sum())

    res = minimize(nll, theta0, method="BFGS")
    res = minimize(nll, res.x, method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10})
    theta = res.x
    # central-difference Hessian of the negative log-likelihood
    h = 1e-4 * np.maximum(1.0, np.abs(theta))
    H = np.zeros((k + 1, k + 1))
    f0 = nll(theta)
    for i in range(k + 1):
        for j in range(i, k + 1):
            ei, ej = np.zeros(k + 1), np.zeros(k + 1)
            ei[i], ej[j] = h[i], h[j]
            if i == j:
                H[i, i] = (nll(theta + ei) - 2 * f0 + nll(theta - ei)) / h[i] ** 2
            else:
                mixed = (nll(theta + ei + ej) - nll(theta + ei - ej)
                         - nll(theta - ei + ej) + nll(theta - ei - ej))
                H[i, j] = H[j, i] = mixed / (4 * h[i] * h[j])
    cov = np.linalg.pinv(H)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    loglik = -f0
    return dict(beta=theta[:k], se=se[:k], sigma=float(np.exp(theta[k])), loglik=float(loglik),
                aic=float(2 * (k + 1) - 2 * loglik), converged=bool(res.success))


def seasonal_kendall(values, years, seasons):
    """Seasonal Kendall test (Hirsch, Slack & Smith 1982; Helsel et al. 2020 section 12.4).
    Within each season compare values across years (Kendall's S with tie correction), sum the
    S and the variances over seasons, apply the continuity correction. Multiple observations in
    the same season-year are reduced to their median first. Sen slope = median of the
    within-season pairwise slopes (units per year). Returns dict(S, var, z, p, tau, sen,
    per_season DataFrame)."""
    import pandas as pd
    df = (pd.DataFrame(dict(v=np.asarray(values, float), y=np.asarray(years),
                            s=np.asarray(seasons)))
          .groupby(["s", "y"], as_index=False).v.median())
    S_tot, var_tot, n_pairs, rows, slopes = 0.0, 0.0, 0, [], []
    for s, g in df.groupby("s"):
        g = g.sort_values("y")
        yv, vv = g.y.values.astype(float), g.v.values
        n = len(vv)
        if n < 2:
            continue
        dv = vv[None, :] - vv[:, None]      # v_j - v_i for i < j
        dy = yv[None, :] - yv[:, None]
        iu = np.triu_indices(n, 1)
        S = float(np.sum(np.sign(dv[iu])))
        _, counts = np.unique(vv, return_counts=True)
        var = (n * (n - 1) * (2 * n + 5) - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18
        with np.errstate(divide="ignore", invalid="ignore"):
            sl = (dv / dy)[iu]
        slopes.extend(sl[np.isfinite(sl)].tolist())
        S_tot += S
        var_tot += var
        n_pairs += n * (n - 1) / 2
        ok = np.isfinite(sl)
        rows.append(dict(season=s, n_years=n, S=S, tau=S / (n * (n - 1) / 2),
                         sen=float(np.median(sl[ok])) if ok.any() else np.nan))
    z = (S_tot - np.sign(S_tot)) / np.sqrt(var_tot) if var_tot > 0 else 0.0
    return dict(S=S_tot, var=var_tot, z=float(z), p=float(2 * sps.norm.sf(abs(z))),
                tau=S_tot / n_pairs if n_pairs else np.nan, sen=float(np.median(slopes)),
                per_season=pd.DataFrame(rows))


def dunn_test(groups, labels=None, adjust: str = "holm"):
    """Dunn's (1964) pairwise rank-sum comparisons after Kruskal-Wallis, with the tie
    correction and Holm or Benjamini-Hochberg adjustment. Returns a DataFrame of pairs."""
    import pandas as pd
    groups = [np.asarray(g, float) for g in groups]
    labels = list(labels) if labels is not None else list(range(len(groups)))
    allv = np.concatenate(groups)
    N = len(allv)
    ranks = sps.rankdata(allv)
    _, counts = np.unique(allv, return_counts=True)
    tie_term = np.sum(counts ** 3 - counts) / (12 * (N - 1))
    idx = np.cumsum([0] + [len(g) for g in groups])
    rbar = [ranks[idx[i]:idx[i + 1]].mean() for i in range(len(groups))]
    rows = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            se = np.sqrt((N * (N + 1) / 12 - tie_term) * (1 / len(groups[i]) + 1 / len(groups[j])))
            z = (rbar[i] - rbar[j]) / se
            rows.append(dict(a=labels[i], b=labels[j], z=z, p=2 * sps.norm.sf(abs(z))))
    out = pd.DataFrame(rows)
    if adjust == "bh":
        out["p_adj"] = sps.false_discovery_control(out.p, method="bh")
    else:  # Holm step-down
        order = np.argsort(out.p.values)
        stepped = np.minimum(1, np.maximum.accumulate(
            out.p.values[order] * np.arange(len(out), 0, -1)))
        out["p_adj"] = stepped[np.argsort(order)]
    return out


def ucl_table(values, censored=None, alpha: float = 0.05, n_boot: int = 2000, rng=None):
    """Upper confidence limits for the mean (the ProUCL menu; EPA 2022 chapter 2).
    Uncensored: Student-t, approximate gamma (Grice & Bain adjusted shape), Chebyshev at 95 /
    97.5 / 99 %, bootstrap-t, BCa bootstrap. Censored: Kaplan-Meier mean and SE, then
    KM-t and KM-Chebyshev. Land's H-UCL is deliberately omitted (ProUCL stopped recommending
    it in 2009). Returns a dict name -> UCL."""
    v = np.asarray(values, float)
    c = np.zeros(len(v), bool) if censored is None else np.asarray(censored, bool)
    rng = np.random.default_rng(rng)
    n = len(v)
    out = {}
    if c.any():
        mean, se = km_mean_se(v, c)
        out["KM mean"] = mean
        out["KM (t) UCL95"] = mean + sps.t.ppf(1 - alpha, n - 1) * se
        for q in (0.95, 0.975, 0.99):
            out[f"KM Chebyshev UCL{q * 100:g}"] = mean + np.sqrt(1 / (1 - q) - 1) * se
        return out
    xbar, s = v.mean(), v.std(ddof=1)
    se = s / np.sqrt(n)
    out["mean"] = xbar
    out["Student-t UCL95"] = xbar + sps.t.ppf(1 - alpha, n - 1) * se
    # gamma: MLE shape with the Grice-Bain small-sample correction (ProUCL 'adjusted gamma')
    k_hat, _, scale_hat = sps.gamma.fit(v, floc=0)
    k_adj = (n - 3) * k_hat / n + 2 / (3 * n)
    out["Gamma UCL95 (adjusted)"] = 2 * n * xbar * k_adj / sps.chi2.ppf(alpha, 2 * n * k_adj)
    for q in (0.95, 0.975, 0.99):
        out[f"Chebyshev UCL{q * 100:g}"] = xbar + np.sqrt(1 / (1 - q) - 1) * se
    boots = rng.choice(v, (n_boot, n), replace=True)
    bm = boots.mean(axis=1)
    bs = boots.std(axis=1, ddof=1) / np.sqrt(n)
    t_star = (bm - xbar) / np.where(bs > 0, bs, np.nan)
    out["Bootstrap-t UCL95"] = xbar - np.nanquantile(t_star, alpha) * se
    bca = sps.bootstrap((v,), np.mean, confidence_level=1 - 2 * alpha, method="BCa",
                        n_resamples=n_boot, random_state=rng)
    out["BCa bootstrap UCL95"] = float(bca.confidence_interval.high)
    return out


def ucl_recommend(values, censored=None):
    """A condensed version of the ProUCL decision tables (EPA 2022, chapter 2): which UCL to
    report, by distribution fit and the skewness of the logs. ProUCL's full tables also move
    the thresholds with n; this keeps the ranking and the rationale, not every cell.
    Returns dict(choice, rationale, sd_log, fits)."""
    v = np.asarray(values, float)
    c = np.zeros(len(v), bool) if censored is None else np.asarray(censored, bool)
    det = v[~c & (v > 0)]
    sd_log = float(np.std(np.log(det), ddof=1)) if len(det) > 2 else np.nan
    fits = {}
    if c.any():
        if sd_log < 0.5:
            return dict(choice="KM (t) UCL95", sd_log=sd_log, fits=fits,
                        rationale="non-detects present, mild skew: KM mean with a t-interval")
        q = {True: "95", False: "97.5"}[sd_log < 1.0] if sd_log < 1.5 else "99"
        return dict(choice=f"KM Chebyshev UCL{q}", sd_log=sd_log, fits=fits,
                    rationale=f"non-detects present, sd of log detects {sd_log:.2f}: KM mean "
                              "with a Chebyshev interval, widened as skew grows")
    fits["normal p"] = float(sps.shapiro(v).pvalue) if len(v) <= 5000 else np.nan
    k_hat, _, scale_hat = sps.gamma.fit(v, floc=0)
    fits["gamma p (KS)"] = float(sps.kstest(v, "gamma", args=(k_hat, 0, scale_hat)).pvalue)
    fits["lognormal p"] = float(sps.shapiro(np.log(v)).pvalue)
    if fits["normal p"] > 0.05:
        return dict(choice="Student-t UCL95", sd_log=sd_log, fits=fits,
                    rationale="normality not rejected: the t-interval is exact")
    if fits["gamma p (KS)"] > 0.05:
        return dict(choice="Gamma UCL95 (adjusted)", sd_log=sd_log, fits=fits,
                    rationale="gamma fit not rejected: the gamma UCL is stable where "
                              "lognormal-based limits explode")
    if sd_log < 0.5:
        choice, why = "Student-t UCL95", "skew is mild (sd of logs < 0.5): t-interval, or BCa"
    elif sd_log < 1.0:
        choice, why = "Chebyshev UCL95", "moderate skew: nonparametric Chebyshev"
    elif sd_log < 1.5:
        choice, why = "Chebyshev UCL97.5", "strong skew: Chebyshev at 97.5 % for nominal coverage"
    else:
        choice, why = "Chebyshev UCL99", "extreme skew: Chebyshev at 99 %; consider bootstrap-t"
    return dict(choice=choice, sd_log=sd_log, fits=fits, rationale=why)


def prediction_limit(background, k_future: int = 1, alpha: float = 0.05):
    """Upper prediction limits for the next k_future samples (EPA Unified Guidance 2009,
    chapters 18-19). Parametric: xbar + t s sqrt(1 + 1/n) (Bonferroni over k). Nonparametric:
    the background maximum, whose confidence for k future values is n / (n + k) per sample
    set. Returns dict(normal, nonparametric, nonparametric_confidence)."""
    x = np.asarray(background, float)
    n = len(x)
    t = sps.t.ppf(1 - alpha / k_future, n - 1)
    normal = x.mean() + t * x.std(ddof=1) * np.sqrt(1 + 1 / n)
    return dict(normal=float(normal), nonparametric=float(x.max()),
                nonparametric_confidence=n / (n + k_future))


def tolerance_limit(background, coverage: float = 0.95, confidence: float = 0.95):
    """One-sided upper normal tolerance limit covering `coverage` of the population with
    `confidence` (noncentral-t factor; Unified Guidance chapter 17)."""
    x = np.asarray(background, float)
    n = len(x)
    k = sps.nct.ppf(confidence, n - 1, sps.norm.ppf(coverage) * np.sqrt(n)) / np.sqrt(n)
    return float(x.mean() + k * x.std(ddof=1))


# ================================================================= reliability
def kaplan_meier(time, event):
    """Right-censored Kaplan-Meier survival curve. Returns (times, survival, at_risk)."""
    t, e = np.asarray(time, float), np.asarray(event, bool)
    order = np.argsort(t)
    t, e = t[order], e[order]
    times, S, risk = [0.0], [1.0], [len(t)]
    s = 1.0
    for ti in np.unique(t[e]):
        n_risk = (t >= ti).sum()
        d = (e & (t == ti)).sum()
        s *= 1 - d / n_risk
        times.append(ti)
        S.append(s)
        risk.append(n_risk)
    return np.array(times), np.array(S), np.array(risk)


@dataclass
class WeibullFit:
    shape: float
    scale: float
    loglik: float

    @property
    def mttf(self):
        from scipy.special import gamma as G
        return self.scale * G(1 + 1 / self.shape)

    def hazard(self, t):
        t = np.asarray(t, float)
        return self.shape / self.scale * (t / self.scale) ** (self.shape - 1)

    def reliability(self, t):
        return np.exp(-(np.asarray(t, float) / self.scale) ** self.shape)


def weibull_mle(time, event) -> WeibullFit:
    """Weibull maximum likelihood with right censoring. beta < 1 infant mortality (decreasing
    hazard), beta = 1 random/exponential, beta > 1 wear-out. Censored units contribute the
    survival probability only."""
    t, e = np.asarray(time, float), np.asarray(event, bool)
    t = np.maximum(t, 1e-9)

    def nll(p):
        beta, eta = np.exp(np.clip(p, -20, 20))
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            logpdf = np.log(beta / eta) + (beta - 1) * np.log(t / eta) - (t / eta) ** beta
            logsf = -(t / eta) ** beta
            val = -(np.sum(logpdf[e]) + np.sum(logsf[~e]))
        return val if np.isfinite(val) else 1e300

    best = None
    for b0 in (0.7, 1.5, 3.0):  # several starts: the likelihood can be flat for beta near 1
        res = minimize(nll, [np.log(b0), np.log(t.mean())], method="Nelder-Mead",
                       options=dict(xatol=1e-8, fatol=1e-10, maxiter=5000))
        if best is None or res.fun < best.fun:
            best = res
    beta, eta = np.exp(best.x)
    return WeibullFit(float(beta), float(eta), float(-best.fun))


# ================================================================ control charts
def shewhart_limits(x, phase1: int | None = None):
    """Individuals (I-MR) chart: centre = mean, sigma from the average moving range / 1.128,
    limits at +-3 sigma. Estimate from the phase-1 window only."""
    x = np.asarray(x, float)
    ref = x if phase1 is None else x[:phase1]
    mr = np.abs(np.diff(ref))
    sigma = mr.mean() / 1.128
    centre = ref.mean()
    return centre, sigma, centre - 3 * sigma, centre + 3 * sigma


def cusum(x, target, sigma, k: float = 0.5, h: float = 5.0, reset: bool = True):
    """Tabular CUSUM: C+_t = max(0, C+_{t-1} + (x_t - target)/sigma - k); signal when > h.
    k = half the shift (in sigma) you want to detect, h ~ 4-5 for ARL0 of several hundred.
    Detects small sustained shifts (~1 sigma) far faster than a Shewhart chart. After a signal
    both statistics restart at zero (reset=True). Returns (C+, C-, signal indices)."""
    z = (np.asarray(x, float) - target) / sigma
    cp, cm = np.zeros(len(z)), np.zeros(len(z))
    signals = []
    for i in range(1, len(z)):
        cp[i] = max(0, cp[i - 1] + z[i] - k)
        cm[i] = max(0, cm[i - 1] - z[i] - k)
        if cp[i] > h or cm[i] > h:
            signals.append(i)
            if reset:
                cp[i] = cm[i] = 0.0
    return cp, cm, np.array(signals, dtype=int)


def ewma(x, lam: float = 0.2, target=None, sigma=None, L: float = 2.7):
    """Exponentially weighted moving average chart with time-varying limits."""
    x = np.asarray(x, float)
    target = x.mean() if target is None else target
    sigma = np.std(x, ddof=1) if sigma is None else sigma
    z = np.empty(len(x))
    z[0] = target
    for i in range(1, len(x)):
        z[i] = lam * x[i] + (1 - lam) * z[i - 1]
    t = np.arange(1, len(x) + 1)
    half = L * sigma * np.sqrt(lam / (2 - lam) * (1 - (1 - lam) ** (2 * t)))
    return z, target - half, target + half


def u_chart_limits(counts, exposure, scale: float = 1.0):
    """u-chart for rates (events per unit exposure, Poisson): centre u = sum c / sum n;
    limits u +- 3 sqrt(u / n_i) (vary with exposure)."""
    c, n = np.asarray(counts, float), np.asarray(exposure, float) / scale
    u = c.sum() / n.sum()
    se = np.sqrt(u / n)
    return u, np.maximum(u - 3 * se, 0), u + 3 * se


# ============================================================== sampling theory
def _liberation(d_cm, liberation, liberation_size_cm, alpha):
    if liberation_size_cm is not None:
        return np.minimum(1.0, (liberation_size_cm / np.asarray(d_cm, float)) ** (3 - alpha))
    return 1.0 if liberation is None else liberation


def gy_fse_variance(d_cm, mass_g, c_g_cm3, f=0.5, g=0.25, liberation=None,
                    liberation_size_cm=None, alpha=3.0):
    """Gy's fundamental sampling error: relative variance s^2 = c f g l d^3 / M_s
    (sample mass only; lot mass assumed >> sample). d = nominal top size (cm), c = mineralogical
    constant (g/cm3), f shape factor (0.5), g granulometric factor (0.25 for crushed material),
    l liberation factor (0-1).
    With liberation_size_cm given, l = (d_l / d)^(3 - alpha) capped at 1, so that
    s^2 = c f g d_l^(3-alpha) d^alpha / M: alpha = 2.5 is Gy's classic sqrt(d_l/d) liberation
    factor, alpha = 1.5 is the Francois-Bongarcon calibration that fits gold deposits.
    Returns relative standard deviation (fraction)."""
    l = _liberation(d_cm, liberation, liberation_size_cm, alpha)
    return np.sqrt(c_g_cm3 * f * g * l * np.asarray(d_cm, float) ** 3 / np.asarray(mass_g, float))


def gy_min_mass(d_cm, target_rsd, c_g_cm3, f=0.5, g=0.25, liberation=None, liberation_size_cm=None,
                alpha=3.0):
    """Minimum sample mass (g) to hold the fundamental error at target_rsd for top size d."""
    l = _liberation(d_cm, liberation, liberation_size_cm, alpha)
    return c_g_cm3 * f * g * l * np.asarray(d_cm, float) ** 3 / target_rsd ** 2


def mineralogical_constant(grade_frac, density_mineral, density_gangue):
    """c = ((1-a)/a) [ (1-a) rho_m + a rho_g ] for critical content a (mass fraction of mineral).
    Low-grade gold: a ~ 1e-6 gives c in the millions, which is why gold needs kilograms."""
    a = grade_frac
    return (1 - a) / a * ((1 - a) * density_mineral + a * density_gangue)


# ================================================================== economics
def gbm_paths(s0, mu, sigma, n_steps, n_paths, dt=1 / 12, rng=None):
    """Geometric Brownian motion S_t+dt = S_t exp((mu - sigma^2/2) dt + sigma sqrt(dt) Z)."""
    rng = np.random.default_rng(rng)
    z = rng.normal(size=(n_paths, n_steps))
    logret = (mu - 0.5 * sigma ** 2) * dt + sigma * np.sqrt(dt) * z
    return s0 * np.exp(np.cumsum(logret, axis=1))


def fit_ou(log_price, dt=1 / 12):
    """Fit a mean-reverting (Ornstein-Uhlenbeck) process to log prices via the AR(1) regression
    x_{t+1} = a + b x_t + e. kappa = -ln(b)/dt, long-run mean = a/(1-b), half-life = ln2/kappa.
    Returns dict(kappa, mean, sigma, half_life_years, b, p_unit_root)."""
    x = np.asarray(log_price, float)
    x = x[~np.isnan(x)]
    b, a = np.polyfit(x[:-1], x[1:], 1)
    resid = x[1:] - (a + b * x[:-1])
    kappa = -np.log(b) / dt if 0 < b < 1 else np.nan
    if 0 < b < 1:
        sigma = resid.std(ddof=2) * np.sqrt(2 * kappa / (1 - b ** 2))
    else:
        sigma = resid.std(ddof=2) / np.sqrt(dt)
    from statsmodels.tsa.stattools import adfuller
    p_ur = adfuller(x, regression="c", autolag="AIC", result_object=False)[1]
    return dict(kappa=kappa, mean=a / (1 - b) if 0 < b < 1 else np.nan, sigma=sigma,
                half_life_years=np.log(2) / kappa if kappa == kappa else np.nan, b=b,
                p_unit_root=p_ur)


def correlated_normals(n, corr, rng=None):
    """Draw n samples of standard normals with the given correlation matrix (Cholesky)."""
    rng = np.random.default_rng(rng)
    L = np.linalg.cholesky(np.asarray(corr, float))
    return rng.normal(size=(n, len(corr))) @ L.T
