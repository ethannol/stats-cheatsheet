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
