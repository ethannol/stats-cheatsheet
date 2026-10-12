# Statistics for earth-science data: the cheat sheet

Condensed from the eleven scripts in `analyses/`. Section numbers match the scripts.

## 0. Which method?

| Question | Method | Where |
|---|---|---|
| What does this variable look like? | histogram, QQ plot, median/IQR/MAD, skew; climatology + anomaly | [01](analyses/01_eda.py) |
| Are two groups different? | Welch t (means), Mann-Whitney (ranks), permutation (anything); report effect size + CI | [02](analyses/02_hypothesis_tests.py) |
| Did a proportion change? | two-proportion z / chi-square; Fisher if small | [02](analyses/02_hypothesis_tests.py) |
| Many tests at once | Bonferroni (FWER) or Benjamini-Hochberg (FDR) | [02](analyses/02_hypothesis_tests.py) |
| How many samples do I need? | power analysis (effect size / SD, alpha, power) | [02](analyses/02_hypothesis_tests.py) |
| How does y depend on x? | OLS + residual diagnostics; HAC/GLS if autocorrelated; GLM if y is counts/positive | [03](analyses/03_regression.py) |
| Which of several predictors matter? | multiple regression, standardised coefficients, VIF, AIC/BIC | [03](analyses/03_regression.py) |
| Trend + seasonality in a series | STL decomposition; ADF + KPSS; difference; SARIMA with a holdout | [04](analyses/04_time_series.py) |
| Hidden periodicities | periodogram / Welch PSD after detrending; test against red noise | [04](analyses/04_time_series.py) |
| Does x lead y? | cross-correlation of *detrended* (or pre-whitened) series | [04](analyses/04_time_series.py) |
| Is there a monotonic trend? | Mann-Kendall + Sen's slope, with Hamed-Rao correction; OLS with HAC SE | [05](analyses/05_trends_changepoints.py) |
| When did it change? | Pettitt (shift), hinge regression (slope change), CUSUM; PELT for many | [05](analyses/05_trends_changepoints.py) |
| How big is the 100-year event? | GEV on annual maxima or GPD on peaks over threshold; bootstrap the CI | [06](analyses/06_extremes.py) |
| Are extremes changing? | nonstationary GEV (location ~ time), likelihood-ratio test | [06](analyses/06_extremes.py) |
| Uncertainty with no formula | bootstrap (BCa); block bootstrap for dependent data; jackknife for bias | [07](analyses/07_resampling.py) |
| Is my test calibrated? | Monte Carlo from a fitted null (e.g. AR(1)) | [05](analyses/05_trends_changepoints.py), [07](analyses/07_resampling.py) |
| Reduce many variables | standardise, PCA (loadings, scree), then cluster (k-means / Ward) | [08](analyses/08_multivariate.py) |
| Do nearby places resemble each other? | Moran's I on the variable and on model residuals | [09](analyses/09_spatial.py) |
| Fill in a map between stations | variogram -> kriging (regression kriging with covariates); validate LOO | [09](analyses/09_spatial.py) |
| Probability statements about parameters | Bayesian: conjugate update or MCMC; credible intervals, posterior predictive | [10](analyses/10_bayesian.py) |
| Predict, and be honest about skill | time-aware CV, baselines, permutation importance, PR-AUC for rare events | [11](analyses/11_ml.py) |

## 1. Descriptive statistics and EDA

- Audit first: coverage, gaps (reindex to a full calendar), units, sentinel values (-9.99, -999999).
  Missingness type: MCAR / MAR / MNAR (a gauge that fails during floods is MNAR).
- Location and spread: mean & SD for symmetric data; median, IQR, scaled MAD (x 1.4826) otherwise.
  Precipitation is zero-inflated and heavy-tailed: summarise dry-day fraction + wet-day
  distribution separately; log or gamma on wet days.
- Climatology: day-of-year mean smoothed over ~31 days (or first few Fourier harmonics).
  Anomaly = value minus climatology. Do everything downstream on anomalies.
- Correlation: Pearson (linear, outlier-sensitive), Spearman (rank), Kendall (small n). Shared
  seasonality manufactures correlation; correlate anomalies.
- Normality tests reject everything at large n. Use QQ plots and ask whether the departure
  matters for the method (means: CLT helps; tails: it does not).
- Daily environmental series have lag-1 autocorrelation 0.6-0.95. Row count is not sample size.

## 2. Hypothesis testing

- p-value = P(statistic at least this extreme | H0). Not P(H0), not effect size, not importance.
  Non-significant = no evidence of difference, not evidence of none. Report the CI.
- Choose the independent unit: seasonal / annual aggregates, not days. Daily tests on the same
  effect gave p = 1e-17 vs p = 5e-4 (02).
- Welch t by default (Student only if variances are shown equal; Levene to test). Mann-Whitney
  tests stochastic dominance, not medians. Permutation needs only exchangeability.
- Paired design: paired t / Wilcoxon signed-rank. Variance change: Levene. Distribution change: KS.
- Effect size: Cohen's d (0.2 / 0.5 / 0.8), difference with CI, odds ratio for 2x2.
- Multiple comparisons: m tests at 0.05 gives 1-(0.95)^m family-wise error (46% for m = 12).
  Bonferroni = alpha/m; Benjamini-Hochberg controls the expected false-discovery fraction.
- Power = 1 - beta, scales with (effect/SD) x sqrt(n). Decide n before, not "observed power" after.
  Detecting a 0.25 C summer shift at one station needs > 150 years per group.

## 3. Regression

- OLS assumptions (LINE): Linear mean, Independent errors, Normal errors (small n only), Equal
  variance. X need not be normal. Check: residuals vs fitted, QQ, residuals vs time/space,
  Durbin-Watson (2 = none), Ljung-Box, Breusch-Pagan.
- Autocorrelated residuals: coefficients unbiased, SEs too small by ~sqrt(n/n_eff) with
  n_eff = n(1-r1)/(1+r1). Fixes: Newey-West HAC SEs, GLS-AR(1) (Cochrane-Orcutt), or aggregate.
- Heteroscedasticity: HC (White) robust SEs, transform y, or a GLM with the right variance
  function (Gamma for positive skewed, Poisson/negative binomial for counts).
- Multicollinearity: VIF = 1/(1-R2_j); > 5-10 is trouble. Predictions fine, coefficients unstable
  (ln CO2 coefficient doubled when 'year' was added; VIF = 91).
- Model comparison: adjusted R2, AIC (prediction-oriented), BIC (parsimony-oriented), on the same y.
  Count estimated breakpoints and tuned hyperparameters as parameters.
- Omitted-variable bias: the ln(CO2) slope absorbs methane, aerosols and everything else that
  trended. Observational regression describes; it does not attribute.
- Partial residual (component + residual) plots check each predictor's linearity.

## 4. Time series

- Decompose: STL (loess; seasonal pattern may evolve) vs classical (fixed seasonal). CO2 seasonal
  amplitude grew 15% since the 1960s; STL sees it, classical averages it away.
- Stationarity: ADF H0 = unit root; KPSS H0 = stationary. Both together: ADF rejects and KPSS
  does not -> stationary. Differencing: d = 1 for trend, D = 1 at lag s for seasonality.
  Over-differencing shows as ACF(1) near -0.5.
- Identification: ACF tails off / PACF cuts at p -> AR(p); ACF cuts at q -> MA(q). Spikes at
  lags 1 and s after (1-B)(1-B^s): the airline model SARIMA(0,1,1)(0,1,1)s.
- Evaluate on a time-ordered holdout and against naive / seasonal naive / climatology baselines.
  Rolling-origin CV for honest error growth with horizon. Check interval coverage.
- Residuals should be white: Ljung-Box p > 0.05 at lags s and 2s.
- Spectra: detrend first (polynomial if the trend is curved). Red noise (AR(1)) has power at low
  frequencies; a 'cycle' must beat the red-noise background to be real. ENSO shows at 3-7 years.
- Cross-correlation of two trending series is meaningless; detrend or pre-whiten both.
- Exponential smoothing (Holt-Winters) is a strong, cheap competitor to SARIMA.

## 5. Trends and change points

- Report slope per decade with a CI that accounts for autocorrelation. Print r1 of the residuals
  next to every trend.
- Mann-Kendall: S = sum of pairwise signs; Var(S) = n(n-1)(2n+5)/18 minus tie terms; z with
  continuity correction. Sen's slope = median of pairwise slopes (29% breakdown point).
- Autocorrelation inflates MK false positives (33% at rho = 0.6 for nominal 5%). Fixes compared
  by simulation: simple pre-whitening (right size, loses power), Yue's trend-free pre-whitening
  (worse than raw under H0!), Hamed-Rao variance correction (size ~9%, best trade-off).
- Pettitt: single shift in the median; on a gradual trend it 'finds' a change mid-series.
  Hinge (broken-stick) regression estimates a slope change (global warming: 1975).
  CUSUM visualises regime shifts. PELT / binary segmentation for multiple change points.
- Endpoint sensitivity: trend CIs explode for windows < 20 years; 1998 start = 'hiatus'.
  Fix the window before looking.
- Hydrology: test several annual statistics (mean, low-flow quantile, maximum). That is multiple
  testing.

## 6. Extremes

- Block maxima -> GEV(mu, sigma, xi). xi > 0 Frechet (floods, rain), xi = 0 Gumbel, xi < 0
  Weibull bounded (temperature). scipy: `genextreme(c)` with c = -xi.
- Return level z_T = mu + sigma/xi [(-ln(1-1/T))^(-xi) - 1]; the T-year event has annual
  exceedance probability 1/T and a 1-(1-1/T)^T (~63%) chance over T years.
- Plotting positions: Weibull i/(n+1), Gringorten (i-0.44)/(n+0.12). US floods: log-Pearson III
  (Bulletin 17C).
- Peaks over threshold -> GPD(sigma_u, xi). Choose u where the mean-residual-life plot is linear
  and the shape / modified scale are stable. Decluster (runs declustering) so events are
  independent. z_T = u + sigma/xi [(lambda T)^xi - 1], lambda = events per year.
- Uncertainty: parametric or nonparametric bootstrap; profile likelihood. With ~100 years the
  100-year level has roughly +-30-40% uncertainty and the sign of xi may not be settled.
- Nonstationarity: mu(t) = mu0 + mu1 t (or a covariate), fit by maximising the log-likelihood
  with an optimiser, compare with a likelihood-ratio test (chi-square, df = extra parameters).
  Return levels then drift; "the 100-year event" needs a reference year.

## 7. Resampling

- Bootstrap: resample rows with replacement, recompute. Percentile CI is simple; basic reflects
  it; BCa corrects bias and skew (default). All under-cover for the mean of heavy-tailed data at
  small n (coverage 84-88% at n = 30 lognormal).
- Dependent data: moving-block bootstrap (block ~ several decorrelation times or n^(1/3)),
  stationary bootstrap, residual bootstrap from a fitted ARMA. iid bootstrap on autocorrelated
  data reproduces the naive OLS SE exactly.
- Regression: case (pairs) bootstrap is robust to heteroscedasticity; residual bootstrap assumes
  the model is right.
- Jackknife: leave-one-out; bias = (n-1)(mean(theta_-i) - theta); SE from the spread. Not for
  medians / quantiles.
- Permutation test: shuffle labels; exact under exchangeability.
- Monte Carlo: simulate from a fitted null (AR(1) with the residual phi) to get the real null
  distribution of your statistic. Null SD of the trend t-statistic was 1.6 instead of 1.0.
- Always set the seed and report the number of resamples.

## 8. Multivariate

- Standardise unless units are shared and variance weighting is intended.
- PCA: orthogonal directions of maximal variance; loadings interpret, scores locate. Keep
  components by scree elbow / cumulative variance / Kaiser (eigenvalue > 1). Signs are
  arbitrary. Linear; outlier-sensitive. On a space x time field this is EOF analysis.
- k-means: needs k, spherical equal-variance clusters, multiple inits. Elbow (inertia) and
  silhouette rarely agree sharply on a continuous gradient: say the choice is a judgement.
- Hierarchical (Ward) gives a dendrogram; compare partitions with the adjusted Rand index.
- Gaussian mixtures for soft membership; DBSCAN for arbitrary shapes and noise.
- Validate clusters against external knowledge (geography, Koppen classes).

## 9. Spatial

- Spatial autocorrelation violates independence exactly as temporal autocorrelation does.
- Moran's I = (n/S0) z'Wz / z'z with W from kNN, distance band or contiguity; E[I] = -1/(n-1);
  permutation p-value. Test the *residuals* of your model, not just the raw field.
- Variogram gamma(h) = half mean squared difference at separation h. Nugget (micro-scale +
  measurement error), sill (variance), range (decorrelation distance). Models: spherical,
  exponential, Gaussian. Fit only to about half the domain width. A variogram that never
  reaches a sill means trend remains: improve the mean model.
- Interpolation: IDW (deterministic, bulls-eyes, no uncertainty); ordinary kriging (BLUE given
  the variogram, gives a variance, assumes a stationary mean); regression / universal kriging
  (covariate trend + kriged residuals). Validate by leave-one-out, never by fit at data points.
- Use great-circle or projected distances. Check anisotropy. Use a DEM for elevation covariates.

## 10. Bayesian

- posterior ∝ likelihood x prior. Weakly informative priors by default; flat is a choice too.
- Credible interval: P(theta in interval | data). Confidence interval: long-run coverage of the
  procedure. With much data and weak priors they coincide numerically.
- Conjugate pairs: Beta-Binomial, Gamma-Poisson, Normal-Normal, Normal-Inverse-Gamma.
  Posterior sampling makes derived quantities (ratios, return levels) and probability statements
  (P(slope > 0.15 C/decade)) trivial.
- Check the likelihood: hot-day counts have variance 4x the mean (heat waves cluster), so a
  Poisson understates uncertainty; use negative binomial.
- MCMC (Metropolis-Hastings): propose, accept with min(1, p_new/p_old) in log space; tune step for
  25-45% acceptance; >= 4 dispersed chains; burn-in; R-hat < 1.01; ESS in the thousands;
  trace plots; posterior predictive checks. Use NUTS (PyMC / Stan / NumPyro) in practice.
- Bayesian advantage: small n, hierarchical / partial pooling across stations, propagating
  uncertainty through nonlinear transforms, encoding physical constraints as priors.

## 11. Machine learning

- Split by time; TimeSeriesSplit with a gap. Shuffled K-fold leaked 0.18 of R2 here.
- Baselines: persistence, climatology. Rivers are integrators: persistence is hard to beat one
  day ahead; the interesting question is skill without it.
- Features are domain knowledge: rolling precipitation sums at several scales, an exponentially
  weighted antecedent precipitation index, sin/cos of day-of-year, temperature for
  evapotranspiration, snow for melt.
- Regularisation: Ridge (L2) shrinks correlated coefficients together; Lasso (L1) selects.
  Standardise first; tune alpha with time-aware CV.
- Trees / boosting: nonlinear, interaction-aware, scale-free, but cannot extrapolate beyond the
  training range (record floods, warming climate). When linear ties boosting, the ceiling is the
  information in the features, not the model.
- Interpretation: permutation importance on held-out data (impurity importance is biased);
  partial dependence for shape; correlated features share importance, so read groups.
- Hydrology metrics: NSE = 1 - SSE/SST on flow (rewards peaks); KGE combines correlation,
  variability ratio and bias; log-space RMSE weights low flows.
- Rare-event classification: accuracy is meaningless (95.6% for 'never'); report PR-AUC vs the
  base rate, calibration (reliability diagram, Brier), and pick the threshold from the costs of
  misses vs false alarms.

## Earth-science gotchas, in one list

1. Autocorrelation in time: row count is not n. Aggregate, adjust SEs (HAC, n_eff), block bootstrap, Hamed-Rao.
2. Autocorrelation in space: same problem; test residuals with Moran's I; krige the residuals.
3. Seasonality: work on anomalies; otherwise correlations and trends are seasonal artefacts.
4. Heavy tails and zeros: precipitation and discharge want logs, gammas, or extreme-value models.
5. Nonstationarity: differencing for ARIMA; time-varying parameters for extremes; trend models for kriging.
6. Unit of analysis and pseudo-replication: days within a season are not independent samples.
7. Multiple testing: months, stations, variables and quantiles multiply your tests.
8. Endpoint and window choice: short windows and cherry-picked starts manufacture or hide trends.
9. Extrapolation: 500-year return levels from 95 years, tree models beyond the training range, ln(CO2) outside the record.
10. Leakage: time-aware splits with gaps; strictly antecedent features (shift by 1).
11. Collinear forcings: observational regression describes; attribution needs physics or experiments.
12. Baselines and uncertainty: a forecast without a naive comparison and an interval is not a result.

## Formulas worth having in your head

- Effective sample size (AR1): n_eff = n (1-r1)/(1+r1); SE_adj = SE sqrt(n/n_eff)
- Mann-Kendall: S = sum_{i<j} sign(x_j - x_i); Var(S) = [n(n-1)(2n+5) - sum_t t(t-1)(2t+5)]/18
- Cohen's d = (mean_a - mean_b) / pooled SD
- Family-wise error for m independent tests: 1 - (1-alpha)^m
- GEV return level: z_T = mu + (sigma/xi)[(-ln(1-1/T))^(-xi) - 1]; Gumbel: mu - sigma ln(-ln(1-1/T))
- GPD return level: z_T = u + (sigma/xi)[(lambda T)^xi - 1]
- Moran's I = (n/S0) (z'Wz)/(z'z); E[I] = -1/(n-1)
- Semivariance: gamma(h) = 1/(2N(h)) sum (z_i - z_j)^2
- Jackknife bias: (n-1)(mean(theta_(-i)) - theta_hat)
- NSE = 1 - sum(obs-sim)^2 / sum(obs-mean)^2; KGE = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)
- Beta-Binomial update: Beta(a+k, b+n-k); Gamma-Poisson update: Gamma(a + sum x, rate b + n)
- Metropolis acceptance: min(1, p(theta')/p(theta)); R-hat = sqrt(((n-1)/n W + B/n) / W)
