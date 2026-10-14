# Statistics for environmental and mining data: the cheat sheet

Condensed from the twenty-one scripts in `analyses/`. Section numbers match the scripts:
01-11 are general methods worked on climate, river and station data, 12-21 the methods specific
to mining practice. For the question-first decision chart see [FLOWCHART.md](FLOWCHART.md).

## 0. Which method?

Each method appears once; the last column lists every section that demonstrates it.

| Question | Method | Where |
|---|---|---|
| What does this variable look like? | histogram, QQ plot, median/IQR/MAD, skew; climatology + anomaly; for grades: log histogram, lognormal probability plot, CV, per-domain stats | [01](analyses/01_eda.py), [12](analyses/12_assay_eda.py) |
| Is the sampling biased toward high values? | cell declustering; declustered-mean vs cell-size curve | [12](analyses/12_assay_eda.py) |
| Should I cap outliers? | probability-plot break, P97.5-P99 table of metal removed, per domain | [12](analyses/12_assay_eda.py) |
| Samples have different lengths | composite to fixed length (length-weighted); check metal is conserved | [12](analyses/12_assay_eda.py) |
| Elements are percentages / ppm | log-ratio transforms (clr for biplots, ilr for models), variation matrix | [13](analyses/13_compositional.py) |
| Are two groups different? | Welch t (means), Mann-Whitney (ranks), permutation (anything); report effect size + CI; with non-detects: Gehan or Peto-Prentice, never substitute then t-test | [02](analyses/02_hypothesis_tests.py), [21](analyses/21_censored_inference.py) |
| More than two groups | Kruskal-Wallis then Dunn with Holm/BH; Welch ANOVA + Tukey on logs | [21](analyses/21_censored_inference.py) |
| Did a proportion change? | two-proportion z / chi-square; Fisher if small; exceedance fraction per period with a Wilson CI | [02](analyses/02_hypothesis_tests.py), [20](analyses/20_environmental.py) |
| Many tests at once | Bonferroni (FWER) or Benjamini-Hochberg (FDR) | [02](analyses/02_hypothesis_tests.py), [21](analyses/21_censored_inference.py) |
| How many samples do I need? | power analysis (effect size / SD, alpha, power); for physical samples Gy's FSE c f g l d^3 / M, and (z CV / E)^2 for a domain mean | [02](analyses/02_hypothesis_tests.py), [16](analyses/16_sampling_theory.py) |
| How does y depend on x? | OLS + residual diagnostics; HAC/GLS if autocorrelated; GLM if y is counts/positive; MLE (Tobit) or Akritas-Theil-Sen if y is censored | [03](analyses/03_regression.py), [21](analyses/21_censored_inference.py) |
| Which of several predictors matter? | multiple regression, standardised coefficients, VIF, AIC/BIC | [03](analyses/03_regression.py) |
| Concentration depends on discharge | C-Q regression on log-log axes (sign tells dilution from mobilisation); flow-adjust before any trend test | [20](analyses/20_environmental.py) |
| Trend + seasonality in a series | STL decomposition; ADF + KPSS; difference; SARIMA with a holdout | [04](analyses/04_time_series.py) |
| Hidden periodicities | periodogram / Welch PSD after detrending; test against red noise | [04](analyses/04_time_series.py) |
| Does x lead y? | cross-correlation of *detrended* (or pre-whitened) series | [04](analyses/04_time_series.py) |
| Is there a monotonic trend? | Mann-Kendall + Sen's slope with Hamed-Rao correction, OLS with HAC SE (annual series); flow-adjusted Kendall/Sen (river chemistry); Seasonal Kendall + Sen slope (monthly or seasonal sampling) | [05](analyses/05_trends_changepoints.py), [20](analyses/20_environmental.py), [21](analyses/21_censored_inference.py) |
| When did it change? | Pettitt (shift), hinge regression (slope change), CUSUM of annual (flow-adjusted) medians; PELT for many | [05](analyses/05_trends_changepoints.py), [20](analyses/20_environmental.py) |
| Has a process / KPI / lab shifted? | I-MR for big shifts; CUSUM or EWMA for small sustained shifts (CRM bias, production); deseasonalise first | [14](analyses/14_qaqc.py), [17](analyses/17_reliability.py) |
| How big is the 100-year event? | GEV on annual maxima or GPD on peaks over threshold; bootstrap the CI | [06](analyses/06_extremes.py) |
| Are extremes changing? | nonstationary GEV (location ~ time), likelihood-ratio test | [06](analyses/06_extremes.py) |
| Uncertainty with no formula | bootstrap (BCa); block bootstrap for dependent data; jackknife for bias | [07](analyses/07_resampling.py) |
| Is my test calibrated? | Monte Carlo from a fitted null (e.g. AR(1)) | [05](analyses/05_trends_changepoints.py), [07](analyses/07_resampling.py) |
| Reduce many variables | standardise, PCA (loadings, scree), then cluster (k-means / Ward) | [08](analyses/08_multivariate.py) |
| Do nearby places resemble each other? | Moran's I on the variable and on model residuals | [09](analyses/09_spatial.py) |
| How continuous is the field / grade? | variogram: nugget, sill, range; directional variograms and the anisotropy ellipse | [09](analyses/09_spatial.py), [15](analyses/15_geostatistics.py) |
| Fill in a map between samples | regression kriging with covariates (stations); ordinary kriging with a moving neighbourhood vs IDW vs NN (deposits); validate by leave-one-out | [09](analyses/09_spatial.py), [15](analyses/15_geostatistics.py) |
| Tonnes and grade above cut-off | grade-tonnage curves at the SMU support; beware point-support curves | [15](analyses/15_geostatistics.py) |
| How uncertain is the resource / the map? | sequential Gaussian simulation; re-block; spread across realisations | [15](analyses/15_geostatistics.py) |
| How much of the nugget is measurement? | duplicate variance vs variogram nugget | [16](analyses/16_sampling_theory.py) |
| Probability statements about parameters | Bayesian: conjugate update or MCMC; credible intervals, posterior predictive | [10](analyses/10_bayesian.py) |
| Rank units by a rate | empirical-Bayes (Gamma-Poisson) shrinkage before ranking | [18](analyses/18_safety.py) |
| Predict, and be honest about skill | time-aware CV, baselines, permutation importance, PR-AUC for rare events | [11](analyses/11_ml.py) |
| Is the lab precise? | duplicate scatter, HARD (90% < 20%), Thompson-Howarth sd = s0 + k c | [14](analyses/14_qaqc.py) |
| Is the lab biased? | CRM control chart +-2SD/3SD, CUSUM for small sustained bias | [14](analyses/14_qaqc.py) |
| Contamination? | blanks vs DL and vs the preceding sample's grade | [14](analyses/14_qaqc.py) |
| Values below detection | Kaplan-Meier (default), ROS (lognormal), never substitution; KM per era for low-level analytes; KM mean + SE for limits | [14](analyses/14_qaqc.py), [20](analyses/20_environmental.py), [21](analyses/21_censored_inference.py) |
| Upper confidence limit on a mean (exposure point) | ProUCL logic: t if normal, gamma if gamma fits, Chebyshev by skew, KM with non-detects | [21](analyses/21_censored_inference.py) |
| Is a new sample consistent with background? | prediction limit for k future samples; tolerance limit for a population fraction | [21](analyses/21_censored_inference.py) |
| Equipment lifetimes with units still running | Kaplan-Meier; Weibull MLE with censoring; beta tells the failure mode | [17](analyses/17_reliability.py) |
| When to replace preventively? | age-replacement cost rate from Weibull; only if beta > 1 | [17](analyses/17_reliability.py) |
| Compare rates across units / years | rates per 200,000 h; Poisson/NB GLM with log(hours) offset; u-chart | [18](analyses/18_safety.py) |
| Price risk | log returns, annualised vol, GBM vs mean-reversion fans | [19](analyses/19_economics.py) |
| Project value under uncertainty | Monte Carlo NPV with correlated inputs; tornado; P(NPV < 0) | [19](analyses/19_economics.py) |
| Is the river getting cleaner? | C-Q regression, flow-adjusted Kendall/Sen trend, exceedance fractions with CIs, CUSUM change point | [20](analyses/20_environmental.py) |

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
- Skewed positive variables (grades, metals): section 12.

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
- Groups that contain non-detects, and more than two groups: section 21.

## 3. Regression

- OLS assumptions (LINE): Linear mean, Independent errors, Normal errors (small n only), Equal
  variance. X need not be normal. Check: residuals vs fitted, QQ, residuals vs time/space,
  Durbin-Watson (2 = none), Ljung-Box, Breusch-Pagan.
- Autocorrelated residuals: coefficients unbiased, SEs too small by ~sqrt(n/n_eff) with
  n_eff = n(1-r1)/(1+r1). Fixes: Newey-West HAC SEs, GLS-AR(1) (Cochrane-Orcutt), or aggregate.
- Heteroscedasticity: HC (White) robust SEs, transform y, or a GLM with the right variance
  function (Gamma for positive skewed, Poisson/negative binomial for counts; rates with an
  exposure offset in section 18).
- Multicollinearity: VIF = 1/(1-R2_j); > 5-10 is trouble. Predictions fine, coefficients unstable
  (ln CO2 coefficient doubled when 'year' was added; VIF = 91).
- Model comparison: adjusted R2, AIC (prediction-oriented), BIC (parsimony-oriented), on the same y.
  Count estimated breakpoints and tuned hyperparameters as parameters.
- Omitted-variable bias: the ln(CO2) slope absorbs methane, aerosols and everything else that
  trended. Observational regression describes; it does not attribute.
- Partial residual (component + residual) plots check each predictor's linearity.
- A censored response (non-detects): section 21.

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
- Seasonal sampling (monthly water quality): Seasonal Kendall, section 21. A covariate such as
  discharge: flow-adjust first, section 20.

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
- Parts of a whole (percentages, ppm) need log-ratios before any of this: section 13.

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
- Preferential sampling, directional variograms, the kriging system, change of support and
  simulation: section 15.

## 10. Bayesian

- posterior ∝ likelihood x prior. Weakly informative priors by default; flat is a choice too.
- Credible interval: P(theta in interval | data). Confidence interval: long-run coverage of the
  procedure. With much data and weak priors they coincide numerically.
- Conjugate pairs: Beta-Binomial, Gamma-Poisson, Normal-Normal, Normal-Inverse-Gamma.
  Posterior sampling makes derived quantities (ratios, return levels) and probability statements
  (P(slope > 0.15 C/decade)) trivial. Gamma-Poisson as empirical Bayes for rates: section 18.
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

## 12. Assay EDA

- Grades are positive and right-skewed; often lognormal within a domain. Quote the CV: < 1 benign,
  1-2 care, > 2 nuggety (a few samples carry the metal). Probability plots on a log axis show
  mixed populations as kinks and guide capping.
- Preferential sampling (drilling where it is rich) biases the naive mean high. Cell declustering:
  weight = 1 / (samples in the cell), averaged over grid offsets; pick the cell size at the minimum
  of the declustered mean (or the nominal spacing). It removes most but not all of the bias
  (Walker Lake: 435 -> 287 vs truth 278).
- Capping: tabulate samples cut, metal removed, CV before/after for candidate caps; decide per domain;
  report it. For low-CV variables it hardly matters; for gold it can be 10-20% of the metal.
- Compositing to fixed lengths length-weights the assays, conserves metal, and lowers variance
  (volume-variance). Estimate on composites, not raw assays.
- Domain first (lithology, alteration, oxidation): stationarity is a decision.

## 13. Compositional data

- Parts of a whole are constrained; independent amounts closed to 100% acquire negative
  correlations (~ -0.5 for three equal parts). Raw correlation matrices of assay percentages are
  partly artefacts; sub-compositions change the correlations; Euclidean distance is meaningless.
- alr: log(x_i / x_D), D-1 coordinates, depends on the reference. clr: log(x_i / g(x)), symmetric,
  rows sum to zero, singular covariance, use for biplots. ilr: orthonormal rotation of clr, D-1
  unconstrained coordinates for regression / clustering / kriging; back-transform results.
- Variation matrix var(ln(x_i/x_j)): near zero means proportional elements (same mineral host).
- Compositional mean = closed geometric mean; distance = Aitchison. Replace zeros / non-detects
  (e.g. 65% of DL, multiplicative) before logs.
- Trace elements (Jura: 0.03% of the soil) suffer little closure to the total, but the ratios still
  carry the geology; log-ratio thinking is about which elements travel together.

## 14. QA/QC and censored data

- Insert blind at ~5% each: field duplicates (precision), pulp duplicates (assay precision),
  CRMs (accuracy), coarse blanks (contamination).
- HARD = |a-b| / (a+b); rule of thumb 90% of field duplicates < 20% (gold), pulp < 10%. Precision is
  concentration-dependent: Thompson-Howarth sd = s0 + k c from median absolute differences x 1.048.
  Quote 'X% at grade Y'.
- Sample swaps between similar grades are invisible to duplicate statistics; chain of custody beats
  statistics for gross errors.
- CRM charts: certified +-2SD warning, +-3SD action. A sustained bias below 3SD hides from the limits;
  CUSUM (mechanics in section 17) or 'run of n points on one side' catches it. Remove mix-ups
  (single wild points) from the bias estimate; investigate them. Bias proportional across CRMs =
  calibration; constant offset = blank/contamination.
- Blanks: plot vs DL and vs the previous sample's grade; failing blanks after high grades =
  carry-over in preparation (fix the prep, re-assay the followers), not a lab-wide problem.
- Non-detects are left-censored. Kaplan-Meier (flip, right-censored KM, flip back) is
  distribution-free and the regulatory default; ROS imputes below the DL from a lognormal fit to the
  detects' probability plot; substitution (0, DL/2, DL) biases both mean and variance and breaks
  when detection limits changed (Clark Fork lead: means 0.23 / 1.01 / 1.79 vs KM 0.28). Above ~50%
  censoring report upper percentiles only. Tests, trends and limits with non-detects: section 21.

## 15. Geostatistics

- Variogram, nugget, sill, range and models as in section 9. Mining specifics: lag ~ sample
  spacing; fit the short lags (they drive the weights) and stop at half the domain. Directional
  variograms (azimuth cone +- 22.5 deg, bandwidth) give anisotropy: ratio of minor to major range;
  model it by stretching distances.
- Ordinary kriging: solve [C 1; 1' 0][w; mu] = [c0; 1] in a moving neighbourhood. Weights reflect
  distance and redundancy (clusters share weight): kriging declusters implicitly, IDW does not
  (Walker Lake IDW bias +57 vs OK -3). Kriging variance depends only on geometry.
- Cross-validation: mean error ~0, RMSE, correlation, slope of truth-on-estimate ~1 (conditional
  unbiasedness), standardised-error variance ~1 (otherwise the nugget/sill scale is off).
- Smoothing: estimates have less variance than reality; NN has more than it should; squared-error
  optimality costs variance.
- Volume-variance (Krige's relation): block variance falls with block size; grade-tonnage curves
  depend on support. Point-support curves overstate high-grade tonnes; smoothed estimates
  understate them at high cut-offs and overstate just above low cut-offs. State the SMU.
- Simulation (SGS): normal-score transform (with declustering weights), simple kriging along a
  random path, draw, condition on the draw. Reproduces histogram and variogram; realisations show
  uncertainty, E-type ~ kriging. It reproduces the *input* histogram, so a biased histogram gives a
  biased uncertainty band. Re-block realisations to SMU before grade-tonnage.

## 16. Sampling theory

- Gy: relative FSE variance = c f g l d^3 / M. c ~ rho_mineral / a for low critical content a:
  gold at 1.5 g/t c ~ 1.3e7, copper at 1% c ~ 140, iron ore ~ 0.5. f = 0.5, g = 0.25 (crushed),
  l liberation (<= 1). Classic l = sqrt(d_l / d) (alpha 2.5); Francois-Bongarcon alpha ~1.5 for gold.
- Protocol = chain of crush and split stages; split variances add. The nomogram makes the worst
  stage obvious: splitting while coarse. Crush fine before the first big mass reduction.
- Nugget decomposition: duplicate-derived measurement variance / variogram nugget = share that a
  better protocol could remove; the rest is geology at sub-sample-spacing scale.
- n for a domain mean to +-E: (z CV / E)^2; CV 2 at 10% needs ~1,500 independent samples. Kish
  effective n = (sum w)^2 / sum w^2 under declustering weights; spatial correlation lowers it more.
  Power analysis for a designed comparison: section 2.

## 17. Reliability and SPC

- Right censoring: units still running are 'at least this long'. Kaplan-Meier S(t) = prod (1 - d/n).
  Weibull MLE: likelihood = pdf for failures x survival for censored.
- beta < 1 infant mortality (replace on failure, fix quality), ~1 random (exponential, PM useless),
  > 1 wear-out (PM pays). eta = 63.2% life. MTBF = eta Gamma(1 + 1/beta), not the mean of observed
  TBFs. Probability plot: ln(-ln S) vs ln t, slope beta.
- Availability = MTBF / (MTBF + MTTR); components in series: unavailabilities roughly add.
- Age replacement interval minimises [c_pm R(T) + c_f (1-R(T))] / expected cycle length.
- Homogeneity across units: chi-square on counts, var/mean of counts, before naming a bad actor.
- SPC: I-MR (sigma = mean moving range / 1.128) for shifts > 2-3 sigma; CUSUM (k = half the target
  shift, h ~ 4-5) and EWMA (lambda 0.1-0.3) for 0.5-1.5 sigma sustained shifts, 3-5x faster at the
  same false-alarm rate. Set limits from phase I; deseasonalise first. The same CUSUM catches
  laboratory bias (section 14) and regime shifts in a river (section 20).

## 18. Safety analytics

- Rate = events x 200,000 / hours (100 worker-years). Fatalities: exact Poisson CI from chi-square
  quantiles because counts are small.
- Overdispersion: Poisson deviance/df >> 1 even nationally; across mines NB alpha ~0.5. Rate ratios
  agree between Poisson and NB; SEs do not (NB ~2x). Use NB or quasi-Poisson for inference.
- GLM with log(hours) offset (GLM families in section 3): exp(beta) = rate ratio per unit of
  covariate (per year, coal vs metal).
- Small numbers: one injury in 2,000 hours = rate 90. Gamma-Poisson (section 10) empirical Bayes
  by maximum likelihood of the NB marginal; shrunk rate = (k + a) / (E + b). Zero-injury tiny
  mines are uninformative, not safe.
- u-chart limits u +- 3 sqrt(u / n_i) widen with less exposure; at very large exposure Poisson limits
  are too tight (overdispersion) and I-chart limits from observed variability are more honest.
- Lost days are heavy-tailed (median 24, mean 57, top 10% of cases = 45% of days). Counts by
  experience band are not rates without hours by band.

## 19. Economics

- Log returns; annualise monthly vol by sqrt(12); fat tails (copper excess kurtosis ~5); base metals
  co-move, gold does not.
- GBM: S exp((mu - sigma^2/2) dt + sigma sqrt(dt) Z). OU/mean reversion from AR(1) on log price:
  kappa = -ln(b)/dt, half-life ln2/kappa. ADF (section 4) rarely rejects a unit root for
  commodities; the reversion assumption must be argued from the cost curve. Fans differ hugely at
  10 years.
- Monte Carlo NPV: distributions per input, correlations via Cholesky, asymmetric capex overruns;
  report mean vs deterministic, P10/P90, P(NPV < 0). Tornado for ranking drivers only.
- Cut-off grade g_c = opex per tonne ore / (price x recovery x payability); it moves with price and
  resizes the ore inventory via the grade-tonnage curve (section 15).
- Option to defer: E[max(NPV_t, 0)] discounted vs commit-now NPV; valuable when marginal and volatile.

## 20. Environmental monitoring

- C-Q regression on log-log axes: slope < 0 dilution (point source), > 0 mobilisation (sediment /
  tailings; Clark Fork +0.9). Seasonality follows the hydrograph.
- Flow-adjust (residuals from C-Q) before trend tests: Mann-Kendall and Sen (section 5) on
  annual medians, Seasonal Kendall (section 21) on monthly sampling; compare raw vs adjusted tau
  and Sen slope (as %/yr from the log slope). Era tables communicate better than p-values.
- Compliance: exceedance fraction with Wilson CI per period; hardness-dependent criteria need
  paired hardness. A chronically exceeded benchmark is a sediment-source problem, not an event one.
- Change points: CUSUM of annual flow-adjusted medians (section 5 for Pettitt and hinge
  alternatives); match the turning point to the intervention record (hypothesis, not proof).
- Censored low-level analytes: Kaplan-Meier per era (section 14); naive means with DL-as-value
  run 2-3x high.

## 21. Censored inference, seasonal trends and compliance limits

- Seasonal Kendall: Mann-Kendall (section 5) within each season (month), S and variances summed
  across seasons, continuity-corrected z. Multiple values in a season-year reduce to their median.
  Seasonality never enters the comparisons, so power is higher than on annual medians and the
  variance is not inflated as in a whole-series test. Flow-adjust first (LOWESS of log C on
  log Q, section 20); a C-Q slope and drifting sampled flows otherwise leak hydrology into the
  trend. WRTDS (Hirsch 2010) is the estimation successor: time, discharge and season in one
  weighted regression, reporting flow-normalised concentration.
- Two groups with non-detects: Gehan's generalised Wilcoxon scores each observation +1 against
  every value it is certainly above and -1 against every value it is certainly below; a
  '< 1' against a detected 0.3 scores 0. W = sum of scores in group 1; permutation variance
  n1 n2 / (n (n-1)) sum(score^2). Without censoring it is Mann-Whitney (section 2). Substituting
  the DL turned a p = 2e-3 decrease into p = 1e-15; substituting zero reversed the sign.
  Peto-Prentice weights pairs more evenly when censoring fractions differ between groups.
- Censored regression: each non-detect contributes log Phi((DL - Xb)/sigma) to the
  likelihood; detects contribute the normal density. SEs from the Hessian. Detects-only OLS
  keeps the early samples that were high enough to measure and exaggerates the decline;
  DL/2 substitution injects the laboratory's detection-limit history into the response.
  Akritas-Theil-Sen: the slope at which censored Kendall's S of the residuals is zero;
  permutation p-value. When MLE and ATS disagree, the residual distribution, a missing
  covariate (flow) or a step rather than a monotone trend is the reason.
- UCL95 on the mean (ProUCL): Student-t if normality is not rejected; adjusted gamma UCL
  2 n xbar k / chi2_alpha(2 n k) if a gamma fits; otherwise Chebyshev xbar + sqrt(1/alpha - 1) SE
  at 95 / 97.5 / 99 % as the SD of the logs passes 0.5 / 1 / 1.5; bootstrap-t (section 7) for
  large n. Land's lognormal H-UCL is unstable and no longer recommended. With non-detects: KM
  mean (section 14) and Greenwood-type SE, then the same t-or-Chebyshev choice. A UCL assumes
  independent samples; seasonal, autocorrelated monitoring data have a smaller effective n.
- Prediction limit for the next k samples: xbar + t_{1-alpha/k, n-1} s sqrt(1 + 1/n)
  (Bonferroni; the Unified Guidance uses exact kappa tables and retest rules).
  Nonparametric: the background maximum, confidence n / (n + k). Tolerance limit covering
  fraction p with confidence 1-a: xbar + k s with k from the noncentral t. These are the
  compliance-monitoring questions; they need independent background, which a river
  record is not.
- k groups: Kruskal-Wallis omnibus, then Dunn's z on pooled mean ranks with the tie term and
  Holm (FWER) or BH (FDR) (section 2). Tukey HSD on logs when groups are roughly lognormal.
  Eras are an observational grouping; method, detection limit and flow regime changed with them.

## Gotchas, in one list

1. Autocorrelation in time: row count is not n. Aggregate, adjust SEs (HAC, n_eff), block bootstrap, Hamed-Rao.
2. Autocorrelation in space: same problem; test residuals with Moran's I; krige the residuals.
3. Seasonality: work on anomalies or within seasons; otherwise correlations and trends are seasonal artefacts.
4. Heavy tails, zeros and skew: precipitation, discharge and grades want medians, logs, CV, gammas or extreme-value models; capping is a documented judgement.
5. Nonstationarity: differencing for ARIMA; time-varying parameters for extremes; trend models for kriging.
6. Unit of analysis and pseudo-replication: days within a season are not independent samples.
7. Multiple testing: months, stations, variables and quantiles multiply your tests.
8. Endpoint and window choice: short windows and cherry-picked starts manufacture or hide trends.
9. Extrapolation: 500-year return levels from 95 years, tree models beyond the training range, ln(CO2) outside the record.
10. Leakage: time-aware splits with gaps; strictly antecedent features (shift by 1).
11. Collinear forcings: observational regression describes; attribution needs physics or experiments.
12. Baselines and uncertainty: a forecast without a naive comparison and an interval is not a result.
13. Preferential sampling: decluster before any global statistic or histogram.
14. Support: assays, composites, blocks and SMUs have different variances; grade-tonnage depends on it.
15. Closure: ratios, not percentages, for multi-element work.
16. Censoring: non-detects (left) and running equipment (right) need censored estimators, for
    inference (Gehan, censored MLE, Akritas-Theil-Sen) as well as for means.
17. QA/QC before geostatistics: precision and bias propagate into every block.
18. Smoothing and conditional bias at the cut-off: validate with slope and CV, not R^2.
19. Overdispersion everywhere counts appear: injuries, failures, exceedances.
20. Small denominators: shrink before ranking.
21. Flow (or any covariate) confounds environmental trends; adjust first.
22. Price process assumptions dominate long-horizon economics; show both GBM and reversion.
23. Synthetic inputs must be labelled and their recovered estimates compared with the truth.

## Formulas worth having in your head

- Effective sample size (AR1): n_eff = n (1-r1)/(1+r1); SE_adj = SE sqrt(n/n_eff)
- Mann-Kendall: S = sum_{i<j} sign(x_j - x_i); Var(S) = [n(n-1)(2n+5) - sum_t t(t-1)(2t+5)]/18
- Seasonal Kendall: S = sum_seasons S_k; Var = sum_seasons Var(S_k); z = (S - sign S) / sqrt(Var)
- Cohen's d = (mean_a - mean_b) / pooled SD
- Family-wise error for m independent tests: 1 - (1-alpha)^m
- GEV return level: z_T = mu + (sigma/xi)[(-ln(1-1/T))^(-xi) - 1]; Gumbel: mu - sigma ln(-ln(1-1/T))
- GPD return level: z_T = u + (sigma/xi)[(lambda T)^xi - 1]
- Moran's I = (n/S0) (z'Wz)/(z'z); E[I] = -1/(n-1)
- Semivariance: gamma(h) = 1/(2N(h)) sum (z_i - z_j)^2; spherical: nugget + (sill-nugget)(1.5r - 0.5r^3), r = h/range
- Kriging system: [C 1; 1' 0][w; mu] = [c0; 1]; variance = sill - w'c0 - mu
- Krige's relation: D^2(point, domain) = D^2(point, block) + D^2(block, domain)
- Declustering weight: w_i = (1 / n_cell(i)) x n / sum_j (1 / n_cell(j))
- clr(x) = ln(x) - mean(ln x); ilr = clr V' with V orthonormal (D-1 x D); Aitchison distance = ||clr(x) - clr(y)||
- Jackknife bias: (n-1)(mean(theta_(-i)) - theta_hat)
- NSE = 1 - sum(obs-sim)^2 / sum(obs-mean)^2; KGE = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)
- Beta-Binomial update: Beta(a+k, b+n-k); Gamma-Poisson update: Gamma(a + sum x, rate b + n)
- Metropolis acceptance: min(1, p(theta')/p(theta)); R-hat = sqrt(((n-1)/n W + B/n) / W)
- HARD = |a-b| / (a+b); Thompson-Howarth sd = 1.048 x median|a-b| per concentration bin
- Kaplan-Meier: S(t) = prod_{t_i <= t} (1 - d_i / n_i); left censoring: flip x -> M - x
- Weibull: R(t) = exp(-(t/eta)^beta); h(t) = (beta/eta)(t/eta)^(beta-1); MTTF = eta Gamma(1 + 1/beta)
- Gy: s^2_FSE = c f g l d^3 / M; c = (1-a)/a [(1-a) rho_m + a rho_g]; l = (d_l/d)^(3-alpha)
- Shewhart I-chart sigma = mean(MR) / 1.128; CUSUM C+ = max(0, C+ + z - k), signal > h
- Incidence rate = injuries x 200,000 / hours; NB variance = mu + alpha mu^2; EB rate = (k + a)/(E + b)
- GBM: S_{t+dt} = S_t exp((mu - sigma^2/2)dt + sigma sqrt(dt) Z); OU half-life = ln 2 / kappa, kappa = -ln(b)/dt
- Cut-off grade: g_c = opex / (price x recovery x payability)
- Wilson interval for a proportion; C-Q: log C = a + b log Q, b < 0 dilution, b > 0 mobilisation
- Gehan: W = sum_{i in group 1} u_i, u_i = #(certainly above) - #(certainly below); Var = n1 n2 sum u^2 / (n (n-1))
- Censored likelihood: sum_det log phi((y - Xb)/s)/s + sum_nd log Phi((DL - Xb)/s)
- Chebyshev UCL: xbar + sqrt(1/alpha - 1) s/sqrt(n); gamma UCL: 2 n xbar k / chi2_alpha(2 n k)
- Prediction limit (k future): xbar + t_{1-alpha/k, n-1} s sqrt(1 + 1/n); nonparametric max has confidence n/(n+k)
