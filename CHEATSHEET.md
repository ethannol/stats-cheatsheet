# Statistics for mining data: the cheat sheet

Condensed from the nine scripts in `analyses/`. Section numbers match the scripts.

## 0. Which method?

| Question | Method | Where |
|---|---|---|
| What do these assays look like? | log histogram, lognormal probability plot, CV, per-domain stats | [01](analyses/01_assay_eda.py) |
| Is the drilling biased toward high grade? | cell declustering; declustered-mean vs cell-size curve | [01](analyses/01_assay_eda.py) |
| Should I cap outliers? | probability-plot break, P97.5-P99 table of metal removed, per domain | [01](analyses/01_assay_eda.py) |
| Samples have different lengths | composite to fixed length (length-weighted); check metal is conserved | [01](analyses/01_assay_eda.py) |
| Elements are percentages / ppm | log-ratio transforms (clr for biplots, ilr for models), variation matrix | [02](analyses/02_compositional.py) |
| Is the lab precise? | duplicate scatter, HARD (90% < 20%), Thompson-Howarth sd = s0 + k c | [03](analyses/03_qaqc.py) |
| Is the lab biased? | CRM control chart +-2SD/3SD, CUSUM for small sustained bias | [03](analyses/03_qaqc.py) |
| Contamination? | blanks vs DL and vs the preceding sample's grade | [03](analyses/03_qaqc.py) |
| Values below detection | Kaplan-Meier (default), ROS (lognormal), never substitution | [03](analyses/03_qaqc.py), [09](analyses/09_environmental.py) |
| How continuous is the grade? | omnidirectional then directional variograms; nugget/sill/range; anisotropy ellipse | [04](analyses/04_geostatistics.py) |
| Estimate blocks | ordinary kriging with a moving neighbourhood; validate by leave-one-out | [04](analyses/04_geostatistics.py) |
| Which estimator? | compare bias, RMSE, slope of truth-on-estimate, smoothing (sd of estimates) | [04](analyses/04_geostatistics.py) |
| Tonnes and grade above cut-off | grade-tonnage curves at the SMU support; beware point-support curves | [04](analyses/04_geostatistics.py) |
| How uncertain is the resource? | sequential Gaussian simulation; re-block; spread across realisations | [04](analyses/04_geostatistics.py) |
| How big a sample? | Gy's FSE: c f g l d^3 / M; nomogram of the preparation protocol | [05](analyses/05_sampling_theory.py) |
| How much of the nugget is measurement? | duplicate variance vs variogram nugget | [05](analyses/05_sampling_theory.py) |
| Equipment lifetimes with units still running | Kaplan-Meier; Weibull MLE with censoring; beta tells the failure mode | [06](analyses/06_reliability.py) |
| When to replace preventively? | age-replacement cost rate from Weibull; only if beta > 1 | [06](analyses/06_reliability.py) |
| Has production / a KPI shifted? | I-MR for big shifts; CUSUM or EWMA for small sustained shifts; deseasonalise first | [06](analyses/06_reliability.py) |
| Compare safety across mines/years | rates per 200,000 h; Poisson/NB GLM with log(hours) offset | [07](analyses/07_safety.py) |
| Rank mines by rate | empirical-Bayes shrinkage (Gamma-Poisson) before ranking | [07](analyses/07_safety.py) |
| Price risk | log returns, annualised vol, GBM vs mean-reversion fans | [08](analyses/08_economics.py) |
| Project value under uncertainty | Monte Carlo NPV with correlated inputs; tornado; P(NPV < 0) | [08](analyses/08_economics.py) |
| Is the river getting cleaner? | C-Q regression, flow-adjusted Kendall/Sen trend, exceedance fractions with CIs | [09](analyses/09_environmental.py) |

## 1. Assay EDA

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

## 2. Compositional data

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

## 3. QA/QC and censored data

- Insert blind at ~5% each: field duplicates (precision), pulp duplicates (assay precision),
  CRMs (accuracy), coarse blanks (contamination).
- HARD = |a-b| / (a+b); rule of thumb 90% of field duplicates < 20% (gold), pulp < 10%. Precision is
  concentration-dependent: Thompson-Howarth sd = s0 + k c from median absolute differences x 1.048.
  Quote 'X% at grade Y'.
- Sample swaps between similar grades are invisible to duplicate statistics; chain of custody beats
  statistics for gross errors.
- CRM charts: certified +-2SD warning, +-3SD action. A sustained bias below 3SD hides from the limits;
  CUSUM or 'run of n points on one side' catches it. Remove mix-ups (single wild points) from the
  bias estimate; investigate them. Bias proportional across CRMs = calibration; constant offset =
  blank/contamination.
- Blanks: plot vs DL and vs the previous sample's grade; failing blanks after high grades =
  carry-over in preparation (fix the prep, re-assay the followers), not a lab-wide problem.
- Non-detects are left-censored. Kaplan-Meier (flip, right-censored KM, flip back) is
  distribution-free and the regulatory default; ROS imputes below the DL from a lognormal fit to the
  detects' probability plot; substitution (0, DL/2, DL) biases both mean and variance and breaks
  when detection limits changed (Clark Fork lead: means 0.23 / 1.01 / 1.79 vs KM 0.28). Above ~50%
  censoring report upper percentiles only.

## 4. Geostatistics

- Variogram gamma(h) = 1/(2N) sum (z_i - z_j)^2. Nugget = sampling error + micro-scale; sill ~
  variance; range = correlation length. Lag ~ sample spacing; fit the short lags (they drive the
  weights) and stop at half the domain. Directional variograms (azimuth cone +- 22.5 deg, bandwidth)
  give anisotropy: ratio of minor to major range; model it by stretching distances.
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

## 5. Sampling theory

- Gy: relative FSE variance = c f g l d^3 / M. c ~ rho_mineral / a for low critical content a:
  gold at 1.5 g/t c ~ 1.3e7, copper at 1% c ~ 140, iron ore ~ 0.5. f = 0.5, g = 0.25 (crushed),
  l liberation (<= 1). Classic l = sqrt(d_l / d) (alpha 2.5); Francois-Bongarcon alpha ~1.5 for gold.
- Protocol = chain of crush and split stages; split variances add. The nomogram makes the worst
  stage obvious: splitting while coarse. Crush fine before the first big mass reduction.
- Nugget decomposition: duplicate-derived measurement variance / variogram nugget = share that a
  better protocol could remove; the rest is geology at sub-sample-spacing scale.
- n for a domain mean to +-E: (z CV / E)^2; CV 2 at 10% needs ~1,500 independent samples. Kish
  effective n = (sum w)^2 / sum w^2 under declustering weights; spatial correlation lowers it more.

## 6. Reliability and SPC

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
  same false-alarm rate. Set limits from phase I; deseasonalise first.

## 7. Safety analytics

- Rate = events x 200,000 / hours (100 worker-years). Fatalities: exact Poisson CI from chi-square
  quantiles because counts are small.
- Overdispersion: Poisson deviance/df >> 1 even nationally; across mines NB alpha ~0.5. Rate ratios
  agree between Poisson and NB; SEs do not (NB ~2x). Use NB or quasi-Poisson for inference.
- GLM with log(hours) offset: exp(beta) = rate ratio per unit of covariate (per year, coal vs metal).
- Small numbers: one injury in 2,000 hours = rate 90. Gamma-Poisson empirical Bayes by maximum
  likelihood of the NB marginal; shrunk rate = (k + a) / (E + b). Zero-injury tiny mines are
  uninformative, not safe.
- u-chart limits u +- 3 sqrt(u / n_i) widen with less exposure; at very large exposure Poisson limits
  are too tight (overdispersion) and I-chart limits from observed variability are more honest.
- Lost days are heavy-tailed (median 24, mean 57, top 10% of cases = 45% of days). Counts by
  experience band are not rates without hours by band.

## 8. Economics

- Log returns; annualise monthly vol by sqrt(12); fat tails (copper excess kurtosis ~5); base metals
  co-move, gold does not.
- GBM: S exp((mu - sigma^2/2) dt + sigma sqrt(dt) Z). OU/mean reversion from AR(1) on log price:
  kappa = -ln(b)/dt, half-life ln2/kappa. ADF rarely rejects a unit root for commodities; the
  reversion assumption must be argued from the cost curve. Fans differ hugely at 10 years.
- Monte Carlo NPV: distributions per input, correlations via Cholesky, asymmetric capex overruns;
  report mean vs deterministic, P10/P90, P(NPV < 0). Tornado for ranking drivers only.
- Cut-off grade g_c = opex per tonne ore / (price x recovery x payability); it moves with price and
  resizes the ore inventory via the grade-tonnage curve.
- Option to defer: E[max(NPV_t, 0)] discounted vs commit-now NPV; valuable when marginal and volatile.

## 9. Environmental monitoring

- C-Q regression on log-log axes: slope < 0 dilution (point source), > 0 mobilisation (sediment /
  tailings; Clark Fork +0.9). Seasonality follows the hydrograph.
- Flow-adjust (residuals from C-Q) before trend tests; compare raw vs adjusted Kendall tau and Sen
  slope (as %/yr from the log slope). Era tables communicate better than p-values.
- Compliance: exceedance fraction with Wilson CI per period; hardness-dependent criteria need
  paired hardness. A chronically exceeded benchmark is a sediment-source problem, not an event one.
- Change points: CUSUM of annual flow-adjusted medians; match the turning point to the
  intervention record (hypothesis, not proof).
- Censored low-level analytes: Kaplan-Meier per era; naive means with DL-as-value run 2-3x high.

## Mining data gotchas, in one list

1. Preferential sampling: decluster before any global statistic or histogram.
2. Support: assays, composites, blocks and SMUs have different variances; grade-tonnage depends on it.
3. Skewness: medians, logs, CV; capping is a documented judgement.
4. Closure: ratios, not percentages, for multi-element work.
5. Censoring: non-detects (left) and running equipment (right) need censored estimators.
6. QA/QC before geostatistics: precision and bias propagate into every block.
7. Smoothing and conditional bias at the cut-off: validate with slope and CV, not R^2.
8. Overdispersion everywhere counts appear: injuries, failures, exceedances.
9. Small denominators: shrink before ranking.
10. Flow (or any covariate) confounds environmental trends; adjust first.
11. Price process assumptions dominate long-horizon economics; show both GBM and reversion.
12. Synthetic inputs must be labelled and their recovered estimates compared with the truth.

## Formulas worth having in your head

- Declustering weight: w_i = (1 / n_cell(i)) x n / sum_j (1 / n_cell(j))
- Semivariance: gamma(h) = 1/(2N(h)) sum (z_i - z_j)^2; spherical: nugget + (sill-nugget)(1.5r - 0.5r^3), r = h/range
- Kriging system: [C 1; 1' 0][w; mu] = [c0; 1]; variance = sill - w'c0 - mu
- Krige's relation: D^2(point, domain) = D^2(point, block) + D^2(block, domain)
- clr(x) = ln(x) - mean(ln x); ilr = clr V' with V orthonormal (D-1 x D); Aitchison distance = ||clr(x) - clr(y)||
- HARD = |a-b| / (a+b); Thompson-Howarth sd = 1.048 x median|a-b| per concentration bin
- Kaplan-Meier: S(t) = prod_{t_i <= t} (1 - d_i / n_i); left censoring: flip x -> M - x
- Weibull: R(t) = exp(-(t/eta)^beta); h(t) = (beta/eta)(t/eta)^(beta-1); MTTF = eta Gamma(1 + 1/beta)
- Gy: s^2_FSE = c f g l d^3 / M; c = (1-a)/a [(1-a) rho_m + a rho_g]; l = (d_l/d)^(3-alpha)
- Shewhart I-chart sigma = mean(MR) / 1.128; CUSUM C+ = max(0, C+ + z - k), signal > h
- Incidence rate = injuries x 200,000 / hours; NB variance = mu + alpha mu^2; EB rate = (k + a)/(E + b)
- GBM: S_{t+dt} = S_t exp((mu - sigma^2/2)dt + sigma sqrt(dt) Z); OU half-life = ln 2 / kappa, kappa = -ln(b)/dt
- Cut-off grade: g_c = opex / (price x recovery x payability)
- Wilson interval for a proportion; C-Q: log C = a + b log Q, b < 0 dilution, b > 0 mobilisation
