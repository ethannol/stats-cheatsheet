# Choosing a statistical method for earth-science and mining data

A decision chart from the question and the shape of the data to a method, and from the
method to a runnable worked example in [CHEATSHEET.md](CHEATSHEET.md). It starts from the
earth-science question rather than from the statistics, because that first layer (is the
response censored, seasonal, autocorrelated, spatial, compositional, preferentially sampled?)
is where generic "which test" charts stop and where environmental data usually break the
assumptions.

Codes in the boxes are section numbers: **01-11** are the general methods worked on
environmental data and **12-21** the methods specific to mining practice. The index at the end
links each code to its script in [`analyses/`](analyses).

The charts are original and drawn in [Mermaid](https://mermaid.js.org), so they render on
GitHub and can be edited as text. The decision logic condenses the sources listed at the
end; the USGS and EPA documents among them are works of the United States government and in
the public domain.

## 0. Before any method

```mermaid
flowchart TD
    A([New dataset]) --> B[Audit: coverage, gaps, units,<br>sentinel values, duplicates<br>01]
    B --> C{Are the rows<br>independent?}
    C -- "repeated in time" --> D[Autocorrelation: row count is not n.<br>Aggregate, adjust SEs, block bootstrap<br>01 03 07]
    C -- "spread over space" --> E[Spatial autocorrelation:<br>Moran's I on residuals<br>09]
    C -- "sampled where it was rich" --> F[Preferential sampling:<br>decluster first<br>12]
    C -- "yes" --> G
    D --> G
    E --> G
    F --> G
    G{Does the response<br>have structure to<br>remove first?}
    G -- "seasonal" --> H[Work on anomalies or<br>within seasons<br>01 21]
    G -- "depends on a covariate<br>such as discharge" --> I[Flow-adjust: residuals from<br>LOWESS of C on Q<br>20 21]
    G -- "parts of a whole<br>(% or ppm)" --> J[Log-ratio transform<br>13]
    G -- "values below a<br>detection limit" --> K[Censored methods, never<br>substitution: chart 5]
    G -- "no" --> L
    H --> L
    I --> L
    J --> L
    K --> L
    L([Now choose the method:<br>charts 1 to 6])
```

## 1. One variable, or a few groups

```mermaid
flowchart TD
    A{What is the question?} -- "describe it" --> B{Symmetric?}
    B -- yes --> B1[Mean, SD, histogram, QQ<br>01]
    B -- "skewed or heavy-tailed" --> B2[Median, IQR, MAD, log scale,<br>CV, probability plot<br>01 12]
    A -- "two groups differ?" --> C{Non-detects?}
    C -- yes --> C1[Gehan or Peto-Prentice<br>score test<br>21]
    C -- no --> C2{Paired?}
    C2 -- yes --> C3[Paired t or<br>Wilcoxon signed-rank<br>02]
    C2 -- no --> C4[Welch t on logs,<br>Mann-Whitney, or<br>permutation test<br>02 07]
    A -- "k groups differ?" --> D[Kruskal-Wallis then Dunn<br>with Holm or BH;<br>Welch ANOVA + Tukey on logs<br>21]
    A -- "a proportion changed?" --> E[Two-proportion z, chi-square,<br>Fisher if small; Wilson CI<br>02 20]
    A -- "many tests at once" --> F[Bonferroni for FWER,<br>Benjamini-Hochberg for FDR<br>02 21]
    A -- "how many samples?" --> G[Power analysis, or Gy's<br>sampling theory for<br>physical samples<br>02 16]
    A -- "upper limit on the mean<br>for a risk assessment" --> H[UCL95 by ProUCL logic:<br>chart 5]
```

## 2. Relationships and regression

```mermaid
flowchart TD
    A{What is the response?} -- "continuous" --> B{Non-detects<br>in y?}
    B -- yes --> B1[Censored MLE regression;<br>Akritas-Theil-Sen for a<br>distribution-free line<br>21]
    B -- no --> C[OLS + residual diagnostics<br>03]
    C --> D{Residuals<br>autocorrelated?}
    D -- yes --> D1[Newey-West HAC SEs,<br>GLS-AR1, or aggregate<br>03]
    D -- no --> E{Variance grows<br>with the mean?}
    E -- yes --> E1[Log y, robust HC SEs,<br>or Gamma GLM<br>03]
    E -- no --> F{Several<br>predictors?}
    F -- yes --> F1[Standardised coefficients, VIF,<br>AIC/BIC, partial residuals<br>03]
    A -- "counts or rates" --> G[Poisson GLM with offset;<br>negative binomial when<br>overdispersed<br>18 10]
    A -- "a composition" --> H[Regress on ilr coordinates,<br>back-transform<br>13]
    A -- "a concentration that<br>depends on discharge" --> I[C-Q regression on log-log<br>axes; sign tells dilution<br>from mobilisation<br>20]
    A -- "predict, with honest skill" --> J[Time-aware or spatial CV,<br>baselines, permutation<br>importance<br>11]
```

## 3. Time series, trends and extremes

```mermaid
flowchart TD
    A{What is the question?} -- "is there a<br>monotonic trend?" --> B{Seasonal<br>sampling?}
    B -- "monthly or seasonal" --> B1[Seasonal Kendall + Sen slope,<br>flow-adjusted if C depends on Q<br>21 20]
    B -- "annual series" --> B2[Mann-Kendall with Hamed-Rao<br>correction, Sen slope;<br>OLS with HAC SEs<br>05]
    B1 --> B3{Non-detects?}
    B3 -- yes --> B4[Akritas-Theil-Sen or<br>censored MLE with time<br>21]
    A -- "when did it change?" --> C[Pettitt for a shift, hinge<br>regression for a slope change,<br>CUSUM, PELT for many<br>05 20]
    A -- "forecast it" --> D[STL, ADF/KPSS, SARIMA or<br>Holt-Winters; time-ordered<br>holdout against naive baselines<br>04]
    A -- "hidden periodicity?" --> E[Periodogram after detrending;<br>test against red noise<br>04]
    A -- "does x lead y?" --> F[Cross-correlation of detrended<br>or pre-whitened series<br>04]
    A -- "how big is the<br>100-year event?" --> G{Block maxima or<br>peaks over threshold?}
    G -- "annual maxima" --> G1[GEV; bootstrap the<br>return level<br>06]
    G -- "all exceedances" --> G2[GPD after declustering<br>06]
    G1 --> G3{Climate or land use<br>changing?}
    G2 --> G3
    G3 -- yes --> G4[Nonstationary GEV,<br>likelihood-ratio test<br>06]
    A -- "has a process or KPI<br>shifted?" --> H[Shewhart I-MR for big shifts;<br>CUSUM or EWMA for small<br>sustained ones<br>17 14]
    A -- "uncertainty on any<br>of these?" --> I[Block bootstrap, Monte Carlo<br>from a fitted AR1 null<br>07]
```

## 4. Spatial data

```mermaid
flowchart TD
    A{What is the question?} -- "do nearby places<br>resemble each other?" --> B[Moran's I on the variable<br>and on model residuals<br>09]
    A -- "fill in a map<br>between samples" --> C{Were samples placed<br>preferentially?}
    C -- yes --> C1[Cell declustering before any<br>global statistic<br>12]
    C -- no --> D
    C1 --> D[Experimental variogram:<br>nugget, sill, range<br>09 15]
    D --> E{Direction<br>matters?}
    E -- yes --> E1[Directional variograms,<br>anisotropy ratio<br>15]
    E -- no --> F
    E1 --> F{Is there a trend or<br>a covariate?}
    F -- yes --> F1[Regression or universal<br>kriging on the residuals<br>09]
    F -- no --> F2[Ordinary kriging, moving<br>neighbourhood; validate by<br>leave-one-out<br>15]
    F1 --> G
    F2 --> G{Need block<br>estimates?}
    G -- yes --> G1[Change of support: volume-variance,<br>grade-tonnage at the SMU<br>15]
    G -- no --> H
    G1 --> H{Need the uncertainty<br>of the whole map?}
    H -- yes --> H1[Sequential Gaussian simulation;<br>re-block before grade-tonnage<br>15]
    A -- "classify or predict<br>across space" --> I[Spatial cross-validation;<br>random folds leak location<br>11]
```

## 5. Censored values and compliance limits

```mermaid
flowchart TD
    A{Values below a<br>detection limit?} -- "more than ~80 %" --> A1[Report the fraction detected<br>and upper percentiles only<br>14]
    A -- "some" --> B{What is needed?}
    B -- "a mean or percentiles" --> B1[Kaplan-Meier, distribution-free;<br>ROS if lognormal and n is small;<br>never substitute<br>14 20]
    B -- "an upper confidence<br>limit on the mean" --> B2[KM mean and SE, then<br>t or Chebyshev by skew<br>21]
    B -- "compare two groups" --> B3[Gehan or Peto-Prentice<br>21]
    B -- "a trend" --> B4[Akritas-Theil-Sen;<br>censored MLE with covariates<br>21]
    B -- "a map or a model" --> B5[Replace below-limit values<br>multiplicatively, then log-ratio<br>or kriging<br>13]
    A -- none --> C{What is needed?}
    C -- "exposure-point<br>concentration" --> C1{Distribution?}
    C1 -- normal --> C2[Student-t UCL95<br>21]
    C1 -- "gamma fits" --> C3[Adjusted gamma UCL95<br>21]
    C1 -- "skewed, no fit" --> C4[Chebyshev at 95, 97.5 or 99 %<br>by the SD of the logs;<br>bootstrap-t for large n<br>21]
    C -- "is a new sample<br>consistent with<br>background?" --> D[Prediction limit for k future<br>samples; tolerance limit for<br>a population fraction<br>21]
    C -- "how often is a<br>criterion exceeded?" --> E[Exceedance fraction per period<br>with a Wilson interval;<br>paired hardness if the criterion<br>depends on it<br>20]
    C -- "did the site get<br>cleaner after a<br>remediation?" --> F[Flow-adjusted seasonal Kendall<br>by era; CUSUM change point<br>against the intervention record<br>20 21]
```

## 6. Compositions, counts, lifetimes, money and samples

```mermaid
flowchart TD
    A{What kind of data?} -- "percentages or ppm<br>that sum to a whole" --> B[Closure makes raw correlations<br>spurious: alr/clr/ilr, variation<br>matrix, Aitchison distance<br>13]
    A -- "assays from a laboratory" --> C[QA/QC first: duplicates for<br>precision, CRMs and CUSUM for<br>bias, blanks for carry-over<br>14]
    C --> C1{Capping?}
    C1 --> C2[Probability-plot break; table<br>of metal removed per cap<br>12]
    A -- "injuries, failures,<br>exceedances per exposure" --> D[Rates with an offset; Poisson vs<br>negative binomial; empirical-Bayes<br>shrinkage before ranking; u-chart<br>18]
    A -- "lifetimes with units<br>still running" --> E[Right-censored Kaplan-Meier;<br>Weibull MLE; beta tells the<br>failure mode<br>17]
    A -- "prices and cash flows" --> F[Log returns, GBM vs mean<br>reversion, Monte Carlo NPV<br>with correlated inputs<br>19]
    A -- "how much rock to<br>sample or crush?" --> G[Gy's fundamental sampling<br>error; nomogram of the<br>preparation protocol<br>16]
    A -- "many correlated<br>variables" --> H[Standardise, PCA, then cluster;<br>validate against geography<br>08]
    A -- "a probability statement<br>about a parameter" --> I[Bayesian: conjugate update or<br>MCMC; hierarchical pooling<br>across stations<br>10]
```

## Script index

| Code | Script |
|---|---|
| 01 | [EDA and climatology](analyses/01_eda.py) |
| 02 | [Hypothesis tests, power, multiplicity](analyses/02_hypothesis_tests.py) |
| 03 | [Regression and diagnostics](analyses/03_regression.py) |
| 04 | [Time series and spectra](analyses/04_time_series.py) |
| 05 | [Trends and change points](analyses/05_trends_changepoints.py) |
| 06 | [Extremes](analyses/06_extremes.py) |
| 07 | [Resampling](analyses/07_resampling.py) |
| 08 | [Multivariate](analyses/08_multivariate.py) |
| 09 | [Spatial statistics](analyses/09_spatial.py) |
| 10 | [Bayesian](analyses/10_bayesian.py) |
| 11 | [Machine learning](analyses/11_ml.py) |
| 12 | [Assay EDA, declustering, capping](analyses/12_assay_eda.py) |
| 13 | [Compositional data](analyses/13_compositional.py) |
| 14 | [QA/QC and non-detects](analyses/14_qaqc.py) |
| 15 | [Geostatistics](analyses/15_geostatistics.py) |
| 16 | [Sampling theory](analyses/16_sampling_theory.py) |
| 17 | [Reliability and SPC](analyses/17_reliability.py) |
| 18 | [Safety rates](analyses/18_safety.py) |
| 19 | [Economics](analyses/19_economics.py) |
| 20 | [Environmental monitoring](analyses/20_environmental.py) |
| 21 | [Censored inference, seasonal trends, compliance limits](analyses/21_censored_inference.py) |

## Sources the logic condenses

- Helsel, D.R., Hirsch, R.M., Ryberg, K.R., Archfield, S.A. and Gilroy, E.J. (2020). *Statistical
  Methods in Water Resources*. USGS Techniques and Methods 4-A3. https://doi.org/10.3133/tm4A3.
  Public domain. The backbone of charts 1 to 3 and 5.
- Helsel, D.R. (2012). *Statistics for Censored Environmental Data Using Minitab and R*, 2nd ed.
  Wiley. Kaplan-Meier, ROS, Gehan and Peto-Prentice, Akritas-Theil-Sen, censored MLE.
- Hirsch, R.M., Slack, J.R. and Smith, R.A. (1982). Techniques of trend analysis for monthly
  water quality data. *Water Resources Research* 18(1), 107-121. The seasonal Kendall test.
- Hirsch, R.M., Moyer, D.L. and Archfield, S.A. (2010). Weighted regressions on time, discharge,
  and season (WRTDS). *JAWRA* 46(5), 857-880.
- U.S. EPA (2022). *ProUCL Version 5.2 Technical Guide*. Office of Research and Development.
  Public domain. The UCL decision tables in chart 5.
- U.S. EPA (2009). *Statistical Analysis of Groundwater Monitoring Data at RCRA Facilities:
  Unified Guidance*, EPA 530/R-09-007. Public domain. Prediction and tolerance limits.
- Isaaks, E.H. and Srivastava, R.M. (1989). *An Introduction to Applied Geostatistics*. Oxford.
  Rossi, M.E. and Deutsch, C.V. (2014). *Mineral Resource Estimation*. Springer. Chart 4.
- Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman and Hall.
  Pawlowsky-Glahn, V., Egozcue, J.J. and Tolosana-Delgado, R. (2015). *Modeling and Analysis
  of Compositional Data*. Wiley. Chart 6.
- Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.
- Gy, P. (1979). *Sampling of Particulate Materials*. Elsevier. Chart 6.
- Generic "which statistical test" decision trees published by university statistics services
  (for example the University of Manchester and Loughborough University charts) cover the
  two-group and k-group branches of chart 1. They are copyrighted and are not reproduced here.
