# Generated worked-example results

These are executed educational examples, not a new scientific study.
Data cutoff: `2026-09-30`. Source: `9122a16fd7f87003eaa0bf1a0e89733a75ce12a8`.
The run manifest records the exact executed source, dependencies and hashes.
Full transcripts retain diagnostic warnings, assumptions and interpretation.
The excerpts below omit diagnostic warnings and read numerical output directly
from those hash-checked transcripts. Synthetic examples are labelled in source.

## 01 eda

[Source](../analyses/01_eda.py) · [Full output](outputs/01_eda.txt)

```text
        count   mean    std   min  25%   50%   75%    max
TMAX  28032.0  17.10  10.23 -12.8  8.9  17.8  26.1   40.0
TMIN  28032.0   8.83   9.20 -18.9  1.7   8.9  16.7   28.9
PRCP  28032.0   3.31   9.35   0.0  0.0   0.0   1.3  192.3
SNOW  28031.0   1.79  16.11   0.0  0.0   0.0   0.0  693.0
 TMAX    0.0
TMIN    0.0
PRCP    0.0
SNOW    0.0
Longest run of missing TMAX days: 0
n            28032.00            28032.00             9385.00               9385.00
mean            17.10                3.31                9.88                  1.39
median          17.80                0.00                4.60                  1.53
sd              10.23                9.35               14.01                  1.47
iqr             17.20                1.30               11.70                  2.30
mad_scaled      12.45                0.00                5.93                  1.66
skew            -0.20                5.54                3.44                 -0.18
excess_kurt     -0.96               49.52               20.33                 -0.91
p05              0.60                0.00                0.30                 -1.20
p95             31.70               19.80               36.30                  3.59
max             40.00              192.30              192.30                  5.26
Fraction of dry days: 0.665
TMAX (n=all)           Shapiro W=0.9711 p=1.34e-30 | D'Agostino K2=5462.5 p=0.00e+00
TMAX (n=50 sample)     Shapiro W=0.9370 p=1.03e-02 | D'Agostino K2=5.2 p=7.39e-02
ln wet PRCP (n=all)    Shapiro W=0.9676 p=3.81e-32 | D'Agostino K2=1409.3 p=9.29e-307
Autocorrelation of daily TMAX, lag 1: 0.912  | of daily anomaly, lag 1: 0.583
    cycle. Anomalies still have 0.58, which is real weather persistence. Either way, daily
TMAX  1.000  0.953 -0.011
TMIN  0.953  1.000  0.046
PRCP -0.011  0.046  1.000
TMAX  1.000  0.953 -0.049
TMIN  0.953  1.000  0.036
PRCP -0.049  0.036  1.000
TMAX_anom      1.000      0.773 -0.058
TMIN_anom      0.773      1.000  0.071
```

## 02 hypothesis tests

[Source](../analyses/02_hypothesis_tests.py) · [Full output](outputs/02_hypothesis_tests.txt)

```text
Summer mean TMIN (nights): 1951-1980 vs 1991-2020 (n=30 summers each)
early: mean=18.91 sd=0.76 | late: mean=19.64 sd=0.79 | diff=+0.74 C
Welch t:        t=3.67  p=0.0005  95% CI diff=(0.33, 1.14)
Student t:      t=3.67  p=0.0005
Mann-Whitney U: U=677  p=0.0008
Permutation:    diff=+0.74  p=0.0007  (20000 shuffles)
Levene (equal var?): W=0.00 p=0.966
KS 2-sample:    D=0.400 p=0.0156
Cohen's d = 0.95
    are 0.7 C warmer with d = 0.95 (large). Levene finds no variance difference so Student
diff=-0.09 C  Welch p=0.737  95% CI=(-0.61, 0.43)  d=-0.09  Levene p=0.022
    with anything from a small cooling to a ~0.4 C warming. 'Not significant' does not
    mean 'no change'. Note also Levene p=0.02: summer TMAX variance increased even though
Daily values: n=2760 vs 2760, t=8.6, p=8.1e-18
Seasonal means: n=30 vs 30, t=3.7, p=5.3e-04
lag-1 autocorrelation of JJA daily TMIN: 0.75
    magnitude smaller because it pretends to have 2,760 independent observations per
Bootstrap 95% CI for difference in means: (0.34, 1.12)  vs Welch CI (0.33, 1.14)
Proportions: warm nights (TMIN >= 21.1 C / 70 F) by period
warm-night fraction: early=0.289 late=0.377
two-proportion z=6.94 p=3.89e-12 | chi2=47.8 dof=1 p=4.76e-12
odds ratio (late vs early) = 1.49
Warm nights per summer: early median=28 late median=36  Mann-Whitney p=0.0031  (n=30 vs 30 summers)
    chi2). Both inherit the daily-independence problem. The cleaner design is one count
Multiple comparisons: monthly TMIN change, 12 tests
     1 1.0432 0.0932        1.0000    0.1243    False     False   False
     2 1.3062 0.0294        0.3524    0.0503     True     False   False
     3 0.6345 0.1381        1.0000    0.1507    False     False   False
     4 0.6826 0.0185        0.2216    0.0424     True     False    True
     5 0.6196 0.0778        0.9342    0.1168    False     False   False
     6 0.6897 0.0033        0.0396    0.0307     True      True    True
     7 0.7891 0.0109        0.1303    0.0326     True     False    True
     8 0.7281 0.0212        0.2545    0.0424     True     False    True
     9 0.9182 0.0088        0.1061    0.0326     True     False    True
    10 0.5749 0.1188        1.0000    0.1426    False     False   False
```

## 03 regression

[Source](../analyses/03_regression.py) · [Full output](outputs/03_regression.txt)

```text
(a) Simple OLS: annual GISTEMP anomaly ~ ln(CO2 / 280)
                 coef    std err          t      P>|t|      [0.025      0.975]
Intercept     -0.5792      0.035    -16.346      0.000      -0.650      -0.508
lnco2          3.9134      0.134     29.228      0.000       3.646       4.181
R^2=0.928  AIC=-118.8  n=68
Implied warming per CO2 doubling (slope * ln 2): 2.71 C (OLS 95% CI 2.53 to 2.90)
    2.7 C here. The IPCC transient climate response is ~1.8 C (likely 1.4-2.2), so this
    single-predictor fit overshoots. Why: methane, N2O, halocarbons and declining aerosol
    cooling all rose alongside CO2, and the coefficient absorbs every one of them (omitted
Durbin-Watson=1.38 (2 = none; 1.38 => positive autocorrelation)
residual lag-1 r=0.28  effective n=38 of 68
Ljung-Box p (lag1, lag5): [0.0208, 0.0796]
Breusch-Pagan p=0.573 (heteroscedasticity)  Shapiro p=0.061
Standard errors under autocorrelation: OLS vs HAC vs GLS-AR(1)
OLS                       3.9134  0.1339  29.2284
OLS + n_eff adj           3.9134  0.1792  21.8353
HAC (Newey-West, 5 lags)  3.9134  0.1921  20.3699
GLSAR(1) Cochrane-Orcutt  3.9921  0.1776  22.4814
GLSAR estimated rho = 0.28
    sqrt(n / n_eff). HAC (Newey-West) and GLS with AR(1) errors both widen the interval.
Model comparison: ln CO2 vs linear time vs both
                 adj R2      AIC      BIC
T~lnco2           0.927 -118.824 -114.385
T~year            0.894  -93.009  -88.570
T~lnco2+year      0.940 -130.724 -124.066
T~lnco2+lnco2^2   0.934 -124.726 -118.068
In T~lnco2+year: coef(lnco2)=8.36 (was 3.91 alone), p(lnco2)=0.000, p(year)=0.000, corr(lnco2, year)=0.995, VIF=91
    so did warming. Adding 'year' still improves AIC, but with r=0.995 between the
    predictors the VIF is 91 and the CO2 coefficient swings from 3.9 to 8.4:
    excellent. Observational data alone cannot separate 'CO2' from 'anything else that
                 coef    std err          t      P>|t|      [0.025      0.975]
Intercept     38.3199      0.639     59.933      0.000      37.063      39.577
lat           -0.8229      0.013    -65.024      0.000      -0.848      -0.798
elev_km       -3.3821      0.146    -23.146      0.000      -3.670      -3.095
lon           -0.0865      0.005    -17.634      0.000      -0.096      -0.077
```

## 04 time series

[Source](../analyses/04_time_series.py) · [Full output](outputs/04_time_series.txt)

```text
Seasonal amplitude (peak-to-trough): 1960s 5.90 ppm, 2010s 6.81 ppm
Remainder SD: classical 0.299, STL 0.266
    since the 1960s (longer growing seasons, more vigorous Northern Hemisphere
Stationarity tests on CO2 and its differences
CO2 level                          ADF p=1.000  KPSS p=0.010  -> non-stationary
first difference                   ADF p=0.000  KPSS p=0.100  -> stationary
seasonal difference (lag 12)       ADF p=0.006  KPSS p=0.010  -> ambiguous
first + seasonal difference        ADF p=0.000  KPSS p=0.100  -> stationary
STL remainder                      ADF p=0.000  KPSS p=0.100  -> stationary
    carries the growth trend; both together give white-ish noise. That is the (d=1, D=1)
    little else: the signature of an MA(1) x seasonal MA(1), i.e. the classic 'airline
    model' SARIMA(0,1,1)(0,1,1)12. Let the data suggest the order, then confirm with AIC.
  (0,1,1)(0,1,1)12 airline     AIC=   386.4  BIC=   400.1
  (1,1,1)(1,1,1)12             AIC=   387.1  BIC=   410.0
  (2,1,0)(0,1,1)12             AIC=   394.8  BIC=   413.2
best by AIC: (0,1,1)(0,1,1)12 airline
                 coef    std err          z      P>|z|      [0.025      0.975]
ma.L1         -0.3804      0.033    -11.523      0.000      -0.445      -0.316
ma.S.L12      -0.8628      0.021    -40.734      0.000      -0.904      -0.821
sigma2         0.0964      0.005     20.384      0.000       0.087       0.106
Ljung-Box on residuals p(12)=0.625 p(24)=0.830
Holdout 2020-2025 (72 months):
SARIMA                  0.747  0.991  0.176  0.600
Holt-Winters            0.785  1.023  0.185  0.660
Seasonal naive + drift  1.554  1.731  0.369 -1.554
SARIMA 95% interval coverage on holdout: 1.00
    truth ran above the forecast: CO2 growth accelerated after 2020 (El Nino 2023-24),
    which a model fit through 2019 cannot know. Interval coverage is 1.00: the intervals
Rolling-origin (expanding window) backtest, 12-month horizon
12-month RMSE by origin year: 2005: 0.26, 2007: 0.44, 2009: 0.72, 2011: 0.32, 2013: 0.34, 2015: 0.81, 2017: 0.40, 2019: 0.31, 2021: 0.42, 2023: 0.82
mean RMSE=0.48 ppm
CO2 dominant frequency: 1.002 cycles/yr (period 1.00 yr); strongest above 1.5 cycles/yr: 2.005
Variance explained by the annual line (f in 0.9-1.1): 0.80
    rising toward low frequencies, as any AR(1)/integrating system does. A bump in the 3-7
    dangerous: a significance test against an AR(1) background spectrum is required.
```

## 05 trends changepoints

[Source](../analyses/05_trends_changepoints.py) · [Full output](outputs/05_trends_changepoints.txt)

```text
GISTEMP 1880-2025 (C)            slope=+0.083/decade  SE naive=0.004  SE AR1-adj=0.012  SE HAC=0.008  r1(resid)=0.83  n_eff=14/146
GISTEMP 1975-2025 (C)            slope=+0.207/decade  SE naive=0.010  SE AR1-adj=0.014  SE HAC=0.014  r1(resid)=0.33  n_eff=25/51
Central Park TMIN 1950-2025 (C)  slope=+0.229/decade  SE naive=0.032  SE AR1-adj=0.036  SE HAC=0.040  r1(resid)=0.12  n_eff=59/76
Potomac mean flow (kcfs)         slope=+0.164/decade  SE naive=0.153  SE AR1-adj=0.167  SE HAC=0.152  r1(resid)=0.08  n_eff=80/95
    1975-2025 the residual autocorrelation is modest. Always print r1 of the residuals
GISTEMP 1975-2025          MK z=+8.51 p=0.0000 | Hamed-Rao z=+2.66 p=0.0079 (r1=0.33) | Sen +0.206/decade (95% +0.182, +0.227)
Central Park TMIN          MK z=+5.71 p=0.0000 | Hamed-Rao z=+2.20 p=0.0279 (r1=0.12) | Sen +0.233/decade (95% +0.165, +0.300)
Potomac mean flow (kcfs)   MK z=+0.56 p=0.5759 | Hamed-Rao z=+0.56 p=0.5759 (r1=0.09) | Sen +0.076/decade (95% -0.201, +0.370)
0.0 size    0.048      0.051       0.050      0.059
    power   0.825      0.732       0.824      0.788
0.3 size    0.124      0.053       0.150      0.090
    power   0.738      0.378       0.740      0.578
0.6 size    0.327      0.043       0.415      0.094
    power   0.585      0.117       0.652      0.273
    even more than raw MK (Bayazit & Onoz 2007; Hamed 2009). The Hamed-Rao variance
GISTEMP annual       Pettitt change after 1937  K=4903  p=2.04e-20
Potomac mean flow    Pettitt change after 1970  K=518  p=3.12e-01
Central Park TMIN    Pettitt change after 1997  K=1002  p=2.62e-06
GISTEMP hinge: breakpoint=1975, slope before=+0.037/decade, after=+0.206/decade
AIC linear=-61.3  AIC hinge=-204.1  (hinge uses 2 extra params incl. breakpoint)
    changed; for global temperature the acceleration lands in the 1970s. The breakpoint is
Endpoint sensitivity: trend to 2025 as a function of start year
1950   0.161  0.137  0.185  76
1975   0.207  0.180  0.234  51
1998   0.253  0.203  0.304  28
2005   0.297  0.214  0.380  21
2010   0.363  0.274  0.452  16
2015   0.299  0.113  0.486  11
    has since filled in: the 1998-2025 trend is now the steepest in the table). For the
    11-year window the HAC interval runs 0.11 to 0.49 C/decade, from well below to nearly
    triple the 1950-2025 rate. Short-window trends are endpoint noise; the honest answer
  annual mean                  MK p=0.576  Sen slope=+0.7% of mean per decade
  annual 10th pct (low flow)   MK p=0.782  Sen slope=+0.4% of mean per decade
  annual max                   MK p=0.956  Sen slope=-0.1% of mean per decade
    multiple testing (section 02).
```

## 06 extremes

[Source](../analyses/06_extremes.py) · [Full output](outputs/06_extremes.txt)

```text
n=95 water years; median peak 109 kcfs; record 426 kcfs in WY1936
GEV: xi=+0.181 loc=89.0 scale=41.6   Gumbel: loc=93.4 scale=45.3
AIC GEV=1033.7  AIC Gumbel=1039.2  LRT p=0.006
parametric bootstrap: 1000 fits (0 discarded as implausible); xi 95% CI = (+0.01, +0.35)
          GEV  Gumbel    LP3  GEV 2.5%  GEV 97.5%
2       105.0   110.0  104.0      94.0      117.0
5       161.0   161.0  163.0     141.0      182.0
10      205.0   195.0  208.0     174.0      239.0
25      269.0   238.0  273.0     215.0      338.0
50      325.0   270.0  327.0     243.0      436.0
100     388.0   302.0  385.0     275.0      563.0
200     458.0   333.0  449.0     304.0      723.0
500     567.0   375.0  544.0     342.0     1008.0
    over Gumbel (p=0.006), yet the bootstrap CI on xi runs +0.01 to +0.35: 95 years
    constrain the tail only loosely. The three models agree to ~10 years and diverge
    beyond 100, where the GEV 100-year level has a roughly +-37% interval. Extrapolating
    to 500 years from 95 years is a modelling assumption, not a measurement. Say so in the
threshold u=30 mm: 684 exceedance days -> 634 declustered events (8.34/yr); GPD xi=+0.110 sigma=16.6
2       85.0        76.0
5      107.0       104.0
10     125.0       125.0
25     151.0       153.0
50     172.0       176.0
100    195.0       200.0
200    220.0       225.0
500    256.0       262.0
GEV on annual maxima: xi=+0.094 (n=76)
    Both approaches put the 100-year daily rainfall near 195-200 mm. Declustering matters:
    without it, a 3-day storm counts as 3 'independent' events and lambda is inflated. The
stationary: mu=25.28 sigma=1.31 xi=-0.22
trend:      mu0=24.44 mu1=+0.241 C/decade sigma=1.26 xi=-0.28
likelihood-ratio stat=11.67  p=0.0006
20-yr return level: 27.0 C in 1950 -> 28.8 C in 2025; the 1950 20-yr night is now a 1-in-2-year event
    worth knowing: it generalises to any covariate. The negative shape (xi=-0.28) is
    shifted by 1.8 C since 1950, which turns a 1950s 20-year event into a roughly
```

## 07 resampling

[Source](../analyses/07_resampling.py) · [Full output](outputs/07_resampling.txt)

```text
Bootstrap CI for the mean and median of wet-day precipitation (2015, n=110)
  mean      est=  9.46  percentile=( 7.37,11.68)  basic=( 7.25,11.56)  BCa=( 7.55,11.90)
  median    est=  4.95  percentile=( 3.15, 7.60)  basic=( 2.30, 6.75)  BCa=( 3.14, 7.60)
  90th pct  est= 28.53  percentile=(18.30,34.88)  basic=(22.18,38.76)  BCa=(18.40,35.60)
  mean: classical t-interval = (7.27, 11.65)
Coverage check: which interval actually has 95% coverage? (lognormal, n=30)
  mean:t             coverage=0.843
  mean:percentile    coverage=0.840
  mean:BCa           coverage=0.875
  median:percentile  coverage=0.930
  median:BCa         coverage=0.935
Trend SE: iid bootstrap vs block bootstrap vs HAC (GISTEMP monthly 1975-2025)
slope = 0.0207 C/yr   residual r1=0.67   n_eff=121 of 612
  SE OLS (naive)             0.00042
  SE iid residual bootstrap  0.00043
  SE iid case bootstrap      0.00045
  SE block bootstrap (24 mo) 0.00112
  SE HAC Newey-West (24 lag) 0.00133
  SE AR(1) n_eff adjusted    0.00094
    autocorrelation that makes the naive SE wrong. Block bootstrap keeps 2-year chunks
    intact and lands in the same place as HAC and the AR(1) adjustment. If the resampling
skewness=0.801  jackknife bias=-0.061  jackknife SE=0.364  bootstrap SE=0.304
bias-corrected skewness = 0.862
Monte Carlo: calibrate a trend t-test against an AR(1) null (GISTEMP annual 2000-2025)
slope=+0.266 C/decade  t=10.16  naive p=0.0000  residual phi=0.41
Monte Carlo AR(1) null: p=0.0000  (null SD of t = 1.61 vs 1.0 for iid)
    ask how often the t-statistic is as large as observed. The null SD of t exceeds 1 by
    roughly sqrt(n/n_eff), the same factor as the analytic adjustment in section 05.
```

## 08 multivariate

[Source](../analyses/08_multivariate.py) · [Full output](outputs/08_multivariate.txt)

```text
361 stations x 15 climate variables
tavg  htdd   -0.986
      tmin    0.972
      tmax    0.971
tmin  htdd   -0.962
dt32  emnt   -0.957
tmax  htdd   -0.955
    dt32...), 'how wet' (prcp, dp01, dp10), 'how snowy / continental'. PCA will find
explained variance by component: [0.637 0.226 0.057 0.025 0.024 0.012] ... cumulative(3) = 0.921
Kaiser rule (eigenvalue > 1): keep 2 components
       PC1   PC2   PC3
tavg -0.32  0.06  0.14
tmax -0.32 -0.04  0.11
tmin -0.30  0.15  0.16
prcp -0.05  0.50  0.18
snow  0.26  0.03  0.31
dp01  0.15  0.46  0.05
dp10  0.06  0.52  0.05
dx90 -0.27 -0.19  0.36
dx32  0.29 -0.05  0.41
dt32  0.29 -0.16 -0.00
dt00  0.24 -0.13  0.51
emxt -0.20 -0.36  0.03
emnt -0.30  0.12 -0.11
htdd  0.32 -0.08  0.02
cldd -0.28  0.00  0.48
    with the other), PC2 is a wetness axis (prcp, dp01, dp10), PC3 is continentality: many
    the large annual range of the interior. Three components carry >90% of the variance in
  k= 2 inertia=  3092.0 silhouette=0.371
  k= 3 inertia=  2332.9 silhouette=0.308
  k= 4 inertia=  1920.5 silhouette=0.329
  k= 5 inertia=  1580.7 silhouette=0.337
  k= 6 inertia=  1369.0 silhouette=0.337
  k= 7 inertia=  1226.0 silhouette=0.334
  k= 8 inertia=  1109.6 silhouette=0.336
```

## 09 spatial

[Source](../analyses/09_spatial.py) · [Full output](outputs/09_spatial.txt)

```text
366 stations; median nearest-neighbour distance 95 km
  tavg (raw)                       I=+0.866  E[I]=-0.003  z=+35.0  p_perm=0.002
  OLS residual                     I=+0.422  E[I]=-0.003  z=+16.9  p_perm=0.002
  shuffled residual (null check)   I=-0.018  E[I]=-0.003  z=-0.7  p_perm=0.536
    regression's standard errors in section 03 were therefore too small, exactly as with
spherical:   nugget=0.552 sill=1.798 range=3132 km
exponential: nugget=0.495 sill=1.594 range=3564 km (practical)
residual variance = 1.576
    differences, microclimate. Semivariance rises toward the sill over ~3132 km, the
regression kriging  0.802  0.535
OK on tavg          1.256  0.770
OLS trend only      1.274  0.911
IDW p=2             2.045  1.478
IDW p=1             3.720  2.952
```

## 10 bayesian

[Source](../analyses/10_bayesian.py) · [Full output](outputs/10_bayesian.txt)

```text
Beta-Binomial: probability that a July day in NYC is wet (>= 1 mm)
all data: 626/2356 wet; posterior mean=0.2659  95% credible=(0.2483, 0.2839)
frequentist Clopper-Pearson 95% CI=(0.2480, 0.2840)
    season they have nearly converged; with 2,300 days the prior is irrelevant and the
Gamma-Poisson: hot days (TMAX >= 35 C) and warm nights (TMIN >= 23.9 C) per summer
hot days TMAX>=35C
  early:  107 events / 30 summers -> rate 3.54 (2.90, 4.24)   var/mean=4.3
  late:   107 events / 30 summers -> rate 3.54 (2.90, 4.24)   var/mean=4.1
  rate ratio late/early: median 1.00, 95% credible (0.76, 1.31);  P(ratio > 1) = 0.499
warm nights TMIN>=23.9C
  early:  164 events / 30 summers -> rate 5.41 (4.62, 6.27)   var/mean=4.1
  late:   245 events / 30 summers -> rate 8.07 (7.09, 9.10)   var/mean=3.7
  rate ratio late/early: median 1.49, 95% credible (1.23, 1.82);  P(ratio > 1) = 1.000
    ratio's credible interval and P(ratio > 1) fall out with no delta method. Hot days
    show no change (consistent with the TMAX null in section 02); warm nights have become
    var/mean = 1, but heat waves cluster events, so the variance is several times the mean
Metropolis-Hastings: Bayesian linear trend for GISTEMP 1975-2025 (annual)
acceptance rates: [0.31 0.31 0.31 0.31]  (target ~0.25-0.45)
  intercept          R-hat=1.000  ESS=5585 of 64000
  slope (C/decade)   R-hat=1.000  ESS=8939 of 64000
  log sigma          R-hat=1.000  ESS=3531 of 64000
slope posterior: mean=0.2074  95% credible=(0.1873, 0.2276)
OLS:             est =0.2074  95% CI      =(0.1871, 0.2277)
P(slope > 0.15 C/decade) = 1.000;  P(slope > 0.2) = 0.767
with a sceptical prior N(0, 0.05) on the slope: posterior mean=0.1987 (shrunk from 0.2074); the data still overwhelm it
    acceptance tuning, R-hat near 1 and a healthy ESS, then a posterior that matches OLS
    0.15 C/decade) is a direct probability statement, which no p-value provides. In
```

## 11 ml

[Source](../analyses/11_ml.py) · [Full output](outputs/11_ml.txt)

```text
27639 days; train through 2010 (22160), test 2011-2025 (5479); 12 precipitation/season features (+1 persistence feature)
Regression: log10(discharge) from antecedent weather. Test period 2011-2025
                                   rmse_log  r2_log  NSE_flow  KGE_flow
climatology (day-of-year mean)        0.346   0.312     0.053     0.034
persistence (yesterday's flow)        0.097   0.946     0.776     0.888
linear regression [weather only]      0.241   0.665     0.416     0.490
ridge (CV alpha) [weather only]       0.241   0.665     0.415     0.489
lasso (CV alpha) [weather only]       0.241   0.666     0.415     0.486
random forest [weather only]          0.242   0.664     0.400     0.444
gradient boosting [weather only]      0.245   0.656     0.383     0.475
gradient boosting [+ persistence]     0.083   0.960     0.846     0.886
lasso kept 11 of 12 features; dropped: ['p30']
    the last few months of weather at one gauge? All five land at R2 0.66-0.67 in log
    the features (one rain gauge for an 11,500 sq mi basin), not model flexibility. Adding
shuffled KFold R2:    0.792 +- 0.007
TimeSeriesSplit R2:   0.612 +- 0.087  (60-day gap between train and test)
true holdout 2011-25: 0.664
    30-120 day rolling features they are almost the same row, so the score is inflated.
sin_doy    0.5506
api        0.2287
tmax30     0.2135
p120       0.1523
p60        0.0373
p30        0.0319
p7         0.0256
p14        0.0241
cos_doy    0.0234
p1         0.0072
p3         0.0071
snow30     0.0048
    temperature (again evapotranspiration) and the 120-day rainfall total. Short windows
    (1-7 days) matter little because a single DC gauge is a poor proxy for a storm over an
    11,500 sq mi basin, while months of accumulated rainfall control baseflow. The partial
    dependence plots show the shape of each effect (flow rises with p120, falls with
    tmax30). Permutation importance on held-out data avoids the impurity-importance bias
```

## 12 assay eda

[Source](../analyses/12_assay_eda.py) · [Full output](outputs/12_assay_eda.txt)

```text
n             470.000
mean          435.299
median        424.000
sd            299.882
CV              0.689
skew            0.459
min             0.000
P90           817.380
P99          1207.589
max          1528.100
frac zero       0.047
exhaustive truth: mean=278.0 median=221.2 CV=0.90
lognormal fit to V>0: sigma_ln=1.09  Shapiro on ln(V>0) p=0.000
1     45   39.54    18.3   51.73  1.31
2    425  477.20   473.3  284.28  0.60
    be. CV of 0.69 is moderate for a base-metal style variable. The two T types have very
naive mean=435.3  declustered mean (cell=20)=287.3  truth=278.0
weights: min=0.29 max=2.41 (1 = average sample)
bias removed: 94% of the naive bias
GSLIB cluster.dat: naive mean=4.35, declustered (ours)=2.59, declustered (shipped weights)=2.53; weight correlation r=0.96
    15).
P95     916.920     24     0.051      435.299     426.516          0.020      0.689     0.656
P97.5  1011.932     12     0.026      435.299     430.167          0.012      0.689     0.666
P99    1207.589      5     0.011      435.299     433.427          0.004      0.689     0.679
P99.5  1346.818      3     0.006      435.299     434.445          0.002      0.689     0.684
Within T=2 only: P99=1213, max=1528, top 1% of samples hold 2.8% of the metal
    from 0.69 to 0.68. For this moderately skewed variable capping barely matters; for
    nuggety gold the same table can remove 10-20% of the metal, and the choice becomes the
Compositing irregular samples to 2-unit lengths (drillhole built from Walker Lake column)
raw samples: n=163, lengths 0.5-3.5, mean length=1.83, grade sd=177.8
2-unit composites: n=149, grade sd=177.8;  4-unit composites: n=75, grade sd=174.6
metal (grade x length) raw=58720 composited=58720 -> conserved
naive mean of raw sample grades=197.2 vs length-weighted mean=196.7
    the arithmetic mean of raw assays is biased high. The 2-unit composite sd (178) barely
    differs from the raw sd (178) because the raw mean length is already ~2, but 4-unit
```

## 13 compositional

[Source](../analyses/13_compositional.py) · [Full output](outputs/13_compositional.txt)

```text
359 samples, 7 metals; rock types: {'Kimmeridgian': 124, 'Sequanian': 89, 'Argovian': 76, 'Quaternary': 64, 'Portlandian': 6}
independent amounts: pairwise r = [-0.03  0.02 -0.04]
after closing to 100%: pairwise r = [-0.5  -0.48 -0.52]
clr of closed data:     r = [-0.51 -0.49 -0.51] (clr correlations are not interpretable as independence either: rows sum to zero)
cd  1.00  0.26  0.58  0.15  0.49  0.22  0.62
co  0.26  1.00  0.48  0.19  0.74  0.16  0.44
cr  0.58  0.48  1.00  0.21  0.71  0.26  0.61
cu  0.15  0.19  0.21  1.00  0.22  0.82  0.66
ni  0.49  0.74  0.71  0.22  1.00  0.27  0.59
pb  0.22  0.16  0.26  0.82  0.27  1.00  0.67
zn  0.62  0.44  0.61  0.66  0.59  0.67  1.00
cd  0.00  0.40  0.29  0.76  0.26  0.45  0.25
co  0.40  0.00  0.13  0.56  0.08  0.32  0.16
cr  0.29  0.13  0.00  0.51  0.09  0.24  0.10
cu  0.76  0.56  0.51  0.00  0.58  0.23  0.30
ni  0.26  0.08  0.09  0.58  0.00  0.31  0.13
pb  0.45  0.32  0.24  0.23  0.31  0.00  0.13
zn  0.25  0.16  0.10  0.30  0.13  0.13  0.00
most proportional pairs: co/ni (0.08), cr/ni (0.09), cr/zn (0.10)
least proportional pairs: co/cu (0.56), cu/ni (0.58), cd/cu (0.76)
explained variance PC1=0.54 PC2=0.24 (clr covariance has rank D-1=6)
         1  1.6680 0.1569
         2 17.3443 0.0000
         3  6.1649 0.0001
         4 36.0841 0.0000
         5  3.4529 0.0087
         6  7.5543 0.0000
Argovian vs Kimmeridgian: Euclidean distance between mean compositions=25.0 ppm (dominated by Zn), Aitchison distance=0.61 (scale-free)
```

## 14 qaqc

[Source](../analyses/14_qaqc.py) · [Full output](outputs/14_qaqc.txt)

```text
400 pairs; HARD median=0.060, 90th percentile=0.164 (passes the 90% < 20% criterion)
pairs with HARD > 0.5: 3  (planted swaps: 6, of which flagged: 3)
Thompson-Howarth: sd = 0.010 + 0.100 * c   -> precision at 1 g/t = 11% (generated with sd = sqrt((0.12 c)^2 + 0.02^2), i.e. 12% at high grade, 0.02 g/t floor)
mean duplicate / mean original = 1.014 (should be ~1; a systematic departure means the duplicate stream is treated differently)
    how one should state precision: 'X% at grade Y', never a single percentage. Only 3 of
    the 6 planted swaps show up as HARD > 0.5: a swap between two samples of similar grade
CRM-A: certified 0.52; 2 points beyond 3SD (planted mix-ups in this CRM: 2); 6 warnings; overall bias (ex-outliers)=+0.6%; CUSUM signals at batch 78
CRM-B: certified 1.85; 4 points beyond 3SD (planted mix-ups in this CRM: 1); 7 warnings; overall bias (ex-outliers)=+1.1%; CUSUM signals at batch 78
CRM-C: certified 4.10; 2 points beyond 3SD (planted mix-ups in this CRM: 0); 10 warnings; overall bias (ex-outliers)=+0.8%; CUSUM signals at batch 71
    only chart passes a lab that is over-reporting grade by 8%. The CUSUM, or simply 'how
    many of the last 10 points are above the line', catches it within a few batches. Bias
180 blanks; 167 below DL; 7 above 5xDL (planted carry-over events: 7, flagged: 7)
mean blank after a >3 g/t sample: 31.7 ppb vs 2.1 ppb otherwise
Lead           Dissolved    121   335  36.0
Cadmium        Dissolved    114   327  35.0
               Total         89   310  29.0
Zinc           Dissolved     16   335   5.0
Lead           Recoverable   13   316   4.0
Copper         Dissolved      3   334   1.0
using dissolved lead: 335 results, 121 non-detects (36%), detection limits used: [np.float64(0.04), np.float64(0.06), np.float64(0.08), np.float64(0.12), np.float64(0.5), np.float64(0.6), np.float64(1.0), np.float64(5.0), np.float64(100.0)]
   substitute 0  0.23    0.06  0.84
substitute DL/2  1.01    0.19  5.47
  substitute DL  1.79    0.19 10.88
   Kaplan-Meier  0.28    0.11   NaN
ROS (lognormal)  0.63    0.19  1.36
    method reported '< 100'), substitution alone moves the mean from 0.23 to 1.79 ug/L
    against a Kaplan-Meier estimate of 0.28. Kaplan-Meier needs no distributional
    when detection limits changed over time, as they did here. Above ~50% censoring, stop
```

## 15 geostatistics

[Source](../analyses/15_geostatistics.py) · [Full output](outputs/15_geostatistics.txt)

```text
470 samples -> 3120 blocks of 5x5; true block mean=278.0, true block sd=228.7, point sd=249.8, sample sd=299.6
omnidirectional spherical: nugget=29165 sill=92634 range=38.2  (sample variance 89738); relative nugget=0.31
  azimuth   0: range=61.7
  azimuth  45: range=27.0
  azimuth  90: range=30.5
  azimuth 135: range=43.2
anisotropy: major axis azimuth 0 (range 62), minor 45 (range 27), ratio=0.44
    Directionally, continuity is longest along azimuth 0 and shortest along 45 (ratio
    0.44): the mineralisation is elongated, as the sample map suggests, and the kriging
Estimating 5x5 blocks: OK vs IDW vs nearest neighbour, scored against truth
OK (block, anisotropic)    275.03  -2.95  114.86  0.87  215.46                0.92
OK (point)                 275.47  -2.51  116.11  0.87  218.07                0.91
IDW p=2                    334.52  56.54  138.24  0.84  181.35                1.05
nearest neighbour          275.82  -2.16  150.17  0.80  245.10                0.75
truth: mean=278.0 sd=228.7;  declustered sample mean=287.3
    essentially unbiased (-3 and -2): both decluster implicitly, kriging through the
    and carries a +57 bias straight from the preferential sampling. OK has the lowest RMSE
    and the truth-on-estimate slope closest to 1 (least conditional bias). Every estimator
    is smoother than the truth (sd 229); NN is the least smooth and the least accurate.
LOO: mean error=+7.1  RMSE=183.2  r=0.79  slope(truth on est)=0.94
standardised errors: mean=+0.03 variance=0.73 (want ~0 and ~1)
    error variance of 0.73 says the kriging variance is mis-scaled for the variogram sill
    used; a value far from 1 means the nugget/sill is wrong. The slope below 1 is
1         78000   62423.2  249.8
2         19500   57511.9  239.8
5          3120   52304.1  228.7
10          780   46753.8  216.2
20          195   37810.4  194.4
50           36   20074.3  141.7
variance reduction factor point -> 5x5: 0.84  (Krige's relation: D2(point,domain) = D2(point,block) + D2(block,domain))
at cut-off 500: tonnes above (blocks) truth=521  OK=437  IDW=474  point-support truth=587
    of 500, point-support truth reports +13% tonnes relative to the block truth. Smoothed
normal-score variogram: nugget=0.30 range=45.9
realisation sd (mean over 8): 285  vs point truth sd 250, E-type sd 222, OK sd 218
realisation mean: 291  (declustered sample mean 287; truth 278)
```

## 16 sampling theory

[Source](../analyses/16_sampling_theory.py) · [Full output](outputs/16_sampling_theory.txt)

```text
                                                       c_g_cm3  mass_10pct_FSE_1cm_g  mass_10pct_FSE_2mm_g  mass_10pct_FSE_75um_g
gold 1.5 g/t (native Au, 50 um liberation, alpha 1.5) 1.29e+07              5.69e+04              5.09e+03                   36.9
gold 1.5 g/t (classic Gy alpha 2.5)                   1.29e+07              1.14e+07              2.03e+05                   55.4
copper 1% as chalcopyrite (100 um)                         140                   175                  3.12               0.000737
iron ore 60% Fe as hematite (500 um)                     0.512                  1.43                0.0256                2.7e-06
    from iron ore. At 1 cm top size a representative gold sample weighs many kilograms to
    tonnes; only after pulverising to 75 microns does a sub-sample of tens of grams carry
    a 10% fundamental error. The liberation factor is the least knowable input: Gy's
    classic sqrt(d_l/d) (alpha 2.5) is far more pessimistic than the Francois-Bongarcon
    alpha 1.5 calibration that practitioners fit to gold duplicate data. Protocols are
A sample preparation protocol on the Gy nomogram (gold, 1.5 g/t)
core as received           2.50    3000        86.56  crush/receive
crush to 1 cm              1.00    3000        43.54  crush/receive
split to 1 kg              1.00    1000        75.41          split
crush to 2 mm              0.20    1000        22.55  crush/receive
split to 250 g             0.20     250        45.10          split
pulverise to 75 um         0.01     250         3.84  crush/receive
30 g fire assay            0.01      30        11.10          split
total FSE of protocol (splits add in variance): 88.6% relative
core as received           2.50    3000        86.56  crush/receive
crush to 2 mm              0.20    3000        13.02  crush/receive
split to 1 kg              0.20    1000        22.55          split
pulverise to 75 um         0.01    1000         1.92  crush/receive
50 g fire assay            0.01      50         8.59          split
total FSE of protocol (splits add in variance): 24.1% relative
    error (89%) is dominated by the 'split to 1 kg' stage: splitting while the gold is
    still locked in centimetre particles. Protocol B crushes everything to 2 mm before the
    first split and pulverises the whole kilogram, cutting the total to 24%. Changing the
    final aliquot from 30 to 50 g is nearly irrelevant. Match the top size to the mass you
variogram nugget=29165, sill=92634 -> relative nugget 31%
duplicate relative sd of a single assay=0.161; at the mean grade 435 that is a measurement variance of 4910
measurement error explains 17% of the nugget; the rest (83%) is geological variability at distances shorter than the closest sample pairs
How many samples to estimate the domain mean to +-10%?
CV=0.69 -> n=182 independent samples for +-10% at 95%; we have 470 samples but an effective n of 289 after declustering weights (Kish), and spatial correlation lowers it further
  CV=0.5: n=96
```

## 17 reliability

[Source](../analyses/17_reliability.py) · [Full output](outputs/17_reliability.txt)

```text
engine         151          12  2.266  2266.936       2.40    2400.0      2008.008    1977.521                    1912.483  20.095         0.990
hydraulics     324          11  0.907   925.575       0.85     900.0       969.618     931.200                     934.069   6.843         0.993
tyres          204          12  3.077  1681.091       3.00    1800.0      1502.881    1491.547                    1454.827   4.215         0.997
    failures where scheduled replacement would waste parts; tyres beta=3: wear-out where
    the ones that get censored. Inherent availability ranges from 0.990 to 0.997; at fleet
optimal age-replacement interval ~790 h (0.47 eta) at failure/PM cost ratio 6; run-to-failure cost rate=0.00399 vs optimum 0.00190 per h
    < 1) the same curve is monotone decreasing: replace only on failure.
failures per truck: {1: 60, 2: 48, 3: 65, 4: 60, 5: 55, 6: 58, 7: 53, 8: 52, 9: 54, 10: 59, 11: 56, 12: 59}
chi-square homogeneity: chi2=4.0 df=11 p=0.971  (dispersion var/mean=0.36)
    real fleet with p < 0.05 and var/mean >> 1 has a bad actor or a maintenance-history
Production SPC: detecting a sustained 7% drop in daily tonnes
planted shift at day 601: -7% = 1.15 sigma of daily noise
Shewhart I-chart: false alarms before shift=4, first signal after shift: day 647
CUSUM (k=0.5, h=5): false alarms before shift=1, first signal after shift: day 616
EWMA (lambda=0.2): false alarms before shift=4, first signal after shift: day 609
    day 601: Shewhart 46 days, CUSUM 15 days, EWMA 8 days. The accumulating charts win on
    delay but pay in false alarms over a 600-day in-control run (ARL0 ~ 465 for these
```

## 18 safety

[Source](../analyses/18_safety.py) · [Full output](outputs/18_safety.txt)

```text
345,316 mine-years with hours, 220,087 injuries, 1,152 fatalities, 6.12 million worker-years
National injury and fatality rates per 200,000 hours
2000    8.15            5.09
2005    5.70            4.07
2010    4.33            2.87
2015    3.66            2.50
2020    3.36            1.92
2024    3.37            1.95
coal            injury rate ratio per year=0.961 (0.960-0.962); deviance/df=28.2 (Poisson, should be ~1)
metal/nonmetal  injury rate ratio per year=0.959 (0.958-0.959); deviance/df=14.0 (Poisson, should be ~1)
    GLM reports a deviance/df far above 1: even aggregated to the national level the year-
2019, mines with >= 50,000 h: n=1674; Poisson deviance/df=2.24; NB alpha=0.46
coal vs metal/nonmetal rate ratio: Poisson 1.66 (SE 0.029) | NB 1.57 (SE 0.054)
    is 1.9x the Poisson one. With alpha=0.46, a mine with an expected 10 injuries has a
    variance of 10 + alpha x 100: mines differ in underlying risk (commodity, method,
Gamma-Poisson prior by maximum likelihood: mean rate 1.89, sd 1.57 per 200k h (a=1.45, b=0.77); a mine's rate is shrunk toward 1.89 with weight b/(b+E) = 0.81 at median exposure
worst 5 by raw rate (2020-24):
4201852       1.0   2186  91.5      3.2
2200460       4.0   9059  88.3      6.7
0503791       1.0   2266  88.3      3.2
0505125       1.0   2286  87.5      3.2
4608549       1.0   2324  86.1      3.2
worst 5 by shrunken rate:
2700083      45.0  411090  21.9     16.5
1103182     118.0 1700575  13.9     12.9
2302273      33.0  390836  16.9     12.7
1103141     159.0 2547870  12.5     11.9
4601437     540.0 9390674  11.5     11.3
95% of mines with < 20,000 h over five years recorded zero injuries (rate 0.0)
mine 0200024 (Freeport-McMoRan Morenci Inc.): centre line 3.38 per 200k h; 11 of 25 years outside Poisson limits: [2001, 2006, 2007, 2008, 2010, 2011, 2012, 2016, 2017, 2022, 2023]
Poisson sigma at median exposure=0.32; observed year-to-year sigma (I-chart)=0.73
    at tens of millions of hours per year the Poisson limits are so tight (sigma 0.32)
    that ordinary year-to-year variation (0.73 by the moving range) breaches them. Rates
days lost (where > 0): n=105,719  mean=56.9  median=24  P90=164  P99=350  max=3616
share of all lost days from the worst 10% of cases: 45%
```

## 19 economics

[Source](../analyses/19_economics.py) · [Full output](outputs/19_economics.txt)

```text
Log returns, volatility, correlation (monthly, 1990-2024)
copper      20.8%
gold        11.9%
iron_ore    27.8%
nickel      27.9%
zinc        21.7%
copper      1.00  0.27      0.33    0.55  0.63
gold        0.27  1.00      0.16    0.27  0.21
iron_ore    0.33  0.16      1.00    0.30  0.27
nickel      0.55  0.27      0.30    1.00  0.51
zinc        0.63  0.21      0.27    0.51  1.00
copper monthly return skew=-0.61, excess kurtosis=4.83 (normal = 0; fat tails)
copper   ADF p (unit root)=0.630  AR(1) b=0.9965  half-life=16.6 yr  long-run mean=10,054
gold     ADF p (unit root)=0.991  AR(1) b=1.0022 >= 1: no mean reversion estimable (pure random walk or trending)
    random walk is low. Yet the two assumptions give very different 10-year price fans and
10-year 90% band width: GBM 9.2x, mean reversion 6.2x
base case NPV = $0.91 billion at copper 8,916 $/t
tornado (+-10% one at a time, $bn):
discount         1.03  0.79   0.25
capex            1.07  0.75   0.32
tonnes_per_year  0.65  1.16   0.50
opex_per_t       1.19  0.62   0.57
grade_pct        0.37  1.44   1.07
recovery         0.37  1.44   1.07
price            0.37  1.44   1.07
Monte Carlo (20,000 runs): mean NPV=$0.91bn, median=$0.73bn, P10=$-0.72bn, P90=$2.75bn, P(NPV<0)=26.9%
    clipped), and asymmetric inputs (capex overruns). The result is a mean NPV of $0.91bn
    against a deterministic $0.91bn and a 27% chance of destroying value, which is the
g_c = opex / (price x recovery x payability) = 0.292% Cu at 8,916 $/t; at -30% price: 0.417%
commit now: expected NPV $0.91bn (P(<0)=32%); wait 2 years and build only if positive: $0.98bn -> option value $+0.06bn
```

## 20 environmental

[Source](../analyses/20_environmental.py) · [Full output](outputs/20_environmental.txt)

```text
309 total-recoverable copper samples 1985-2023 with matched daily flow; 0 non-detects
      median_cu   n    p90
1985       51.5   8   95.0
1991       64.5   4  324.0
1996       96.0   7  528.0
2001       25.0   8   79.3
2006       30.5   7   72.6
2011       92.6  12  222.4
2016       41.9   8   60.4
2021       25.5   8   36.9
log10 C = -0.59 + 0.90 log10 Q   (R2=0.53, n=309)
a doubling of flow multiplies copper by 1.87
monthly median Cu (ug/L): 1:38, 2:62, 3:36, 4:38, 5:76, 6:79, 7:19, 8:14, 9:38, 10:16, 11:16, 12:24
raw log10 Cu               Kendall tau=-0.24 p=0.0343  Sen slope=-0.77 log10-units/century = -1.8%/yr
flow-adjusted residual     Kendall tau=-0.51 p=0.0000  Sen slope=-1.00 log10-units/century = -2.3%/yr
1969-90   35      54.00        271.0        1.59
1991-05  116      32.05        241.0        1.01
2006-15   83      29.90        289.0        0.78
2016-23   75      29.60        339.0        0.79
1969-90   35      34  0.971  0.855  0.995
1991-05  116     107  0.922  0.859  0.959
2006-15   83      74  0.892  0.807  0.942
2016-23   75      64  0.853  0.756  0.916
fraction of exceedances occurring on top-quartile flow days: 28% (exceedance rate on those days 100% vs 87% otherwise)
    of exceedances this year' is a count without a denominator. The decline from 97% to
    85% is real but the river still exceeds this benchmark in most samples at any flow: a
    when it is worst, not when it is compliant. The criterion (13.0 ug/L) is the
    hardness-100 chronic value for illustration; the real standard moves with each
CUSUM of flow-adjusted annual medians peaks in 1996: regime shift from above-average to below-average residuals
  cadmium 1969-90: n=87, non-detects=69%, KM median=0.00 ug/L, KM mean=0.58, naive mean (DL as value)=1.77
  cadmium 1991-05: n=234, non-detects=61%, KM median=0.07 ug/L, KM mean=0.16, naive mean (DL as value)=0.45
  cadmium 2006-23: n=316, non-detects=0%, KM median=0.10 ug/L, KM mean=0.18, naive mean (DL as value)=0.18
    expanded in the early- to-mid 1990s; the 2010s floodplain tailings removals are a
    detect fraction was 60-70% in the early eras with changing detection limits, so only
    'DL as value' means are 2-3x too high.
```

## 21 censored inference

[Source](../analyses/21_censored_inference.py) · [Full output](outputs/21_censored_inference.txt)

```text
recoverable copper 309 samples (0 non-detects), dissolved lead 311 (117 non-detects), total arsenic 310 (0 non-detects)
years sampled per calendar month: 1:2, 2:7, 3:35, 4:34, 5:39, 6:34, 7:30, 8:30, 9:6, 10:16, 11:15, 12:3
raw log10 Cu    seasonal Kendall tau=-0.122 z=-2.60 p=9.19e-03 Sen=-1.04%/yr   |   annual-median Mann-Kendall tau=-0.236 p=3.43e-02 Sen=-1.76%/yr
flow-adjusted   seasonal Kendall tau=-0.350 z=-7.49 p=7.14e-14 Sen=-1.98%/yr   |   annual-median Mann-Kendall tau=-0.555 p=6.63e-07 Sen=-2.23%/yr
did the sampled flows drift? Kendall tau of log10 Q on time = +0.089 (p=0.019)
1             2   -1.0 -1.000 -0.022   -5.043
2             7   -9.0 -0.429 -0.035   -7.734
3            35 -169.0 -0.284 -0.008   -1.798
4            34 -201.0 -0.358 -0.008   -1.834
5            39 -285.0 -0.385 -0.008   -1.821
6            34 -129.0 -0.230 -0.007   -1.505
7            30 -207.0 -0.476 -0.015   -3.343
8            30 -191.0 -0.439 -0.010   -2.213
9             6   -3.0 -0.200 -0.007   -1.680
10           16  -46.0 -0.383 -0.010   -2.231
11           15  -15.0 -0.143 -0.006   -1.290
12            3   -1.0 -0.333 -0.008   -1.751
    0.009 to p = 7e-14: copper rides on discharge (section 20), so the flow the sampler
    test. The sampled flows themselves drifted over the record (tau +0.09), so part of the
    raw trend is a flow artefact. Report the flow-adjusted slope (-1.98%/yr) with the raw
    slope (-1.04%/yr) beside it. Modern USGS practice (WRTDS, Hirsch et al. 2010) fits the
  1969-90: n=35, non-detects=69%, detection limits [np.float64(1.0), np.float64(5.0)]
  1991-05: n=117, non-detects=74%, detection limits [np.float64(0.08), np.float64(0.5), np.float64(0.6), np.float64(1.0)]
  2006-15: n=84, non-detects=7%, detection limits [np.float64(0.04), np.float64(0.08), np.float64(0.12)]
  2016-23: n=75, non-detects=1%, detection limits [np.float64(0.06)]
                Gehan (censored)  3.13  0.00177
     Mann-Whitney, substitute DL   NaN 1.48e-15
  Welch t on logs, substitute DL  15.9 1.25e-26
   Mann-Whitney, substitute DL/2   NaN 5.96e-15
Welch t on logs, substitute DL/2  12.4 7.71e-20
      Mann-Whitney, substitute 0   NaN  0.00132
   Welch t on logs, substitute 0 -2.21    0.033
Gehan: 1969-90 is the higher group (positive scores favour the first label)
KM medians by era: 1969-90: 1.00 ug/L, 1991-05: 0.07 ug/L, 2006-15: 0.11 ug/L, 2016-23: 0.14 ug/L
    detects almost everything, with a median of 0.14 ug/L. Substituting the detection
```
