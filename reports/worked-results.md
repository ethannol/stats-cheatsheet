# Generated worked-example results

These are executed educational examples, not a new scientific study.
Data cutoff: `2026-09-30`. Source: `d618e0a4de4f330d1acd7cd0206423136d6f5d50`.
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
ma.S.L12      -0.8628      0.021    -40.736      0.000      -0.904      -0.821
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
spherical:   nugget=0.522 sill=1.746 range=2889 km
exponential: nugget=0.470 sill=1.750 range=4140 km (practical)
residual variance = 1.576
    differences, microclimate. Semivariance rises toward the sill over ~2889 km, the
regression kriging  0.803  0.537
OK on tavg          1.256  0.769
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
