# Generated worked-example results

These are executed educational examples, not a new scientific study.
Data cutoff: `2026-09-30`. Source: `515607f9b2b954e517d839c0886811f1db3f962a`.
The run manifest records the exact executed source, dependencies and hashes.
Full transcripts retain diagnostic warnings, assumptions and interpretation.
The excerpts below omit diagnostic warnings and read numerical output directly
from those hash-checked transcripts. Synthetic examples are labelled in source.

## 01 assay eda

[Source](../analyses/01_assay_eda.py) · [Full output](outputs/01_assay_eda.txt)

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
    04).
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

## 02 compositional

[Source](../analyses/02_compositional.py) · [Full output](outputs/02_compositional.txt)

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

## 03 qaqc

[Source](../analyses/03_qaqc.py) · [Full output](outputs/03_qaqc.txt)

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

## 04 geostatistics

[Source](../analyses/04_geostatistics.py) · [Full output](outputs/04_geostatistics.txt)

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

## 05 sampling theory

[Source](../analyses/05_sampling_theory.py) · [Full output](outputs/05_sampling_theory.txt)

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

## 06 reliability

[Source](../analyses/06_reliability.py) · [Full output](outputs/06_reliability.txt)

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

## 07 safety

[Source](../analyses/07_safety.py) · [Full output](outputs/07_safety.txt)

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

## 08 economics

[Source](../analyses/08_economics.py) · [Full output](outputs/08_economics.txt)

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

## 09 environmental

[Source](../analyses/09_environmental.py) · [Full output](outputs/09_environmental.txt)

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

## 10 censored inference

[Source](../analyses/10_censored_inference.py) · [Full output](outputs/10_censored_inference.txt)

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
    0.009 to p = 7e-14: copper rides on discharge (section 09), so the flow the sampler
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
