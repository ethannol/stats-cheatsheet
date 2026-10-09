# %% [markdown]
# # 03 Regression: OLS, diagnostics, robust SEs, multiple regression
#
# (a) Global temperature anomaly vs ln(CO2): the simplest climate sensitivity regression.
#     Physics: radiative forcing is ~ 5.35 ln(C/C0) W/m2, so T should be linear in ln(CO2).
# (b) Station mean temperature ~ latitude + elevation + longitude across ~370 US stations:
#     recovers the environmental lapse rate (~ -6.5 C per km) and the latitudinal gradient.
#
# **Cheat sheet**
# - OLS assumptions (LINE): Linearity, Independence of errors, Normality of errors
#   (only matters for small-n inference), Equal variance. Normality of X is NOT required.
# - Always check residual plots: vs fitted (nonlinearity, heteroscedasticity), QQ, vs time /
#   space (autocorrelation). Durbin-Watson ~2 means no lag-1 autocorrelation; <1 is bad.
# - Autocorrelated residuals leave coefficients unbiased but make SEs too small. Fixes:
#   Newey-West / HAC standard errors, GLS with AR(1) errors, or Cochrane-Orcutt.
# - Heteroscedasticity: Breusch-Pagan test; fix with HC (White) robust SEs or transform y.
# - Multicollinearity: VIF > 5-10 is a warning. Coefficients become unstable, predictions fine.
# - R^2 is not a measure of model validity. Adjusted R^2 / AIC / BIC for model comparison;
#   AIC penalises less, BIC prefers simpler models. Compare on the same y.
# - Standardised coefficients compare importance across predictors with different units.
# - Spatial data: residuals are spatially autocorrelated -> same SE problem (see 09).

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats as sps
from statsmodels.stats.diagnostic import acorr_ljungbox, het_breuschpagan
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.stats.stattools import durbin_watson

from earthstats import data
from earthstats.plotting import banner, note, savefig, setup
from earthstats.stats import ar1_effective_n

setup()
SEC = "03_regression"

# %% ----------------------------------------------------------- (a) temperature vs ln CO2
banner("(a) Simple OLS: annual GISTEMP anomaly ~ ln(CO2 / 280)")
co2 = data.co2().co2.resample("YE").mean()
gt = data.gistemp().anomaly.resample("YE").mean()
df = pd.DataFrame({"co2": co2, "T": gt}).dropna()
df = df[df.index.year <= 2025]  # full years only
df["lnco2"] = np.log(df.co2 / 280.0)  # 280 ppm ~ pre-industrial
df["year"] = df.index.year

m1 = smf.ols("T ~ lnco2", data=df).fit()
print(m1.summary().tables[1])
print(f"R^2={m1.rsquared:.3f}  AIC={m1.aic:.1f}  n={int(m1.nobs)}")
slope = m1.params["lnco2"]
tcr_like = slope * np.log(2)
ci_lo, ci_hi = m1.conf_int().loc["lnco2"].values * np.log(2)
print(f"\nImplied warming per CO2 doubling (slope * ln 2): {tcr_like:.2f} C "
      f"(OLS 95% CI {ci_lo:.2f} to {ci_hi:.2f})")
note(f"""The slope is the warming per e-fold of CO2; times ln(2) it is warming per doubling:
{tcr_like:.1f} C here. The IPCC transient climate response is ~1.8 C (likely 1.4-2.2), so this
single-predictor fit overshoots. Why: methane, N2O, halocarbons and declining aerosol
cooling all rose alongside CO2, and the coefficient absorbs every one of them (omitted
variable bias). The regression is a useful descriptive fit; it is not an attribution
study. Physics plus a model of all forcings is what makes the IPCC number defensible.""")

# %% ----------------------------------------------------------- diagnostics
banner("Residual diagnostics")
resid = m1.resid
dw = durbin_watson(resid)
bp = het_breuschpagan(resid, m1.model.exog)
lb = acorr_ljungbox(resid, lags=[1, 5], return_df=True)
sw = sps.shapiro(resid)
n_eff, r1 = ar1_effective_n(resid)
print(f"Durbin-Watson={dw:.2f} (2 = none; {dw:.2f} => positive autocorrelation)")
print(f"residual lag-1 r={r1:.2f}  effective n={n_eff:.0f} of {len(resid)}")
print(f"Ljung-Box p (lag1, lag5): {lb.lb_pvalue.round(4).tolist()}")
print(f"Breusch-Pagan p={bp[1]:.3f} (heteroscedasticity)  Shapiro p={sw.pvalue:.3f}")

fig, axes = plt.subplots(2, 2, figsize=(11, 7))
axes[0, 0].scatter(df.lnco2, df["T"], s=12, alpha=.7)
xs = np.linspace(df.lnco2.min(), df.lnco2.max(), 50)
axes[0, 0].plot(xs, m1.params.iloc[0] + slope * xs, color="C3")
axes[0, 0].set_xlabel("ln(CO2/280)")
axes[0, 0].set_ylabel("global anomaly (C)")
axes[0, 0].set_title(f"Fit: slope={slope:.2f} C per e-fold, R2={m1.rsquared:.2f}")
axes[0, 1].scatter(m1.fittedvalues, resid, s=12)
axes[0, 1].axhline(0, color="k", lw=.8)
axes[0, 1].set_xlabel("fitted")
axes[0, 1].set_ylabel("residual")
axes[0, 1].set_title("Residuals vs fitted (look for curvature / fanning)")
axes[1, 0].plot(df.index, resid, lw=1)
axes[1, 0].axhline(0, color="k", lw=.8)
axes[1, 0].set_title(f"Residuals vs time: runs => autocorrelation (DW={dw:.2f})")
sm.qqplot(resid, line="s", ax=axes[1, 1])
axes[1, 1].set_title("QQ plot of residuals")
fig.tight_layout()
savefig(fig, "ols_diagnostics", SEC)

# %% ----------------------------------------------------------- fixing the SEs
banner("Standard errors under autocorrelation: OLS vs HAC vs GLS-AR(1)")
hac = smf.ols("T ~ lnco2", data=df).fit(cov_type="HAC", cov_kwds={"maxlags": 5})
gls = sm.GLSAR(df["T"], sm.add_constant(df.lnco2), rho=1).iterative_fit(maxiter=20)
rows = pd.DataFrame({
    "OLS": [m1.params["lnco2"], m1.bse["lnco2"]],
    "OLS + n_eff adj": [m1.params["lnco2"], m1.bse["lnco2"] * np.sqrt(len(resid) / n_eff)],
    "HAC (Newey-West, 5 lags)": [hac.params["lnco2"], hac.bse["lnco2"]],
    "GLSAR(1) Cochrane-Orcutt": [gls.params["lnco2"], gls.bse["lnco2"]],
}, index=["slope", "SE"]).T
rows["t"] = rows.slope / rows.SE
print(rows.round(4).to_string())
print(f"GLSAR estimated rho = {gls.model.rho[0]:.2f}")
note("""Autocorrelation does not bias the slope but the naive OLS SE is too small by roughly
sqrt(n / n_eff). HAC (Newey-West) and GLS with AR(1) errors both widen the interval. The
conclusion survives here because the signal is enormous; for marginal trends it often
does not. Interviewers like to hear 'coefficients unbiased, SEs wrong' in one breath.""")

# %% ----------------------------------------------------------- model comparison
banner("Model comparison: ln CO2 vs linear time vs both")
m_time = smf.ols("T ~ year", data=df).fit()
m_both = smf.ols("T ~ lnco2 + year", data=df).fit()
m_quad = smf.ols("T ~ lnco2 + I(lnco2**2)", data=df).fit()
cmp = pd.DataFrame({m: [f.rsquared_adj, f.aic, f.bic] for m, f in
                    [("T~lnco2", m1), ("T~year", m_time), ("T~lnco2+year", m_both),
                     ("T~lnco2+lnco2^2", m_quad)]},
                   index=["adj R2", "AIC", "BIC"]).T
print(cmp.round(3).to_string())
X2 = sm.add_constant(df[["lnco2", "year"]])
vif2 = variance_inflation_factor(X2.values, 1)
print(f"\nIn T~lnco2+year: coef(lnco2)={m_both.params['lnco2']:.2f} (was {m1.params['lnco2']:.2f} alone), "
      f"p(lnco2)={m_both.pvalues['lnco2']:.3f}, p(year)={m_both.pvalues['year']:.3f}, "
      f"corr(lnco2, year)={df.lnco2.corr(df.year):.3f}, VIF={vif2:.0f}")
note(f"""ln(CO2) beats a linear time trend (lower AIC/BIC) because CO2 growth accelerated and so
did warming. Adding 'year' still improves AIC, but with r=0.995 between the predictors the
VIF is {vif2:.0f} and the CO2 coefficient swings from {m1.params['lnco2']:.1f} to
{m_both.params['lnco2']:.1f}: multicollinearity makes individual coefficients unstable even
when the joint fit is excellent. Observational data alone cannot separate 'CO2' from
'anything else that rose smoothly'.""")

# %% ----------------------------------------------------------- (b) multiple regression
banner("(b) Multiple regression: station mean temperature ~ lat + elev + lon (US airports)")
st = data.station_normals()
st["elev_km"] = st.elev / 1000
m2 = smf.ols("tavg ~ lat + elev_km + lon", data=st).fit()
print(m2.summary().tables[1])
print(f"R^2={m2.rsquared:.3f}  n={int(m2.nobs)}  RMSE={np.sqrt(m2.mse_resid):.2f} C")
X = sm.add_constant(st[["lat", "elev_km", "lon"]])
vif = pd.Series([variance_inflation_factor(X.values, i) for i in range(1, X.shape[1])],
                index=X.columns[1:])
print("\nVIF:\n" + vif.round(2).to_string())
std_coef = m2.params[1:] * st[["lat", "elev_km", "lon"]].std() / st.tavg.std()
print("\nStandardised coefficients (SD of tavg per SD of x):\n" + std_coef.round(3).to_string())
note(f"""Elevation coefficient {m2.params['elev_km']:.2f} C/km is a surface lapse rate. It is
shallower than the free-air 6.5 C/km because high stations sit in the dry, sunny
interior West where surface heating is strong and because longitude (r with elevation is
nonzero) soaks up part of the effect; drop lon and the lapse rate steepens. The latitude
coefficient {m2.params['lat']:.2f} C/deg is the meridional gradient and dominates the
standardised coefficients. VIFs are small so the coefficients are interpretable. Residual
RMSE {np.sqrt(m2.mse_resid):.1f} C is coasts, aspect, urban heat; those residuals are spatially
autocorrelated, which section 09 tests.""")

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
for ax, var in zip(axes, ["lat", "elev_km", "lon"]):
    others = [v for v in ["lat", "elev_km", "lon"] if v != var]
    # partial residual (component + residual) plot
    pr = m2.resid + m2.params[var] * st[var]
    ax.scatter(st[var], pr, s=10, alpha=.6)
    xs = np.linspace(st[var].min(), st[var].max(), 20)
    ax.plot(xs, m2.params[var] * xs, color="C3")
    ax.set_xlabel(var)
    ax.set_ylabel("partial residual")
    ax.set_title(f"{var}: {m2.params[var]:.2f} C per unit")
fig.suptitle("Component + residual plots (check each predictor's linearity)", y=1.02)
fig.tight_layout()
savefig(fig, "multiple_regression_partials", SEC)

# %% ----------------------------------------------------------- transforms & GLM teaser
banner("Precipitation ~ predictors: when y is positive and skewed")
m3 = smf.ols("prcp ~ lat + elev_km + lon", data=st).fit()
m3log = smf.ols("np.log(prcp) ~ lat + elev_km + lon", data=st).fit()
m3gam = smf.glm("prcp ~ lat + elev_km + lon", data=st,
                family=sm.families.Gamma(link=sm.families.links.Log())).fit()
print(f"OLS on prcp:      R2={m3.rsquared:.3f}  Breusch-Pagan p="
      f"{het_breuschpagan(m3.resid, m3.model.exog)[1]:.1e}")
print(f"OLS on log(prcp): R2={m3log.rsquared:.3f}  Breusch-Pagan p="
      f"{het_breuschpagan(m3log.resid, m3log.model.exog)[1]:.3f}")
print("coefficients (multiplicative scale):\n" + pd.DataFrame(
    {"log-OLS": m3log.params, "Gamma GLM log-link": m3gam.params}).round(4).to_string())
note("""Annual precipitation is strictly positive and heteroscedastic (wet places vary more).
Log-transforming y or using a Gamma GLM with a log link gives multiplicative effects
(coefficient = fractional change per unit of x). The two sets of coefficients agree closely.
Breusch-Pagan still rejects after the log because the variance also depends on longitude
(a wet East / dry West split that lat+lon+elev cannot fully describe); the GLM handles
mean-variance coupling by construction and predicts on the natural scale without the
back-transform bias of exp(E[log y]).""")
