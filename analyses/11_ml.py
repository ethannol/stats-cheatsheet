# %% [markdown]
# # 11 Machine learning on environmental data: rainfall-runoff
#
# Predict daily Potomac discharge from antecedent precipitation and temperature at Washington
# National (DCA). Then classify high-flow days. Emphasis on evaluation that respects time.
#
# **Cheat sheet**
# - Split by time. Random K-fold on autocorrelated data leaks the future into training and
#   produces fantasy scores. Use TimeSeriesSplit / blocked CV with a gap.
# - Baselines first: persistence (yesterday's value), climatology (day-of-year mean). A model
#   that does not beat them is not a model.
# - Feature engineering is where domain knowledge enters: antecedent precipitation index,
#   rolling sums at several time scales, season encodings (sin/cos), snowmelt proxies.
# - Linear + regularisation: Ridge (L2, shrinks correlated coefficients together), Lasso (L1,
#   sparse selection), ElasticNet. Tune alpha by CV; standardise features first.
# - Trees / boosting: nonlinear, interaction-aware, no scaling needed, but cannot extrapolate
#   beyond the training range (a problem for record floods and a warming climate).
# - Importance: impurity importance is biased toward high-cardinality features; permutation
#   importance on the test set is more honest; partial dependence shows the shape.
# - Hydrology metrics: NSE (= R^2 vs mean), KGE, log-space RMSE (floods vs low flows).
# - Classification with rare events: accuracy is useless; use PR-AUC, recall at fixed
#   precision, calibration. Choose the threshold from the cost of misses vs false alarms.

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.inspection import PartialDependenceDisplay, permutation_importance
from sklearn.linear_model import LassoCV, LinearRegression, LogisticRegression, RidgeCV
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    mean_squared_error,
    precision_recall_curve,
    r2_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import KFold, TimeSeriesSplit, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from earthstats import data
from earthstats.plotting import banner, note, savefig, setup

setup()
SEC = "11_ml"

# %% ----------------------------------------------------------- features
q = data.usgs_daily().flow_cfs
w = data.ghcn_daily(data.REAGAN_NATIONAL)
df = pd.DataFrame({"q": q}).join(w[["PRCP", "TMAX", "TMIN", "SNOW"]]).loc["1950":"2025"]
df["PRCP"] = df.PRCP.fillna(0)
df[["TMAX", "TMIN"]] = df[["TMAX", "TMIN"]].interpolate()
df["SNOW"] = df.SNOW.fillna(0)
feat = pd.DataFrame(index=df.index)
for win in [1, 3, 7, 14, 30, 60, 120]:
    feat[f"p{win}"] = df.PRCP.rolling(win).sum().shift(1)  # strictly antecedent
feat["api"] = df.PRCP.shift(1).ewm(halflife=10).mean()  # antecedent precipitation index
feat["tmax30"] = df.TMAX.rolling(30).mean().shift(1)
feat["snow30"] = df.SNOW.rolling(30).sum().shift(1)
doy = df.index.dayofyear
feat["sin_doy"], feat["cos_doy"] = np.sin(2 * np.pi * doy / 365.25), np.cos(2 * np.pi * doy / 365.25)
feat["logq_lag1"] = np.log10(df.q).shift(1)
y = np.log10(df.q)
ok = feat.notna().all(axis=1) & y.notna()
feat, y = feat[ok], y[ok]
train = feat.index.year <= 2010
test = ~train
no_persist = [c for c in feat.columns if c != "logq_lag1"]
print(f"{ok.sum()} days; train through 2010 ({train.sum()}), test 2011-2025 ({test.sum()}); "
      f"{len(no_persist)} precipitation/season features (+1 persistence feature)")


def nse(obs, sim):
    return 1 - np.sum((obs - sim) ** 2) / np.sum((obs - obs.mean()) ** 2)


def kge(obs, sim):
    r = np.corrcoef(obs, sim)[0, 1]
    alpha = sim.std() / obs.std()
    beta = sim.mean() / obs.mean()
    return 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)


def report(name, yhat_log):
    obs, sim = 10 ** y[test].values, 10 ** yhat_log
    return dict(model=name, rmse_log=np.sqrt(mean_squared_error(y[test], yhat_log)),
                r2_log=r2_score(y[test], yhat_log), NSE_flow=nse(obs, sim), KGE_flow=kge(obs, sim))


# %% ----------------------------------------------------------- baselines and models
banner("Regression: log10(discharge) from antecedent weather. Test period 2011-2025")
rows = []
clim = y[train].groupby(doy[ok][train]).mean()
rows.append(report("climatology (day-of-year mean)", clim.reindex(doy[ok][test]).values))
rows.append(report("persistence (yesterday's flow)", feat.logq_lag1[test].values))

Xtr, Xte = feat.loc[train, no_persist], feat.loc[test, no_persist]
models = {
    "linear regression": make_pipeline(StandardScaler(), LinearRegression()),
    "ridge (CV alpha)": make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 3, 25))),
    "lasso (CV alpha)": make_pipeline(StandardScaler(), LassoCV(cv=TimeSeriesSplit(5), random_state=0)),
    "random forest": RandomForestRegressor(300, min_samples_leaf=20, n_jobs=-1, random_state=0),
    "gradient boosting": HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=0),
}
fitted = {}
for name, m in models.items():
    m.fit(Xtr, y[train])
    fitted[name] = m
    rows.append(report(f"{name} [weather only]", m.predict(Xte)))
gb_p = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=0).fit(
    feat.loc[train], y[train])
rows.append(report("gradient boosting [+ persistence]", gb_p.predict(feat.loc[test])))
res = pd.DataFrame(rows).set_index("model")
print(res.round(3).to_string())
lasso = fitted["lasso (CV alpha)"][-1]
print(f"\nlasso kept {np.sum(lasso.coef_ != 0)} of {len(no_persist)} features; dropped: "
      f"{[f for f, c in zip(no_persist, lasso.coef_) if c == 0]}")
wo = res[res.index.str.contains("weather only")]
note(f"""Read the table top-down. Climatology explains the seasonal cycle; persistence is nearly
unbeatable one day ahead because rivers are integrators. Weather-only models answer a
different and more interesting question: how much of the flow is explained by the last few
months of weather at one gauge? All five land at R2 {wo.r2_log.min():.2f}-{wo.r2_log.max():.2f}
in log space: when a linear model ties gradient boosting, the ceiling is the information in
the features (one rain gauge for an 11,500 sq mi basin), not model flexibility. Adding
persistence to the boosted model gives the best one-day forecast. NSE rewards matching big
floods; KGE also penalises bias and variability errors.""")

# %% ----------------------------------------------------------- leakage demo
banner("Cross-validation that leaks: shuffled K-fold vs TimeSeriesSplit (random forest, train period)")
rf = RandomForestRegressor(200, min_samples_leaf=20, n_jobs=-1, random_state=0)
kf = cross_val_score(rf, Xtr, y[train], cv=KFold(5, shuffle=True, random_state=0), scoring="r2")
ts = cross_val_score(rf, Xtr, y[train], cv=TimeSeriesSplit(5, gap=60), scoring="r2")
print(f"shuffled KFold R2:    {kf.mean():.3f} +- {kf.std():.3f}")
print(f"TimeSeriesSplit R2:   {ts.mean():.3f} +- {ts.std():.3f}  (60-day gap between train and test)")
print(f"true holdout 2011-25: {res.loc['random forest [weather only]', 'r2_log']:.3f}")
note("""Shuffled K-fold puts Tuesday in the training set and Wednesday in the test set; with
30-120 day rolling features they are almost the same row, so the score is inflated. The
time-aware split with a gap matches the real holdout. This is the single most common
evaluation mistake with environmental data.""")

# %% ----------------------------------------------------------- interpretation
banner("Interpretation: permutation importance (test set) and partial dependence")
gb = fitted["gradient boosting"]
pi = permutation_importance(gb, Xte, y[test], n_repeats=10, random_state=0, n_jobs=-1)
imp = pd.Series(pi.importances_mean, index=no_persist).sort_values(ascending=False)
print(imp.round(4).to_string())
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
imp.iloc[::-1].plot.barh(ax=axes[0], xerr=pd.Series(pi.importances_std, index=no_persist)[imp.index[::-1]])
axes[0].set_title("Permutation importance (drop in R2 on test set)")
PartialDependenceDisplay.from_estimator(gb, Xte, ["p30", "p120", "tmax30"], ax=axes[1], kind="average")
axes[1].set_title("Partial dependence: 30-day, 120-day precipitation; 30-day TMAX")
fig.tight_layout()
savefig(fig, "importance_pdp", SEC)
top = ", ".join(imp.index[:4])
note(f"""The top features are {top}: the seasonal cycle of evapotranspiration, the exponentially
weighted antecedent precipitation index, recent temperature (again evapotranspiration) and
the 120-day rainfall total. Short windows (1-7 days) matter little because a single DC
gauge is a poor proxy for a storm over an 11,500 sq mi basin, while months of accumulated
rainfall control baseflow. The partial dependence plots show the shape of each effect
(flow rises with p120, falls with tmax30). Permutation importance on held-out data avoids
the impurity-importance bias toward continuous, high-cardinality features; correlated
features (p60/p120/api) share importance, so read groups, not individual bars.""")

# predictions vs observed
fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
seg = slice("2018-01-01", "2019-12-31")
axes[0].plot(10 ** y[seg] / 1000, color="k", lw=1, label="observed")
axes[0].plot(10 ** pd.Series(gb.predict(feat.loc[seg, no_persist]), index=y[seg].index) / 1000,
             lw=1, label="gradient boosting [weather only]")
axes[0].set_yscale("log")
axes[0].set_ylabel("discharge (kcfs)")
axes[0].set_title("2018-2019: timing right, peaks under-predicted (tree models cannot extrapolate)")
axes[0].legend()
axes[1].scatter(y[test], gb.predict(Xte), s=3, alpha=.3)
lims = [y.min(), y.max()]
axes[1].plot(lims, lims, "k--", lw=.8)
axes[1].set_xlabel("observed log10 Q")
axes[1].set_ylabel("predicted log10 Q")
axes[1].set_title(f"Test set: NSE(flow)={res.loc['gradient boosting [weather only]', 'NSE_flow']:.2f}")
savefig(fig, "predictions", SEC)

# %% ----------------------------------------------------------- classification of high-flow days
banner("Classification: will tomorrow be a high-flow day (> 95th percentile)?")
thr = np.quantile(10 ** y[train], .95)
lab = (10 ** y > thr).astype(int)
# predict tomorrow's exceedance from today's features including persistence
Xc = feat.copy()
yc = lab.shift(-1).dropna().astype(int)
Xc = Xc.loc[yc.index]
trc, tec = Xc.index.year <= 2010, Xc.index.year > 2010
print(f"high-flow threshold {thr / 1000:.0f} kcfs; positive rate train={yc[trc].mean():.3f} test={yc[tec].mean():.3f}")
clfs = {"logistic regression": make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
        "gradient boosting": HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=0)}
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
for name, clf in clfs.items():
    clf.fit(Xc[trc], yc[trc])
    p = clf.predict_proba(Xc[tec])[:, 1]
    auc, ap, brier = roc_auc_score(yc[tec], p), average_precision_score(yc[tec], p), brier_score_loss(yc[tec], p)
    acc_naive = 1 - yc[tec].mean()
    print(f"  {name:20s} ROC-AUC={auc:.3f}  PR-AUC={ap:.3f}  Brier={brier:.4f}  "
          f"(accuracy of 'never' = {acc_naive:.3f})")
    fpr, tpr, _ = roc_curve(yc[tec], p)
    prec, rec, th = precision_recall_curve(yc[tec], p)
    axes[0].plot(fpr, tpr, label=f"{name} AUC={auc:.2f}")
    axes[1].plot(rec, prec, label=f"{name} AP={ap:.2f}")
    frac_pos, mean_pred = calibration_curve(yc[tec], p, n_bins=10, strategy="quantile")
    axes[2].plot(mean_pred, frac_pos, "o-", label=name)
    if name == "gradient boosting":
        # threshold choice: maximise F1, then show the trade-off
        f1 = 2 * prec[:-1] * rec[:-1] / (prec[:-1] + rec[:-1] + 1e-12)
        k = np.argmax(f1)
        print(f"    best-F1 threshold={th[k]:.2f}: precision={prec[k]:.2f} recall={rec[k]:.2f}; "
              f"at threshold 0.5: precision={prec[np.searchsorted(th, 0.5)]:.2f} "
              f"recall={rec[np.searchsorted(th, 0.5)]:.2f}")
axes[0].plot([0, 1], [0, 1], "k--", lw=.8)
axes[0].set_title("ROC")
axes[0].set_xlabel("false positive rate")
axes[0].set_ylabel("true positive rate")
axes[1].axhline(yc[tec].mean(), color="k", ls="--", lw=.8, label="base rate")
axes[1].set_title("Precision-recall (the honest curve for rare events)")
axes[1].set_xlabel("recall")
axes[1].set_ylabel("precision")
axes[2].plot([0, 1], [0, 1], "k--", lw=.8)
axes[2].set_title("Calibration (reliability diagram)")
axes[2].set_xlabel("predicted probability")
axes[2].set_ylabel("observed frequency")
for ax in axes:
    ax.legend(fontsize=8)
fig.tight_layout()
savefig(fig, "classification", SEC)
note("""ROC-AUC looks superb because the negatives are easy (dry August days). PR-AUC is the
number to quote for a 5% event; compare it with the base rate. Calibration tells you
whether 'p = 0.7' means 70%; it matters more than accuracy if the probability feeds a
decision. The threshold is a policy choice (cost of a missed flood vs a false alarm),
not a statistic.""")
