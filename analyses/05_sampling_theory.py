# %% [markdown]
# # 05 Sampling theory: Gy's fundamental sampling error and the nugget effect
#
# **Cheat sheet**
# - Every assay is a tiny sub-sample of a lot. Gy's theory decomposes sampling error into:
#   fundamental sampling error (FSE, from particle heterogeneity; irreducible except by mass or
#   comminution), grouping & segregation error, delimitation/extraction errors (bias from bad
#   equipment), preparation errors, and analytical error.
# - FSE relative variance s^2 = c f g l d^3 / M. Doubling the top size d costs 8x the mass.
#   c = mineralogical constant ~ rho_mineral / a for low grade a: gold at 1 g/t has c ~ 1e7 g/cm3,
#   copper at 1% has c ~ 1e2. That is why gold samples are kilograms and copper samples grams.
# - A sampling protocol is a chain of (crush -> split) stages; the FSE variances add. The
#   nomogram (log mass vs log variance, slope -1 lines per top size) shows which stage breaks
#   the budget; fix that stage, not the whole protocol.
# - Variogram nugget = geological micro-variability + sampling + preparation + assay error.
#   Duplicate programmes let you say how much of the nugget is measurement error and hence how
#   much better a bigger sample or finer crush would make the resource model.
# - For the domain mean, n = (z CV / E)^2 samples for relative error E; CV 1 and 10% error needs
#   ~400 independent samples. Clustered samples count for less (section 01).

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from miningstats import data
from miningstats.plotting import banner, note, savefig, setup
from miningstats.stats import (
    cell_declustering,
    experimental_variogram,
    fit_variogram,
    gy_fse_variance,
    gy_min_mass,
    mineralogical_constant,
    spherical,
)

setup()
SEC = "05_sampling_theory"

# %% ----------------------------------------------------------- mineralogical constant
banner("Mineralogical constant: why gold is hard")
cases = {
    "gold 1.5 g/t (native Au, 50 um liberation, alpha 1.5)": dict(a=1.5e-6, rho_m=19.3, rho_g=2.7, dl=0.005, alpha=1.5),
    "gold 1.5 g/t (classic Gy alpha 2.5)": dict(a=1.5e-6, rho_m=19.3, rho_g=2.7, dl=0.005, alpha=2.5),
    "copper 1% as chalcopyrite (100 um)": dict(a=0.01 / 0.346, rho_m=4.2, rho_g=2.7, dl=0.01, alpha=2.5),
    "iron ore 60% Fe as hematite (500 um)": dict(a=0.60 / 0.70, rho_m=5.3, rho_g=2.7, dl=0.05, alpha=2.5),
}
rows = []
for name, k in cases.items():
    c = mineralogical_constant(k["a"], k["rho_m"], k["rho_g"])
    rows.append(dict(case=name, c_g_cm3=c,
                     mass_10pct_FSE_1cm_g=gy_min_mass(1.0, 0.10, c, liberation_size_cm=k["dl"], alpha=k["alpha"]),
                     mass_10pct_FSE_2mm_g=gy_min_mass(0.2, 0.10, c, liberation_size_cm=k["dl"], alpha=k["alpha"]),
                     mass_10pct_FSE_75um_g=gy_min_mass(0.0075, 0.10, c, liberation_size_cm=k["dl"], alpha=k["alpha"])))
tab = pd.DataFrame(rows).set_index("case")
with pd.option_context("display.float_format", "{:,.3g}".format):
    print(tab.to_string())
note("""The mineralogical constant scales as 1/grade: five orders of magnitude separate gold
from iron ore. At 1 cm top size a representative gold sample weighs many kilograms to
tonnes; only after pulverising to 75 microns does a sub-sample of tens of grams carry a
10% fundamental error. The liberation factor is the least knowable input: Gy's classic
sqrt(d_l/d) (alpha 2.5) is far more pessimistic than the Francois-Bongarcon alpha 1.5
calibration that practitioners fit to gold duplicate data. Protocols are calibrated against
duplicates, not derived from theory alone.""")

# %% ----------------------------------------------------------- protocol nomogram
banner("A sample preparation protocol on the Gy nomogram (gold, 1.5 g/t)")
c_au = mineralogical_constant(1.5e-6, 19.3, 2.7)
GY = dict(liberation_size_cm=0.005, alpha=1.5)
# stage: (top size cm, mass g) after each step; a split keeps top size, a crush keeps mass
protocols = {
    "A: split coarse (common shortcut)": [("core as received", 2.5, 3000), ("crush to 1 cm", 1.0, 3000),
                                          ("split to 1 kg", 1.0, 1000), ("crush to 2 mm", 0.2, 1000),
                                          ("split to 250 g", 0.2, 250), ("pulverise to 75 um", 0.0075, 250),
                                          ("30 g fire assay", 0.0075, 30)],
    "B: crush fine before splitting": [("core as received", 2.5, 3000), ("crush to 2 mm", 0.2, 3000),
                                       ("split to 1 kg", 0.2, 1000), ("pulverise to 75 um", 0.0075, 1000),
                                       ("50 g fire assay", 0.0075, 50)],
}
paths = {}
for pname, protocol in protocols.items():
    rows, total_var = [], 0.0
    for i, (stage, d, m) in enumerate(protocol):
        rsd = gy_fse_variance(d, m, c_au, **GY)
        is_split = i > 0 and protocol[i - 1][1] == d
        if is_split:  # error is incurred when mass is reduced at a given top size
            total_var += rsd ** 2
        rows.append(dict(stage=stage, top_size_cm=d, mass_g=m, FSE_rsd_pct=rsd * 100,
                         counted=("split" if is_split else "crush/receive")))
    prot = pd.DataFrame(rows).set_index("stage")
    paths[pname] = (prot, np.sqrt(total_var))
    print(f"\n{pname}\n" + prot.round(2).to_string())
    print(f"total FSE of protocol (splits add in variance): {np.sqrt(total_var) * 100:.1f}% relative")
worst = paths["A: split coarse (common shortcut)"][0].query("counted == 'split'").FSE_rsd_pct.idxmax()
note(f"""On the nomogram each split is a vertical drop in mass at fixed top size (error grows) and
each crush a horizontal move to a finer top size (error shrinks). Protocol A's error
({paths['A: split coarse (common shortcut)'][1] * 100:.0f}%) is dominated by the '{worst}' stage: splitting
while the gold is still locked in centimetre particles. Protocol B crushes everything to 2 mm
before the first split and pulverises the whole kilogram, cutting the total to
{paths['B: crush fine before splitting'][1] * 100:.0f}%. Changing the final aliquot from 30 to 50 g is nearly
irrelevant. Match the top size to the mass you are about to take.""")

fig, ax = plt.subplots(figsize=(9, 6))
masses = np.logspace(0, 4, 100)
for d, lab in [(2.5, "2.5 cm"), (1.0, "1 cm"), (0.2, "2 mm"), (0.0075, "75 um")]:
    ax.loglog(masses, gy_fse_variance(d, masses, c_au, **GY) ** 2, lw=1, label=f"top size {lab}")
for (pname, (prot, tot)), col in zip(paths.items(), ["k", "C3"]):
    ax.loglog(prot.mass_g, (prot.FSE_rsd_pct / 100) ** 2, "-o", color=col, lw=2, ms=6,
              label=f"{pname}: total {tot * 100:.0f}%")
for stage, r in paths["A: split coarse (common shortcut)"][0].iterrows():
    ax.annotate(stage, (r.mass_g, (r.FSE_rsd_pct / 100) ** 2), fontsize=7, xytext=(4, 4), textcoords="offset points")
ax.axhline(0.1 ** 2, color="C1", ls="--", label="10% RSD safety line")
ax.set_xlabel("sample mass (g)")
ax.set_ylabel("relative variance of FSE")
ax.set_title("Gy sampling nomogram, gold at 1.5 g/t (liberation 50 um, alpha 1.5)")
ax.legend(fontsize=8)
savefig(fig, "nomogram", SEC)

# %% ----------------------------------------------------------- nugget decomposition
banner("Decomposing the variogram nugget (Walker Lake V + the synthetic duplicate precision)")
w = data.walker()
V = w.V.values
h, g, n = experimental_variogram(w.x, w.y, V, 5.0, 16)
p = fit_variogram(h, g, spherical, counts=n, nugget_max=0.6 * V.var())
qa = data.assay_qaqc()
dup = qa["duplicates"]
# precision as relative sd from duplicates (pairs): sd_single = sd(diff)/sqrt(2); use relative form
rel = (dup.original - dup.duplicate) / (0.5 * (dup.original + dup.duplicate))
rel_sd_single = rel.std(ddof=1) / np.sqrt(2)
meas_var = (rel_sd_single * V.mean()) ** 2
print(f"variogram nugget={p[0]:.0f}, sill={p[1]:.0f} -> relative nugget {p[0] / p[1]:.0%}")
print(f"duplicate relative sd of a single assay={rel_sd_single:.3f}; at the mean grade {V.mean():.0f} that is a "
      f"measurement variance of {meas_var:.0f}")
print(f"measurement error explains {meas_var / p[0]:.0%} of the nugget; the rest ({1 - meas_var / p[0]:.0%}) is "
      f"geological variability at distances shorter than the closest sample pairs")
note("""(The duplicate precision here is the synthetic gold programme's, applied to Walker Lake
for illustration.) The decomposition tells management what a better protocol can buy: if
measurement error is a minor part of the nugget, spending on assay precision will not
sharpen the resource model, but closer-spaced drilling might. Reported precision from
duplicates must be of the same sample type as the variogram data (field vs pulp duplicates
capture different error stages).""")

# %% ----------------------------------------------------------- how many samples
banner("How many samples to estimate the domain mean to +-10%?")
cv = V.std(ddof=1) / V.mean()
n_needed = (1.96 * cv / 0.10) ** 2
wts = cell_declustering(w.x.values, w.y.values, 20)
n_eff = wts.sum() ** 2 / np.sum(wts ** 2)  # Kish effective sample size under unequal weights
print(f"CV={cv:.2f} -> n={n_needed:.0f} independent samples for +-10% at 95%; we have {len(V)} samples but an "
      f"effective n of {n_eff:.0f} after declustering weights (Kish), and spatial correlation lowers it further")
for cv_ in [0.5, 1.0, 2.0, 3.0]:
    print(f"  CV={cv_:.1f}: n={(1.96 * cv_ / 0.10) ** 2:.0f}")
note("""Nuggety gold with CV 2-3 needs thousands of independent samples for a +-10% domain
mean; no drilling programme delivers that, which is why gold resources carry wide
classification bands and why grade control relies on dense blast-hole sampling. Clustering
and spatial correlation reduce the effective sample size below the row count.""")
