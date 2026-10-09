"""Dataset loaders. Each loader downloads once, caches under data/, returns a tidy pandas object.

Real data
---------
walker(), walker_exhaustive()  Walker Lake (Isaaks & Srivastava 1989): 470 samples of V, U, T and
                               the 78,000-point exhaustive grid. Built from a real Nevada DEM; used
                               for 35 years as the standard mining-geostatistics teaching deposit.
jura()                         Swiss Jura soil geochemistry (Goovaerts 1997): 359 samples, 7 metals,
                               rock type and land use. Compositional / multivariate geostatistics.
gslib_cluster()                GSLIB cluster.dat: 140 clustered samples with reference declustering
                               weights.
msha_mines(), msha_production(), msha_accidents()
                               US Mine Safety and Health Administration open data: every US mine,
                               annual hours/employees/coal production per mine, and every reportable
                               accident/injury since 2000.
clark_fork_wq(), clark_fork_flow()
                               Clark Fork River at Deer Lodge, Montana (USGS 12324200), downstream
                               of the Butte/Anaconda copper district: metals and pH from the Water
                               Quality Portal (with non-detects), daily discharge from NWIS.
commodity_prices()             World Bank "Pink Sheet" monthly metal/energy prices, 1960-2024.

Synthetic data (clearly labelled; generated with fixed seeds from stated parameters)
---------------
assay_qaqc()                   field duplicates, certified reference materials, blanks with planted
                               problems. No public assay QA/QC database exists.
equipment_failures()           haul-truck component failures with censoring, repair times, daily
                               production. Mine maintenance records are never public.
drillhole_from_walker()        a 1-D 'drillhole' with irregular sample intervals built from a
                               Walker Lake exhaustive-grid column, for compositing.

Run `python -m miningstats.data` to fetch everything up front.
"""

from __future__ import annotations

import io
import json
import sys
import zipfile
from datetime import date

import numpy as np
import pandas as pd
import requests

from . import DATA_DIR

RAW = DATA_DIR / "raw"
UA = {"User-Agent": "miningstats-cheatsheet/0.1 (educational)"}
GSTAT = "https://raw.githubusercontent.com/r-spatial/gstat/master/"
MSHA = "https://arlweb.msha.gov/OpenGovernmentData/DataSets/"
WQP = "https://www.waterqualitydata.us/data/Result/search"
CLARK_FORK = "12324200"
END_DATE = "2026-09-30"


def _get(url, params=None, timeout=900) -> bytes:
    r = requests.get(url, params=params, headers=UA, timeout=timeout)
    r.raise_for_status()
    return r.content


def _cached_bytes(name, url, params=None) -> bytes:
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / name
    if not path.exists():
        print(f"  downloading {name} ...", file=sys.stderr)
        path.write_bytes(_get(url, params))
    return path.read_bytes()


def _parquet_cached(name, builder):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{name}.parquet"
    if path.exists():
        return pd.read_parquet(path)
    df = builder()
    df.to_parquet(path)
    return df


def _rda(name):
    import rdata

    _cached_bytes(name, GSTAT + "data/" + name)
    return rdata.read_rda(RAW / name)


# ------------------------------------------------------------------ Walker Lake
def walker() -> pd.DataFrame:
    """470 Walker Lake samples: x, y (m-ish grid units), V (ppm-like grade, all samples),
    U (second variable, 275 samples), T (indicator 1/2 for geological 'type')."""

    def build():
        w = _rda("walker.rda")["walker"]
        xy = np.asarray(w.coords)
        df = w.data.reset_index(drop=True)
        df.insert(0, "y", xy[:, 1])
        df.insert(0, "x", xy[:, 0])
        df["T"] = df["T"].astype(int)
        return df.drop(columns=["Id"])

    return _parquet_cached("walker", build)


def walker_exhaustive() -> pd.DataFrame:
    """78,000-node exhaustive Walker Lake grid (260 x 300, unit spacing): x, y, U, V.
    sp stores grid cells x-fastest with y decreasing from the top row."""

    def build():
        e = _rda("walker.rda")["walker.exh"]
        g = vars(e.grid)
        off, cs, dim = (np.asarray(g["cellcentre.offset"]), np.asarray(g["cellsize"]),
                        np.asarray(g["cells.dim"]).astype(int))
        xs = off[0] + cs[0] * np.arange(dim[0])
        ys = off[1] + cs[1] * np.arange(dim[1])
        X, Y = np.meshgrid(xs, ys[::-1])  # top row first
        df = e.data.reset_index(drop=True)
        df.insert(0, "y", Y.ravel())
        df.insert(0, "x", X.ravel())
        return df

    return _parquet_cached("walker_exhaustive", build)


# ------------------------------------------------------------------------ Jura
def jura() -> pd.DataFrame:
    """359 Jura soil samples: x, y (km), landuse, rock, cd co cr cu ni pb zn (ppm), set."""

    def build():
        d = _rda("jura.rda")
        pred = d["jura.pred"].assign(set="prediction")
        val = d["jura.val"].assign(set="validation")
        df = pd.concat([pred, val], ignore_index=True)
        df.columns = [str(c).lower() for c in df.columns]
        df = df.rename(columns={"xloc": "x", "yloc": "y"})
        for c in ["landuse", "rock", "set"]:
            df[c] = df[c].astype(str)
        return df[["x", "y", "long", "lat", "landuse", "rock", "cd", "co", "cr", "cu", "ni", "pb",
                   "zn", "set"]]

    return _parquet_cached("jura", build)


# --------------------------------------------------------------- GSLIB cluster
def gslib_cluster() -> pd.DataFrame:
    """GSLIB cluster.dat: x, y, primary, secondary, declus_weight (reference weights)."""

    def build():
        raw = _cached_bytes("cluster.txt", GSTAT + "inst/external/cluster.txt")
        df = pd.read_csv(io.BytesIO(raw), sep=r"\s+")
        df.columns = ["x", "y", "primary", "secondary", "declus_weight"]
        return df

    return _parquet_cached("gslib_cluster", build)


# ------------------------------------------------------------------------ MSHA
def _msha_txt(zip_name: str, usecols=None) -> pd.DataFrame:
    raw = _cached_bytes(zip_name, MSHA + zip_name)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        member = z.namelist()[0]
        with z.open(member) as f:
            return pd.read_csv(f, sep="|", encoding="latin-1", usecols=usecols, dtype=str,
                               quotechar='"', na_values=["", "?", "NO VALUE FOUND"],
                               keep_default_na=True, low_memory=False)


def msha_mines() -> pd.DataFrame:
    """One row per US mine: id, name, coal/metal indicator, mine type, status, commodity (SIC),
    state, lat/lon, employees."""

    def build():
        cols = ["MINE_ID", "CURRENT_MINE_NAME", "COAL_METAL_IND", "CURRENT_MINE_TYPE",
                "CURRENT_MINE_STATUS", "CURRENT_STATUS_DT", "STATE", "PRIMARY_SIC",
                "PRIMARY_CANVASS",
                "LONGITUDE", "LATITUDE", "NO_EMPLOYEES", "CURRENT_CONTROLLER_NAME"]
        df = _msha_txt("Mines.zip", usecols=cols)
        df.columns = [c.lower() for c in df.columns]
        for c in ["longitude", "latitude", "no_employees"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df["current_status_dt"] = pd.to_datetime(df.current_status_dt, errors="coerce")
        return df

    return _parquet_cached("msha_mines", build)


def msha_production() -> pd.DataFrame:
    """Mine-year rows: mine_id, subunit, year, annual_hrs, annual_coal_prod (short tons, coal
    only), avg_annual_empl, coal_metal_ind."""

    def build():
        df = _msha_txt("MinesProdYearly.zip")
        df.columns = [c.lower() for c in df.columns]
        for c in ["calendar_yr", "annual_hrs", "annual_coal_prod", "avg_annual_empl",
                  "avg_employee_hours"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        df = df.rename(columns={"calendar_yr": "year", "c_m_ind": "coal_metal_ind"})
        return df[["mine_id", "curr_mine_nm", "state_abbr", "subunit_cd", "subunit_desc", "year",
                   "annual_hrs", "annual_coal_prod", "avg_annual_empl", "coal_metal_ind"]]

    return _parquet_cached("msha_production", build)


def msha_accidents() -> pd.DataFrame:
    """One row per reportable accident/injury/illness (~275k since 2000), selected columns."""

    def build():
        cols = ["MINE_ID", "SUBUNIT_CD", "SUBUNIT", "ACCIDENT_DT", "CAL_YR", "DEGREE_INJURY_CD",
                "DEGREE_INJURY", "MINING_EQUIP", "CLASSIFICATION", "ACCIDENT_TYPE", "NO_INJURIES",
                "TOT_EXPER", "MINE_EXPER", "JOB_EXPER", "OCCUPATION", "ACTIVITY", "INJURY_SOURCE",
                "NATURE_INJURY", "INJ_BODY_PART", "DAYS_RESTRICT", "DAYS_LOST", "COAL_METAL_IND"]
        df = _msha_txt("Accidents.zip", usecols=cols)
        df.columns = [c.lower() for c in df.columns]
        df["accident_dt"] = pd.to_datetime(df.accident_dt, format="%m/%d/%Y", errors="coerce")
        for c in ["cal_yr", "degree_injury_cd", "no_injuries", "tot_exper", "mine_exper",
                  "job_exper", "days_restrict", "days_lost"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        return df.rename(columns={"cal_yr": "year"})

    return _parquet_cached("msha_accidents", build)


# ------------------------------------------------------------------ Clark Fork
def clark_fork_wq() -> pd.DataFrame:
    """Discrete water-quality results: date, characteristic, fraction, value, unit, nondetect flag,
    detection limit. Metals in ug/L; pH in standard units."""

    def build():
        raw = _cached_bytes("wqp_clark_fork.csv", WQP, dict(
            siteid=f"USGS-{CLARK_FORK}", mimeType="csv", zip="no",
            characteristicName=["Copper", "Zinc", "Arsenic", "Cadmium", "Lead", "pH"]))
        df = pd.read_csv(io.BytesIO(raw), low_memory=False)
        keep = {"ActivityStartDate": "date", "CharacteristicName": "characteristic",
                "ResultSampleFractionText": "fraction", "ResultMeasureValue": "value",
                "ResultMeasure/MeasureUnitCode": "unit",
                "ResultDetectionConditionText": "detect_cond",
                "DetectionQuantitationLimitMeasure/MeasureValue": "dl",
                "DetectionQuantitationLimitMeasure/MeasureUnitCode": "dl_unit",
                "ActivityMediaName": "media"}
        df = df[list(keep)].rename(columns=keep)
        df = df[df.media == "Water"].drop(columns="media")
        df["date"] = pd.to_datetime(df.date)
        df["value"] = pd.to_numeric(df.value, errors="coerce")
        df["dl"] = pd.to_numeric(df.dl, errors="coerce")
        df["nondetect"] = df.detect_cond.notna()
        # for nondetects the reported value (and its unit) is empty; the detection limit is the
        # censoring level
        df.loc[df.nondetect & df.value.isna(), "value"] = df.dl
        df["unit"] = df.unit.fillna(df.dl_unit).str.lower().str.replace("ug/l", "ug/l")
        df = df.drop(columns="dl_unit")
        df["fraction"] = df.fraction.fillna("Total")
        return df.dropna(subset=["value"]).sort_values("date").reset_index(drop=True)

    return _parquet_cached("clark_fork_wq", build)


def clark_fork_flow() -> pd.DataFrame:
    """Daily mean discharge (cfs) at Clark Fork at Deer Lodge."""

    def build():
        raw = _cached_bytes("nwis_clark_fork.json", "https://waterservices.usgs.gov/nwis/dv/", dict(
            format="json", sites=CLARK_FORK, parameterCd="00060", statCd="00003",
            startDT="1969-01-01", endDT=END_DATE))
        ts = json.loads(raw)["value"]["timeSeries"][0]
        df = pd.DataFrame(ts["values"][0]["value"])
        df["date"] = pd.to_datetime(df.dateTime.str[:10])
        df["flow_cfs"] = pd.to_numeric(df["value"]).replace(-999999.0, np.nan)
        return df.set_index("date")[["flow_cfs"]]

    return _parquet_cached("clark_fork_flow", build)


# ------------------------------------------------------------ commodity prices
CMO_URL = ("https://thedocs.worldbank.org/en/doc/5d903e848db1d1b83e0ec8f744e55570-0350012021/"
           "related/CMO-Historical-Data-Monthly.xlsx")
METALS = {"gold": "Gold", "silver": "Silver", "platinum": "Platinum", "copper": "Copper",
          "aluminum": "Aluminum", "nickel": "Nickel", "zinc": "Zinc", "lead": "Lead", "tin": "Tin",
          "iron_ore": "Iron ore", "coal_aus": "Coal, Australian", "crude_oil": "Crude oil, average"}


def commodity_prices() -> pd.DataFrame:
    """Monthly nominal USD prices (World Bank Pink Sheet) for the metals in METALS plus coal and
    oil. Units: gold/silver/platinum $/troy oz, base metals $/mt, iron ore $/dmt, coal $/mt,
    oil $/bbl."""

    def build():
        raw = _cached_bytes("cmo_monthly.xlsx", CMO_URL)
        wide = pd.read_excel(io.BytesIO(raw), sheet_name="Monthly Prices", header=None)
        names = wide.iloc[4].astype(str).str.strip()
        units = wide.iloc[5].astype(str).str.strip()
        body = wide.iloc[6:].copy()
        body.index = pd.to_datetime(body.iloc[:, 0].astype(str), format="%YM%m")
        body = body.iloc[:, 1:]
        body.columns = names.iloc[1:]
        out = pd.DataFrame(index=body.index)
        attrs = {}
        for key, label in METALS.items():
            matches = [c for c in body.columns if c.lower().startswith(label.lower())]
            col = matches[0]
            out[key] = pd.to_numeric(body[col].replace("…", np.nan), errors="coerce")
            attrs[key] = units.iloc[1:][names.iloc[1:] == col].iloc[0]
        out.index.name = "date"
        out.attrs["units"] = attrs
        return out

    return _parquet_cached("commodity_prices", build)


# ======================================================================= synthetic
def assay_qaqc(seed: int = 7) -> dict[str, pd.DataFrame]:
    """SYNTHETIC gold assay QA/QC for a drilling programme (g/t Au), 120 lab batches.

    duplicates : 400 field-duplicate pairs. Each member has sd = sqrt((0.12 c)^2 + 0.02^2) g/t; 6
        pairs are sample swaps (planted gross errors).
    crm        : 3 certified reference materials inserted 60 times each. Lab A is unbiased;
                 batches 70-90 carry a planted +8% bias; 3 planted CRM mix-ups.
    blanks     : 180 coarse blanks; 7 planted carry-over contaminations after high-grade samples.
    All grades are floored at the detection limit DL = 0.005 g/t.
    """
    rng = np.random.default_rng(seed)
    DL = 0.005
    n_dup = 400
    true = np.exp(rng.normal(np.log(0.8), 1.0, n_dup))  # lognormal, median 0.8 g/t
    rel_sd = 0.12
    abs_sd = 0.02
    sd = np.sqrt((rel_sd * true) ** 2 + abs_sd ** 2)  # each member carries its own error
    orig = np.maximum(true + rng.normal(0, 1, n_dup) * sd, DL)
    dup = np.maximum(true + rng.normal(0, 1, n_dup) * sd, DL)
    batch = rng.integers(1, 121, n_dup)
    swap = rng.choice(n_dup, 6, replace=False)
    dup[swap] = np.exp(rng.normal(np.log(0.8), 1.0, 6))  # unrelated sample
    duplicates = pd.DataFrame({"batch": batch, "original": orig, "duplicate": dup,
                               "planted_swap": np.isin(np.arange(n_dup),
                                   swap)}).sort_values("batch")

    crms = {"CRM-A": (0.52, 0.025), "CRM-B": (1.85, 0.07), "CRM-C": (4.10, 0.14)}
    rows = []
    for name, (cert, sd) in crms.items():
        b = np.sort(rng.integers(1, 121, 60))
        val = rng.normal(cert, sd * 0.9, 60)
        bias = np.where((b >= 70) & (b <= 90), 1.08, 1.0)
        rows.append(pd.DataFrame({"crm": name, "batch": b, "certified": cert, "cert_sd": sd,
                                  "assay": val * bias, "planted_bias": bias > 1}))
    crm = pd.concat(rows, ignore_index=True)
    mix = rng.choice(len(crm), 3, replace=False)
    crm["planted_mixup"] = False
    crm.loc[mix, "assay"] = crm.loc[mix, "assay"].values * rng.choice([0.3, 2.5], 3)
    crm.loc[mix, "planted_mixup"] = True

    n_bl = 180
    bl_batch = np.sort(rng.integers(1, 121, n_bl))
    blank = np.maximum(rng.normal(0.002, 0.0015, n_bl), 0)
    prev_grade = np.exp(rng.normal(np.log(0.8), 1.0, n_bl))
    contaminated = rng.choice(np.where(prev_grade > 3)[0], 7, replace=False)
    blank[contaminated] += prev_grade[contaminated] * rng.uniform(0.004, 0.02, 7)
    blanks = pd.DataFrame({"batch": bl_batch, "assay": np.maximum(blank, 0),
        "prev_grade": prev_grade,
                           "below_dl": blank < DL, "planted_contamination":
                           np.isin(np.arange(n_bl), contaminated)})
    blanks["assay_reported"] = np.where(blanks.below_dl, DL / 2, blanks.assay)
    return {"duplicates": duplicates, "crm": crm, "blanks": blanks, "dl": DL}


def equipment_failures(seed: int = 11) -> dict[str, pd.DataFrame]:
    """SYNTHETIC haul-truck reliability for a fleet of 12 trucks over 3 years (26,280 h).

    failures   : one row per component failure or censoring event: truck, component, operating
                 hours between failures (TBF), repair hours (TTR), event (1 = failure, 0 = right-
                 censored at study end).  Components: 'engine' Weibull(beta=2.4, eta=2400 h)
                 wear-out; 'hydraulics' Weibull(beta=0.85, eta=900 h) infant-mortality / random;
                 'tyres' Weibull(beta=3.0, eta=1800 h).
    production : daily tonnes hauled by the fleet with a planted step change (-7%) from day 600.
    """
    rng = np.random.default_rng(seed)
    horizon = 3 * 365 * 24.0
    comps = {"engine": (2.4, 2400.0, (np.log(18), 0.5)),
             "hydraulics": (0.85, 900.0, (np.log(6), 0.6)),
             "tyres": (3.0, 1800.0, (np.log(4), 0.3))}
    rows = []
    for truck in range(1, 13):
        for comp, (beta, eta, (mu, s)) in comps.items():
            t = 0.0
            while t < horizon:
                tbf = eta * rng.weibull(beta)
                if t + tbf > horizon:
                    rows.append(dict(truck=truck, component=comp, tbf=horizon - t, ttr=np.nan,
                        event=0))
                    break
                ttr = np.exp(rng.normal(mu, s))
                rows.append(dict(truck=truck, component=comp, tbf=tbf, ttr=ttr, event=1))
                t += tbf + ttr
    failures = pd.DataFrame(rows)
    days = 3 * 365
    base = 48000 + 2500 * np.sin(2 * np.pi * np.arange(days) / 365)  # seasonal (weather) effect
    tonnes = rng.normal(base, 3000)
    tonnes[600:] *= 0.93
    production = pd.DataFrame({"day": np.arange(1, days + 1), "tonnes": tonnes,
                               "planted_shift": np.arange(days) >= 600})
    return {"failures": failures, "production": production}


def drillhole_from_walker(column_x: int = 120, seed: int = 3) -> pd.DataFrame:
    """A 1-D 'drillhole' built from Walker Lake exhaustive column x=column_x: V along y is treated
    as grade down-hole, sampled over irregular intervals (0.5-3.5 units) as a lab would receive
    them. Columns: from, to, length, grade (length-weighted mean of the exhaustive V)."""
    rng = np.random.default_rng(seed)
    ex = walker_exhaustive()
    col = ex[ex.x == column_x].sort_values("y")
    depth = col.y.values - col.y.min()
    v = col.V.values
    rows, start = [], 0.0
    while start < depth.max() - 0.5:
        length = float(rng.choice([0.5, 1.0, 1.0, 1.5, 2.0, 2.0, 3.0, 3.5]))
        end = min(start + length, depth.max())
        m = (depth >= start) & (depth < end)
        if m.sum() == 0:
            m = np.argmin(np.abs(depth - start))
        rows.append(dict(**{"from": start, "to": end}, length=end - start,
            grade=float(np.mean(v[m]))))
        start = end
    return pd.DataFrame(rows)


def fetch_all():
    print("Fetching / verifying cached datasets", file=sys.stderr)
    for fn in (walker, walker_exhaustive, jura, gslib_cluster, msha_mines, msha_production,
               msha_accidents, clark_fork_wq, clark_fork_flow, commodity_prices):
        df = fn()
        print(f"  {fn.__name__:18s} {df.shape}", file=sys.stderr)
    print(f"Done ({date.today()}).", file=sys.stderr)


if __name__ == "__main__":
    fetch_all()
