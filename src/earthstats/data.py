"""Dataset loaders. Every loader downloads once, caches under data/, and returns a tidy
pandas object with a DatetimeIndex (or lat/lon columns for the station table).

Datasets
--------
co2()             NOAA GML Mauna Loa monthly mean CO2 (ppm), 1958-present.
                  Textbook trend + seasonality series.
gistemp()         NASA GISTEMP v4 global land-ocean temperature anomaly (deg C vs 1951-80),
                  monthly 1880-present. Trend-detection and regression target.
ghcn_daily(id)    NOAA GHCN-Daily station record: TMAX/TMIN (deg C), PRCP (mm), SNOW (mm).
                  Default station is Central Park NYC (USW00094728), 1950-present.
                  Daily data -> distributions, extremes, hypothesis tests.
usgs_daily(site)  USGS NWIS daily mean discharge (cfs). Default is the Potomac at Little
                  Falls (01646500), 1930-present. Floods, return periods, hydrology.
station_normals() NOAA GSOY annual summaries for ~1,700 first-order US airport stations,
                  1991-2020, collapsed to one climatology row per station.
                  Spatial statistics, multivariate analysis, multiple regression.

Run `python -m earthstats.data` to fetch everything up front.
"""

from __future__ import annotations

import io
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from . import DATA_DIR
from .provenance import cached_bytes, cached_frame, sha256

RAW = DATA_DIR / "raw"
END_DATE = "2026-09-30"  # fixed so cached results are reproducible; bump + `make refresh`
NCEI = "https://www.ncei.noaa.gov/access/services/data/v1"
UA = {"User-Agent": "earthstats-cheatsheet/0.1 (educational; contact via github)"}


def _get(url: str, params: dict | None = None, timeout: int = 180) -> bytes:
    r = requests.get(url, params=params, headers=UA, timeout=timeout)
    r.raise_for_status()
    return r.content


def _cached_bytes(name: str, url: str, params: dict | None = None) -> bytes:
    return cached_bytes(RAW / name, url, params, _get, pinned="raw.githubusercontent.com" in url)


def _parquet_cached(name: str, builder):
    # A loader edit or a changed date cutoff always rebuilds old parsed caches.
    return cached_frame(DATA_DIR / f"{name}.parquet", builder,
                        identity={"loader_sha256": sha256(Path(__file__)),
                                  "end_date": END_DATE}, data_dir=DATA_DIR, pd=pd)


# --------------------------------------------------------------------------- CO2
def co2() -> pd.DataFrame:
    """Monthly Mauna Loa CO2. Columns: co2 (ppm), deseasonalized, ndays, sdev, unc."""

    def build():
        raw = _cached_bytes("co2_mm_mlo.csv",
                            "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.csv")
        df = pd.read_csv(io.BytesIO(raw), comment="#")
        df.index = pd.to_datetime(dict(year=df.year, month=df.month, day=1))
        df.index.name = "date"
        df = df.rename(columns={"average": "co2", "decimal date": "decimal_date"})
        df = df.replace({-9.99: np.nan, -0.99: np.nan, -1: np.nan})
        return df.loc[:END_DATE, ["decimal_date", "co2", "deseasonalized", "ndays", "sdev", "unc"]]

    return _parquet_cached("co2", build)


# ----------------------------------------------------------------------- GISTEMP
def gistemp() -> pd.DataFrame:
    """Monthly global land-ocean anomaly (deg C relative to 1951-1980). Column: anomaly."""

    def build():
        raw = _cached_bytes("gistemp_glb.csv",
                            "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv")
        wide = pd.read_csv(io.BytesIO(raw), skiprows=1, na_values="***")
        months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        long = wide.melt(id_vars="Year", value_vars=months, var_name="mon", value_name="anomaly")
        long["date"] = pd.to_datetime(long.Year.astype(str) + "-" + long.mon, format="%Y-%b")
        out = long.set_index("date").sort_index()[["anomaly"]].dropna()
        return out.loc[:END_DATE]

    return _parquet_cached("gistemp", build)


# ------------------------------------------------------------------- GHCN daily
CENTRAL_PARK = "USW00094728"
REAGAN_NATIONAL = "USW00013743"  # Washington DC; precipitation driver for the Potomac


def ghcn_daily(station: str = CENTRAL_PARK, start: str = "1950-01-01") -> pd.DataFrame:
    """Daily TMAX, TMIN (deg C), PRCP, SNOW (mm) for one GHCN-Daily station."""

    def build():
        raw = _cached_bytes(f"ghcnd_{station}.csv" if start == "1950-01-01" else f"ghcnd_{station}_{start}.csv", NCEI, dict(
            dataset="daily-summaries", stations=station, startDate=start, endDate=END_DATE,
            dataTypes="TMAX,TMIN,PRCP,SNOW", format="csv", units="metric",
            includeStationLocation="1"))
        df = pd.read_csv(io.BytesIO(raw), parse_dates=["DATE"]).set_index("DATE")
        df.index.name = "date"
        df.attrs = dict(station=station, lat=float(df.LATITUDE.iloc[0]),
                        lon=float(df.LONGITUDE.iloc[0]), elev=float(df.ELEVATION.iloc[0]))
        df = df.loc[:END_DATE, ["TMAX", "TMIN", "PRCP", "SNOW"]].astype(float)
        # reindex to a complete daily calendar so gaps show up as NaN, not silently missing
        full = pd.date_range(df.index.min(), df.index.max(), freq="D", name="date")
        return df.reindex(full)

    return _parquet_cached(f"ghcnd_{station}_{start}_{END_DATE}", build)


# ------------------------------------------------------------------ USGS daily Q
POTOMAC = "01646500"


def usgs_daily(site: str = POTOMAC, start: str = "1930-01-01") -> pd.DataFrame:
    """Daily mean discharge (cfs) for a USGS gauge. Columns: flow_cfs, qualifier."""

    def build():
        raw = _cached_bytes(f"usgs_{site}.json" if start == "1930-01-01" else f"usgs_{site}_{start}.json", "https://waterservices.usgs.gov/nwis/dv/", dict(
            format="json", sites=site, parameterCd="00060", statCd="00003",
            startDT=start, endDT=END_DATE))
        ts = json.loads(raw)["value"]["timeSeries"][0]
        vals = ts["values"][0]["value"]
        df = pd.DataFrame(vals)
        df["date"] = pd.to_datetime(df.dateTime.str[:10])
        df["flow_cfs"] = pd.to_numeric(df["value"]).replace(-999999.0, np.nan)
        df["qualifier"] = df.qualifiers.str.join(",")
        df = df.set_index("date").loc[:END_DATE, ["flow_cfs", "qualifier"]]
        df.attrs = dict(site=site, name=ts["sourceInfo"]["siteName"])
        full = pd.date_range(df.index.min(), df.index.max(), freq="D", name="date")
        return df.reindex(full)

    return _parquet_cached(f"usgs_{site}_{start}_{END_DATE}", build)


# ------------------------------------------------------------ station normals
GSOY_TYPES = ["TAVG", "TMAX", "TMIN", "PRCP", "SNOW", "DP01", "DP10", "DX90", "DX32",
              "DT32", "DT00", "EMXT", "EMNT", "HTDD", "CLDD"]
CONUS = dict(lat_min=24.5, lat_max=49.5, lon_min=-125.0, lon_max=-66.5)


def _conus_airport_stations() -> pd.DataFrame:
    raw = _cached_bytes("ghcnd-stations.txt",
                        "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt")
    rows = []
    for line in raw.decode("latin-1").splitlines():
        if not line.startswith("USW"):
            continue
        rows.append(dict(station=line[0:11], lat=float(line[12:20]), lon=float(line[21:30]),
                         elev=float(line[31:37]), name=line[41:71].strip()))
    st = pd.DataFrame(rows)
    m = st.lat.between(CONUS["lat_min"], CONUS["lat_max"]) & \
        st.lon.between(CONUS["lon_min"], CONUS["lon_max"])
    return st[m].reset_index(drop=True)


def station_annual(start: str = "1991-01-01", end: str = "2020-12-31") -> pd.DataFrame:
    """Raw GSOY annual rows (one per station-year) for CONUS airport stations."""

    def build():
        stations = _conus_airport_stations()
        frames = []
        chunk = 100
        for i in range(0, len(stations), chunk):
            ids = ",".join(stations.station.iloc[i:i + chunk])
            raw = _cached_bytes(f"gsoy_{i:04d}.csv" if (start, end) == ("1991-01-01", "2020-12-31") else f"gsoy_{i:04d}_{start}_{end}.csv", NCEI, dict(
                dataset="global-summary-of-the-year", stations=ids, startDate=start,
                endDate=end, dataTypes=",".join(GSOY_TYPES), format="csv", units="metric",
                includeStationLocation="1", includeStationName="1"))
            if raw.strip():
                frames.append(pd.read_csv(io.BytesIO(raw)))
        df = pd.concat(frames, ignore_index=True)
        df.columns = [c.lower() for c in df.columns]
        df = df.rename(columns={"date": "year", "latitude": "lat", "longitude": "lon",
                                "elevation": "elev"})
        return df

    return _parquet_cached(f"gsoy_annual_{start}_{end}", build)


def station_normals(min_years: int = 25) -> pd.DataFrame:
    """One row per station: 1991-2020 means of the GSOY variables, lat/lon/elev/name.

    Stations are kept only if they have >= min_years of TAVG and PRCP, so the table is
    dominated by long-running first-order stations rather than new ASOS sites.
    Variables (all annual): tavg/tmax/tmin (deg C), prcp/snow (mm), dp01/dp10 (days with
    >= 0.25 mm / >= 2.5 mm precip), dx90/dx32 (days TMAX >= 32.2 C / <= 0 C), dt32/dt00
    (days TMIN <= 0 C / <= -17.8 C), emxt/emnt (annual extreme max/min T), htdd/cldd
    (heating/cooling degree days, base 18.3 C).
    """
    ann = station_annual()
    vars_ = [v.lower() for v in GSOY_TYPES]
    complete = ann.dropna(subset=["tavg", "prcp"])
    counts = complete.groupby("station").size()
    keep = counts[counts >= min_years].index
    sub = ann[ann.station.isin(keep)]
    normals = sub.groupby("station")[vars_].mean()
    meta = sub.groupby("station")[["name", "lat", "lon", "elev"]].first()
    out = meta.join(normals).join(counts.rename("n_years"))
    return out.reset_index()


def fetch_all():
    print("Fetching / verifying cached datasets", file=sys.stderr)
    for fn in (co2, gistemp, ghcn_daily, lambda: ghcn_daily(REAGAN_NATIONAL), usgs_daily,
               station_normals):
        df = fn()
        print(f"  {getattr(fn, '__name__', 'ghcn_daily(DC)'):18s} {df.shape}", file=sys.stderr)
    print(f"Done ({date.today()}).", file=sys.stderr)


if __name__ == "__main__":
    fetch_all()
