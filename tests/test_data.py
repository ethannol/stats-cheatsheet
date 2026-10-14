"""Offline regression fixtures for mutable-provider monthly snapshots."""
import pandas as pd

from statsheet import data


def test_monthly_loaders_clip_newer_rows_even_in_cached_downloads(tmp_path, monkeypatch):
    monkeypatch.setattr(data, "DATA_DIR", tmp_path)
    monkeypatch.setattr(data, "RAW", tmp_path / "raw")
    monkeypatch.setattr(data, "END_DATE", "2026-09-30")
    co2 = (b"year,month,decimal date,average,deseasonalized,ndays,sdev,unc\n"
           b"2026,9,2026.7,425,424,30,1,0.1\n"
           b"2026,10,2026.8,426,425,31,1,0.1\n")
    gistemp = ("provider comment\nYear,Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec\n"
               "2026," + ",".join(str(i) for i in range(12)) + "\n").encode()
    monkeypatch.setattr(data, "_get", lambda url, params: co2 if "co2" in url else gistemp)
    for load in (data.co2, data.gistemp):
        result = load()
        assert result.index.max() == pd.Timestamp("2026-09-01")
        assert load().equals(result)
