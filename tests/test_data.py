"""Offline regression fixtures for mutable-provider snapshots (monthly series, water quality)."""
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


def test_water_quality_loader_clips_provider_rows_after_cutoff(tmp_path, monkeypatch):
    monkeypatch.setattr(data, "DATA_DIR", tmp_path)
    monkeypatch.setattr(data, "RAW", tmp_path / "raw")
    monkeypatch.setattr(data, "END_DATE", "2026-09-30")
    raw = (b"ActivityStartDate,CharacteristicName,ResultSampleFractionText,ResultMeasureValue,"
           b"ResultMeasure/MeasureUnitCode,ResultDetectionConditionText,"
           b"DetectionQuantitationLimitMeasure/MeasureValue,"
           b"DetectionQuantitationLimitMeasure/MeasureUnitCode,ActivityMediaName\n"
           b"2026-09-30,Copper,Recoverable,12,ug/l,,,,Water\n"
           b"2026-10-01,Copper,Recoverable,15,ug/l,,,,Water\n")
    monkeypatch.setattr(data, "_get", lambda url, params: raw)
    result = data.clark_fork_wq()
    assert result.date.max() == pd.Timestamp("2026-09-30")
    assert len(result) == 1
    assert data.clark_fork_wq().equals(result)
