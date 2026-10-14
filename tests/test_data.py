"""Offline regression fixture for a mutable water-quality response."""
import pandas as pd

from statsheet import data


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
