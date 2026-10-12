"""Small, offline tests of cache integrity and snapshot provenance (no rendering)."""
import json

import pandas as pd
import pytest

from miningstats.provenance import cached_bytes, cached_frame, data_manifest


def test_raw_integrity_and_request_identity(tmp_path):
    path = tmp_path / "raw" / "series.csv"
    calls = []

    def fetch(url, params):
        calls.append((url, params))
        return str(params["end"]).encode()

    assert cached_bytes(path, "https://example.test", {"end": 1}, fetch) == b"1"
    assert cached_bytes(path, "https://example.test", {"end": 1}, fetch) == b"1"
    assert len(calls) == 1
    assert cached_bytes(path, "https://example.test", {"end": 2}, fetch) == b"2"
    path.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="changed"):
        cached_bytes(path, "https://example.test", {"end": 2}, fetch)


def test_legacy_history_is_unknown_and_pinned_inputs_refetched(tmp_path):
    dynamic = tmp_path / "raw" / "legacy.csv"
    dynamic.parent.mkdir()
    dynamic.write_bytes(b"old")
    assert cached_bytes(dynamic, "https://example.test", {}, lambda *_: b"new") == b"old"
    meta = json.loads(dynamic.with_name(dynamic.name + ".metadata.json").read_text())
    assert meta["retrieved_at"] is None
    assert meta["origin"] == "legacy-cache"
    pinned = dynamic.with_name("pinned.rda")
    pinned.write_bytes(b"old")
    assert cached_bytes(pinned, "https://example.test/commit", {}, lambda *_: b"pinned",
                        pinned=True) == b"pinned"
    snapshot = data_manifest(tmp_path)
    assert snapshot["files"]["raw/pinned.rda"]["retrieved_at"]


def test_parsed_cache_rebuilds_on_cutoff_and_rejects_mutation(tmp_path):
    path = tmp_path / "series.parquet"
    calls = []

    def build():
        calls.append(1)
        return pd.DataFrame({"value": [len(calls)]})

    kwargs = {"data_dir": tmp_path, "pd": pd}
    assert cached_frame(path, build, identity={"cutoff": 1}, **kwargs).value.iloc[0] == 1
    assert cached_frame(path, build, identity={"cutoff": 1}, **kwargs).value.iloc[0] == 1
    assert cached_frame(path, build, identity={"cutoff": 2}, **kwargs).value.iloc[0] == 2
    path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="changed"):
        cached_frame(path, build, identity={"cutoff": 2}, **kwargs)
