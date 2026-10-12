"""Hash-checked download and parsed caches; metadata never invents retrieval history."""
from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cached_bytes(path: Path, url: str, params: dict | None, fetch, *, pinned=False) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    sidecar = path.with_name(path.name + ".metadata.json")
    request = {"url": url, "params": params or {}}
    meta = json.loads(sidecar.read_text()) if sidecar.exists() else None
    if meta and path.exists() and sha256(path) != meta["sha256"]:
        raise ValueError(f"cached download changed: {path}; remove it explicitly to refresh")
    # Legacy dynamic downloads may be reused, but their retrieval time is unknown.
    # Pinned teaching inputs must be fetched from the pinned revision at least once.
    refresh = not path.exists() or (meta and meta["request"] != request) or (pinned and not meta)
    if refresh:
        payload = fetch(url, params)
        temporary = path.with_name(path.name + ".download")
        temporary.write_bytes(payload)
        temporary.replace(path)
        now = datetime.now(UTC).isoformat()
        meta = {"request": request, "retrieved_at": now, "observed_at": now,
                "origin": "downloaded", "sha256": sha256(path), "bytes": path.stat().st_size}
    elif meta is None:
        meta = {"request": request, "retrieved_at": None,
                "observed_at": datetime.now(UTC).isoformat(), "origin": "legacy-cache",
                "note": "Request describes the loader; original retrieval history is unavailable.",
                "sha256": sha256(path), "bytes": path.stat().st_size}
    sidecar.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return path.read_bytes()


def cached_frame(path: Path, builder, *, identity: dict, data_dir: Path, pd):
    path.parent.mkdir(parents=True, exist_ok=True)
    sidecar = path.with_name(path.name + ".metadata.json")
    meta = json.loads(sidecar.read_text()) if sidecar.exists() else None
    raw = data_dir / "raw"
    inputs = {p.name: sha256(p) for p in sorted(raw.glob("*"))
              if p.is_file() and not p.name.endswith((".metadata.json", ".download"))}
    if path.exists() and meta and meta["identity"] == identity and meta["raw_inputs"] == inputs:
        if sha256(path) != meta["sha256"]:
            raise ValueError(f"parsed cache changed: {path}")
        return pd.read_parquet(path)
    frame = builder()
    frame.to_parquet(path)
    inputs = {p.name: sha256(p) for p in sorted(raw.glob("*"))
              if p.is_file() and not p.name.endswith((".metadata.json", ".download"))}
    meta = {"identity": identity, "raw_inputs": inputs, "sha256": sha256(path),
            "rows": len(frame), "columns": list(frame.columns)}
    sidecar.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return frame


def data_manifest(data_dir: Path) -> dict:
    files = {}
    for p in sorted(data_dir.rglob("*")):
        if p.is_file() and not p.name.endswith((".metadata.json", ".download")):
            meta_path = p.with_name(p.name + ".metadata.json")
            if not meta_path.exists():
                continue
            metadata = json.loads(meta_path.read_text())
            if sha256(p) != metadata["sha256"]:
                raise ValueError(f"manifest input changed: {p}")
            files[p.relative_to(data_dir).as_posix()] = metadata
    return {"schema_version": 1, "files": files}
