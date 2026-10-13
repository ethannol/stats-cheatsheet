"""Run worked examples and publish hash-addressed outputs (workstation/CI only)."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from miningstats import data
from miningstats.provenance import data_manifest, sha256

ROOT = Path(__file__).resolve().parents[1]
PREVIEWS = ['04_geostatistics/estimation_maps.png', '03_qaqc/crm_control_charts.png',
            '10_censored_inference/seasonal_kendall.png']


def verify_raw_snapshot(manifest_path: Path) -> None:
    snapshot = json.loads(manifest_path.read_text())
    for name, record in snapshot["files"].items():
        if not name.startswith("raw/"):
            continue
        path = data.DATA_DIR / name
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"exact input snapshot missing or changed: {name}")


def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-inputs", type=Path)
    args = parser.parse_args()
    if args.verify_inputs:
        verify_raw_snapshot(args.verify_inputs)
    env = dict(os.environ, MPLBACKEND="Agg", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")
    sources = {p.relative_to(ROOT).as_posix(): sha256(p)
               for folder in ("src", "analyses", "scripts")
               for p in sorted((ROOT / folder).rglob("*.py"))}
    sources["uv.lock"] = sha256(ROOT / "uv.lock")
    patch = subprocess.check_output(["git", "-C", str(ROOT), "diff", "HEAD", "--",
                                     "src", "analyses", "scripts", "uv.lock"])
    manifest = {"schema_version": 1, "started_at": datetime.now(UTC).isoformat(),
                "source_commit": git("rev-parse", "HEAD"),
                "source_patch_sha256": hashlib.sha256(patch).hexdigest(),
                "source_files": sources, "end_date": data.END_DATE,
                "python": sys.version, "packages": {d.metadata["Name"]: d.version
                    for d in importlib.metadata.distributions()},
                "threads": {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                               "MKL_NUM_THREADS")}, "analyses": {}}
    if hasattr(data, "GSTAT_REVISION"):
        manifest["upstream_gstat_commit"] = data.GSTAT_REVISION
    with tempfile.TemporaryDirectory(prefix="worked-results-", dir=ROOT) as temp:
        stage = Path(temp)
        (stage / "outputs").mkdir()
        (stage / "figures").mkdir()
        if patch:
            (stage / "source.patch").write_bytes(patch)
        chapters = ["# Generated worked-example results", "",
                    "These are executed educational examples, not a new scientific study.",
                    f"Data cutoff: `{data.END_DATE}`. Source: `{manifest['source_commit']}`.",
                    "The run manifest records any source patch and its hash, dependencies,",
                    "input hashes, script hashes and output hashes. Full transcripts are linked",
                    "for assumptions and interpretation; the numeric excerpts below are emitted",
                    "directly by each script and may include illustrative synthetic examples.", ""]
        for script in sorted((ROOT / "analyses").glob("*.py")):
            print(f"Running {script.name}", flush=True)
            output = stage / "outputs" / f"{script.stem}.txt"
            try:
                with output.open("w") as stream:
                    subprocess.run([sys.executable, str(script)], cwd=ROOT, env=env,
                                   stdout=stream, stderr=subprocess.STDOUT, check=True,
                                   timeout=14400)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                manifest["failed_script"] = script.name
                manifest["failure"] = str(error)
                manifest["failed_at"] = datetime.now(UTC).isoformat()
                (stage / "results.json").write_text(
                    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
                )
                stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
                failure_dir = ROOT / "reports" / "failed-runs" / stamp
                shutil.copytree(stage, failure_dir)
                print(f"Failure evidence preserved at {failure_dir}", file=sys.stderr)
                print("\n".join(output.read_text().splitlines()[-60:]), file=sys.stderr)
                raise
            transcript = output.read_text().replace(str(ROOT) + "/", "")
            output.write_text(transcript)
            # Select numerical output, not a separately maintained copy of headline values.
            lines = [line for line in transcript.splitlines()
                     if any(c.isdigit() for c in line) and "[fig]" not in line
                     and not line.lstrip().startswith((">", "File ", "warnings.", "return "))]
            chapters += [f"## {script.stem.replace('_', ' ')}", "",
                         f"[Source](../analyses/{script.name}) · "
                         f"[Full output](outputs/{script.stem}.txt)", "", "```text",
                         *lines[:35], "```", ""]
            manifest["analyses"][script.name] = {
                "script_sha256": sha256(script), "output": f"outputs/{output.name}",
                "output_sha256": sha256(output), "exit_code": 0}
        for relative in PREVIEWS:
            source = ROOT / "figures" / relative
            shutil.copyfile(source, stage / "figures" / source.name)
        manifest["figures"] = {p.name: sha256(p) for p in sorted((stage / "figures").glob("*.png"))}
        manifest["finished_at"] = datetime.now(UTC).isoformat()
        inputs = data_manifest(data.DATA_DIR)
        (stage / "data-manifest.json").write_text(
            json.dumps(inputs, indent=2, sort_keys=True) + "\n"
        )
        (stage / "worked-results.md").write_text("\n".join(chapters))
        (stage / "results.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        # Publish only after every analysis and preview succeeded; leave old results on failure.
        destination = ROOT / "reports"
        destination.mkdir(exist_ok=True)
        for item in stage.iterdir():
            target = destination / item.name
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()
            shutil.move(str(item), target)
    print("Published reports/worked-results.md and reproducibility manifests", flush=True)


if __name__ == "__main__":
    main()
