"""Format recorded transcripts without rerunning analyses or rendering figures."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    reports = ROOT / "reports"
    path = reports / "results.json"
    manifest = json.loads(path.read_text())
    chapters = [
        "# Generated worked-example results",
        "",
        "These are executed educational examples, not a new scientific study.",
        f"Data cutoff: `{manifest['end_date']}`. Source: `{manifest['source_commit']}`.",
        "The run manifest records the exact executed source, dependencies and hashes.",
        "Full transcripts retain diagnostic warnings, assumptions and interpretation.",
        "The excerpts below omit diagnostic warnings and read numerical output directly",
        "from those hash-checked transcripts. Synthetic examples are labelled in source.",
        "",
    ]
    for name, analysis in sorted(manifest["analyses"].items()):
        output = (reports / analysis["output"]).resolve()
        if (
            not output.is_relative_to(reports.resolve())
            or digest(output) != analysis["output_sha256"]
        ):
            raise ValueError(f"recorded output changed: {name}")
        lines = [
            line
            for line in output.read_text().splitlines()
            if re.search(r"\d", line)
            and not re.search(r"Warning:|/site-packages/|\[fig\]", line)
            and not line.lstrip().startswith((">", "File ", "warnings.", "return "))
        ]
        stem = Path(name).stem
        chapters += [
            f"## {stem.replace('_', ' ')}",
            "",
            f"[Source](../analyses/{name}) · [Full output]({analysis['output']})",
            "",
            "```text",
            *lines[:35],
            "```",
            "",
        ]
    target = reports / "worked-results.md"
    target.write_text("\n".join(chapters))
    manifest["presentation"] = {
        "formatted_at": datetime.now(UTC).isoformat(),
        "formatter": "scripts/format_results.py",
        "formatter_sha256": digest(Path(__file__)),
        "output": "worked-results.md",
        "output_sha256": digest(target),
        "note": "Text-only formatting of the original hash-checked numerical transcripts.",
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
