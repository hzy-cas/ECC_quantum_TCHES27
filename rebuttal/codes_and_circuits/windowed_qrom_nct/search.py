"""Resumable, audit-safe exact ``(s,p)`` search.

No surrogate row is ever written to an exact search file.  A row is appended
only after :func:`synthesize_design` has completed and the complete per-line
ASAP result has been serialized.  The companion coverage manifest enumerates
the entire requested grid and therefore makes interrupted searches obvious.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

from .curves import CURVES
from .synthesis import ExactDesign, synthesize_design, write_design


SEARCH_SCHEMA = "windowed-exact-grid-v1"
SEARCH_FIELDS = (
    "schema",
    "model_version",
    "n",
    "s",
    "p",
    "arithmetic",
    "addends",
    "toffoli",
    "cnot",
    "width",
    "toffoli_depth",
    "nct_depth",
    "current_depth",
    "dw",
    "tdw",
    "table_seed",
    "scheduler",
    "artifact",
)


def grid_points(n: int, s_min: int = 2, s_max: int = 18):
    if n not in CURVES:
        raise ValueError(f"unsupported field size {n}")
    if not 2 <= s_min <= s_max <= 18:
        raise ValueError("require 2 <= s_min <= s_max <= 18")
    scalar_bits = n + 1
    for s in range(s_min, s_max + 1):
        addends = 2 * ((scalar_bits + s - 1) // s)
        for p in range(1, addends + 1):
            yield s, p


def _key(row: dict[str, object]) -> tuple[int, int]:
    return int(row["s"]), int(row["p"])


def _read_completed(path: Path, n: int, arithmetic: str) -> set[tuple[int, int]]:
    if not path.is_file():
        return set()
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    completed: set[tuple[int, int]] = set()
    for row in rows:
        if row.get("schema") != SEARCH_SCHEMA:
            raise ValueError(f"foreign search schema in {path}")
        if int(row["n"]) != n or row["arithmetic"] != arithmetic:
            raise ValueError(f"mixed search scope in {path}")
        key = _key(row)
        if key in completed:
            raise ValueError(f"duplicate exact point {key} in {path}")
        completed.add(key)
    return completed


def _row(design: ExactDesign, artifact: Path) -> dict[str, object]:
    payload = asdict(design)
    return {
        "schema": SEARCH_SCHEMA,
        **{field: payload[field] for field in SEARCH_FIELDS if field in payload},
        "artifact": str(artifact),
    }


def run_exact_grid(
    *,
    n: int,
    arithmetic: str,
    output: Path,
    artifact_dir: Path,
    s_min: int = 2,
    s_max: int = 18,
) -> None:
    """Run or resume an exact grid, flushing every completed design point."""
    if arithmetic != "ac":
        raise ValueError("this submission artifact supports only AC arithmetic")
    completed = _read_completed(output, n, arithmetic)
    output.parent.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    write_header = not output.is_file() or output.stat().st_size == 0
    with output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SEARCH_FIELDS)
        if write_header:
            writer.writeheader()
            handle.flush()
        for s, p in grid_points(n, s_min, s_max):
            if (s, p) in completed:
                continue
            design = synthesize_design(n, s, p, arithmetic)
            artifact = artifact_dir / f"n{n}_s{s}_p{p}_{arithmetic}.json"
            write_design(artifact, design)
            writer.writerow(_row(design, artifact))
            handle.flush()


def index_artifacts(
    *, n: int, arithmetic: str, artifact_dir: Path, output: Path
) -> int:
    """Build a deterministic exact-row index from already completed artifacts."""
    rows: dict[tuple[int, int], tuple[ExactDesign, Path]] = {}
    for artifact in sorted(artifact_dir.glob(f"n{n}_s*_p*_{arithmetic}.json")):
        design = ExactDesign(**json.loads(artifact.read_text(encoding="utf-8")))
        if design.n != n or design.arithmetic != arithmetic:
            raise ValueError(f"artifact scope mismatch: {artifact}")
        key = (design.s, design.p)
        if key in rows:
            raise ValueError(f"duplicate exact artifact for {key}")
        rows[key] = (design, artifact)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SEARCH_FIELDS)
        writer.writeheader()
        for key in sorted(rows):
            design, artifact = rows[key]
            writer.writerow(_row(design, artifact))
    return len(rows)


def write_coverage(
    path: Path,
    *,
    n: int,
    arithmetic: str,
    search_csv: Path,
    s_min: int = 2,
    s_max: int = 18,
) -> dict[str, object]:
    completed = _read_completed(search_csv, n, arithmetic)
    points = list(grid_points(n, s_min, s_max))
    missing = [{"s": s, "p": p} for s, p in points if (s, p) not in completed]
    payload: dict[str, object] = {
        "schema": "windowed-exact-grid-coverage-v1",
        "search_schema": SEARCH_SCHEMA,
        "n": n,
        "arithmetic": arithmetic,
        "s_min": s_min,
        "s_max": s_max,
        "expected_points": len(points),
        "exact_points": len(completed),
        "complete": not missing,
        "missing_points": missing,
        "search_csv": str(search_csv),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload
