#!/usr/bin/env python3
"""Fast module-level Window-QROM scan from the bundled AC layer caches.

This command never emits a quantum gate.  It enumerates the complete
``s=2..18, p=1..L_s`` grid, composes every point whose required uncontrolled
``reduction(k)`` records already exist, and leaves unavailable points in the
coverage CSV with their missing lane counts.

The QROM data words are represented by the paper's declared closed-form
upper-weight model: every word weight is set to ``2n``.  This is a
table-independent upper bound for the coherent unary-iteration and
balanced-fanout QROM modules.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .paths import AC_RESULTS, PACKAGE_ROOT


DEFAULT_RESULTS = AC_RESULTS
DEFAULT_OUTPUT_DIR = PACKAGE_ROOT / "artifacts" / "analytic"
FIELD_SIZES = (163, 233, 283, 571)
DEFAULT_S_MIN = 2
DEFAULT_S_MAX = 18
RESOURCE_FIELDS = (
    "qubits",
    "toffoli",
    "cnot",
    "full_depth",
    "current_depth",
    "toffoli_depth",
)


@dataclass(frozen=True)
class Stats:
    qubits: int
    toffoli: int
    cnot: int
    full_depth: int
    current_depth: int
    toffoli_depth: int


@dataclass(frozen=True)
class TableShape:
    window_bits: int

    @property
    def entries(self) -> int:
        return 1 << self.window_bits


def _ceil_log2(value: int) -> int:
    if value <= 0:
        raise ValueError("ceil(log2(h)) requires h > 0")
    return (value - 1).bit_length()


def scheduled_tables(n: int, s: int) -> list[TableShape]:
    scalar_bits = n + 1
    windows = math.ceil(scalar_bits / s)
    top_bits = scalar_bits - s * (windows - 1)
    result: list[TableShape] = []
    for index in range(windows):
        bits = top_bits if index == windows - 1 else s
        result.extend((TableShape(bits), TableShape(bits)))
    return result


def qrom_stats(n: int, table: TableShape, word_weight: int) -> tuple[Stats, Stats]:
    """Closed-form lookup and lookup/precompute/unlookup resources."""
    b = table.window_bits
    k = table.entries
    if not 1 <= word_weight <= 2 * n:
        raise ValueError("QROM word weight outside 1..2n")
    total_hamming = word_weight * k
    log_sum = _ceil_log2(word_weight) * k
    last_log = _ceil_log2(word_weight)

    lookup_t = 2 * (k - 2)
    lookup_c = (k - 2) + 3 * total_hamming - 2 * k
    lookup_depth = 3 * k - 4 + 2 * log_sum + k
    lookup = Stats(
        qubits=4 * n + 2 * b - 2,
        toffoli=lookup_t,
        cnot=lookup_c,
        full_depth=lookup_depth,
        current_depth=lookup_depth,
        toffoli_depth=lookup_t,
    )

    pair_t = 4 * k - 2 * b - 6
    pair_depth = 2 * lookup_depth - 2 * b - 2 * last_log + 3
    pair = Stats(
        qubits=8 * n + 2 * b - 2,
        toffoli=pair_t,
        cnot=2 * lookup_c + 4 * n - 2 * (word_weight - 1),
        full_depth=pair_depth,
        current_depth=pair_depth,
        toffoli_depth=pair_t,
    )
    return lookup, pair


def load_reductions(results: Path, n: int) -> tuple[dict[int, Stats], list[str]]:
    records: dict[int, Stats] = {}
    sources: dict[int, str] = {}
    paths = sorted(results.glob(f"layers_n{n}*.csv"))
    if not paths:
        raise FileNotFoundError(f"no layer cache for n={n} under {results}")
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("mode") != "reduction" or int(row.get("n", 0)) != n:
                    continue
                lanes = int(row["lanes"])
                stats = Stats(**{field: int(row[field]) for field in RESOURCE_FIELDS})
                previous = records.get(lanes)
                if previous is not None and previous != stats:
                    raise ValueError(
                        f"conflicting reduction({lanes}) rows in {sources[lanes]} and {path}"
                    )
                records[lanes] = stats
                sources[lanes] = str(path)
    return records, [str(path) for path in paths]


def required_reduction_lanes(addends: int, p: int) -> list[int]:
    required: set[int] = set()
    for start in range(p, addends, p):
        required.add(min(p, addends - start))
    items = p
    while items > 1:
        required.add(items // 2)
        items = items // 2 + items % 2
    return sorted(required)


def required_grid_reduction_lanes(
    n: int,
    s_min: int = DEFAULT_S_MIN,
    s_max: int = DEFAULT_S_MAX,
) -> list[int]:
    """Return every ``reduction(k)`` module needed by the full ``(s,p)`` grid.

    The windowed architecture uses only uncontrolled quantum--quantum point
    additions.  Consequently these are Phase-II ``reduction`` lane counts;
    the legacy controlled ``accumulation`` cache is intentionally excluded.
    """
    if not 1 <= s_min <= s_max:
        raise ValueError("require 1 <= s_min <= s_max")
    required: set[int] = set()
    for s in range(s_min, s_max + 1):
        addends = len(scheduled_tables(n, s))
        for p in range(1, addends + 1):
            required.update(required_reduction_lanes(addends, p))
    return sorted(required)


def _early_tail(n: int, reduction: Stats, lanes: int) -> Stats:
    return Stats(
        qubits=reduction.qubits - n * lanes,
        toffoli=reduction.toffoli,
        cnot=reduction.cnot - 3 * n * lanes,
        full_depth=reduction.full_depth - 1,
        current_depth=reduction.current_depth - 1,
        toffoli_depth=reduction.toffoli_depth,
    )


def _parallel_initialization(n: int, tables: list[TableShape], word_weight: int) -> Stats:
    parts = [qrom_stats(n, table, word_weight)[0] for table in tables]
    return Stats(
        qubits=sum(part.qubits for part in parts),
        toffoli=sum(part.toffoli for part in parts),
        cnot=sum(part.cnot for part in parts),
        full_depth=max(part.full_depth for part in parts),
        current_depth=max(part.current_depth for part in parts),
        toffoli_depth=max(part.toffoli_depth for part in parts),
    )


def _window_layer(
    n: int,
    tables: list[TableShape],
    word_weight: int,
    reduction: Stats,
) -> Stats:
    lanes = len(tables)
    pairs = [qrom_stats(n, table, word_weight)[1] for table in tables]
    address_bits = sum(table.window_bits for table in tables)
    qrom_scratch = sum(
        part.qubits - table.window_bits - 4 * n
        for table, part in zip(tables, pairs)
    )
    tail = _early_tail(n, reduction, lanes)
    tail_scratch = tail.qubits - 4 * n * lanes
    return Stats(
        qubits=address_bits + 4 * n * lanes + max(qrom_scratch, tail_scratch),
        toffoli=sum(part.toffoli for part in pairs) + tail.toffoli,
        cnot=sum(part.cnot for part in pairs) + tail.cnot,
        full_depth=max(part.full_depth for part in pairs) + tail.full_depth,
        current_depth=max(part.current_depth for part in pairs) + tail.current_depth,
        toffoli_depth=max(part.toffoli_depth for part in pairs) + tail.toffoli_depth,
    )


def evaluate_design(
    n: int,
    s: int,
    p: int,
    reductions: dict[int, Stats],
    word_weight: int,
) -> dict[str, object]:
    tables = scheduled_tables(n, s)
    addends = len(tables)
    classical_entries = sum(table.entries for table in tables)
    required = required_reduction_lanes(addends, p)
    missing = [lanes for lanes in required if lanes not in reductions]
    base: dict[str, object] = {
        "n": n,
        "scalar_bits": n + 1,
        "s": s,
        "p": p,
        "addends": addends,
        "point_additions": addends - 1,
        "required_reduction_lanes": ";".join(map(str, required)),
        "missing_reduction_lanes": ";".join(map(str, missing)),
        "available": "false" if missing else "true",
        "qrom_word_weight": word_weight,
        "classical_table_entries": classical_entries,
        "classical_table_logical_bits": 2 * n * classical_entries,
        "classical_table_storage_bytes": 2 * math.ceil(n / 8) * classical_entries,
        "classical_table_progression_additions": classical_entries - addends,
    }
    if missing:
        return base

    init_tables = tables[:p]
    initialization = _parallel_initialization(n, init_tables, word_weight)
    active_address = sum(table.window_bits for table in init_tables)
    peak = 2 * (n + 1) + initialization.qubits - active_address
    clear = initialization.qubits - active_address - 2 * n * p
    forward_t = initialization.toffoli
    forward_c = initialization.cnot
    forward_d = initialization.full_depth
    forward_cd = initialization.current_depth
    forward_td = initialization.toffoli_depth

    accumulation_layers = 0
    for start in range(p, addends, p):
        batch = tables[start : start + p]
        lanes = len(batch)
        layer = _window_layer(n, batch, word_weight, reductions[lanes])
        address_bits = sum(table.window_bits for table in batch)
        demand = layer.qubits - address_bits - 2 * n * lanes
        if demand > clear:
            peak += demand - clear
            clear = layer.qubits - address_bits - 7 * n * lanes
        else:
            clear = layer.qubits - address_bits - 7 * n * lanes + (clear - demand)
        forward_t += layer.toffoli
        forward_c += layer.cnot
        forward_d += layer.full_depth
        forward_cd += layer.current_depth
        forward_td += layer.toffoli_depth
        accumulation_layers += 1

    tree_layers = 0
    items = p
    while items > 1:
        lanes = items // 2
        layer = reductions[lanes]
        demand = layer.qubits - 4 * n * lanes
        if demand > clear:
            peak += demand - clear
            clear = layer.qubits - 8 * n * lanes
        else:
            clear = layer.qubits - 8 * n * lanes + (clear - demand)
        forward_t += layer.toffoli
        forward_c += layer.cnot
        forward_d += layer.full_depth
        forward_cd += layer.current_depth
        forward_td += layer.toffoli_depth
        tree_layers += 1
        items = lanes + items % 2

    width = peak + 2 * n
    toffoli = 2 * forward_t
    cnot = 2 * forward_c + 2 * n
    toffoli_depth = 2 * forward_td
    nct_depth = 2 * forward_d + 2
    current_depth = 2 * forward_cd + 2
    base.update(
        {
            "window_accumulation_layers": accumulation_layers,
            "tree_layers": tree_layers,
            "toffoli": toffoli,
            "cnot": cnot,
            "width": width,
            "toffoli_depth": toffoli_depth,
            "nct_depth_barrier": nct_depth,
            "current_depth_diagnostic": current_depth,
            "dw_barrier": nct_depth * width,
            "tdw": toffoli_depth * width,
        }
    )
    return base


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    preferred = [
        "n", "scalar_bits", "s", "p", "addends", "point_additions",
        "available", "required_reduction_lanes", "missing_reduction_lanes",
        "qrom_word_weight", "classical_table_entries",
        "classical_table_logical_bits", "classical_table_storage_bytes",
        "classical_table_progression_additions",
        "window_accumulation_layers", "tree_layers",
        "toffoli", "cnot", "width", "toffoli_depth", "nct_depth_barrier",
        "current_depth_diagnostic", "dw_barrier", "tdw",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=preferred, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _objective(rows: list[dict[str, object]], field: str) -> dict[str, object]:
    row = min(rows, key=lambda item: int(item[field]))
    return {key: row[key] for key in (
        "n", "s", "p", "addends", "toffoli", "cnot", "width",
        "toffoli_depth", "nct_depth_barrier", "dw_barrier", "tdw",
        "classical_table_entries", "classical_table_logical_bits",
        "classical_table_storage_bytes", "classical_table_progression_additions",
    )}


def _write_minima_csv(
    path: Path,
    model: str,
    per_field: dict[str, object],
    field_sizes: Iterable[int],
) -> None:
    rows: list[dict[str, object]] = []
    objective_names = (
        "minimum_dw_barrier",
        "minimum_toffoli",
        "minimum_tdw",
    )
    for n in field_sizes:
        field = per_field[str(n)]
        objectives = field["objectives_over_available_points"]  # type: ignore[index]
        for objective in objective_names:
            rows.append(
                {
                    "qrom_model": model,
                    "objective": objective,
                    **objectives[objective],  # type: ignore[index]
                }
            )
    fields = (
        "qrom_model", "n", "objective", "s", "p", "addends",
        "toffoli", "cnot", "width", "toffoli_depth", "nct_depth_barrier",
        "dw_barrier", "tdw", "classical_table_entries",
        "classical_table_logical_bits", "classical_table_storage_bytes",
        "classical_table_progression_additions",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(
    results: Path,
    output_dir: Path,
    models: Iterable[str] = ("upper",),
    *,
    field_sizes: Iterable[int] = FIELD_SIZES,
    s_min: int = DEFAULT_S_MIN,
    s_max: int = DEFAULT_S_MAX,
) -> dict[str, object]:
    selected_fields = tuple(sorted(set(field_sizes)))
    if not selected_fields or any(n not in FIELD_SIZES for n in selected_fields):
        raise ValueError(f"field sizes must be selected from {FIELD_SIZES}")
    if not DEFAULT_S_MIN <= s_min <= s_max <= DEFAULT_S_MAX:
        raise ValueError(
            f"window range must satisfy {DEFAULT_S_MIN} <= s_min <= s_max <= {DEFAULT_S_MAX}"
        )
    started = time.perf_counter()
    reductions_by_n: dict[int, dict[int, Stats]] = {}
    sources_by_n: dict[int, list[str]] = {}
    for n in selected_fields:
        reductions_by_n[n], sources_by_n[n] = load_reductions(results, n)

    model_summaries: dict[str, object] = {}
    for model in models:
        if model == "upper":
            weight = lambda n: 2 * n
            description = "table-independent unary-QROM upper bound: h_u=2n for every row"
        else:
            raise ValueError(f"unknown QROM model: {model}")

        all_rows: list[dict[str, object]] = []
        per_field: dict[str, object] = {}
        for n in selected_fields:
            required_lanes = required_grid_reduction_lanes(n, s_min, s_max)
            missing_lanes = [
                lanes for lanes in required_lanes if lanes not in reductions_by_n[n]
            ]
            field_rows = [
                evaluate_design(n, s, p, reductions_by_n[n], weight(n))
                for s in range(s_min, s_max + 1)
                for p in range(1, len(scheduled_tables(n, s)) + 1)
            ]
            all_rows.extend(field_rows)
            available = [row for row in field_rows if row["available"] == "true"]
            per_field[str(n)] = {
                "grid_points": len(field_rows),
                "evaluated_points": len(available),
                "unavailable_points": len(field_rows) - len(available),
                "required_reduction_lanes": required_lanes,
                "available_reduction_lanes": sorted(reductions_by_n[n]),
                "missing_reduction_lanes": missing_lanes,
                "source_files": sources_by_n[n],
                "objectives_over_available_points": {
                    "minimum_dw_barrier": _objective(available, "dw_barrier"),
                    "minimum_tdw": _objective(available, "tdw"),
                    "minimum_toffoli": _objective(available, "toffoli"),
                    "minimum_nct_depth_barrier": _objective(available, "nct_depth_barrier"),
                    "minimum_width": _objective(available, "width"),
                },
            }
        csv_path = output_dir / f"windowed_theory_{model}.csv"
        _write_csv(csv_path, all_rows)
        minima_path = output_dir / f"windowed_theory_{model}_minima.csv"
        _write_minima_csv(minima_path, model, per_field, selected_fields)
        model_summaries[model] = {
            "qrom_model": description,
            "csv": str(csv_path.resolve()),
            "minima_csv": str(minima_path.resolve()),
            "rows": len(all_rows),
            "fields": per_field,
        }

    elapsed = time.perf_counter() - started
    summary: dict[str, object] = {
        "schema": "ac-window-existing-layers-theory-scan-v1",
        "elapsed_seconds": elapsed,
        "gate_generation": False,
        "complete_grid_enumeration": True,
        "field_sizes": list(selected_fields),
        "window_range": [s_min, s_max],
        "resource_rows_only_when_all_required_reduction_layers_exist": True,
        "nct_depth_status": "module-barrier composition; not a monolithic per-line ASAP DAG",
        "fixed_table_status": (
            "QROM CNOT/depth are theoretical hamming-weight models, not fixed-seed table-exact values"
        ),
        "models": model_summaries,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.json"
    summary["summary"] = str(summary_path.resolve())
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--model",
        action="append",
        choices=("upper",),
        help="the coherent upper-weight model reported in the paper",
    )
    parser.add_argument(
        "--n", type=int, choices=FIELD_SIZES, action="append",
        help="field size to scan; repeat as needed; default scans all four",
    )
    parser.add_argument("--s-min", type=int, default=DEFAULT_S_MIN)
    parser.add_argument("--s-max", type=int, default=DEFAULT_S_MAX)
    args = parser.parse_args(argv)
    models = args.model or ["upper"]
    summary = run(
        args.results.resolve(),
        args.output_dir.resolve(),
        models,
        field_sizes=args.n or FIELD_SIZES,
        s_min=args.s_min,
        s_max=args.s_max,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
