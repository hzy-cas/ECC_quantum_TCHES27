#!/usr/bin/env python3
"""Estimate one in-place point addition with the bundled analyzer."""

from __future__ import annotations

import argparse
import csv
import gc
import json
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
PACKAGE_PARENT = PACKAGE_ROOT.parent
if not __package__ and str(PACKAGE_PARENT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_PARENT))

if __package__:
    from ..circuit.config import (
        BINARY_ECC,
        CONFIGS,
        INVERSION_MODES,
        SUPPORTED_SIZES,
    )
    from .stats_utils import get_exact_resources_optimized
    from ..verification.verify_point_addition import build_test_circuit
else:
    from algorithm1_inplace_point_add.circuit.config import (
        BINARY_ECC,
        CONFIGS,
        INVERSION_MODES,
        SUPPORTED_SIZES,
    )
    from algorithm1_inplace_point_add.resources.stats_utils import (
        get_exact_resources_optimized,
    )
    from algorithm1_inplace_point_add.verification.verify_point_addition import (
        build_test_circuit,
    )


def low_level_multiplication_calls(n):
    """Two clean divisions and two clean products emit ``4I+8`` multiplications."""

    return 4 * CONFIGS[n]["inversion_multiplications"] + 8


def raw_toffoli_closed_form(n):
    return (
        low_level_multiplication_calls(n)
        * CONFIGS[n]["multiplication_targets"]
        + 2 * n
    )


def estimate_one(n, inversion_mode=BINARY_ECC):
    start = time.perf_counter()
    circuit = build_test_circuit(n, inversion_mode)
    build_seconds = time.perf_counter() - start
    raw_toffoli, raw_cnot, raw_x = circuit["gm"].get_stats()
    expected_toffoli = raw_toffoli_closed_form(n)
    if raw_toffoli != expected_toffoli:
        raise AssertionError(
            f"GF(2^{n}) Toffoli={raw_toffoli}; closed-form expectation={expected_toffoli}"
        )

    start = time.perf_counter()
    stats, full_depth, current_depth, toffoli_depth = (
        get_exact_resources_optimized(circuit["gm"], circuit["width"])
    )
    estimation_seconds = time.perf_counter() - start
    row = {
        "n": n,
        "inversion_mode": inversion_mode,
        "qubits": circuit["width"],
        "raw_cnot": raw_cnot,
        "raw_toffoli": raw_toffoli,
        "cnot": stats["CNOT_count"],
        "toffoli": stats["Toffoli_count"],
        "full_depth": full_depth,
        "current_depth": current_depth,
        "toffoli_depth": toffoli_depth,
        "DW": circuit["width"] * full_depth,
        "TDW": circuit["width"] * toffoli_depth,
        "low_level_multiplications": low_level_multiplication_calls(n),
        "build_seconds": round(build_seconds, 6),
        "estimation_seconds": round(estimation_seconds, 6),
        "estimator": "circuit_backend/python (nct-v2)",
        "multiplication": "AC-based",
        "inversion_schedule": (
            "Binary_ECC sequential Itoh-Tsujii"
            if inversion_mode == BINARY_ECC
            else "parallel AC Toffoli-depth-optimized inversion"
        ),
    }
    circuit["gm"].clear()
    del circuit
    gc.collect()
    return row


def write_json(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        writer.writeheader()
        writer.writerows(rows)


def print_table(rows):
    print(
        "| n | inversion | qubits | Toffoli | CNOT | full depth | "
        "T-depth | DW | TDW |"
    )
    print("|---:|:---|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        print(
            f"| {row['n']} | {row['inversion_mode']} | {row['qubits']} | "
            f"{row['toffoli']} | "
            f"{row['cnot']} | {row['full_depth']} | {row['toffoli_depth']} | "
            f"{row['DW']} | {row['TDW']} |"
        )


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate one AC-based in-place point addition with either inversion backend."
    )
    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        choices=SUPPORTED_SIZES,
        default=[163],
        help="run n=163 by default; large parameter sets should be run separately",
    )
    parser.add_argument(
        "--inversion",
        choices=INVERSION_MODES + ("both",),
        default=BINARY_ECC,
        help="select low width, optimal Toffoli depth, or estimate both in sequence",
    )
    parser.add_argument("--json", type=Path)
    parser.add_argument("--csv", type=Path)
    return parser.parse_args()


def main():
    args = _parse_args()
    rows = []
    modes = INVERSION_MODES if args.inversion == "both" else (args.inversion,)
    for n in args.sizes:
        for inversion_mode in modes:
            print(
                f"Building and analyzing GF(2^{n}), inversion={inversion_mode} ...",
                flush=True,
            )
            rows.append(estimate_one(n, inversion_mode))
    print_table(rows)
    if args.json:
        write_json(args.json, rows)
    if args.csv:
        write_csv(args.csv, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
