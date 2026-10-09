"""Exact Toffoli-count search for the maximally parallel windowed design.

For the Toffoli objective, ``p=L_s`` is optimal: moving a table out of the
initialization set replaces one QROM call by a lookup/unlookup pair, while no
choice of fewer lanes can use fewer point-addition layers than the maximally
parallel reduction tree.  The scan over ``s=2..18`` is therefore exhaustive in
both parameters without materializing data-dependent Clifford gates.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .synthesis import scheduled_descriptors, synthesize_tree_layer


@dataclass(frozen=True)
class ToffoliSearchRow:
    n: int
    arithmetic: str
    s: int
    p: int
    addends: int
    qrom_toffoli: int
    arithmetic_toffoli: int
    total_toffoli: int


def _primitive_toffoli(n: int, arithmetic: str) -> tuple[int, int]:
    """Recover exact multiplication/inversion counts from k=1,2 PA streams."""
    one = synthesize_tree_layer(n, arithmetic, 1).toffoli
    two = synthesize_tree_layer(n, arithmetic, 2).toffoli
    multiplication, remainder = divmod(two - one, 10)
    if remainder or (one - 4 * multiplication) % 2:
        raise AssertionError("point-addition streams violate 2I+(10k-6)M")
    inversion = (one - 4 * multiplication) // 2
    return multiplication, inversion


def _tree_toffoli(items: int, multiplication: int, inversion: int) -> int:
    layers = (items - 1).bit_length() if items > 1 else 0
    return 2 * layers * inversion + (10 * (items - 1) - 6 * layers) * multiplication


def scan_toffoli(n: int, arithmetic_name: str) -> list[ToffoliSearchRow]:
    if arithmetic_name != "ac":
        raise ValueError("this submission artifact supports only AC arithmetic")
    rows: list[ToffoliSearchRow] = []
    multiplication, inversion = _primitive_toffoli(n, arithmetic_name)
    for s in range(2, 19):
        tables = scheduled_descriptors(n, s)
        # One forward initialization lookup and its outer-map inverse.
        qrom = 2 * sum(2 * (table.entries - 2) for table in tables)
        arithmetic_count = 2 * _tree_toffoli(len(tables), multiplication, inversion)
        rows.append(
            ToffoliSearchRow(
                n=n,
                arithmetic=arithmetic_name,
                s=s,
                p=len(tables),
                addends=len(tables),
                qrom_toffoli=qrom,
                arithmetic_toffoli=arithmetic_count,
                total_toffoli=qrom + arithmetic_count,
            )
        )
    return rows


def write_toffoli_scan(path: Path, n: int, arithmetic_name: str) -> list[ToffoliSearchRow]:
    rows = scan_toffoli(n, arithmetic_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ToffoliSearchRow.__dataclass_fields__)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)
    return rows
