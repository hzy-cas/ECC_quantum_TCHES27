"""Resources of one paper-reordered in-place or out-of-place point addition."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from init import basic_gates as gates
from init.config import SUPPORTED_FIELDS
from init.stats_utils import get_exact_resources_optimized
from point_addition import (
    IN_PLACE,
    OUT_OF_PLACE,
    CircuitBuild,
    build_point_addition,
    parse_mode,
)


@dataclass(frozen=True)
class StageResources:
    """Same resource tuple used by ``recursive_karatsuba`` Python code."""

    width: int
    toffoli_count: int
    cnot_count: int
    full_depth: int
    current_depth: int
    toffoli_depth: int


def _measure_current_build(build: CircuitBuild) -> StageResources:
    stats, full_depth, current_depth, toffoli_depth = (
        get_exact_resources_optimized(gates.gm, build.width)
    )
    return StageResources(
        width=build.width,
        toffoli_count=stats["Toffoli_count"],
        cnot_count=stats["CNOT_count"],
        full_depth=full_depth,
        current_depth=current_depth,
        toffoli_depth=toffoli_depth,
    )


def inplace_point_addition_resources(target_n: int) -> StageResources:
    return _measure_current_build(build_point_addition(target_n, IN_PLACE))


def outofplace_point_addition_resources(target_n: int) -> StageResources:
    return _measure_current_build(build_point_addition(target_n, OUT_OF_PLACE))


def point_addition_resources(target_n: int, mode: str) -> StageResources:
    selected_mode = parse_mode(mode)
    if selected_mode == IN_PLACE:
        return inplace_point_addition_resources(target_n)
    return outofplace_point_addition_resources(target_n)


def _resource_row(
    target_n: int, mode: str, resource: StageResources
) -> Dict[str, object]:
    return {"n": target_n, "mode": mode, **asdict(resource)}


def _print_resources(
    target_n: int, mode: str, resource: StageResources
) -> str:
    return "\n".join(
        (
            f"GF(2^{target_n}) {mode}",
            f"  Toffoli Gates : {resource.toffoli_count:,}",
            f"  CNOT Gates     : {resource.cnot_count:,}",
            f"  Full Depth     : {resource.full_depth:,}",
            f"  Current Depth  : {resource.current_depth:,}",
            f"  Toffoli Depth  : {resource.toffoli_depth:,}",
            f"  Qubits         : {resource.width:,}",
        )
    )


def _write_csv(rows: Iterable[Dict[str, object]], stream) -> None:
    entries = list(rows)
    if not entries:
        return
    writer = csv.DictWriter(stream, fieldnames=list(entries[0]))
    writer.writeheader()
    writer.writerows(entries)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--n", type=int, choices=SUPPORTED_FIELDS, default=163)
    selection.add_argument("--all", action="store_true")
    parser.add_argument(
        "--mode", choices=(IN_PLACE, OUT_OF_PLACE, "both"), default="both"
    )
    parser.add_argument(
        "--format", choices=("text", "json", "csv"), default="text"
    )
    parser.add_argument("--output")
    args = parser.parse_args(argv)

    sizes = SUPPORTED_FIELDS if args.all else (args.n,)
    modes = (
        (IN_PLACE, OUT_OF_PLACE)
        if args.mode == "both"
        else (parse_mode(args.mode),)
    )
    entries: List[tuple[int, str, StageResources]] = []
    for target_n in sizes:
        for mode in modes:
            entries.append(
                (target_n, mode, point_addition_resources(target_n, mode))
            )

    stream = (
        open(args.output, "w", encoding="utf-8", newline="")
        if args.output
        else sys.stdout
    )
    try:
        if args.format == "json":
            json.dump(
                [_resource_row(*entry) for entry in entries],
                stream,
                indent=2,
                ensure_ascii=False,
            )
            stream.write("\n")
        elif args.format == "csv":
            _write_csv((_resource_row(*entry) for entry in entries), stream)
        else:
            stream.write(
                "\n\n".join(_print_resources(*entry) for entry in entries)
            )
            stream.write("\n")
    finally:
        if args.output:
            stream.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
