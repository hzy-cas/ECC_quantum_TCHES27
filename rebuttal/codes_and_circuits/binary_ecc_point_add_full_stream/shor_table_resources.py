"""Complete-Shor estimate composed from the paper-reordered point stream."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from init.config import SUPPORTED_FIELDS, get_config
from point_addition import IN_PLACE, OUT_OF_PLACE, parse_mode
from point_addition_resources import StageResources, point_addition_resources


SHOR_WIDTHS = {
    (163, IN_PLACE): 53_134,
    (163, OUT_OF_PLACE): 15_847_322,
    (233, IN_PLACE): 82_899,
    (233, OUT_OF_PLACE): 35_588_718,
    (283, IN_PLACE): 144_672,
    (283, OUT_OF_PLACE): 75_930_486,
    (571, IN_PLACE): 500_450,
    (571, OUT_OF_PLACE): 535_201_622,
}


@dataclass(frozen=True)
class ShorEstimate:
    field_degree: int
    mode: str
    point_additions: int
    resource_multiplier: int
    result_copy_cnot: int
    toffoli_count: int
    cnot_count: int
    full_depth: int
    current_depth: int
    toffoli_depth: int
    qubits: int
    dw_cost: int
    dw_power: int
    dw_coefficient: float
    cleanup: str
    width_model: str

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


def _coherent_width(target_n: int, point: StageResources, mode: str) -> int:
    if mode == IN_PLACE:
        return point.width
    config = get_config(target_n)
    point_additions = 2 * target_n + 2
    persistent_per_addition = target_n + (
        config.inversion_multiplications + 2
    ) * config.multiplication_targets
    return (
        point_additions
        + 6 * target_n
        + config.karatsuba_ancillas
        + point_additions * persistent_per_addition
    )


def estimate_shor_resources(
    target_n: int,
    mode: str,
    width_model: str = "table6",
) -> ShorEstimate:
    """Match ``recursive_karatsuba`` resource-composition conventions."""

    if width_model not in ("table6", "coherent"):
        raise ValueError("width_model must be table6 or coherent")
    selected_mode = parse_mode(mode)
    point = point_addition_resources(target_n, selected_mode)
    point_additions = 2 * target_n + 2
    resource_multiplier = (
        point_additions if selected_mode == IN_PLACE else 2 * point_additions
    )
    result_copy_cnot = 0 if selected_mode == IN_PLACE else 2 * target_n
    copy_depth = 1 if result_copy_cnot else 0

    toffoli_count = resource_multiplier * point.toffoli_count
    cnot_count = resource_multiplier * point.cnot_count + result_copy_cnot
    full_depth = resource_multiplier * point.full_depth + copy_depth
    current_depth = resource_multiplier * point.current_depth + copy_depth
    toffoli_depth = resource_multiplier * point.toffoli_depth
    qubits = (
        SHOR_WIDTHS[(target_n, selected_mode)]
        if width_model == "table6"
        else _coherent_width(target_n, point, selected_mode)
    )
    dw_cost = qubits * full_depth
    dw_power = int(math.floor(math.log2(dw_cost)))

    return ShorEstimate(
        field_degree=target_n,
        mode=selected_mode,
        point_additions=point_additions,
        resource_multiplier=resource_multiplier,
        result_copy_cnot=result_copy_cnot,
        toffoli_count=toffoli_count,
        cnot_count=cnot_count,
        full_depth=full_depth,
        current_depth=current_depth,
        toffoli_depth=toffoli_depth,
        qubits=qubits,
        dw_cost=dw_cost,
        dw_power=dw_power,
        dw_coefficient=dw_cost / (1 << dw_power),
        cleanup=("none" if selected_mode == IN_PLACE else "compute-copy-uncompute"),
        width_model=width_model,
    )


def _print_estimate(estimate: ShorEstimate) -> str:
    label = "In-place" if estimate.mode == IN_PLACE else "Out-of-place"
    return "\n".join(
        (
            f"GF(2^{estimate.field_degree}) {label}",
            f"  Toffoli Gates : {estimate.toffoli_count:,}",
            f"  CNOT Gates     : {estimate.cnot_count:,}",
            f"  Full Depth     : {estimate.full_depth:,}",
            f"  Current Depth  : {estimate.current_depth:,}",
            f"  Toffoli Depth  : {estimate.toffoli_depth:,}",
            f"  Qubits         : {estimate.qubits:,}",
            f"  DW-cost        : "
            f"{estimate.dw_coefficient:.6f} * 2^{estimate.dw_power}",
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
        "--width-model", choices=("table6", "coherent"), default="table6"
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
    estimates: List[ShorEstimate] = []
    for target_n in sizes:
        for mode in modes:
            estimates.append(
                estimate_shor_resources(target_n, mode, args.width_model)
            )

    stream = (
        open(args.output, "w", encoding="utf-8", newline="")
        if args.output
        else sys.stdout
    )
    try:
        if args.format == "json":
            json.dump(
                [estimate.to_dict() for estimate in estimates],
                stream,
                indent=2,
                ensure_ascii=False,
            )
            stream.write("\n")
        elif args.format == "csv":
            _write_csv((estimate.to_dict() for estimate in estimates), stream)
        else:
            stream.write("\n\n".join(_print_estimate(item) for item in estimates))
            stream.write("\n")
    finally:
        if args.output:
            stream.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
