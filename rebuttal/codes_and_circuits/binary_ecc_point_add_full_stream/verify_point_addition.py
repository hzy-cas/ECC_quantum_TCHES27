"""Computational-basis verification of both JSB+25 point-addition streams.

"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, MutableSequence, Sequence, Tuple

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


def gf_multiply(left: int, right: int, n: int, polynomial_tail: int) -> int:
    """Polynomial-basis multiplication modulo x^n+polynomial_tail."""

    mask = (1 << n) - 1
    result = 0
    multiplicand = left & mask
    multiplier = right & mask
    while multiplier:
        if multiplier & 1:
            result ^= multiplicand
        multiplier >>= 1
        carry = (multiplicand >> (n - 1)) & 1
        multiplicand = (multiplicand << 1) & mask
        if carry:
            multiplicand ^= polynomial_tail
    return result & mask


def gf_inverse(value: int, n: int, polynomial_tail: int) -> int:
    if value == 0:
        raise ZeroDivisionError("zero has no field inverse")
    exponent = (1 << n) - 2
    result = 1
    base = value
    while exponent:
        if exponent & 1:
            result = gf_multiply(result, base, n, polynomial_tail)
        exponent >>= 1
        if exponent:
            base = gf_multiply(base, base, n, polynomial_tail)
    return result


@dataclass(frozen=True)
class ClassicalPointAdd:
    x1: int
    y1: int
    denominator: int
    numerator: int
    inverse: int
    lambda_value: int
    x3: int
    y3: int


def classical_point_add(build: CircuitBuild) -> ClassicalPointAdd:
    """Create one deterministic non-exceptional test vector."""

    cfg = build.config
    tail = cfg.polynomial_tail
    # Fixed differences exercise several low polynomial coefficients while
    # remaining valid for all four supported field sizes.
    for denominator, numerator in ((0x53, 0xA7), (0xD3, 0x6D), (0x1D, 0xB5)):
        x1 = cfg.x2 ^ denominator
        y1 = cfg.y2 ^ numerator
        inverse = gf_inverse(denominator, cfg.n, tail)
        lambda_value = gf_multiply(numerator, inverse, cfg.n, tail)
        lambda_square = gf_multiply(lambda_value, lambda_value, cfg.n, tail)
        x3 = lambda_square ^ lambda_value ^ x1 ^ cfg.x2 ^ 3
        y3 = (
            gf_multiply(lambda_value, cfg.x2 ^ x3, cfg.n, tail)
            ^ cfg.y2
            ^ x3
        )
        if x3 != cfg.x2:
            return ClassicalPointAdd(
                x1,
                y1,
                denominator,
                numerator,
                inverse,
                lambda_value,
                x3,
                y3,
            )
    raise AssertionError("failed to construct a non-exceptional test vector")


def write_register(
    state: MutableSequence[int], wires: Sequence[int], value: int
) -> None:
    for bit, wire in enumerate(wires):
        state[wire] = (value >> bit) & 1


def read_register(state: Sequence[int], wires: Sequence[int]) -> int:
    result = 0
    for bit, wire in enumerate(wires):
        if state[wire]:
            result |= 1 << bit
    return result


def _register_bits(state: Sequence[int], wires: Sequence[int]) -> bytes:
    return bytes(state[wire] for wire in wires)


def _flatten(registers: Iterable[Sequence[int]]) -> List[int]:
    return [wire for register in registers for wire in register]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _require_value(
    state: Sequence[int], wires: Sequence[int], expected: int, name: str
) -> None:
    actual = read_register(state, wires)
    _require(
        actual == expected,
        f"{name}: expected {hex(expected)}, got {hex(actual)}",
    )


def _require_clean(
    state: Sequence[int], wires: Iterable[int], name: str
) -> None:
    dirty = [wire for wire in wires if state[wire]]
    _require(not dirty, f"{name}: {len(dirty)} qubits remain nonzero")


def _require_snapshot(
    state: Sequence[int], snapshots: Mapping[str, Tuple[Sequence[int], bytes]]
) -> None:
    for name, (wires, expected) in snapshots.items():
        actual = _register_bits(state, wires)
        _require(actual == expected, f"{name} changed after it was retained")


def _new_input_state(
    build: CircuitBuild, sample: ClassicalPointAdd, q: int, extra: int = 0
) -> bytearray:
    state = bytearray(build.width + extra)
    state[build.layout["q"]] = q  # type: ignore[index]
    write_register(state, build.layout["x"], sample.x1)  # type: ignore[arg-type]
    write_register(state, build.layout["y"], sample.y1)  # type: ignore[arg-type]
    return state


def _in_place_workspace(build: CircuitBuild) -> List[int]:
    layout = build.layout
    return _flatten(
        [
            layout["karatsuba_ancillas"],  # type: ignore[list-item]
            layout["square0"],  # type: ignore[list-item]
            layout["square1"],  # type: ignore[list-item]
            *layout["inversion_products"],  # type: ignore[misc]
            layout["arithmetic_product"],  # type: ignore[list-item]
        ]
    )


def _verify_in_place_stages(
    build: CircuitBuild, sample: ClassicalPointAdd, q: int
) -> bytearray:
    layout = build.layout
    checkpoints = build.checkpoints
    state = _new_input_state(build, sample, q)
    previous = checkpoints["start"]

    gates.gm.simulate_range(state, previous, checkpoints["first_division_clean"])
    previous = checkpoints["first_division_clean"]
    numerator = sample.y1 ^ (build.config.y2 if q else 0)
    expected_lambda = gf_multiply(
        numerator,
        gf_inverse(
            sample.denominator, build.n, build.config.polynomial_tail
        ),
        build.n,
        build.config.polynomial_tail,
    )
    _require_value(state, layout["lambda"], expected_lambda, "FLT-in lambda")  # type: ignore[arg-type]
    _require_clean(state, _in_place_workspace(build), "first DIV workspace")

    gates.gm.simulate_range(
        state, previous, checkpoints["first_multiplication_clean"]
    )
    previous = checkpoints["first_multiplication_clean"]
    _require_clean(state, _in_place_workspace(build), "first MUL workspace")

    gates.gm.simulate_range(state, previous, checkpoints["controlled_x_update"])
    previous = checkpoints["controlled_x_update"]
    _require_clean(state, _in_place_workspace(build), "controlled update workspace")

    gates.gm.simulate_range(
        state, previous, checkpoints["second_multiplication_clean"]
    )
    previous = checkpoints["second_multiplication_clean"]
    _require_clean(state, _in_place_workspace(build), "second MUL workspace")

    gates.gm.simulate_range(state, previous, checkpoints["second_division_clean"])
    previous = checkpoints["second_division_clean"]
    _require_clean(state, _in_place_workspace(build), "second DIV workspace")
    _require_value(state, layout["lambda"], 0, "final FLT-in lambda")  # type: ignore[arg-type]

    gates.gm.simulate_range(state, previous, checkpoints["final"])
    expected_x = sample.x3 if q else sample.x1
    expected_y = sample.y3 if q else sample.y1
    _require_value(state, layout["x"], expected_x, f"FLT-in x(q={q})")  # type: ignore[arg-type]
    _require_value(state, layout["y"], expected_y, f"FLT-in y(q={q})")  # type: ignore[arg-type]
    _require(state[layout["q"]] == q, f"FLT-in changed q={q}")  # type: ignore[index]

    logical = {layout["q"], *layout["x"], *layout["y"]}  # type: ignore[misc]
    _require_clean(
        state,
        (wire for wire in range(build.width) if wire not in logical),
        f"FLT-in final ancillas(q={q})",
    )
    return state


def _snapshot_registers(
    state: Sequence[int], named_registers: Mapping[str, Sequence[int]]
) -> Dict[str, Tuple[Sequence[int], bytes]]:
    return {
        name: (wires, _register_bits(state, wires))
        for name, wires in named_registers.items()
    }


def _verify_out_of_place_stages(
    build: CircuitBuild, sample: ClassicalPointAdd
) -> Tuple[bytearray, Dict[str, Tuple[Sequence[int], bytes]]]:
    layout = build.layout
    checkpoints = build.checkpoints
    state = _new_input_state(build, sample, 0)
    previous = checkpoints["start"]

    gates.gm.simulate_range(state, previous, checkpoints["input_prepared"])
    previous = checkpoints["input_prepared"]
    _require_value(state, layout["x"], sample.x1, "FLT-out prepared x")  # type: ignore[arg-type]
    _require_value(state, layout["y"], sample.numerator, "FLT-out numerator")  # type: ignore[arg-type]
    _require_value(
        state, layout["candidate_x"], sample.denominator, "FLT-out denominator"  # type: ignore[arg-type]
    )
    initially_clean = _flatten(
        [
            layout["karatsuba_ancillas"],  # type: ignore[list-item]
            layout["square0"],  # type: ignore[list-item]
            layout["square1"],  # type: ignore[list-item]
            *layout["inversion_products"],  # type: ignore[misc]
            layout["arithmetic_product"],  # type: ignore[list-item]
            layout["candidate_y_product"],  # type: ignore[list-item]
        ]
    )
    _require_clean(state, initially_clean, "FLT-out pre-inversion workspace")

    gates.gm.simulate_range(state, previous, checkpoints["inversion_ready"])
    previous = checkpoints["inversion_ready"]
    _require_value(state, layout["inverse"], sample.inverse, "FLT-out inverse")  # type: ignore[arg-type]
    inverse_wires = set(layout["inverse"])  # type: ignore[arg-type]
    unused_square_wires = [
        wire
        for wire in (*layout["square0"], *layout["square1"])  # type: ignore[misc]
        if wire not in inverse_wires
    ]
    _require_clean(state, unused_square_wires, "unused inversion square workspace")
    _require_clean(
        state,
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
        "inversion Karatsuba workspace",
    )
    for index, product in enumerate(layout["inversion_products"]):  # type: ignore[union-attr]
        _require(
            any(state[wire] for wire in product),
            f"inversion product block {index} unexpectedly remains all-zero",
        )
    last_product = read_register(state, layout["inverse_last_product"])  # type: ignore[arg-type]
    _require(
        gf_multiply(
            last_product,
            last_product,
            build.n,
            build.config.polynomial_tail,
        )
        == sample.inverse,
        "the final retained inversion product does not square to the inverse",
    )
    inversion_snapshots = _snapshot_registers(
        state,
        {
            f"inversion_product_{index}": product
            for index, product in enumerate(layout["inversion_products"])  # type: ignore[union-attr]
        },
    )
    _require_clean(
        state,
        [
            *layout["arithmetic_product"],  # type: ignore[misc]
            *layout["candidate_y_product"],  # type: ignore[misc]
        ],
        "post-inversion unused products",
    )

    gates.gm.simulate_range(state, previous, checkpoints["lambda_ready"])
    previous = checkpoints["lambda_ready"]
    _require_value(state, layout["lambda"], sample.lambda_value, "FLT-out lambda")  # type: ignore[arg-type]
    _require_snapshot(state, inversion_snapshots)
    _require_clean(
        state,
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
        "lambda multiplication Karatsuba workspace",
    )
    arithmetic_snapshot = _snapshot_registers(
        state,
        {"arithmetic_product": layout["arithmetic_product"]},  # type: ignore[dict-item]
    )

    gates.gm.simulate_range(state, previous, checkpoints["inverse_cleared"])
    previous = checkpoints["inverse_cleared"]
    _require_value(state, layout["y"], sample.y1, "restored FLT-out y")  # type: ignore[arg-type]
    _require_clean(
        state,
        [*layout["square0"], *layout["square1"]],  # type: ignore[misc]
        "cleared FLT-out inverse/square workspace",
    )
    _require_snapshot(state, {**inversion_snapshots, **arithmetic_snapshot})

    gates.gm.simulate_range(state, previous, checkpoints["x_relation_ready"])
    previous = checkpoints["x_relation_ready"]
    _require_value(
        state,
        layout["candidate_x"],  # type: ignore[arg-type]
        build.config.x2 ^ sample.x3,
        "candidate x2+x3",
    )

    gates.gm.simulate_range(state, previous, checkpoints["y_product_ready"])
    previous = checkpoints["y_product_ready"]
    expected_y_product = gf_multiply(
        sample.lambda_value,
        build.config.x2 ^ sample.x3,
        build.n,
        build.config.polynomial_tail,
    )
    _require_value(
        state, layout["candidate_y"], expected_y_product, "candidate y product"  # type: ignore[arg-type]
    )
    _require_clean(
        state,
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
        "candidate-y Karatsuba workspace",
    )
    candidate_y_wires = set(layout["candidate_y"])  # type: ignore[arg-type]
    candidate_y_tail = [
        wire
        for wire in layout["candidate_y_product"]  # type: ignore[union-attr]
        if wire not in candidate_y_wires
    ]
    retained_snapshots = {
        **inversion_snapshots,
        **arithmetic_snapshot,
        **_snapshot_registers(
            state,
            {"candidate_y_product_tail": candidate_y_tail},
        ),
    }

    gates.gm.simulate_range(state, previous, checkpoints["candidate_ready"])
    _require_value(state, layout["candidate_x"], sample.x3, "candidate x3")  # type: ignore[arg-type]
    _require_value(state, layout["candidate_y"], sample.y3, "candidate y3")  # type: ignore[arg-type]
    _require_snapshot(state, retained_snapshots)
    return state, retained_snapshots


def _verify_out_of_place_control(
    build: CircuitBuild,
    sample: ClassicalPointAdd,
    candidate_state: Sequence[int],
    retained_snapshots: Mapping[str, Tuple[Sequence[int], bytes]],
) -> Dict[int, bytearray]:
    layout = build.layout
    result = {}
    for q in (0, 1):
        state = bytearray(candidate_state)
        state[layout["q"]] = q  # type: ignore[index]
        gates.gm.simulate_range(
            state,
            build.checkpoints["candidate_ready"],
            build.checkpoints["final"],
        )
        expected_main = (sample.x3, sample.y3) if q else (sample.x1, sample.y1)
        expected_candidate = (
            (sample.x1, sample.y1) if q else (sample.x3, sample.y3)
        )
        _require_value(state, layout["x"], expected_main[0], f"FLT-out x(q={q})")  # type: ignore[arg-type]
        _require_value(state, layout["y"], expected_main[1], f"FLT-out y(q={q})")  # type: ignore[arg-type]
        _require_value(
            state,
            layout["candidate_x"],  # type: ignore[arg-type]
            expected_candidate[0],
            f"FLT-out candidate_x(q={q})",
        )
        _require_value(
            state,
            layout["candidate_y"],  # type: ignore[arg-type]
            expected_candidate[1],
            f"FLT-out candidate_y(q={q})",
        )
        _require(state[layout["q"]] == q, f"FLT-out changed q={q}")  # type: ignore[index]
        _require_clean(
            state,
            [
                *layout["karatsuba_ancillas"],  # type: ignore[misc]
                *layout["square0"],  # type: ignore[misc]
                *layout["square1"],  # type: ignore[misc]
            ],
            f"FLT-out clean temporary workspace(q={q})",
        )
        _require_snapshot(state, retained_snapshots)
        result[q] = state
    return result


def _verify_compute_copy_uncompute(
    build: CircuitBuild, sample: ClassicalPointAdd, q: int
) -> None:
    n = build.n
    layout = build.layout
    state = _new_input_state(build, sample, q, extra=2 * n)
    gates.gm.simulate_range(state, 0, build.checkpoints["final"])

    result_x = list(range(build.width, build.width + n))
    result_y = list(range(build.width + n, build.width + 2 * n))
    for source, target in zip(layout["x"], result_x):  # type: ignore[union-attr]
        state[target] ^= state[source]
    for source, target in zip(layout["y"], result_y):  # type: ignore[union-attr]
        state[target] ^= state[source]

    gates.gm.simulate_range(
        state, 0, build.checkpoints["final"], reverse=True
    )
    _require_value(state, layout["x"], sample.x1, "round-trip input x")  # type: ignore[arg-type]
    _require_value(state, layout["y"], sample.y1, "round-trip input y")  # type: ignore[arg-type]
    _require(state[layout["q"]] == q, "round-trip changed q")  # type: ignore[index]
    logical_input = {layout["q"], *layout["x"], *layout["y"]}  # type: ignore[misc]
    _require_clean(
        state,
        (wire for wire in range(build.width) if wire not in logical_input),
        f"{build.mode} compute-copy-uncompute workspace(q={q})",
    )
    expected_x = sample.x3 if q else sample.x1
    expected_y = sample.y3 if q else sample.y1
    _require_value(state, result_x, expected_x, "copied result x")
    _require_value(state, result_y, expected_y, "copied result y")


def _verify_optimized_action(
    build: CircuitBuild,
    sample: ClassicalPointAdd,
    raw_states: Mapping[int, Sequence[int]],
) -> None:
    get_exact_resources_optimized(gates.gm, build.width)
    for q, expected in raw_states.items():
        state = _new_input_state(build, sample, q)
        gates.gm.simulate(state)
        _require(
            state == expected,
            f"global cancellation changed the full {build.mode} state for q={q}",
        )


def verify_mode(
    n: int,
    mode: str,
    *,
    check_optimized: bool = True,
) -> None:
    selected_mode = parse_mode(mode)
    build = build_point_addition(n, selected_mode)
    sample = classical_point_add(build)

    if selected_mode == IN_PLACE:
        raw_states = {
            q: _verify_in_place_stages(build, sample, q) for q in (0, 1)
        }
    else:
        candidate_state, retained = _verify_out_of_place_stages(build, sample)
        raw_states = _verify_out_of_place_control(
            build, sample, candidate_state, retained
        )

    for q in (0, 1):
        _verify_compute_copy_uncompute(build, sample, q)
    if check_optimized:
        _verify_optimized_action(build, sample, raw_states)

    toffoli, cnot, x_count = gates.gm.get_stats()
    suffix = ", optimized action preserved" if check_optimized else ""
    print(
        f"PASS n={n} {selected_mode}: q=0/1, stage states, cleanup, "
        f"compute-copy-uncompute{suffix}; "
        f"live stream X={x_count:,}, CNOT={cnot:,}, T={toffoli:,}"
    )


def run_verification(
    n: int,
    mode: str = "both",
    *,
    check_optimized: bool = True,
) -> None:
    if n not in SUPPORTED_FIELDS:
        raise ValueError(f"unsupported n={n}; choose one of {SUPPORTED_FIELDS}")
    modes = (IN_PLACE, OUT_OF_PLACE) if mode == "both" else (parse_mode(mode),)
    for selected_mode in modes:
        verify_mode(
            n,
            selected_mode,
            check_optimized=check_optimized,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=SUPPORTED_FIELDS, default=163)
    parser.add_argument(
        "--mode",
        choices=(IN_PLACE, OUT_OF_PLACE, "both"),
        default="both",
    )
    parser.add_argument(
        "--skip-optimized",
        action="store_true",
        help="skip the expensive global-cancellation equivalence check",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_verification(
        args.n,
        args.mode,
        check_optimized=not args.skip_optimized,
    )


if __name__ == "__main__":
    main()
