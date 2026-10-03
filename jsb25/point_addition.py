"""Full NCT streams for JSB+25 in-place and out-of-place point addition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

from init import basic_gates as gates
from init.config import FieldConfig, SUPPORTED_FIELDS, get_config
from init.inversion import Inverison_Itoh_Tsujii_based


IN_PLACE = "inplace"
OUT_OF_PLACE = "outofplace"
SUPPORTED_MODES = (IN_PLACE, OUT_OF_PLACE)


@dataclass(frozen=True)
class CircuitBuild:
    n: int
    mode: str
    width: int
    config: FieldConfig
    layout: Dict[str, object]
    checkpoints: Dict[str, int]


class QubitAllocator:
    """Allocate monotonically increasing global wire identifiers."""

    def __init__(self) -> None:
        self.offset = 0

    def take(self, size: int) -> List[int]:
        if size < 0:
            raise ValueError("register size must be nonnegative")
        result = list(range(self.offset, self.offset + size))
        self.offset += size
        return result

    def take_one(self) -> int:
        return self.take(1)[0]


def parse_mode(value: str) -> str:
    normalized = value.lower().replace("-", "")
    if normalized in ("in", "inplace"):
        return IN_PLACE
    if normalized in ("out", "outofplace"):
        return OUT_OF_PLACE
    raise ValueError(f"unsupported mode {value!r}; choose {SUPPORTED_MODES}")


def _flatten(registers: Sequence[Sequence[int]]) -> List[int]:
    return [wire for register in registers for wire in register]


def _set_product_targets(registers: Sequence[Sequence[int]]) -> None:
    gates.Toffoli_qubits = _flatten(registers)
    gates.cnt = 0


def _scratch_registers(layout: Dict[str, object]) -> List[List[int]]:
    return [
        layout["karatsuba_ancillas"],  # type: ignore[list-item]
        layout["square0"],  # type: ignore[list-item]
        layout["square1"],  # type: ignore[list-item]
        *_inversion_products(layout),
        layout["arithmetic_product"],  # type: ignore[list-item]
    ]


def _inversion_products(layout: Dict[str, object]) -> List[List[int]]:
    return layout["inversion_products"]  # type: ignore[return-value]


def _division_xor(
    numerator: Sequence[int],
    denominator: Sequence[int],
    destination: Sequence[int],
    layout: Dict[str, object],
) -> List[int]:
    product_targets = _inversion_products(layout) + [
        layout["arithmetic_product"]  # type: ignore[list-item]
    ]
    _set_product_targets(product_targets)
    start = gates.gm.current_pointer()
    inverse = Inverison_Itoh_Tsujii_based(
        denominator,
        layout["square0"],  # type: ignore[arg-type]
        layout["square1"],  # type: ignore[arg-type]
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
        clear_final_square=False,
    )
    quotient = gates.multiply(
        inverse.inverse,
        numerator,
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
    )
    end = gates.gm.current_pointer()
    gates.CNOT_n(quotient, destination)
    gates.gm.replay_reverse(start, end)
    gates.cnt = 0
    return quotient


def _multiplication_xor(
    left: Sequence[int],
    right: Sequence[int],
    destination: Sequence[int],
    layout: Dict[str, object],
) -> None:
    _set_product_targets(
        [layout["arithmetic_product"]]  # type: ignore[list-item]
    )
    start = gates.gm.current_pointer()
    product = gates.multiply(
        left,
        right,
        layout["karatsuba_ancillas"],  # type: ignore[arg-type]
    )
    end = gates.gm.current_pointer()
    gates.CNOT_n(product, destination)
    gates.gm.replay_reverse(start, end)
    gates.cnt = 0


def _controlled_constant_xor(
    control: int,
    target: Sequence[int],
    constant: int,
    fanout_workspace: Sequence[int],
) -> None:
    active_bits = [bit for bit in range(len(target)) if (constant >> bit) & 1]
    if not active_bits:
        return
    start = gates.gm.current_pointer()
    copies = gates.copy_parallel(control, fanout_workspace, len(active_bits))
    end = gates.gm.current_pointer()
    for copied_control, bit in zip(copies, active_bits):
        gates.CNOT(copied_control, target[bit])
    gates.gm.replay_reverse(start, end)


def _build_in_place(layout: Dict[str, object], config: FieldConfig) -> Dict[str, int]:
    q = layout["q"]  # type: ignore[assignment]
    x = layout["x"]  # type: ignore[assignment]
    y = layout["y"]  # type: ignore[assignment]
    lambda_register = layout["lambda"]  # type: ignore[assignment]
    karatsuba = layout["karatsuba_ancillas"]  # type: ignore[assignment]
    checkpoints: Dict[str, int] = {"start": gates.gm.current_pointer()}

    # Algorithm 1, Steps 1--4.
    gates.CONST_ADD_n(x, config.x2)
    _controlled_constant_xor(q, y, config.y2, karatsuba)
    _division_xor(y, x, lambda_register, layout)
    checkpoints["first_division_clean"] = gates.gm.current_pointer()

    # Steps 5--11.
    _multiplication_xor(x, lambda_register, y, layout)
    checkpoints["first_multiplication_clean"] = gates.gm.current_pointer()
    affine_start = gates.gm.current_pointer()
    gates.square_plus_input_xor(lambda_register, y)
    gates.CONST_ADD_n(y, 3)
    gates.CONST_ADD_n(y, config.x2)
    affine_end = gates.gm.current_pointer()

    fanout_start = gates.gm.current_pointer()
    copies = gates.copy_parallel(q, karatsuba, config.n)
    fanout_end = gates.gm.current_pointer()
    for bit in range(config.n):
        gates.Toffoli_gate(copies[bit], y[bit], x[bit])
    gates.gm.replay_reverse(fanout_start, fanout_end)
    gates.gm.replay_reverse(affine_start, affine_end)
    checkpoints["controlled_x_update"] = gates.gm.current_pointer()

    # Steps 12--14.
    _multiplication_xor(x, lambda_register, y, layout)
    checkpoints["second_multiplication_clean"] = gates.gm.current_pointer()
    _division_xor(y, x, lambda_register, layout)
    checkpoints["second_division_clean"] = gates.gm.current_pointer()
    gates.CONST_ADD_n(x, config.x2)

    # Steps 15--17: one n-way control fan-out services both operations.
    final_fanout_start = gates.gm.current_pointer()
    final_copies = gates.copy_parallel(q, karatsuba, config.n)
    final_fanout_end = gates.gm.current_pointer()
    active_y_bits = [bit for bit in range(config.n) if (config.y2 >> bit) & 1]
    for copy_index, bit in enumerate(active_y_bits):
        gates.CNOT(final_copies[copy_index], y[bit])
    for bit in range(config.n):
        gates.Toffoli_gate(final_copies[bit], x[bit], y[bit])
    gates.gm.replay_reverse(final_fanout_start, final_fanout_end)
    checkpoints["final"] = gates.gm.current_pointer()
    return checkpoints


def _build_out_of_place(
    layout: Dict[str, object], config: FieldConfig
) -> Dict[str, int]:
    q = layout["q"]  # type: ignore[assignment]
    x = layout["x"]  # type: ignore[assignment]
    y = layout["y"]  # type: ignore[assignment]
    candidate_x = layout["candidate_x"]  # type: ignore[assignment]
    karatsuba = layout["karatsuba_ancillas"]  # type: ignore[assignment]
    checkpoints: Dict[str, int] = {"start": gates.gm.current_pointer()}

    product_targets = _inversion_products(layout) + [
        layout["arithmetic_product"],  # type: ignore[list-item]
        layout["candidate_y_product"],  # type: ignore[list-item]
    ]
    _set_product_targets(product_targets)

    gates.CNOT_n(x, candidate_x)
    gates.CONST_ADD_n(candidate_x, config.x2)
    gates.CONST_ADD_n(y, config.y2)
    checkpoints["input_prepared"] = gates.gm.current_pointer()

    inverse = Inverison_Itoh_Tsujii_based(
        candidate_x,
        layout["square0"],  # type: ignore[arg-type]
        layout["square1"],  # type: ignore[arg-type]
        karatsuba,
        clear_final_square=True,
    )
    layout["inverse"] = inverse.inverse
    layout["inverse_last_product"] = inverse.last_product
    checkpoints["inversion_ready"] = gates.gm.current_pointer()

    lambda_register = gates.multiply(inverse.inverse, y, karatsuba)
    layout["lambda"] = lambda_register
    checkpoints["lambda_ready"] = gates.gm.current_pointer()

    gates.CONST_ADD_n(y, config.y2)
    # The square of the last inversion product is exactly the inverse.  FLT-out
    # clears this n-bit result while retaining all multiplication products.
    gates.square_xor(inverse.last_product, inverse.inverse, 1)
    checkpoints["inverse_cleared"] = gates.gm.current_pointer()

    gates.CONST_ADD_n(candidate_x, config.x2)
    gates.CONST_ADD_n(candidate_x, 3)
    gates.square_plus_input_xor(lambda_register, candidate_x)
    checkpoints["x_relation_ready"] = gates.gm.current_pointer()

    candidate_y = gates.multiply(lambda_register, candidate_x, karatsuba)
    layout["candidate_y"] = candidate_y
    checkpoints["y_product_ready"] = gates.gm.current_pointer()

    gates.CONST_ADD_n(candidate_y, config.y2)
    gates.CONST_ADD_n(candidate_x, config.x2)
    gates.CNOT_n(candidate_x, candidate_y)
    checkpoints["candidate_ready"] = gates.gm.current_pointer()

    fanout_start = gates.gm.current_pointer()
    copies = gates.copy_parallel(q, karatsuba, 2 * config.n)
    fanout_end = gates.gm.current_pointer()
    for bit in range(config.n):
        gates.CSWAP(copies[bit], x[bit], candidate_x[bit])
        gates.CSWAP(copies[config.n + bit], y[bit], candidate_y[bit])
    gates.gm.replay_reverse(fanout_start, fanout_end)
    checkpoints["final"] = gates.gm.current_pointer()

    expected_targets = (config.inversion_multiplications + 2) * (
        config.multiplication_targets
    )
    if gates.cnt != expected_targets:
        raise AssertionError(
            f"FLT-out consumed {gates.cnt} product targets; expected {expected_targets}"
        )
    return checkpoints


def build_point_addition(n: int, mode: str) -> CircuitBuild:
    """Construct one complete point-addition stream and return its wire map."""

    if n not in SUPPORTED_FIELDS:
        raise ValueError(f"unsupported n={n}; choose one of {SUPPORTED_FIELDS}")
    selected_mode = parse_mode(mode)
    gates.configure_backend(n)
    config = get_config(n)
    allocator = QubitAllocator()

    layout: Dict[str, object] = {
        "q": allocator.take_one(),
        "x": allocator.take(n),
        "y": allocator.take(n),
    }
    if selected_mode == IN_PLACE:
        layout["lambda"] = allocator.take(n)
    else:
        layout["candidate_x"] = allocator.take(n)

    layout["karatsuba_ancillas"] = allocator.take(config.karatsuba_ancillas)
    layout["square0"] = allocator.take(n)
    layout["square1"] = allocator.take(n)
    layout["inversion_products"] = [
        allocator.take(config.multiplication_targets)
        for _ in range(config.inversion_multiplications)
    ]
    layout["arithmetic_product"] = allocator.take(
        config.multiplication_targets
    )

    if selected_mode == OUT_OF_PLACE:
        layout["candidate_y_product"] = allocator.take(
            config.multiplication_targets
        )
        checkpoints = _build_out_of_place(layout, config)
    else:
        checkpoints = _build_in_place(layout, config)

    width = allocator.offset
    expected_width = (
        1
        + 5 * n
        + config.karatsuba_ancillas
        + (
            config.inversion_multiplications
            + (2 if selected_mode == OUT_OF_PLACE else 1)
        )
        * config.multiplication_targets
    )
    if width != expected_width:
        raise AssertionError(f"allocated width {width}; expected {expected_width}")
    return CircuitBuild(n, selected_mode, width, config, layout, checkpoints)


def scratch_registers(build: CircuitBuild) -> List[List[int]]:
    """Return every non-logical register allocated by the selected circuit."""

    registers = _scratch_registers(build.layout)
    if build.mode == OUT_OF_PLACE:
        registers.append(build.layout["candidate_y_product"])  # type: ignore[arg-type]
    return registers
