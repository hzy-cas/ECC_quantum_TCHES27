#!/usr/bin/env python3
"""Classical reference, full NCT simulation, and ancilla checks for Algorithm 1."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
import gc
import json
from pathlib import Path
import random
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
        DATA_ROOT,
        INVERSION_MODES,
        OPTIMAL_DEPTH,
        SUPPORTED_SIZES,
    )
    from ..circuit.point_addition import Point_addition
    from ..backend.qasm import FlatGateManager, MatrixCache
else:
    from algorithm1_inplace_point_add.circuit.config import (
        BINARY_ECC,
        CONFIGS,
        DATA_ROOT,
        INVERSION_MODES,
        OPTIMAL_DEPTH,
        SUPPORTED_SIZES,
    )
    from algorithm1_inplace_point_add.circuit.point_addition import Point_addition
    from algorithm1_inplace_point_add.backend.qasm import (
        FlatGateManager,
        MatrixCache,
    )


POLYNOMIAL_LOW = {
    163: (1 << 7) | (1 << 6) | (1 << 3) | 1,
    233: (1 << 74) | 1,
    283: (1 << 12) | (1 << 7) | (1 << 5) | 1,
    571: (1 << 10) | (1 << 5) | (1 << 2) | 1,
}


@dataclass(frozen=True)
class AffineFixture:
    n: int
    curve_a: int
    curve_b: int
    x1: int
    y1: int
    x2: int
    y2: int
    slope: int
    x3: int
    y3: int


class BasisTransformer:
    """Classical converter between polynomial basis and the AC multiplier L-basis."""

    def __init__(self, n, l2x_path, x2l_path):
        self.n = n
        self.matrix_l2x = self._load_matrix(l2x_path)
        self.matrix_x2l = self._load_matrix(x2l_path)

    def _load_matrix(self, path):
        matrix = []
        with open(path, "r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line:
                    continue
                if len(line) != self.n or set(line) - {"0", "1"}:
                    raise ValueError(f"{path} contains an invalid matrix row")
                matrix.append([int(bit) for bit in line])
        if len(matrix) != self.n:
            raise ValueError(
                f"{path} must contain {self.n} rows; found {len(matrix)}"
            )
        return matrix

    def apply(self, vec_bits, direction="x2l"):
        matrix = self.matrix_x2l if direction == "x2l" else self.matrix_l2x
        result = [0] * self.n
        for row_index, row in enumerate(matrix):
            value = 0
            for column, coefficient in enumerate(row):
                if coefficient and vec_bits[column]:
                    value ^= 1
            result[row_index] = value
        return result


def int_to_bits(value, n):
    return [(value >> bit) & 1 for bit in range(n)]


def bits_to_int(bits):
    value = 0
    for bit, state in enumerate(bits):
        if state:
            value |= 1 << bit
    return value


@lru_cache(maxsize=None)
def basis_transformer(n):
    base_dir = DATA_ROOT / f"quantum_{n}"
    return BasisTransformer(
        n,
        base_dir / f"{n}_L2X.txt",
        base_dir / f"{n}_X2L.txt",
    )


def poly_to_l(value, n):
    return bits_to_int(basis_transformer(n).apply(int_to_bits(value, n), "x2l"))


def l_to_poly(value, n):
    return bits_to_int(basis_transformer(n).apply(int_to_bits(value, n), "l2x"))


def gf_reduce(value, n):
    low = POLYNOMIAL_LOW[n]
    for degree in range(value.bit_length() - 1, n - 1, -1):
        if (value >> degree) & 1:
            value ^= (1 << degree) ^ (low << (degree - n))
    return value


def gf_mul(left, right, n):
    product = 0
    while right:
        if right & 1:
            product ^= left
        left <<= 1
        right >>= 1
    return gf_reduce(product, n)


def gf_square(value, n):
    expanded = 0
    bit = 0
    while value:
        if value & 1:
            expanded |= 1 << (2 * bit)
        value >>= 1
        bit += 1
    return gf_reduce(expanded, n)


def gf_inv(value, n):
    if value == 0:
        raise ZeroDivisionError("zero has no multiplicative inverse")
    result = 1
    base = value
    exponent = (1 << n) - 2
    while exponent:
        if exponent & 1:
            result = gf_mul(result, base, n)
        exponent >>= 1
        if exponent:
            base = gf_square(base, n)
    return result


def affine_add(x1, y1, x2, y2, curve_a, n):
    denominator = x1 ^ x2
    if denominator == 0:
        raise ZeroDivisionError("Algorithm 1 requires x1 != x2")
    slope = gf_mul(y1 ^ y2, gf_inv(denominator, n), n)
    x3 = gf_square(slope, n) ^ slope ^ x1 ^ x2 ^ curve_a
    y3 = gf_mul(slope, x1 ^ x3, n) ^ x3 ^ y1
    return x3, y3


def _curve_residual(x, y, curve_a, n):
    x_squared = gf_square(x, n)
    return (
        gf_square(y, n)
        ^ gf_mul(x, y, n)
        ^ gf_mul(x_squared, x, n)
        ^ gf_mul(curve_a, x_squared, n)
    )


@lru_cache(maxsize=None)
def deterministic_fixture(n, curve_a=1):
    """Construct fixed same-curve test points that avoid both zero divisions."""

    rng = random.Random(0xA11CE000 + n)
    for _ in range(128):
        x1 = rng.getrandbits(n) or 1
        x2 = rng.getrandbits(n) or 2
        if x1 == x2:
            continue
        slope = rng.getrandbits(n) or 1
        delta_x = x1 ^ x2
        y1 = (
            gf_mul(gf_square(slope, n), delta_x, n)
            ^ gf_mul(x2, slope, n)
            ^ gf_square(x1, n)
            ^ gf_mul(x1, x2, n)
            ^ gf_square(x2, n)
            ^ gf_mul(curve_a, delta_x, n)
        )
        y2 = y1 ^ gf_mul(slope, delta_x, n)
        x3, y3 = affine_add(x1, y1, x2, y2, curve_a, n)
        if x2 == x3:
            continue
        residuals = (
            _curve_residual(x1, y1, curve_a, n),
            _curve_residual(x2, y2, curve_a, n),
            _curve_residual(x3, y3, curve_a, n),
        )
        if residuals[0] == residuals[1] == residuals[2]:
            return AffineFixture(
                n,
                curve_a,
                residuals[0],
                x1,
                y1,
                x2,
                y2,
                slope,
                x3,
                y3,
            )
    raise RuntimeError(f"could not construct a non-exceptional test point over GF(2^{n})")


def algorithm1_classical(q, x1, y1, x2, y2, curve_a, n):
    """Execute the 17 assignments of paper Algorithm 1 one by one."""

    if q not in (0, 1):
        raise ValueError("q must be either 0 or 1")
    x = x1 ^ x2
    y = y1 ^ (y2 if q else 0)
    slope = gf_mul(y, gf_inv(x, n), n)
    y ^= gf_mul(x, slope, n)
    y ^= gf_square(slope, n) ^ slope
    y ^= curve_a ^ x2
    if q:
        x ^= y
    y ^= gf_square(slope, n) ^ slope
    y ^= curve_a ^ x2
    y ^= gf_mul(x, slope, n)
    slope ^= gf_mul(y, gf_inv(x, n), n)
    x ^= x2
    if q:
        y ^= y2
        y ^= x
    return x, y, slope


def build_test_circuit(n, inversion_mode=BINARY_ECC):
    """Allocate wires for tests/resources; ``Point_addition`` emits every quantum gate."""

    if n not in CONFIGS:
        raise ValueError(f"unsupported field degree n={n}")
    if inversion_mode not in INVERSION_MODES:
        raise ValueError(f"unsupported inversion mode {inversion_mode!r}")
    cfg = CONFIGS[n]
    base_dir = DATA_ROOT / f"quantum_{n}"
    path_config = {
        "TD": str(base_dir / "CNOT_mul" / f"result_{n}_TD_seq.txt"),
        "A": str(base_dir / "CNOT_mul" / f"result_{n}_A.txt"),
        "C": str(base_dir / "CNOT_mul" / f"result_{n}_C.txt"),
        "Inv": str(base_dir / "CNOT_mul" / f"result_{n}_CT2D_inv_seq.txt"),
        "mul_dim_1": cfg["mul_dims"][0],
        "mul_dim_2": cfg["mul_dims"][1],
    }
    missing = [
        path
        for path in path_config.values()
        if isinstance(path, str) and not Path(path).is_file()
    ]
    if missing:
        raise FileNotFoundError("missing multiplication matrices: " + ", ".join(missing))

    cursor = 0
    q = cursor
    cursor += 1
    x = list(range(cursor, cursor + n))
    cursor += n
    y = list(range(cursor, cursor + n))
    cursor += n
    lambda_register = list(range(cursor, cursor + n))
    cursor += n
    ancilla = list(range(cursor, cursor + 2 * cfg["block_size"]))
    cursor += len(ancilla)
    side_ancilla = []
    if inversion_mode == OPTIMAL_DEPTH:
        side_ancilla = list(
            range(cursor, cursor + 2 * cfg["block_size"])
        )
        cursor += len(side_ancilla)
    target_count = (
        cfg["inversion_multiplications"] + 1
    ) * cfg["multiplication_targets"]
    toffoli_qubits = list(range(cursor, cursor + target_count))
    cursor += len(toffoli_qubits)
    sqr1 = list(range(cursor, cursor + n))
    cursor += n
    sqr2 = list(range(cursor, cursor + n))
    cursor += n
    sqr3 = []
    sqr4 = []
    if inversion_mode == OPTIMAL_DEPTH:
        sqr3 = list(range(cursor, cursor + n))
        cursor += n
        sqr4 = list(range(cursor, cursor + n))
        cursor += n

    fixture = deterministic_fixture(n)
    gm = FlatGateManager()
    matrix_loader = MatrixCache()
    Point_addition(
        gm,
        matrix_loader,
        q,
        x,
        y,
        lambda_register,
        sqr1,
        sqr2,
        ancilla,
        toffoli_qubits,
        n,
        poly_to_l(fixture.curve_a, n),
        poly_to_l(fixture.x2, n),
        poly_to_l(fixture.y2, n),
        path_config,
        str(base_dir),
        inversion_mode=inversion_mode,
        sqr3=sqr3,
        sqr4=sqr4,
        side_ancilla=side_ancilla,
    )
    return {
        "n": n,
        "inversion_mode": inversion_mode,
        "gm": gm,
        "width": cursor,
        "q": q,
        "x": x,
        "y": y,
        "lambda": lambda_register,
        "ancilla": ancilla,
        "side_ancilla": side_ancilla,
        "toffoli_qubits": toffoli_qubits,
        "sqr1": sqr1,
        "sqr2": sqr2,
        "sqr3": sqr3,
        "sqr4": sqr4,
        "fixture": fixture,
    }


def _set_packed_register(state, wires, values):
    for case, value in enumerate(values):
        mask = 1 << case
        for bit, wire in enumerate(wires):
            if (value >> bit) & 1:
                state[wire] |= mask


def _get_packed_register(state, wires, case):
    value = 0
    mask = 1 << case
    for bit, wire in enumerate(wires):
        if state[wire] & mask:
            value |= 1 << bit
    return value


def simulate_packed(gate_manager, initial_state, case_count):
    """Simulate the q=0 and q=1 basis inputs together as integer bit planes."""

    state = list(initial_state)
    all_case_mask = (1 << case_count) - 1
    mask_arg = (1 << gate_manager.SHIFT_ARG) - 1
    waiting = False
    control_one = control_two = 0
    for chunk in gate_manager.chunks + [gate_manager.current_chunk]:
        for encoded in chunk:
            operation = encoded >> gate_manager.SHIFT_OP
            if waiting:
                state[encoded & mask_arg] ^= state[control_one] & state[control_two]
                waiting = False
            elif operation == gate_manager.OP_NOP:
                continue
            elif operation == gate_manager.OP_X:
                state[encoded & mask_arg] ^= all_case_mask
            elif operation == gate_manager.OP_CNOT:
                control = (encoded >> gate_manager.SHIFT_ARG) & mask_arg
                state[encoded & mask_arg] ^= state[control]
            elif operation == gate_manager.OP_TOFFOLI:
                control_one = (encoded >> gate_manager.SHIFT_ARG) & mask_arg
                control_two = encoded & mask_arg
                waiting = True
            elif operation != gate_manager.OP_DATA:
                raise ValueError(f"unknown gate opcode {operation}")
    if waiting:
        raise ValueError("the gate stream ends with an incomplete Toffoli gate")
    return state


def verify_one(n, inversion_mode=BINARY_ECC, algebra_only=False):
    fixture = deterministic_fixture(n)
    q0 = algorithm1_classical(
        0, fixture.x1, fixture.y1, fixture.x2, fixture.y2, fixture.curve_a, n
    )
    q1 = algorithm1_classical(
        1, fixture.x1, fixture.y1, fixture.x2, fixture.y2, fixture.curve_a, n
    )
    algebra_passed = q0 == (fixture.x1, fixture.y1, 0) and q1 == (
        fixture.x3,
        fixture.y3,
        0,
    )
    result = {
        "n": n,
        "inversion_mode": inversion_mode,
        "algebra_passed": algebra_passed,
        "gate_stream_passed": None,
        "dirty_ancilla_count": None,
        "width": None,
        "raw_toffoli": None,
        "raw_cnot": None,
        "build_seconds": 0.0,
        "simulate_seconds": 0.0,
    }
    if algebra_only:
        return result

    start = time.perf_counter()
    circuit = build_test_circuit(n, inversion_mode)
    result["build_seconds"] = time.perf_counter() - start
    result["width"] = circuit["width"]
    raw_toffoli, raw_cnot = circuit["gm"].get_stats()
    result["raw_toffoli"] = raw_toffoli
    result["raw_cnot"] = raw_cnot

    state = [0] * circuit["width"]
    state[circuit["q"]] = 0b10
    x_l = poly_to_l(fixture.x1, n)
    y_l = poly_to_l(fixture.y1, n)
    _set_packed_register(state, circuit["x"], (x_l, x_l))
    _set_packed_register(state, circuit["y"], (y_l, y_l))
    start = time.perf_counter()
    final = simulate_packed(circuit["gm"], state, 2)
    result["simulate_seconds"] = time.perf_counter() - start

    outputs = []
    for case in range(2):
        x_out = _get_packed_register(final, circuit["x"], case)
        y_out = _get_packed_register(final, circuit["y"], case)
        outputs.append((l_to_poly(x_out, n), l_to_poly(y_out, n)))
    expected = ((fixture.x1, fixture.y1), (fixture.x3, fixture.y3))
    clean_wires = (
        circuit["lambda"]
        + circuit["ancilla"]
        + circuit["side_ancilla"]
        + circuit["toffoli_qubits"]
        + circuit["sqr1"]
        + circuit["sqr2"]
        + circuit["sqr3"]
        + circuit["sqr4"]
    )
    dirty = [wire for wire in clean_wires if final[wire] != 0]
    result["outputs"] = tuple(outputs)
    result["expected"] = expected
    result["dirty_ancilla_count"] = len(dirty)
    result["first_dirty_ancillas"] = dirty[:16]
    result["gate_stream_passed"] = (
        tuple(outputs) == expected and not dirty and final[circuit["q"]] == 0b10
    )
    circuit["gm"].clear()
    del circuit
    gc.collect()
    return result


def _parse_args():
    parser = argparse.ArgumentParser(description="Verify one Algorithm 1 in-place point addition.")
    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        choices=SUPPORTED_SIZES,
        default=[163],
        help="build n=163 by default; select other parameter sets explicitly",
    )
    parser.add_argument("--algebra-only", action="store_true")
    parser.add_argument(
        "--inversion",
        choices=INVERSION_MODES + ("both",),
        default=BINARY_ECC,
        help="select low width, optimal Toffoli depth, or verify both in sequence",
    )
    parser.add_argument("--json", type=Path)
    return parser.parse_args()


def main():
    args = _parse_args()
    results = []
    modes = INVERSION_MODES if args.inversion == "both" else (args.inversion,)
    for n in args.sizes:
        for inversion_mode in modes:
            row = verify_one(n, inversion_mode, args.algebra_only)
            results.append(row)
            gate_status = (
                "SKIP"
                if row["gate_stream_passed"] is None
                else "PASS"
                if row["gate_stream_passed"]
                else "FAIL"
            )
            print(
                f"n={n}, inversion={inversion_mode}: "
                f"algebra={'PASS' if row['algebra_passed'] else 'FAIL'}, "
                f"gate={gate_status}, width={row['width']}, "
                f"dirty={row['dirty_ancilla_count']}, "
                f"build={row['build_seconds']:.3f}s, "
                f"simulate={row['simulate_seconds']:.3f}s"
            )
    if args.json:
        args.json.write_text(
            json.dumps(results, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
    return 0 if all(
        row["algebra_passed"] and row["gate_stream_passed"] is not False
        for row in results
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
