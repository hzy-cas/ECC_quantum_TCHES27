"""Basic gates and low-depth Karatsuba arithmetic.

"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from .config import FieldConfig, get_config
from .qasm import FlatGateManager


gm = FlatGateManager()
config: FieldConfig = get_config(163)
GF_N = config.n
KARATSUBA_ANC_BASE = config.karatsuba_ancillas
TOFFOLI_BASE = config.multiplication_targets
TOFFOLI_OFFSET = config.inversion_multiplications

Toffoli_qubits: List[int] = []
cnt = 0

_square_columns: Dict[int, Tuple[int, ...]] = {}
_square_byte_tables: List[Tuple[int, ...]] = []


def configure_backend(target_n: int = 163) -> None:
    """Select one field and reset all gate-emission state."""

    global config, GF_N, KARATSUBA_ANC_BASE, TOFFOLI_BASE, TOFFOLI_OFFSET
    global Toffoli_qubits, cnt
    global _square_columns, _square_byte_tables

    config = get_config(target_n)
    GF_N = config.n
    KARATSUBA_ANC_BASE = config.karatsuba_ancillas
    TOFFOLI_BASE = config.multiplication_targets
    TOFFOLI_OFFSET = config.inversion_multiplications
    Toffoli_qubits = []
    cnt = 0
    _square_columns = {}
    _square_byte_tables = []
    gm.clear()


def X(target: int) -> List[int]:
    gm.add_X(target)
    return [target]


def CNOT(control: int, target: int) -> List[int]:
    gm.add_CNOT(control, target)
    return [control, target]


def Toffoli_gate(control1: int, control2: int, target: int) -> List[int]:
    gm.add_Toffoli(control1, control2, target)
    return [control1, control2, target]


def CNOT_n(source: Sequence[int], target: Sequence[int]) -> None:
    if len(source) != len(target):
        raise ValueError("CNOT_n requires equal-sized registers")
    for control, output in zip(source, target):
        CNOT(control, output)


def CONST_ADD_n(target: Sequence[int], value: int) -> None:
    for bit, wire in enumerate(target):
        if (value >> bit) & 1:
            X(wire)


def CSWAP(control: int, left: int, right: int) -> None:
    CNOT(right, left)
    Toffoli_gate(control, left, right)
    CNOT(right, left)


def copy_parallel(
    source: int, clean_ancillas: Sequence[int], copies: int
) -> List[int]:
    if copies < 0:
        raise ValueError("copy count must be nonnegative")
    if copies == 0:
        return []
    if len(clean_ancillas) + 1 < copies:
        raise ValueError("copy_parallel has insufficient clean ancillas")

    result = [source] + list(clean_ancillas[: copies - 1])
    largest_power = 1
    while largest_power * 2 <= copies:
        largest_power *= 2
    layer_size = 1
    while layer_size < largest_power:
        for index in range(layer_size):
            CNOT(result[index], result[layer_size + index])
        layer_size *= 2
    for index in range(copies - largest_power):
        CNOT(result[index], result[largest_power + index])
    return result


def _build_square_tables() -> None:
    """Generate Frobenius matrices with Python big-integer arithmetic."""

    global _square_columns, _square_byte_tables
    if _square_columns:
        return

    n = GF_N
    mask = (1 << n) - 1
    tail = config.polynomial_tail

    # remainder[d] = x^d mod f(x), for every degree needed by squaring.
    remainders = [0] * (2 * n - 1)
    value = 1
    for degree in range(2 * n - 1):
        remainders[degree] = value
        carry = (value >> (n - 1)) & 1
        value = (value << 1) & mask
        if carry:
            value ^= tail

    byte_count = (n + 7) // 8
    byte_tables: List[Tuple[int, ...]] = []
    for byte_index in range(byte_count):
        basis = []
        for offset in range(8):
            input_bit = 8 * byte_index + offset
            basis.append(remainders[2 * input_bit] if input_bit < n else 0)
        table = [0] * 256
        for byte_value in range(1, 256):
            least_bit = byte_value & -byte_value
            bit_index = least_bit.bit_length() - 1
            table[byte_value] = table[byte_value ^ least_bit] ^ basis[bit_index]
        byte_tables.append(tuple(table))
    _square_byte_tables = byte_tables

    def square_integer(element: int) -> int:
        result = 0
        value_left = element
        for table in byte_tables:
            result ^= table[value_left & 0xFF]
            value_left >>= 8
        return result

    requested = set(config.square_powers)
    maximum_power = max(requested)
    columns = [1 << column for column in range(n)]
    for power in range(1, maximum_power + 1):
        columns = [square_integer(column) for column in columns]
        if power in requested:
            _square_columns[power] = tuple(columns)


def square_xor(
    source: Sequence[int], target: Sequence[int], power: int
) -> None:
    if len(source) < GF_N or len(target) < GF_N:
        raise ValueError("square registers are shorter than GF_N")
    _build_square_tables()
    try:
        columns = _square_columns[power]
    except KeyError as error:
        raise ValueError(f"uncached square power {power}") from error

    n = GF_N
    # Preserve Binary_ECC's original nested-loop gate order.
    for diagonal in range(n):
        for output_bit in range(n):
            input_bit = (diagonal + output_bit) % n
            if (columns[input_bit] >> output_bit) & 1:
                CNOT(source[input_bit], target[output_bit])


def square_plus_input_xor(
    source: Sequence[int], target: Sequence[int]
) -> None:
    if len(source) < GF_N or len(target) < GF_N:
        raise ValueError("square-plus-input registers are too short")
    _build_square_tables()
    columns = _square_columns[1]
    n = GF_N
    for diagonal in range(n):
        for output_bit in range(n):
            input_bit = (diagonal + output_bit) % n
            enabled = ((columns[input_bit] >> output_bit) & 1) ^ (
                input_bit == output_bit
            )
            if enabled:
                CNOT(source[input_bit], target[output_bit])


def Square(
    combined: Sequence[int], n: int, power: int
) -> Tuple[List[int], List[int]]:
    if n != GF_N or len(combined) < 2 * n:
        raise ValueError("Square expects two GF_N-sized registers")
    source = list(combined[:n])
    target = list(combined[-n:])
    square_xor(source, target, power)
    return source, target


def Squaring_plus_input(
    combined: Sequence[int], n: int
) -> Tuple[List[int], List[int]]:
    if n != GF_N or len(combined) < 2 * n:
        raise ValueError("Squaring_plus_input expects two n-bit registers")
    source = list(combined[:n])
    target = list(combined[-n:])
    square_plus_input_xor(source, target)
    return source, target


def Modular_small(
    source: Sequence[int], result: Sequence[int], size: int
) -> None:
    for tap in config.reduction_taps:
        for index in range(size):
            CNOT(source[index], result[tap + index])


def Reduction(product: Sequence[int]) -> List[int]:
    n = GF_N
    if len(product) != 2 * n - 1:
        raise ValueError("Reduction expects a 2n-1 bit product")
    for tap in config.reduction_taps:
        for index in range(n - 1):
            if index + tap < n:
                CNOT(product[n + index], product[index + tap])
    for tap in config.reduction_taps:
        if tap <= 1:
            continue
        size = tap - 1
        start = 2 * n - tap
        Modular_small(product[start : start + size], product, size)
    return list(product[:n])


def combine(
    low: Sequence[int], high: Sequence[int], cross: Sequence[int], n: int
) -> List[int]:
    if n % 2:
        for index in range(n):
            CNOT(low[index], cross[index])
        for index in range(n - 2):
            CNOT(high[index], cross[index])
        for index in range(n // 2):
            CNOT(low[n // 2 + 1 + index], cross[index])
        for index in range(n // 2):
            CNOT(high[index], cross[n // 2 + 1 + index])
        high_tail = (2 * n - 1) - (n // 2 + 1) - n
        return (
            list(low[: n // 2 + 1])
            + list(cross[:n])
            + list(high[n // 2 : n // 2 + high_tail])
        )

    half = n // 2
    for index in range(n - 1):
        CNOT(low[index], cross[index])
        CNOT(high[index], cross[index])
    for index in range(half - 1):
        CNOT(low[half + index], cross[index])
        CNOT(high[index], cross[half + index])
    return (
        list(low[:half])
        + list(cross[: n - 1])
        + list(high[half - 1 : 2 * half - 1])
    )


def recursive_karatsuba(
    left: Sequence[int],
    right: Sequence[int],
    n: int,
    count: int,
    clean_ancillas: Sequence[int],
):
    global cnt

    if len(left) != n or len(right) != n:
        raise ValueError("recursive_karatsuba operand size mismatch")
    if n == 1:
        if cnt >= len(Toffoli_qubits):
            raise IndexError("Karatsuba multiplication target exhausted")
        target = Toffoli_qubits[cnt]
        cnt += 1
        Toffoli_gate(left[0], right[0], target)
        return [target], count, list(clean_ancillas)

    low_size = (n + 1) // 2
    if count + 2 * low_size > len(clean_ancillas):
        raise IndexError("Karatsuba preparation ancilla exhausted")
    left_sum = list(clean_ancillas[count : count + low_size])
    count += low_size
    right_sum = list(clean_ancillas[count : count + low_size])
    count += low_size

    preparation_start = gm.current_pointer()
    for index in range(low_size):
        CNOT(left[index], left_sum[index])
    for index in range(n // 2):
        CNOT(left[low_size + index], left_sum[index])
    for index in range(low_size):
        CNOT(right[index], right_sum[index])
    for index in range(n // 2):
        CNOT(right[low_size + index], right_sum[index])
    preparation_end = gm.current_pointer()

    if low_size == 1:
        if cnt + 3 > len(Toffoli_qubits):
            raise IndexError("Karatsuba base targets exhausted")
        result = list(Toffoli_qubits[cnt : cnt + 3])
        cnt += 3
        # This is the paper ordering already used in
        # recursive_karatsuba/inv/init/basic_gates_*.py: all three nonlinear
        # product terms are emitted before the two target-XOR CNOTs.  Their
        # controls and targets are disjoint, so the three Toffolis form one
        # Toffoli layer.
        Toffoli_gate(left[0], right[0], result[0])
        Toffoli_gate(left[1], right[1], result[2])
        Toffoli_gate(left_sum[0], right_sum[0], result[1])
        CNOT(result[0], result[1])
        CNOT(result[2], result[1])
        gm.replay_reverse(preparation_start, preparation_end)
        return result, count, list(clean_ancillas)

    low, count, _ = recursive_karatsuba(
        left[:low_size], right[:low_size], low_size, count, clean_ancillas
    )
    high, count, _ = recursive_karatsuba(
        left[low_size:], right[low_size:], n // 2, count, clean_ancillas
    )
    cross, count, _ = recursive_karatsuba(
        left_sum, right_sum, low_size, count, clean_ancillas
    )
    gm.replay_reverse(preparation_start, preparation_end)
    return combine(low, high, cross, n), count, list(clean_ancillas)


def multiply(
    left: Sequence[int], right: Sequence[int], clean_ancillas: Sequence[int]
) -> List[int]:
    product, _, _ = recursive_karatsuba(
        list(left[:GF_N]), list(right[:GF_N]), GF_N, 0, clean_ancillas
    )
    return Reduction(product)
