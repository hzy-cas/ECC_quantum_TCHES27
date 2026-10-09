"""Itoh--Tsujii inversion chains used by the Binary_ECC point additions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

from . import basic_gates as gates


@dataclass(frozen=True)
class InversionResult:
    inverse: List[int]
    last_product: List[int]


def Inverison_Itoh_Tsujii_based(
    a: Sequence[int],
    sqr1: Sequence[int],
    sqr2: Sequence[int],
    ancilla: Sequence[int],
    *,
    clear_final_square: bool,
) -> InversionResult:
    """Emit the exact four NIST-field inversion chains.
    """

    n = gates.GF_N
    if not (
        len(a) == len(sqr1) == len(sqr2) == n
        and len(ancilla) >= gates.KARATSUBA_ANC_BASE
    ):
        raise ValueError("invalid inversion workspace")

    start_target = gates.cnt

    def square(source, target, power):
        gates.square_xor(source, target, power)
        return list(target)

    def multiply(left, right):
        return gates.multiply(left, right, ancilla)

    s0 = list(sqr1)
    s1 = list(sqr2)

    a2 = square(a, s0, 1)
    a_2_0 = multiply(a, a2)
    a2_2 = square(a_2_0, s1, 2)
    square(a, s0, 1)
    b = multiply(a2_2, a_2_0)
    square(a_2_0, s1, 2)
    b2_4 = square(b, s0, 4)
    c = multiply(b, b2_4)
    square(b, s0, 4)
    c2_8 = square(c, s1, 8)
    d = multiply(c, c2_8)
    square(c, s1, 8)
    d2_16 = square(d, s0, 16)

    if n == 163:
        e = multiply(d, d2_16)
        square(d, s0, 16)
        e2_32 = square(e, s1, 32)
        f = multiply(e, e2_32)
        square(e, s1, 32)
        e2_2 = square(e, s1, 2)
        f2_64 = square(f, s0, 64)
        g = multiply(f, f2_64)
        square(f, s0, 64)
        g2_34 = square(g, s0, 34)
        mul1 = multiply(a_2_0, e2_2)
        square(e, s1, 2)
        last_product = multiply(mul1, g2_34)
        if clear_final_square:
            square(g, s0, 34)
        result = square(last_product, s1, 1)

    elif n == 233:
        e = multiply(d, d2_16)
        square(d, s0, 16)
        e2_32 = square(e, s1, 32)
        f = multiply(e, e2_32)
        square(e, s1, 32)
        f2_64 = square(f, s0, 64)
        e2_8 = square(e, s1, 8)
        g = multiply(f, f2_64)
        square(f, s0, 64)
        f2_40 = square(f, s0, 40)
        mul1 = multiply(c, e2_8)
        square(e, s1, 8)
        g2_104 = square(g, s1, 104)
        mul2 = multiply(mul1, f2_40)
        square(f, s0, 40)
        last_product = multiply(mul2, g2_104)
        if clear_final_square:
            square(g, s1, 104)
        result = square(last_product, s0, 1)

    elif n == 283:
        c2_2 = square(c, s1, 2)
        e = multiply(d, d2_16)
        square(d, s0, 16)
        e2_32 = square(e, s0, 32)
        mul1 = multiply(a_2_0, c2_2)
        square(c, s1, 2)
        d2_10 = square(d, s1, 10)
        f = multiply(e, e2_32)
        square(e, s0, 32)
        f2_64 = square(f, s0, 64)
        mul2 = multiply(mul1, d2_10)
        square(d, s1, 10)
        g = multiply(f, f2_64)
        square(f, s0, 64)
        g2_128 = square(g, s1, 128)
        h = multiply(g, g2_128)
        square(g, s1, 128)
        h2_26 = square(h, s0, 26)
        last_product = multiply(mul2, h2_26)
        if clear_final_square:
            square(h, s0, 26)
        result = square(last_product, s1, 1)

    elif n == 571:
        c2_2 = square(c, s1, 2)
        e = multiply(d, d2_16)
        square(d, s0, 16)
        e2_32 = square(e, s0, 32)
        mul1 = multiply(a_2_0, c2_2)
        square(c, s1, 2)
        d2_10 = square(d, s1, 10)
        f = multiply(e, e2_32)
        square(e, s0, 32)
        f2_64 = square(f, s0, 64)
        mul2 = multiply(mul1, d2_10)
        square(d, s1, 10)
        e2_26 = square(e, s1, 26)
        g = multiply(f, f2_64)
        square(f, s0, 64)
        g2_128 = square(g, s0, 128)
        mul3 = multiply(mul2, e2_26)
        square(e, s1, 26)
        h = multiply(g, g2_128)
        square(g, s0, 128)
        h2_256 = square(h, s1, 256)
        i_value = multiply(h, h2_256)
        square(h, s1, 256)
        i2_58 = square(i_value, s0, 58)
        last_product = multiply(mul3, i2_58)
        if clear_final_square:
            square(i_value, s0, 58)
        result = square(last_product, s1, 1)

    else:  # ``configure_backend`` prevents this branch.
        raise AssertionError(f"unsupported inversion chain n={n}")

    consumed = gates.cnt - start_target
    expected = gates.TOFFOLI_OFFSET * gates.TOFFOLI_BASE
    if consumed != expected:
        raise AssertionError(
            f"inversion consumed {consumed} product targets; expected {expected}"
        )
    return InversionResult(result, last_product)
