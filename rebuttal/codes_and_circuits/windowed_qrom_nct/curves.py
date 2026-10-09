"""Frozen binary-curve instances used by the windowed resource study.

Coordinates are in the polynomial basis defined by ``modulus_low``.  The full
field modulus is ``x**n + modulus_low``.  The constants below are the explicit
parameters printed by OpenSSL for the standardized NIST/SECG B-curves.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BinaryCurve:
    n: int
    openssl_name: str
    modulus_low: int
    a: int
    b: int
    gx: int
    gy: int
    order: int
    cofactor: int = 2

    @property
    def coordinate_bytes(self) -> int:
        return (self.n + 7) // 8

    @property
    def modulus(self) -> int:
        return (1 << self.n) | self.modulus_low


CURVES: dict[int, BinaryCurve] = {
    163: BinaryCurve(
        n=163,
        openssl_name="sect163r2",
        modulus_low=0xC9,
        a=1,
        b=int("020a601907b8c953ca1481eb10512f78744a3205fd", 16),
        gx=int("03f0eba16286a2d57ea0991168d4994637e8343e36", 16),
        gy=int("00d51fbc6c71a0094fa2cdd545b11c5c0c797324f1", 16),
        order=int("040000000000000000000292fe77e70c12a4234c33", 16),
    ),
    233: BinaryCurve(
        n=233,
        openssl_name="sect233r1",
        modulus_low=(1 << 74) | 1,
        a=1,
        b=int("66647ede6c332c7f8c0923bb58213b333b20e9ce4281fe115f7d8f90ad", 16),
        gx=int("00fac9dfcbac8313bb2139f1bb755fef65bc391f8b36f8f8eb7371fd558b", 16),
        gy=int("01006a08a41903350678e58528bebf8a0beff867a7ca36716f7e01f81052", 16),
        order=int("01000000000000000000000000000013e974e72f8a6922031d2603cfe0d7", 16),
    ),
    283: BinaryCurve(
        n=283,
        openssl_name="sect283r1",
        modulus_low=0x10A1,
        a=1,
        b=int("027b680ac8b8596da5a4af8a19a0303fca97fd7645309fa2a581485af6263e313b79a2f5", 16),
        gx=int("05f939258db7dd90e1934f8c70b0dfec2eed25b8557eac9c80e2e198f8cdbecd86b12053", 16),
        gy=int("03676854fe24141cb98fe6d4b20d02b4516ff702350eddb0826779c813f0df45be8112f4", 16),
        order=int("03ffffffffffffffffffffffffffffffffffef90399660fc938a90165b042a7cefadb307", 16),
    ),
    571: BinaryCurve(
        n=571,
        openssl_name="sect571r1",
        modulus_low=0x425,
        a=1,
        b=int(
            "02f40e7e2221f295de297117b7f3d62f5c6a97ffcb8ceff1cd6ba8ce4a9a18ad"
            "84ffabbd8efa59332be7ad6756a66e294afd185a78ff12aa520e4de739baca0c7f"
            "feff7f2955727a",
            16,
        ),
        gx=int(
            "0303001d34b856296c16c0d40d3cd7750a93d1d2955fa80aa5f40fc8db7b2abd"
            "bde53950f4c0d293cdd711a35b67fb1499ae60038614f1394abfa3b4c850d927e"
            "1e7769c8eec2d19",
            16,
        ),
        gy=int(
            "037bf27342da639b6dccfffeb73d69d78c6c27a6009cbbca1980f8533921e8a68"
            "4423e43bab08a576291af8f461bb2a8b3531d2f0485c19b16e2f1516e23dd3c1"
            "a4827af1b8ac15b",
            16,
        ),
        order=int(
            "03ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
            "ffe661ce18ff55987308059b186823851ec7dd9ca1161de93d5174d66e8382e9bb"
            "2fe84e47",
            16,
        ),
    ),
}


def gf_reduce(value: int, curve: BinaryCurve) -> int:
    """Reduce a polynomial integer modulo the curve field polynomial."""
    while value.bit_length() > curve.n:
        shift = value.bit_length() - curve.n - 1
        value ^= curve.modulus << shift
    return value


def gf_multiply(left: int, right: int, curve: BinaryCurve) -> int:
    """Polynomial-basis multiplication, used only by verification tests."""
    product = 0
    while right:
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
    return gf_reduce(product, curve)


def gf_square(value: int, curve: BinaryCurve) -> int:
    expanded = 0
    bit = 0
    while value:
        if value & 1:
            expanded |= 1 << (2 * bit)
        value >>= 1
        bit += 1
    return gf_reduce(expanded, curve)


def gf_inverse(value: int, curve: BinaryCurve) -> int:
    if value == 0:
        raise ZeroDivisionError("zero has no field inverse")
    # Fermat: value^(2^n-2).  Repeated square-and-multiply keeps this compact
    # and is used only for classical correctness fixtures, not resource costs.
    result = 1
    base = value
    exponent = (1 << curve.n) - 2
    while exponent:
        if exponent & 1:
            result = gf_multiply(result, base, curve)
        exponent >>= 1
        if exponent:
            base = gf_square(base, curve)
    return result


def affine_add(
    left: tuple[int, int], right: tuple[int, int], curve: BinaryCurve
) -> tuple[int, int]:
    x1, y1 = left
    x2, y2 = right
    denominator = x1 ^ x2
    if denominator == 0:
        raise ValueError("ordinary affine-add branch requires x1 != x2")
    slope = gf_multiply(y1 ^ y2, gf_inverse(denominator, curve), curve)
    x3 = gf_square(slope, curve) ^ slope ^ x1 ^ x2 ^ curve.a
    y3 = gf_multiply(slope, x1 ^ x3, curve) ^ x3 ^ y1
    return x3, y3


def on_curve(x: int, y: int, curve: BinaryCurve) -> bool:
    """Return whether ``(x,y)`` satisfies y^2+xy=x^3+a*x^2+b."""
    lhs = gf_square(y, curve) ^ gf_multiply(x, y, curve)
    x2 = gf_square(x, curve)
    rhs = gf_multiply(x2, x, curve) ^ gf_multiply(curve.a, x2, curve) ^ curve.b
    return lhs == rhs
