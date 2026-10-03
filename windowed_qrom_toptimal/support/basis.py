"""Exact polynomial-basis to AC L-basis conversion."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from .paths import AC_DATA_ROOT


DEFAULT_DATA_ROOT = AC_DATA_ROOT


@lru_cache(maxsize=None)
def _matrix_rows(n: int, direction: str, data_root: str) -> tuple[int, ...]:
    if direction not in {"X2L", "L2X"}:
        raise ValueError("direction must be X2L or L2X")
    path = Path(data_root) / f"quantum_{n}" / f"{n}_{direction}.txt"
    rows: list[int] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        bits = raw.strip()
        if not bits:
            continue
        if len(bits) != n or set(bits) - {"0", "1"}:
            raise ValueError(f"invalid basis matrix row in {path}")
        # Matrix columns correspond to little-endian field-coordinate bits.
        mask = 0
        for column, value in enumerate(bits):
            if value == "1":
                mask |= 1 << column
        rows.append(mask)
    if len(rows) != n:
        raise ValueError(f"{path} has {len(rows)} rows; expected {n}")
    return tuple(rows)


@lru_cache(maxsize=None)
def _byte_tables(n: int, direction: str, data_root: str) -> tuple[tuple[int, ...], ...]:
    """Return 8-bit lookup tables for the selected GF(2) linear map."""
    rows = _matrix_rows(n, direction, data_root)
    columns: list[int] = []
    for input_bit in range(n):
        column = 0
        mask = 1 << input_bit
        for output_bit, row in enumerate(rows):
            if row & mask:
                column |= 1 << output_bit
        columns.append(column)

    result: list[tuple[int, ...]] = []
    for byte_offset in range((n + 7) // 8):
        block = columns[8 * byte_offset : 8 * byte_offset + 8]
        values = [0] * 256
        for byte in range(1, 256):
            lowest = byte & -byte
            bit = lowest.bit_length() - 1
            values[byte] = values[byte ^ lowest] ^ (block[bit] if bit < len(block) else 0)
        result.append(tuple(values))
    return tuple(result)


def transform(value: int, n: int, direction: str, data_root: Path = DEFAULT_DATA_ROOT) -> int:
    tables = _byte_tables(n, direction, str(data_root.resolve()))
    result = 0
    for table in tables:
        result ^= table[value & 0xFF]
        value >>= 8
    return result


def poly_to_l(value: int, n: int, data_root: Path = DEFAULT_DATA_ROOT) -> int:
    return transform(value, n, "X2L", data_root)


def l_to_poly(value: int, n: int, data_root: Path = DEFAULT_DATA_ROOT) -> int:
    return transform(value, n, "L2X", data_root)
