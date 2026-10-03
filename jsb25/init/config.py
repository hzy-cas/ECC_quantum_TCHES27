"""GF(2^n) arithmetic parameters used by the local circuit builders."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple


SUPPORTED_FIELDS: Tuple[int, ...] = (163, 233, 283, 571)


@dataclass(frozen=True)
class FieldConfig:
    n: int
    reduction_taps: Tuple[int, ...]
    karatsuba_ancillas: int
    multiplication_targets: int
    inversion_multiplications: int
    square_powers: Tuple[int, ...]
    x2: int
    y2: int

    @property
    def polynomial_tail(self) -> int:
        """Return f(x)-x^n as a bit polynomial."""

        return sum(1 << tap for tap in self.reduction_taps)


FIELD_CONFIGS: Dict[int, FieldConfig] = {
    163: FieldConfig(
        n=163,
        reduction_taps=(0, 3, 6, 7),
        karatsuba_ancillas=8448,
        multiplication_targets=4387,
        inversion_multiplications=9,
        square_powers=(1, 2, 4, 8, 16, 32, 34, 64),
        x2=0xFF00,
        y2=0x00FF,
    ),
    233: FieldConfig(
        n=233,
        reduction_taps=(0, 74),
        karatsuba_ancillas=12180,
        multiplication_targets=6323,
        inversion_multiplications=10,
        square_powers=(1, 2, 4, 8, 16, 32, 40, 64, 104),
        x2=0xFF00,
        y2=0x00FF,
    ),
    283: FieldConfig(
        n=283,
        reduction_taps=(0, 5, 7, 12),
        karatsuba_ancillas=19980,
        multiplication_targets=10273,
        inversion_multiplications=11,
        square_powers=(1, 2, 4, 8, 10, 16, 26, 32, 64, 128),
        x2=0xFF00,
        y2=0x00FF,
    ),
    571: FieldConfig(
        n=571,
        reduction_taps=(0, 2, 5, 10),
        karatsuba_ancillas=61200,
        multiplication_targets=31171,
        inversion_multiplications=13,
        square_powers=(1, 2, 4, 8, 10, 16, 26, 32, 58, 64, 128, 256),
        x2=0xFF00,
        y2=0x00FF,
    ),
}


def get_config(n: int) -> FieldConfig:
    try:
        return FIELD_CONFIGS[n]
    except KeyError as error:
        choices = ", ".join(map(str, SUPPORTED_FIELDS))
        raise ValueError(f"unsupported n={n}; choose one of {choices}") from error
