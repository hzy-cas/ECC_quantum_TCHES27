"""Static parameters and repository paths; emits no gates."""

import os
from pathlib import Path


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent


def _discover_repository_root() -> Path:
    override = os.environ.get("ECC_SHOR_REPOSITORY_ROOT")
    if override:
        repository_root = Path(override).expanduser().resolve()
        if (repository_root / "mul_line" / "C++" / "data").is_dir():
            return repository_root
        raise RuntimeError(
            "ECC_SHOR_REPOSITORY_ROOT does not contain mul_line/C++/data"
        )
    for candidate in PACKAGE_ROOT.parents:
        if (candidate / "mul_line" / "C++" / "data").is_dir():
            return candidate
    raise RuntimeError(
        "could not locate the repository root; set ECC_SHOR_REPOSITORY_ROOT"
    )


REPOSITORY_ROOT = _discover_repository_root()
DATA_ROOT = REPOSITORY_ROOT / "mul_line" / "C++" / "data"

CONFIGS = {
    163: {
        "block_size": 743,
        "inversion_multiplications": 9,
        "multiplication_targets": 906,
        "mul_dims": (325, 906),
    },
    233: {
        "block_size": 1108,
        "inversion_multiplications": 10,
        "multiplication_targets": 1341,
        "mul_dims": (465, 1341),
    },
    283: {
        "block_size": 1385,
        "inversion_multiplications": 11,
        "multiplication_targets": 1668,
        "mul_dims": (565, 1668),
    },
    571: {
        "block_size": 2998,
        "inversion_multiplications": 13,
        "multiplication_targets": 3569,
        "mul_dims": (1141, 3569),
    },
}

SUPPORTED_SIZES = tuple(CONFIGS)

BINARY_ECC = "binary_ecc"
OPTIMAL_DEPTH = "optimal_depth"
INVERSION_MODES = (BINARY_ECC, OPTIMAL_DEPTH)
