"""Repository path discovery for the portable submission package."""

from __future__ import annotations

import os
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent


def _discover_repository_root() -> Path:
    override = os.environ.get("WINDOWED_QROM_REPO_ROOT")
    if override:
        return Path(override).expanduser().resolve()
    for candidate in PACKAGE_ROOT.parents:
        if (candidate / "mul_line" / "C++" / "data").is_dir():
            return candidate
    raise RuntimeError(
        "could not locate the repository root; set WINDOWED_QROM_REPO_ROOT"
    )


REPOSITORY_ROOT = _discover_repository_root()
AC_ROOT = PACKAGE_ROOT / "ac_backend"
AC_RESULTS = AC_ROOT / "results"
AC_EXECUTABLE = AC_ROOT / "build" / "ac_balanced_estimator"
AC_DATA_ROOT = REPOSITORY_ROOT / "mul_line" / "C++" / "data"
