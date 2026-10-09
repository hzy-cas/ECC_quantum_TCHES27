"""Compatibility imports for the repository's shared NCT backend."""
from pathlib import Path
import sys

for _parent in Path(__file__).resolve().parents:
    if (_parent / "circuit_backend" / "python" / "qasm.py").is_file():
        if str(_parent) not in sys.path:
            sys.path.insert(0, str(_parent))
        break
else:
    raise ImportError("circuit_backend is missing; keep the repository directory structure")

from circuit_backend.python.stats_utils import get_exact_resources_optimized
