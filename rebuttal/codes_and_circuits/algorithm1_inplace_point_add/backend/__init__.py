"""Local ProjectQ-free backend for complete gate streams."""

from .qasm import FlatGateManager, MatrixCache

__all__ = ["FlatGateManager", "MatrixCache"]
