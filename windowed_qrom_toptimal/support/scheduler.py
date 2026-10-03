"""Per-line ASAP scheduler matching the full/Toffoli depth rules of Table 6."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GateStats:
    width: int
    x: int
    cnot: int
    toffoli: int
    full_depth: int
    toffoli_depth: int


class PerLineAsap:
    """Incremental scheduler for coherent X, CNOT, and Toffoli gates.

    The NCT rules are exactly those in the repository's ``GateManager``:
    every gate increments full depth on its operand lines; Clifford gates
    propagate (but do not increment) Toffoli depth; a Toffoli increments it.
    """

    def __init__(self, width: int):
        if width < 0:
            raise ValueError("width must be nonnegative")
        self.width = width
        self.full = [0] * width
        self.toffoli_layer = [0] * width
        self.x_count = 0
        self.cnot_count = 0
        self.toffoli_count = 0

    def _validate(self, *wires: int) -> None:
        if len(set(wires)) != len(wires):
            raise ValueError(f"gate has repeated wire: {wires}")
        if any(wire < 0 or wire >= self.width for wire in wires):
            raise IndexError(f"wire outside width {self.width}: {wires}")

    def x(self, target: int) -> None:
        self._validate(target)
        self.x_count += 1
        self.full[target] += 1

    def cnot(self, control: int, target: int) -> None:
        self._validate(control, target)
        self.cnot_count += 1
        depth = max(self.full[control], self.full[target]) + 1
        self.full[control] = self.full[target] = depth
        layer = max(self.toffoli_layer[control], self.toffoli_layer[target])
        self.toffoli_layer[control] = self.toffoli_layer[target] = layer

    def toffoli(self, control1: int, control2: int, target: int) -> None:
        self._validate(control1, control2, target)
        self.toffoli_count += 1
        depth = max(self.full[control1], self.full[control2], self.full[target]) + 1
        self.full[control1] = self.full[control2] = self.full[target] = depth
        layer = (
            max(
                self.toffoli_layer[control1],
                self.toffoli_layer[control2],
                self.toffoli_layer[target],
            )
            + 1
        )
        self.toffoli_layer[control1] = self.toffoli_layer[control2] = self.toffoli_layer[target] = layer

    def stats(self) -> GateStats:
        return GateStats(
            width=self.width,
            x=self.x_count,
            cnot=self.cnot_count,
            toffoli=self.toffoli_count,
            full_depth=max(self.full, default=0),
            toffoli_depth=max(self.toffoli_layer, default=0),
        )


class ClassicalGateSimulator:
    """Computational-basis simulator used for small reversible QROM tests."""

    def __init__(self, width: int, state: int = 0):
        self.width = width
        self.state = state

    def bit(self, wire: int) -> int:
        return (self.state >> wire) & 1

    def x(self, target: int) -> None:
        self.state ^= 1 << target

    def cnot(self, control: int, target: int) -> None:
        if self.bit(control):
            self.state ^= 1 << target

    def toffoli(self, control1: int, control2: int, target: int) -> None:
        if self.bit(control1) and self.bit(control2):
            self.state ^= 1 << target
