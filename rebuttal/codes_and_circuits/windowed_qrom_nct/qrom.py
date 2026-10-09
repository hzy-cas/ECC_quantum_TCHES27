"""Gate-exact low-width unary-iteration QROM.

The coherent implementation follows the segment-tree traversal used by
Qualtran's unary iteration.  A negative-control AND is decomposed as
X--Toffoli--X.  Its final positive-control inverse is a genuine Toffoli in the
coherent model; it is never replaced by a forward call on dummy registers.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol


class GateSink(Protocol):
    def x(self, target: int) -> None: ...
    def cnot(self, control: int, target: int) -> None: ...
    def toffoli(self, control1: int, control2: int, target: int) -> None: ...


@dataclass(frozen=True)
class QromRegisters:
    # Address is little-endian: address[bit] stores the coefficient of 2**bit.
    address: tuple[int, ...]
    selector_ancillas: tuple[int, ...]
    bus: tuple[int, ...]
    data_fanout: tuple[int, ...] = ()

    def validate(self, entries: int) -> None:
        if entries <= 1 or entries & (entries - 1):
            raise ValueError("unary QROM currently requires a power-of-two table with at least 2 entries")
        bits = entries.bit_length() - 1
        if len(self.address) != bits:
            raise ValueError(f"table has {entries} entries but address has {len(self.address)} bits")
        if len(self.selector_ancillas) != max(0, bits - 1):
            raise ValueError("unary iteration needs address_bits-1 clean selector ancillas")
        wires = self.address + self.selector_ancillas + self.bus + self.data_fanout
        if len(set(wires)) != len(wires):
            raise ValueError("QROM registers overlap")


def _negative_and(sink: GateSink, control: int, selection: int, target: int) -> None:
    sink.x(selection)
    sink.toffoli(control, selection, target)
    sink.x(selection)


def _copy_parallel(
    sink: GateSink, control: int, ancillas: tuple[int, ...], count: int
) -> tuple[list[int], list[tuple[int, int]]]:
    controls = [control]
    gates: list[tuple[int, int]] = []
    ancilla_index = 0
    while len(controls) < count:
        current = tuple(controls)
        for source in current:
            if len(controls) >= count:
                break
            target = ancillas[ancilla_index]
            ancilla_index += 1
            sink.cnot(source, target)
            gates.append((source, target))
            controls.append(target)
    return controls, gates


def _emit_word(
    sink: GateSink,
    control: int,
    bus: tuple[int, ...],
    fanout: tuple[int, ...],
    word: int,
    inverse: bool,
) -> None:
    indices = [bit for bit in range(len(bus)) if (word >> bit) & 1]
    if not indices:
        return
    if not fanout:
        if inverse:
            indices.reverse()
        for bit in indices:
            sink.cnot(control, bus[bit])
        return
    if len(fanout) + 1 < len(bus):
        raise ValueError("fanout register must have at least bus_bits-1 wires")
    controls, copy_gates = _copy_parallel(sink, control, fanout, len(indices))
    order = range(len(indices) - 1, -1, -1) if inverse else range(len(indices))
    for i in order:
        sink.cnot(controls[i], bus[indices[i]])
    for source, target in reversed(copy_gates):
        sink.cnot(source, target)


def emit_unary_qrom(
    sink: GateSink,
    registers: QromRegisters,
    words: Sequence[int],
    *,
    inverse: bool = False,
) -> None:
    """Emit ``bus ^= words[address]`` or its literal reversed gate stream.

    ``inverse=True`` does not call the forward oracle on substitute registers:
    every gate is emitted in the exact reverse order on the same address,
    selector, fanout, and bus wires.
    """
    registers.validate(len(words))
    if any(word < 0 or word.bit_length() > len(registers.bus) for word in words):
        raise ValueError("table word does not fit QROM bus")

    # Qualtran orders the selection sequence from the most-significant tree
    # level down.  Externally we expose the conventional little-endian window.
    selection = tuple(reversed(registers.address))
    ancillas = registers.selector_ancillas
    entries = len(words)

    def subtree(
        control: int,
        level: int,
        left: int,
        right: int,
        reverse: bool,
    ) -> None:
        if right - left == 1:
            _emit_word(
                sink,
                control,
                registers.bus,
                registers.data_fanout,
                words[left],
                reverse,
            )
            return
        middle = (left + right) // 2
        select = selection[level + 1]
        ancilla = ancillas[level]
        if not reverse:
            _negative_and(sink, control, select, ancilla)
            subtree(ancilla, level + 1, left, middle, False)
            sink.cnot(control, ancilla)
            subtree(ancilla, level + 1, middle, right, False)
            sink.toffoli(control, select, ancilla)
        else:
            # Exact inverse of the five forward blocks above.
            sink.toffoli(control, select, ancilla)
            subtree(ancilla, level + 1, middle, right, True)
            sink.cnot(control, ancilla)
            subtree(ancilla, level + 1, left, middle, True)
            _negative_and(sink, control, select, ancilla)

    root = selection[0]
    middle = entries // 2
    if not inverse:
        sink.x(root)
        subtree(root, 0, 0, middle, False)
        sink.x(root)
        subtree(root, 0, middle, entries, False)
    else:
        subtree(root, 0, middle, entries, True)
        sink.x(root)
        subtree(root, 0, 0, middle, True)
        sink.x(root)


def pack_point_word(x: int, y: int, n: int) -> int:
    """Pack two n-bit little-endian field registers into one QROM word."""
    if min(x, y) < 0 or x.bit_length() > n or y.bit_length() > n:
        raise ValueError("coordinate outside field-register width")
    return x | (y << n)
