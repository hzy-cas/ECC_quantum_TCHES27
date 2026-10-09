"""Compact, non-ProjectQ NCT gate-stream backend.
"""

from __future__ import annotations

import array
from dataclasses import dataclass
from typing import Iterator, MutableSequence, Tuple


@dataclass(frozen=True)
class StreamStats:
    words: int = 0
    gates: int = 0
    x: int = 0
    cnot: int = 0
    toffoli: int = 0
    full_depth: int = 0
    ordered_toffoli_depth: int = 0


class FlatGateManager:
    """Store one complete NCT gate stream in compact integer chunks."""

    CHUNK_LIMIT = 1_000_000
    SHIFT_OP = 60
    SHIFT_ARG = 30
    MASK_ARG = (1 << SHIFT_ARG) - 1

    OP_DATA = 0
    OP_X = 1
    OP_CNOT = 2
    OP_TOFFOLI = 3
    OP_NOP = 15

    def __init__(self) -> None:
        self.chunks = []
        self.current_chunk = array.array("Q")

    def _append(self, value: int) -> None:
        if len(self.current_chunk) >= self.CHUNK_LIMIT:
            self.chunks.append(self.current_chunk)
            self.current_chunk = array.array("Q")
        self.current_chunk.append(value)

    def add_X(self, target: int) -> None:
        self._append((self.OP_X << self.SHIFT_OP) | target)

    def add_CNOT(self, control: int, target: int) -> None:
        self._append(
            (self.OP_CNOT << self.SHIFT_OP)
            | (control << self.SHIFT_ARG)
            | target
        )

    def add_Toffoli(self, control1: int, control2: int, target: int) -> None:
        self._append(
            (self.OP_TOFFOLI << self.SHIFT_OP)
            | (control1 << self.SHIFT_ARG)
            | control2
        )
        self._append(target)

    def current_pointer(self) -> int:
        return len(self.chunks) * self.CHUNK_LIMIT + len(self.current_chunk)

    def _word_at(self, absolute_index: int) -> int:
        if absolute_index < 0 or absolute_index >= self.current_pointer():
            raise IndexError("gate-stream word index is out of range")
        chunk_index, offset = divmod(absolute_index, self.CHUNK_LIMIT)
        if chunk_index < len(self.chunks):
            return self.chunks[chunk_index][offset]
        return self.current_chunk[offset]

    def iter_words(self) -> Iterator[int]:
        for chunk in self.chunks:
            yield from chunk
        yield from self.current_chunk

    def get_stats(self) -> Tuple[int, int, int]:
        """Return ``(Toffoli, CNOT, X)`` without changing the stream."""

        toffoli = cnot = x = 0
        for value in self.iter_words():
            opcode = value >> self.SHIFT_OP
            if opcode == self.OP_TOFFOLI:
                toffoli += 1
            elif opcode == self.OP_CNOT:
                cnot += 1
            elif opcode == self.OP_X:
                x += 1
        return toffoli, cnot, x

    def replay_reverse(self, start_idx: int, end_idx: int) -> None:
        """Append the inverse of ``[start_idx, end_idx)`` to the stream."""

        if not 0 <= start_idx <= end_idx <= self.current_pointer():
            raise ValueError("reverse-replay range exceeds the gate stream")
        absolute = end_idx
        while absolute > start_idx:
            value = self._word_at(absolute - 1)
            opcode = value >> self.SHIFT_OP
            if opcode == self.OP_DATA:
                if absolute - start_idx < 2:
                    raise ValueError("reverse range starts inside a Toffoli")
                header = self._word_at(absolute - 2)
                if header >> self.SHIFT_OP != self.OP_TOFFOLI:
                    raise ValueError("Toffoli target has no matching header")
                control1 = (header >> self.SHIFT_ARG) & self.MASK_ARG
                control2 = header & self.MASK_ARG
                self.add_Toffoli(control1, control2, value & self.MASK_ARG)
                absolute -= 2
            elif opcode == self.OP_CNOT:
                control = (value >> self.SHIFT_ARG) & self.MASK_ARG
                self.add_CNOT(control, value & self.MASK_ARG)
                absolute -= 1
            elif opcode == self.OP_X:
                self.add_X(value & self.MASK_ARG)
                absolute -= 1
            elif opcode == self.OP_NOP:
                absolute -= 1
            elif opcode == self.OP_TOFFOLI:
                raise ValueError("reverse range ends inside a Toffoli")
            else:
                raise ValueError(f"unknown opcode {opcode}")

    def _check_range(self, start_idx: int, end_idx: int) -> None:
        if not 0 <= start_idx <= end_idx <= self.current_pointer():
            raise ValueError("simulation range exceeds the gate stream")

    def simulate_range(
        self,
        state: MutableSequence[int],
        start_idx: int,
        end_idx: int,
        *,
        reverse: bool = False,
    ) -> None:
        """Execute a gate interval on a computational-basis state in place."""

        self._check_range(start_idx, end_idx)

        def require_wire(wire: int) -> None:
            if wire < 0 or wire >= len(state):
                raise IndexError(f"gate references wire {wire} outside state")

        def apply_x(target: int) -> None:
            require_wire(target)
            state[target] ^= 1

        def apply_cnot(control: int, target: int) -> None:
            require_wire(control)
            require_wire(target)
            state[target] ^= state[control]

        def apply_toffoli(control1: int, control2: int, target: int) -> None:
            require_wire(control1)
            require_wire(control2)
            require_wire(target)
            state[target] ^= state[control1] & state[control2]

        if reverse:
            absolute = end_idx
            while absolute > start_idx:
                value = self._word_at(absolute - 1)
                opcode = value >> self.SHIFT_OP
                if opcode == self.OP_DATA:
                    if absolute - start_idx < 2:
                        raise ValueError("simulation starts inside a Toffoli")
                    header = self._word_at(absolute - 2)
                    if header >> self.SHIFT_OP != self.OP_TOFFOLI:
                        raise ValueError("Toffoli target has no matching header")
                    apply_toffoli(
                        (header >> self.SHIFT_ARG) & self.MASK_ARG,
                        header & self.MASK_ARG,
                        value & self.MASK_ARG,
                    )
                    absolute -= 2
                elif opcode == self.OP_CNOT:
                    apply_cnot(
                        (value >> self.SHIFT_ARG) & self.MASK_ARG,
                        value & self.MASK_ARG,
                    )
                    absolute -= 1
                elif opcode == self.OP_X:
                    apply_x(value & self.MASK_ARG)
                    absolute -= 1
                elif opcode == self.OP_NOP:
                    absolute -= 1
                elif opcode == self.OP_TOFFOLI:
                    raise ValueError("simulation ends inside a Toffoli")
                else:
                    raise ValueError(f"unknown opcode {opcode}")
            return

        waiting_for_target = False
        control1 = control2 = 0
        for absolute in range(start_idx, end_idx):
            value = self._word_at(absolute)
            opcode = value >> self.SHIFT_OP
            if waiting_for_target:
                if opcode == self.OP_NOP:
                    raise ValueError("only half of a Toffoli was cancelled")
                apply_toffoli(control1, control2, value & self.MASK_ARG)
                waiting_for_target = False
            elif opcode == self.OP_TOFFOLI:
                control1 = (value >> self.SHIFT_ARG) & self.MASK_ARG
                control2 = value & self.MASK_ARG
                waiting_for_target = True
            elif opcode == self.OP_CNOT:
                apply_cnot(
                    (value >> self.SHIFT_ARG) & self.MASK_ARG,
                    value & self.MASK_ARG,
                )
            elif opcode == self.OP_X:
                apply_x(value & self.MASK_ARG)
            elif opcode not in (self.OP_NOP, self.OP_DATA):
                raise ValueError(f"unknown opcode {opcode}")
        if waiting_for_target:
            raise ValueError("simulation range ends inside a Toffoli")

    def simulate(self, state: MutableSequence[int]) -> None:
        self.simulate_range(state, 0, self.current_pointer())

    def count_stream(self, width: int) -> StreamStats:
        """Count the literal stream and perform order-preserving ASAP."""

        if width < 0:
            raise ValueError("width must be nonnegative")
        full_depth = [0] * width
        toffoli_depth = [0] * width
        x = cnot = toffoli = gates = 0
        waiting_for_target = False
        control1 = control2 = 0

        def require_wire(wire: int) -> None:
            if wire < 0 or wire >= width:
                raise IndexError(f"gate references wire {wire} outside width")

        for value in self.iter_words():
            opcode = value >> self.SHIFT_OP
            if waiting_for_target:
                target = value & self.MASK_ARG
                require_wire(control1)
                require_wire(control2)
                require_wire(target)
                depth = max(
                    full_depth[control1], full_depth[control2], full_depth[target]
                ) + 1
                full_depth[control1] = full_depth[control2] = full_depth[target] = depth
                t_depth = max(
                    toffoli_depth[control1],
                    toffoli_depth[control2],
                    toffoli_depth[target],
                ) + 1
                toffoli_depth[control1] = toffoli_depth[control2] = (
                    toffoli_depth[target]
                ) = t_depth
                toffoli += 1
                gates += 1
                waiting_for_target = False
            elif opcode == self.OP_TOFFOLI:
                control1 = (value >> self.SHIFT_ARG) & self.MASK_ARG
                control2 = value & self.MASK_ARG
                waiting_for_target = True
            elif opcode == self.OP_CNOT:
                control = (value >> self.SHIFT_ARG) & self.MASK_ARG
                target = value & self.MASK_ARG
                require_wire(control)
                require_wire(target)
                depth = max(full_depth[control], full_depth[target]) + 1
                full_depth[control] = full_depth[target] = depth
                t_depth = max(toffoli_depth[control], toffoli_depth[target])
                toffoli_depth[control] = toffoli_depth[target] = t_depth
                cnot += 1
                gates += 1
            elif opcode == self.OP_X:
                target = value & self.MASK_ARG
                require_wire(target)
                full_depth[target] += 1
                x += 1
                gates += 1
            elif opcode not in (self.OP_NOP, self.OP_DATA):
                raise ValueError(f"unknown opcode {opcode}")
        if waiting_for_target:
            raise ValueError("truncated Toffoli in gate stream")
        return StreamStats(
            words=self.current_pointer(),
            gates=gates,
            x=x,
            cnot=cnot,
            toffoli=toffoli,
            full_depth=max(full_depth, default=0),
            ordered_toffoli_depth=max(toffoli_depth, default=0),
        )

    def clear(self) -> None:
        self.chunks = []
        self.current_chunk = array.array("Q")

