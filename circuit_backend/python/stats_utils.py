"""Global cancellation and resource estimation for a complete NCT stream.
"""

from __future__ import annotations

import array


def get_exact_resources_optimized(gm, width):
    """Optimize ``gm`` in place and return resource statistics.

    Returns ``(stats, full_depth, current_depth, ordered_toffoli_depth)``.
    The last value is the weighted dependency depth of the emitted gate order;
    it is deliberately not named a globally minimized Toffoli depth.
    """

    if width < 0:
        raise ValueError("width must be nonnegative")

    shift_op = gm.SHIFT_OP
    shift_arg = gm.SHIFT_ARG
    mask_arg = gm.MASK_ARG
    op_x = gm.OP_X
    op_cnot = gm.OP_CNOT
    op_toffoli = gm.OP_TOFFOLI
    op_nop = gm.OP_NOP
    val_nop = (op_nop << shift_op) | mask_arg

    def word_at(absolute):
        return gm._word_at(absolute)  # pylint: disable=protected-access

    def set_word(absolute, value):
        chunk_index, offset = divmod(absolute, gm.CHUNK_LIMIT)
        if chunk_index < len(gm.chunks):
            gm.chunks[chunk_index][offset] = value
        else:
            gm.current_chunk[offset] = value

    total_words = gm.current_pointer()

    # ------------------------------------------------------------------
    # Phase 1: exact cancellation in the literal gate stream.
    # ------------------------------------------------------------------
    wire_stacks = [[] for _ in range(width)]
    pending = None

    for absolute in range(total_words):
        value = word_at(absolute)
        opcode = value >> shift_op
        if pending is not None:
            control1, control2, header_ptr = pending
            pending = None
            target = value & mask_arg
            stacks = (
                wire_stacks[control1],
                wire_stacks[control2],
                wire_stacks[target],
            )
            if all(stacks) and stacks[0][-1] == stacks[1][-1] == stacks[2][-1]:
                old_target_ptr = stacks[0][-1]
                old_header_ptr = old_target_ptr - 1
                old_header = word_at(old_header_ptr)
                old_target = word_at(old_target_ptr)
                old_controls = (
                    (old_header >> shift_arg) & mask_arg,
                    old_header & mask_arg,
                )
                if (
                    old_header >> shift_op == op_toffoli
                    and sorted(old_controls) == sorted((control1, control2))
                    and (old_target & mask_arg) == target
                ):
                    set_word(old_header_ptr, val_nop)
                    set_word(old_target_ptr, val_nop)
                    set_word(header_ptr, val_nop)
                    set_word(absolute, val_nop)
                    stacks[0].pop()
                    stacks[1].pop()
                    stacks[2].pop()
                    continue
            stacks[0].append(absolute)
            stacks[1].append(absolute)
            stacks[2].append(absolute)
        elif opcode == op_cnot:
            control = (value >> shift_arg) & mask_arg
            target = value & mask_arg
            control_stack = wire_stacks[control]
            target_stack = wire_stacks[target]
            if (
                control_stack
                and target_stack
                and control_stack[-1] == target_stack[-1]
            ):
                old_ptr = control_stack[-1]
                old_value = word_at(old_ptr)
                if (
                    old_value >> shift_op == op_cnot
                    and ((old_value >> shift_arg) & mask_arg) == control
                    and (old_value & mask_arg) == target
                ):
                    set_word(old_ptr, val_nop)
                    set_word(absolute, val_nop)
                    control_stack.pop()
                    target_stack.pop()
                    continue
            control_stack.append(absolute)
            target_stack.append(absolute)
        elif opcode == op_toffoli:
            pending = (
                (value >> shift_arg) & mask_arg,
                value & mask_arg,
                absolute,
            )
        elif opcode == op_x:
            target = value & mask_arg
            target_stack = wire_stacks[target]
            if target_stack and word_at(target_stack[-1]) >> shift_op == op_x:
                set_word(target_stack[-1], val_nop)
                set_word(absolute, val_nop)
                target_stack.pop()
                continue
            target_stack.append(absolute)
        elif opcode != op_nop:
            raise ValueError(f"unknown opcode {opcode} in cancellation")
    if pending is not None:
        raise ValueError("truncated Toffoli in cancellation")

    # Release the dominant cancellation workspace before allocating buckets.
    del wire_stacks

    # ------------------------------------------------------------------
    # Phase 2: literal full depth and weighted Toffoli dependency depth.
    # ------------------------------------------------------------------
    wire_tf_depths = array.array("Q", [0]) * width
    wire_full_depths = array.array("Q", [0]) * width
    stats = {"X_count": 0, "CNOT_count": 0, "Toffoli_count": 0}
    clifford_buckets = []
    mt_qubit_buckets = []

    def ensure_buckets(depth):
        while len(clifford_buckets) <= depth:
            clifford_buckets.append(array.array("Q"))
        while len(mt_qubit_buckets) <= depth + 1:
            mt_qubit_buckets.append(set())

    pending_controls = None
    for value in gm.iter_words():
        opcode = value >> shift_op
        if opcode == op_nop:
            continue
        if pending_controls is not None:
            control1, control2 = pending_controls
            pending_controls = None
            target = value & mask_arg
            stats["Toffoli_count"] += 1
            tf_depth = max(
                wire_tf_depths[control1],
                wire_tf_depths[control2],
                wire_tf_depths[target],
            ) + 1
            wire_tf_depths[control1] = tf_depth
            wire_tf_depths[control2] = tf_depth
            wire_tf_depths[target] = tf_depth
            ensure_buckets(tf_depth)
            mt_qubit_buckets[tf_depth].update((control1, control2, target))

            full_depth = max(
                wire_full_depths[control1],
                wire_full_depths[control2],
                wire_full_depths[target],
            ) + 1
            wire_full_depths[control1] = full_depth
            wire_full_depths[control2] = full_depth
            wire_full_depths[target] = full_depth
        elif opcode == op_toffoli:
            pending_controls = (
                (value >> shift_arg) & mask_arg,
                value & mask_arg,
            )
        elif opcode == op_cnot:
            control = (value >> shift_arg) & mask_arg
            target = value & mask_arg
            stats["CNOT_count"] += 1
            tf_depth = max(wire_tf_depths[control], wire_tf_depths[target])
            wire_tf_depths[control] = tf_depth
            wire_tf_depths[target] = tf_depth
            ensure_buckets(tf_depth)
            clifford_buckets[tf_depth].append(value)

            full_depth = max(
                wire_full_depths[control], wire_full_depths[target]
            ) + 1
            wire_full_depths[control] = full_depth
            wire_full_depths[target] = full_depth
        elif opcode == op_x:
            target = value & mask_arg
            stats["X_count"] += 1
            ensure_buckets(wire_tf_depths[target])
            clifford_buckets[wire_tf_depths[target]].append(value)
            wire_full_depths[target] += 1
        elif opcode != gm.OP_DATA:
            raise ValueError(f"unknown opcode {opcode} in depth calculation")
    if pending_controls is not None:
        raise ValueError("truncated Toffoli in depth calculation")

    full_depth = max(wire_full_depths, default=0)
    ordered_toffoli_depth = max(wire_tf_depths, default=0)
    del wire_full_depths
    del wire_tf_depths

    # ------------------------------------------------------------------
    # Phase 3: Toffoli-layer reconstructed depth.
    # ------------------------------------------------------------------
    wire_current_depths = array.array("Q", [0]) * width
    for depth in range(max(len(clifford_buckets), len(mt_qubit_buckets))):
        if depth < len(clifford_buckets):
            for value in clifford_buckets[depth]:
                opcode = value >> shift_op
                target = value & mask_arg
                if opcode == op_cnot:
                    control = (value >> shift_arg) & mask_arg
                    new_depth = max(
                        wire_current_depths[control],
                        wire_current_depths[target],
                    ) + 1
                    wire_current_depths[control] = new_depth
                    wire_current_depths[target] = new_depth
                elif opcode == op_x:
                    wire_current_depths[target] += 1

        target_layer = depth + 1
        if (
            target_layer < len(mt_qubit_buckets)
            and mt_qubit_buckets[target_layer]
        ):
            wires = mt_qubit_buckets[target_layer]
            new_depth = max(wire_current_depths[wire] for wire in wires) + 1
            for wire in wires:
                wire_current_depths[wire] = new_depth

    current_depth = max(wire_current_depths, default=0)
    return stats, full_depth, current_depth, ordered_toffoli_depth

