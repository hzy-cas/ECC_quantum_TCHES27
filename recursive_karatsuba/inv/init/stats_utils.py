import array
import collections

def get_exact_resources_optimized(gm, width):

    
    SHIFT_OP = gm.SHIFT_OP
    SHIFT_ARG = gm.SHIFT_ARG
    MASK_ARG = (1 << SHIFT_ARG) - 1
    OP_X, OP_CNOT, OP_TOFFOLI, OP_NOP = gm.OP_X, gm.OP_CNOT, gm.OP_TOFFOLI, gm.OP_NOP
    VAL_NOP = (OP_NOP << SHIFT_OP) | MASK_ARG 

    all_chunks = gm.chunks + [gm.current_chunk] if gm.current_chunk else gm.chunks

    # ==========================
    # Phase 1: Gate cancellation
    # ==========================
    wire_stacks = [[] for _ in range(width)]
    pending_ptr = None
    
    for c_idx, chunk in enumerate(all_chunks):
        c_idx_shifted = c_idx << 32
        for i, val in enumerate(chunk):
            op = val >> SHIFT_OP
            if pending_ptr is not None:
                h_c1, h_c2, h_ptr = pending_ptr
                pending_ptr = None
                t = val & MASK_ARG
                t_ptr = c_idx_shifted | i
                s_c1, s_c2, s_t = wire_stacks[h_c1], wire_stacks[h_c2], wire_stacks[t]
                if s_c1 and s_c2 and s_t and s_c1[-1] == s_c2[-1] == s_t[-1]:
                    old_ptr = s_c1[-1]
                    oc, oo = old_ptr >> 32, old_ptr & 0xFFFFFFFF
                    if (all_chunks[oc][oo] & MASK_ARG) == t:
                        all_chunks[oc][oo] = all_chunks[oc][oo-1] = VAL_NOP
                        chunk[h_ptr & 0xFFFFFFFF] = chunk[i] = VAL_NOP
                        s_c1.pop(); s_c2.pop(); s_t.pop()
                        continue
                s_c1.append(t_ptr); s_c2.append(t_ptr); s_t.append(t_ptr)
            elif op == OP_CNOT:
                t, c = val & MASK_ARG, (val >> SHIFT_ARG) & MASK_ARG
                s_c, s_t = wire_stacks[c], wire_stacks[t]
                if s_c and s_t and s_c[-1] == s_t[-1]:
                    old_ptr = s_c[-1]
                    oc, oo = old_ptr >> 32, old_ptr & 0xFFFFFFFF
                    old_val = all_chunks[oc][oo]
        
                    if (old_val >> SHIFT_OP) == OP_CNOT and (old_val & MASK_ARG) == t:
                        all_chunks[oc][oo] = chunk[i] = VAL_NOP
                        s_c.pop(); s_t.pop()
                        continue
                wire_stacks[c].append(c_idx_shifted | i)
                wire_stacks[t].append(c_idx_shifted | i)
            elif op == OP_TOFFOLI:
                pending_ptr = ((val >> SHIFT_ARG) & MASK_ARG, val & MASK_ARG, c_idx_shifted | i)
            elif op == OP_X:
                t = val & MASK_ARG
                s_t = wire_stacks[t]
                if s_t and (all_chunks[s_t[-1]>>32][s_t[-1]&0xFFFFFFFF] >> SHIFT_OP) == OP_X:
                    all_chunks[s_t[-1]>>32][s_t[-1]&0xFFFFFFFF] = chunk[i] = VAL_NOP
                    s_t.pop()
                    continue
                wire_stacks[t].append(c_idx_shifted | i)

    # ==========================
    # Phase 2: depth
    # ==========================
    wire_tf_depths = array.array('I', [0] * width)
    wires_f_depth = array.array('I', [0] * width)
    stats = {'CNOT_count': 0, 'Toffoli_count': 0}
    
    clifford_buckets = [] 
    mt_qubit_buckets = [] 
    
    def ensure_buckets(d):
        while len(clifford_buckets) <= d: clifford_buckets.append(array.array('Q'))
        while len(mt_qubit_buckets) <= d + 1: mt_qubit_buckets.append(set())

    pending_tf = None
    for chunk in all_chunks:
        for val in chunk:
            op = val >> SHIFT_OP
            if op == OP_NOP: continue
            
            if pending_tf:
                c1, c2 = pending_tf
                t = val & MASK_ARG
                pending_tf = None
                stats['Toffoli_count'] += 1

                tf_d = max(wire_tf_depths[c1], wire_tf_depths[c2], wire_tf_depths[t]) + 1
                wire_tf_depths[c1] = wire_tf_depths[c2] = wire_tf_depths[t] = tf_d
                ensure_buckets(tf_d)
                mt_qubit_buckets[tf_d].update([c1, c2, t])
                f_mx = max(wires_f_depth[c1], wires_f_depth[c2], wires_f_depth[t]) + 1
                wires_f_depth[c1] = wires_f_depth[c2] = wires_f_depth[t] = f_mx
            elif op == OP_CNOT:
                stats['CNOT_count'] += 1
                c, t = (val >> SHIFT_ARG) & MASK_ARG, val & MASK_ARG
                tf_d = max(wire_tf_depths[c], wire_tf_depths[t])
                wire_tf_depths[c] = wire_tf_depths[t] = tf_d
                ensure_buckets(tf_d)
                clifford_buckets[tf_d].append(val)
                f_mx = max(wires_f_depth[c], wires_f_depth[t]) + 1
                wires_f_depth[c] = wires_f_depth[t] = f_mx
            elif op == OP_TOFFOLI:
                pending_tf = ((val >> SHIFT_ARG) & MASK_ARG, val & MASK_ARG)
            elif op == OP_X:
                t = val & MASK_ARG
                tf_d = wire_tf_depths[t]
                ensure_buckets(tf_d)
                clifford_buckets[tf_d].append(val)
                wires_f_depth[t] += 1

    wires_c_depth = array.array('I', [0] * width)
    limit = max(len(clifford_buckets), len(mt_qubit_buckets))
    
    for d in range(limit):
        if d < len(clifford_buckets):
            for val in clifford_buckets[d]:
                op = val >> SHIFT_OP
                if op == OP_CNOT:
                    c, t = (val >> SHIFT_ARG) & MASK_ARG, val & MASK_ARG
                    new_v = max(wires_c_depth[c], wires_c_depth[t]) + 1
                    wires_c_depth[c] = wires_c_depth[t] = new_v
                elif op == OP_X:
                    wires_c_depth[val & MASK_ARG] += 1
                    
        target_mt = d + 1
        if target_mt < len(mt_qubit_buckets) and mt_qubit_buckets[target_mt]:
            qubits = mt_qubit_buckets[target_mt]
            d_max = max(wires_c_depth[q] for q in qubits) + 1
            for q in qubits:
                wires_c_depth[q] = d_max

    toffoli_depth = max(wire_tf_depths) if width > 0 else 0
    return stats, max(wires_f_depth) if width > 0 else 0, max(wires_c_depth) if width > 0 else 0, toffoli_depth