import random
import sys
import importlib
import shor_out_Montgomery as lib 

# ===========================================
# 0. Global Configuration (Modify TARGET_N here)
# ===========================================
# Optional values: 163, 233, 283, 571
# Note that you should also modify TARGET_N in shor_out_Montgomery at the same time.
TARGET_N = 163
parallel_num =5 # Test the number of parallel branches
POLY_CONFIG = {
    163: (163, (1<<7) | (1<<6) | (1<<3) | 1),
    233: (233, (1<<74) | 1),
    283: (283, (1<<12) | (1<<7) | (1<<5) | 1),
    571: (571, (1<<10) | (1<<5) | (1<<2) | 1)
}

def reload_backend(target_n):
    if target_n not in POLY_CONFIG:
        raise ValueError(f"Unsupported TARGET_N: {target_n}")
        
    print(f"[{'='*10} Loading Backend for GF(2^{target_n}) {'='*10}]")
    
    gates_mod_name = f"init.basic_gates_{target_n}"
    inv_mod_name = f"init.inversion_{target_n}"
    
    try:
        gates_mod = importlib.import_module(gates_mod_name)
        inv_mod = importlib.import_module(inv_mod_name)
    except ImportError as e:
        print(f"Error: Could not import backend modules for {target_n}.")
        print(f"Ensure {gates_mod_name}.py and {inv_mod_name}.py exist.")
        raise e

    lib.GF_N = gates_mod.GF_N
    lib.KARATSUBA_ANC_BASE = gates_mod.KARATSUBA_ANC_BASE
    lib.TOFFOLI_BASE = getattr(gates_mod, 'TOFFOLI_BASE', 0) 
    lib.gm = gates_mod.gm
    lib.gates_mod = gates_mod 

    functions_to_sync = [
        'CNOT_n', 'CONST_ADD_n', 'Reduction', 'recursive_karatsuba', 
        'Squaring_plus_input', 'Square', 'copy_parallel', 'CSWAP', 'combine'
    ]
    for func in functions_to_sync:
        if hasattr(gates_mod, func):
            setattr(lib, func, getattr(gates_mod, func))

    if hasattr(inv_mod, 'Inverison_Itoh_Tsujii_based'):
        lib.Inverison_Itoh_Tsujii_based = inv_mod.Inverison_Itoh_Tsujii_based
        
    print(f"Backend loaded. N={lib.GF_N}, Karatsuba_Base={lib.KARATSUBA_ANC_BASE}")

def list_to_int(bit_list):
    val = 0
    for i, bit in enumerate(bit_list):
        if bit: val |= (1 << i)
    return val

def int_to_list(val, length):
    return [(val >> i) & 1 for i in range(length)]

def gf_mult(a, b):
    p = 0
    while b > 0:
        if b & 1: p ^= a
        a <<= 1
        b >>= 1
    return p

def gf_sq(a):
    return gf_mult(a, a)

def gf_mod(val):
    N, poly_mask = POLY_CONFIG[TARGET_N]
    for i in range(val.bit_length() - 1, N - 1, -1):
        if (val >> i) & 1:
            val ^= (1 << i) 
            val ^= (poly_mask << (i - N)) 
    return val

def gf_inv(val):
    if val == 0: return 0
    res = 1
    base = val
    exp = (1 << TARGET_N) - 2
    while exp > 0:
        if exp % 2 == 1:
            res = gf_mod(gf_mult(res, base))
        base = gf_mod(gf_mult(base, base))
        exp //= 2
    return res

def ec_add(x1, y1, x2, y2, a_param):
    if x1 == x2 and y1 == y2: return 0, 0 

    top = y1 ^ y2
    bot = x1 ^ x2
    if bot == 0: return 0, 0 
    
    inv_bot = gf_inv(bot)
    lam = gf_mod(gf_mult(top, inv_bot))

    lam2 = gf_mod(gf_sq(lam))
    x3 = lam2 ^ lam ^ x1 ^ x2 ^ a_param

    term = gf_mod(gf_mult(lam, x2 ^ x3))
    y3 = term ^ x3 ^ y2
    
    return x3, y3

def shor_reduction_level_sim(parallel_num, x1_list, y1_list, x2_list, y2_list, ancilla, x1x2_list, ancilla2, sqr, lambda_res, mult_res, inv_res):

    a = 3 
    num = lib.KARATSUBA_ANC_BASE 
    n = lib.GF_N                  
    
    # === Step 1: Initialize Inputs ===
    for i in range(parallel_num):
        lib.CNOT_n(x1_list[i], x1x2_list[i])
        lib.CNOT_n(x2_list[i], x1x2_list[i])
        lib.CNOT_n(y2_list[i], y1_list[i])
    
    # === Step 2: Compute Lambda ===
    count = 0
    lib.gates_mod.cnt = 0 
    lamba_flat, _ = lib.Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, lambda_res, inv_res)
    
    lamba_list = []
    for k in range(parallel_num):
        lamba_list.append(lamba_flat[k*n : (k+1)*n])
    lib.gates_mod.cnt = 0 
    
    # === Step 3: Compute Coordinates ===
    for i in range(parallel_num):
        lib.CONST_ADD_n(x1x2_list[i], a)
        lamba_list[i], x1x2_list[i] = lib.Squaring_plus_input(lamba_list[i] + x1x2_list[i], n) 
        
        lib.CNOT_n(x2_list[i], x1x2_list[i]) 
        
        start_mul = lib.gm.current_pointer()
        count = 0
        result2, count, ancilla[num*i:num*(i+1)] = lib.recursive_karatsuba(x1x2_list[i], lamba_list[i], n, count, ancilla[num*i:num*(i+1)])  
        result2 = lib.Reduction(result2)
        end_mul = lib.gm.current_pointer()

        target_y = mult_res[i*n : (i+1)*n] 
        lib.CNOT_n(result2, target_y) 
        lib.gm.replay_reverse(start_mul, end_mul)
        
        lib.CNOT_n(y2_list[i], target_y)
        lib.CNOT_n(x2_list[i], x1x2_list[i]) 
        lib.CNOT_n(x1x2_list[i], target_y)

# ==========================================
#Simulator
# ==========================================
def simulate_circuit(gate_manager, initial_state, total_qubits):
    state = bytearray(initial_state)
    if len(state) < total_qubits:
        state.extend([0] * (total_qubits - len(state)))
    
    all_chunks = gate_manager.chunks + [gate_manager.current_chunk]
    SHIFT_OP, SHIFT_ARG = gate_manager.SHIFT_OP, gate_manager.SHIFT_ARG
    MASK_ARG = (1 << SHIFT_ARG) - 1
    OP_X, OP_CNOT, OP_TOFFOLI, OP_NOP = gate_manager.OP_X, gate_manager.OP_CNOT, gate_manager.OP_TOFFOLI, gate_manager.OP_NOP
    
    pending_c1, pending_c2 = None, None
    waiting_data = False
    
    print(f"Running simulation on {len(all_chunks)} chunks...")
    
    for chunk in all_chunks:
        for val in chunk:
            op = val >> SHIFT_OP
            if op == OP_NOP: continue
                
            if waiting_data:
                t = val & MASK_ARG
                if state[pending_c1] and state[pending_c2]:
                    state[t] ^= 1
                waiting_data = False
            elif op == OP_CNOT:
                t, c = val & MASK_ARG, (val >> SHIFT_ARG) & MASK_ARG
                if state[c]: state[t] ^= 1
            elif op == OP_TOFFOLI:
                pending_c2, pending_c1 = val & MASK_ARG, (val >> SHIFT_ARG) & MASK_ARG
                waiting_data = True
            elif op == OP_X:
                state[val & MASK_ARG] ^= 1
                
    return state

def test_phase3_reduction(parallel_num):
    reload_backend(TARGET_N)
    
    n = lib.GF_N
    print(f"\n{'='*20} Testing Phase 3 (Reduction) Logic for GF(2^{n}) {'='*20}")

    lib.n = n
    lib.cnt = 0
    lib.gates_mod.cnt = 0
    lib.nums1 = lib.KARATSUBA_ANC_BASE * parallel_num + 500 
    lib.gm.clear()

    idx = 0
    
    x1_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    x2_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    x1x2_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    y1_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    y2_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    # Ancillas
    num_ancilla = lib.nums1
    ancilla_indices = list(range(idx, idx + num_ancilla))
    idx += num_ancilla
    # Toffoli Qubits (Safe estimate)
    num_toffoli = 16000000 * parallel_num 
    lib.gates_mod.Toffoli_qubits = list(range(idx, idx + num_toffoli))
    idx += num_toffoli
    
    sqr_indices = list(range(idx, idx + 4*n))
    idx += 4*n
    
    ancilla2_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num

    # Results 
    lambda_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    mult_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    inv_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    
    total_qubits = idx
    print(f"Total Qubits Allocated: {total_qubits}")

    inputs = []
    print("\nGenerating Random Inputs (P1 + P2):")
    for i in range(parallel_num):
        x1 = random.randint(1, (1<<n)-1)
        y1 = random.randint(1, (1<<n)-1)
        x2 = random.randint(1, (1<<n)-1)
        y2 = random.randint(1, (1<<n)-1)
        
        inputs.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
        print(f"  Branch {i}: P1=({hex(x1)[:6]}...), P2=({hex(x2)[:6]}...)")

    print("\nBuilding Circuit (shor_reduction_level)...")
    shor_reduction_level_sim(
        parallel_num,
        x1_indices, y1_indices,
        x2_indices, y2_indices,
        ancilla_indices,
        x1x2_indices,
        ancilla2_indices,
        sqr_indices,
        lambda_res_indices,
        mult_res_indices,
        inv_res_indices
    )
    
    stats = lib.gm.get_stats()
    print(f"Circuit Built. Total Toffoli: {stats[0]}, Total CNOT: {stats[1]}")

    print("\n" + "="*40)
    print("RUNNING SIMULATION")
    print("="*40)
    
    initial_state = [0] * total_qubits

    for i in range(parallel_num):
        for reg_indices, val in [
            (x1_indices[i], inputs[i]['x1']),
            (y1_indices[i], inputs[i]['y1']),
            (x2_indices[i], inputs[i]['x2']),
            (y2_indices[i], inputs[i]['y2'])
        ]:
            for b, bit in enumerate(int_to_list(val, n)):
                if bit: initial_state[reg_indices[b]] = 1
    
    final_state = simulate_circuit(lib.gm, initial_state, total_qubits)

    pass_count = 0
    CURVE_A = 3
    
    print("\nVerifying Results:")
    for i in range(parallel_num):
        exp_x3, exp_y3 = ec_add(inputs[i]['x1'], inputs[i]['y1'], 
                                inputs[i]['x2'], inputs[i]['y2'], CURVE_A)

        res_x = list_to_int([final_state[q] for q in x1x2_indices[i]])

        current_y_indices = mult_res_indices[i*n : (i+1)*n]
        res_y = list_to_int([final_state[q] for q in current_y_indices])
        
        match = (res_x == exp_x3) and (res_y == exp_y3)
        status = "✅ PASS" if match else "❌ FAIL"
        if match: pass_count += 1
        
        print(f"  Branch {i}: {status}")
        if not match:
            print(f"    P1: ({hex(inputs[i]['x1'])}, {hex(inputs[i]['y1'])})")
            print(f"    P2: ({hex(inputs[i]['x2'])}, {hex(inputs[i]['y2'])})")
            print(f"    Exp: ({hex(exp_x3)}, {hex(exp_y3)})")
            print(f"    Got: ({hex(res_x)}, {hex(res_y)})")
            print(f"    Diff: ({hex(res_x ^ exp_x3)}, {hex(res_y ^ exp_y3)})")
            
    if pass_count == parallel_num:
        print(f"\nAll {parallel_num} branches Passed Phase 3 Logic for GF(2^{n})!")
    else:
        print(f"\nWarning: Only {pass_count}/{parallel_num} branches Passed.")

if __name__ == "__main__":
    test_phase3_reduction(parallel_num)