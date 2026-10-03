import random
import sys
import importlib
import shor_out_Montgomery as lib 

# ===========================================
# 0. Global Configuration (Modify TARGET_N here)
# ===========================================
# Optional values: 163, 233, 283, 571
# Note that you should also modify TARGET_N in shor_out_Montgomery at the same time.
TARGET_N = 283
parallel_num =1 # Test the number of parallel branches
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
        
    print(f"Backend loaded. N={lib.GF_N}")

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
    # a^(2^N - 2)
    exp = (1 << TARGET_N) - 2
    while exp > 0:
        if exp % 2 == 1:
            res = gf_mod(gf_mult(res, base))
        base = gf_mod(gf_mult(base, base))
        exp //= 2
    return res

def ec_add(x1, y1, x2, y2, a_param):
    top = y1 ^ y2
    bot = x1 ^ x2
    if bot == 0: return 0, 0
    inv_bot = gf_inv(bot)
    lam = gf_mod(gf_mult(top, inv_bot))

    lam2 = gf_mod(gf_sq(lam))
    x3 = lam2 ^ lam ^ x1 ^ x2 ^ a_param

    term = gf_mod(gf_mult(lam, x1 ^ x3))
    y3 = term ^ x3 ^ y1
    return x3, y3

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

def test_multi_branch(parallel_num):
    reload_backend(TARGET_N)

    n = lib.GF_N
    print(f"\n{'='*20} Testing {parallel_num} Parallel Branches for GF(2^{n}) {'='*20}")

    lib.n = n
    lib.cnt = 0
    lib.gates_mod.cnt = 0
    lib.nums1 = lib.KARATSUBA_ANC_BASE * parallel_num + 500 
    lib.gm.clear()

    idx = 0

    q_indices = list(range(idx, idx + parallel_num))
    idx += parallel_num

    x1_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    
    x1x2_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    
    y1_indices = [list(range(idx + i*n, idx + (i+1)*n)) for i in range(parallel_num)]
    idx += parallel_num * n
    
    # Ancilla (Karatsuba)
    num_ancilla = lib.nums1
    ancilla_indices = list(range(idx, idx + num_ancilla))
    idx += num_ancilla
    
    # Toffoli Ancilla
    num_toffoli = 1600000 * parallel_num 
    lib.gates_mod.Toffoli_qubits = list(range(idx, idx + num_toffoli))
    idx += num_toffoli
    
    # Sqr, Ancilla2
    sqr_indices = list(range(idx, idx + 4*n))
    idx += 4*n
    ancilla2_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    
    # Result buffers (Defined as FLAT LISTS)
    lambda_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    mult_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    inv_res_indices = list(range(idx, idx + n * parallel_num))
    idx += n * parallel_num
    
    # CSWAP Ancilla
    cswap_size = (2 * n - 1) * parallel_num
    cswap_ancilla_indices = list(range(idx, idx + cswap_size))
    idx += cswap_size
    
    total_qubits = idx
    print(f"Total Qubits Allocated: {total_qubits}")

    inputs = []
    x2_list_val = []
    y2_list_val = []
    
    print("\nGenerating Random Inputs:")
    for i in range(parallel_num):
        x1 = random.randint(1, (1<<n)-1)
        y1 = random.randint(1, (1<<n)-1)
        x2 = random.randint(1, (1<<n)-1)
        y2 = random.randint(1, (1<<n)-1)
        
        inputs.append({'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
        x2_list_val.append(x2)
        y2_list_val.append(y2)
        
        print(f"  Branch {i}: P1=({hex(x1)[:6]}...), P2=({hex(x2)[:6]}...)")

    print("\nBuilding Circuit (shor_level)...")
    lib.shor_level(
        parallel_num, 
        q_indices,          # q control bits
        x1_indices,         # x1 registers (List of lists)
        y1_indices,         # y1 registers (List of lists)
        x2_list_val,        # x2 constants (List)
        y2_list_val,        # y2 constants (List)
        ancilla_indices, 
        x1x2_indices,       
        ancilla2_indices, 
        sqr_indices, 
        0, 
        lambda_res_indices, 
        mult_res_indices,   
        inv_res_indices,
        cswap_ancilla_indices
    )
    
    stats = lib.gm.get_stats()
    print(f"Circuit Built. Total Toffoli: {stats[0]}, Total CNOT: {stats[1]}")

    CURVE_A = 3

    # Case: q = 1 
    print("\n" + "="*60)
    print("TEST CASE 1: q = 1 (Expect P1 <- P1 + P2)")
    print("="*60)
    
    initial_state = [0] * total_qubits
    
    for i in range(parallel_num):
        initial_state[q_indices[i]] = 1
        
        val_x1 = inputs[i]['x1']
        val_y1 = inputs[i]['y1']
        for b, bit in enumerate(int_to_list(val_x1, n)): 
            if bit: initial_state[x1_indices[i][b]] = 1
        for b, bit in enumerate(int_to_list(val_y1, n)): 
            if bit: initial_state[y1_indices[i][b]] = 1

    final_state = simulate_circuit(lib.gm, initial_state, total_qubits)
    pass_count = 0
    print("\nVerifying Result (q=1):")
    for i in range(parallel_num):
        exp_x3, exp_y3 = ec_add(inputs[i]['x1'], inputs[i]['y1'], 
                                inputs[i]['x2'], inputs[i]['y2'], CURVE_A)
        
        res_x = list_to_int([final_state[q] for q in x1_indices[i]])
        res_y = list_to_int([final_state[q] for q in y1_indices[i]])
        
        match = (res_x == exp_x3) and (res_y == exp_y3)
        status = "✅ PASS" if match else "❌ FAIL"
        if match: pass_count += 1
        
        print(f"  Branch {i}: {status}")
        if not match:
            print(f"    Exp: ({hex(exp_x3)}, {hex(exp_y3)})")
            print(f"    Got: ({hex(res_x)}, {hex(res_y)})")

    # Case: q = 0 
    print("\n" + "="*60)
    print("TEST CASE 2: q = 0 (Expect P1 Unchanged, Sum in Temp)")
    print("="*60)
    
    initial_state_q0 = [0] * total_qubits
    
    for i in range(parallel_num):
        initial_state_q0[q_indices[i]] = 0  # <--- Set Control to 0
        
        val_x1 = inputs[i]['x1']
        val_y1 = inputs[i]['y1']
        for b, bit in enumerate(int_to_list(val_x1, n)): 
            if bit: initial_state_q0[x1_indices[i][b]] = 1
        for b, bit in enumerate(int_to_list(val_y1, n)): 
            if bit: initial_state_q0[y1_indices[i][b]] = 1

    final_state_q0 = simulate_circuit(lib.gm, initial_state_q0, total_qubits)
    
    pass_count_q0 = 0
    print("\nVerifying Result (q=0):")
    for i in range(parallel_num):
        input_x1 = inputs[i]['x1']
        input_y1 = inputs[i]['y1']
        sum_x3, sum_y3 = ec_add(inputs[i]['x1'], inputs[i]['y1'], 
                                inputs[i]['x2'], inputs[i]['y2'], CURVE_A)

        res_x1_reg = list_to_int([final_state_q0[q] for q in x1_indices[i]])
        res_y1_reg = list_to_int([final_state_q0[q] for q in y1_indices[i]])
        
        res_sum_x = list_to_int([final_state_q0[q] for q in x1x2_indices[i]])
        
        res_sum_y = list_to_int([final_state_q0[q] for q in mult_res_indices[i*n : (i+1)*n]])

        check_input_kept = (res_x1_reg == input_x1) and (res_y1_reg == input_y1)
        check_sum_calc   = (res_sum_x == sum_x3) and (res_sum_y == sum_y3)
        
        if check_input_kept and check_sum_calc:
            status = "✅ PASS"
            pass_count_q0 += 1
        else:
            status = "❌ FAIL"
            
        print(f"  Branch {i}: {status}")
        if not check_input_kept:
            print(f"    [Error] Input NOT kept.")
            print(f"    Exp: ({hex(input_x1)}, {hex(input_y1)})")
            print(f"    Got: ({hex(res_x1_reg)}, {hex(res_y1_reg)})")
        if not check_sum_calc:
            print(f"    [Error] Sum NOT found in temp regs.")
            print(f"    Exp Sum: ({hex(sum_x3)}, {hex(sum_y3)})")
            print(f"    Got Sum: ({hex(res_sum_x)}, {hex(res_sum_y)})")

    if pass_count == parallel_num and pass_count_q0 == parallel_num:
        print(f"\nAll tests passed successfully for GF(2^{n})!")
    else:
        print(f"\nWarning: Some tests failed.")

if __name__ == "__main__":
    test_multi_branch(parallel_num)