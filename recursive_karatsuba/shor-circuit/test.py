import random
import sys
import importlib
import shor_out_Montgomery as lib

# ==========================================

# ==========================================
TARGET_N = 163
parallel_num = 1
POLY_CONFIG = {
    163: (163, (1<<7) | (1<<6) | (1<<3) | 1),
    233: (233, (1<<74) | 1),
    283: (283, (1<<12) | (1<<7) | (1<<5) | 1),
    571: (571, (1<<10) | (1<<5) | (1<<2) | 1)
}

# ==========================================

# ==========================================
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

# ==========================================

# ==========================================
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

# ==========================================

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

# ==========================================

# ==========================================
def debug_internal_steps(parallel_num):
    reload_backend(TARGET_N)
    n = lib.GF_N
    print(f"\n{'='*20} DEBUGGING INTERNAL STEPS (N={n}) {'='*20}")


    x1 = random.randint(1, (1<<n)-1)
    y1 = random.randint(1, (1<<n)-1)
    x2 = random.randint(1, (1<<n)-1)
    y2 = random.randint(1, (1<<n)-1)

    print(f"Inputs:\n  x1={hex(x1)}\n  x2={hex(x2)}\n  y1={hex(y1)}\n  y2={hex(y2)}")

    # ----------------------------------------------------


    # ----------------------------------------------------
    print(f"\n[Step A] Verifying Lambda Calculation...")
    lib.gm.clear()
    lib.gates_mod.cnt = 0


    idx = 0
    x1x2_reg = list(range(idx, idx + n)); idx += n
    y1_reg   = list(range(idx, idx + n)); idx += n

    lib.nums1 = lib.KARATSUBA_ANC_BASE * parallel_num + 1000
    ancilla = list(range(idx, idx + lib.nums1)); idx += lib.nums1
    sqr = list(range(idx, idx + 4*n)); idx += 4*n
    ancilla2 = list(range(idx, idx + n)); idx += n
    lambda_res = list(range(idx, idx + n)); idx += n
    inv_res = list(range(idx, idx + n)); idx += n # Temporary inv
    # Toffoli
    num_toffoli = 500000
    lib.gates_mod.Toffoli_qubits = list(range(idx, idx + num_toffoli)); idx += num_toffoli
    total_qubits = idx



    val_Xsum = x1 ^ x2
    val_Ysum = y1 ^ y2

    initial_state = [0] * total_qubits
    for b, bit in enumerate(int_to_list(val_Xsum, n)):
        if bit: initial_state[x1x2_reg[b]] = 1
    for b, bit in enumerate(int_to_list(val_Ysum, n)):
        if bit: initial_state[y1_reg[b]] = 1




    lib.Montgomerytrick([x1x2_reg], [y1_reg], 0, ancilla, ancilla2, sqr, lambda_res, inv_res)


    final_state = simulate_circuit(lib.gm, initial_state, total_qubits)


    res_lambda = list_to_int([final_state[q] for q in lambda_res])


    term_inv = gf_inv(val_Xsum)
    exp_lambda = gf_mod(gf_mult(val_Ysum, term_inv))

    print(f"  Exp Lambda: {hex(exp_lambda)}")
    print(f"  Got Lambda: {hex(res_lambda)}")

    if res_lambda == exp_lambda:
        print("  ✅ Lambda Calculation: PASS")
    else:
        print("  ❌ Lambda Calculation: FAIL")
        return

    # ----------------------------------------------------



    # ----------------------------------------------------
    print(f"\n[Step B] Verifying Squaring_plus_input...")
    lib.gm.clear()


    idx = 0
    reg_lambda = list(range(idx, idx + n)); idx += n
    reg_x      = list(range(idx, idx + n)); idx += n
    total_qubits_sq = idx


    test_lambda = random.randint(1, (1<<n)-1)
    test_x      = random.randint(1, (1<<n)-1)

    initial_state_sq = [0] * total_qubits_sq
    for b, bit in enumerate(int_to_list(test_lambda, n)):
        if bit: initial_state_sq[reg_lambda[b]] = 1
    for b, bit in enumerate(int_to_list(test_x, n)):
        if bit: initial_state_sq[reg_x[b]] = 1



    lib.Squaring_plus_input(reg_lambda + reg_x, n)


    final_state_sq = simulate_circuit(lib.gm, initial_state_sq, total_qubits_sq)


    res_l_out = list_to_int([final_state_sq[q] for q in reg_lambda])
    res_x_out = list_to_int([final_state_sq[q] for q in reg_x])



    lam_sq = gf_mod(gf_sq(test_lambda))   # lambda^2
    term = lam_sq ^ test_lambda           # lambda^2 + lambda
    exp_x_out = test_x ^ term             # x + (lambda^2 + lambda)

    print(f"  Input Lambda: {hex(test_lambda)}")
    print(f"  Input X:      {hex(test_x)}")
    print(f"  Exp Output X: {hex(exp_x_out)} (Input X + Lambda^2 + Lambda)")
    print(f"  Got Output X: {hex(res_x_out)}")


    print(f"  Got Output L: {hex(res_l_out)} (Should be same as Input Lambda)")

    if res_x_out == exp_x_out:
        print("  ✅ Squaring_plus_input: PASS (Computes x + lambda^2 + lambda)")
    else:
        print("  ❌ Squaring_plus_input: FAIL")

        exp_alt = test_x ^ lam_sq
        if res_x_out == exp_alt:
             print("     -> Note: It seems to compute x + lambda^2 (missing +lambda)")


if __name__ == "__main__":

    debug_internal_steps(parallel_num)
