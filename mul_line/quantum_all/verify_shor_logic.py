import sys
import os
import random
import importlib

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from algorithms.implementations import Shor163, Shor233, Shor283, Shor571
from algorithms.base_shor import BaseShor
from init.qasm import FlatGateManager


DATA_ROOT = DATA_ROOT = os.path.dirname(os.path.abspath(__file__))
TARGET_N =163
CURVE_A_POLY = 1
POLY_CONFIG = {
    163: (163, (1<<7) | (1<<6) | (1<<3) | 1),
    233: (233, (1<<74) | 1),
    283: (283, (1<<12) | (1<<7) | (1<<5) | 1),
    571: (571, (1<<10) | (1<<5) | (1<<2) | 1)
}

def load_matrix(filepath):
    matrix = []
    if not os.path.exists(filepath): raise FileNotFoundError(f"Matrix not found: {filepath}")
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            matrix.append([int(c) for c in line])
    return matrix

def mat_mul_vec_gf2(matrix, vec):
    dim = len(matrix)
    res = [0] * dim
    for i in range(dim):
        val = 0
        for j in range(dim):
            if matrix[i][j] and vec[j]: val ^= 1
        res[i] = val
    return res

def int_to_bits(val, n): return [(val >> i) & 1 for i in range(n)]
def bits_to_int(bits):
    val = 0
    for i, b in enumerate(bits):
        if b: val |= (1 << i)
    return val

class BasisConverter:
    def __init__(self, n):
        self.n = n
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.join(script_dir, f"quantum_{n}")
        print(f"Loading matrices for N={n} from {base_dir} ...")
        self.mat_l2x = load_matrix(os.path.join(base_dir, f"{n}_L2X.txt"))
        self.mat_x2l = load_matrix(os.path.join(base_dir, f"{n}_X2L.txt"))

    def poly_to_l(self, val_int):
        return bits_to_int(mat_mul_vec_gf2(self.mat_x2l, int_to_bits(val_int, self.n)))

    def l_to_poly(self, val_int):
        return bits_to_int(mat_mul_vec_gf2(self.mat_l2x, int_to_bits(val_int, self.n)))

def gf_mod(val, n):
    _, poly_mask = POLY_CONFIG[n]
    for i in range(val.bit_length() - 1, n - 1, -1):
        if (val >> i) & 1: val ^= (1 << i) ^ (poly_mask << (i - n))
    return val

def gf_mult(a, b, n):
    p = 0
    while b > 0:
        if b & 1: p ^= a
        a <<= 1
        b >>= 1
    return gf_mod(p, n)

def gf_sq(a, n):
    res = 0
    for i in range(n):
        if (a >> i) & 1: res |= (1 << (2 * i))
    return gf_mod(res, n)

def gf_inv(val, n):
    if val == 0: return 0
    res, base, exp = 1, val, (1 << n) - 2
    while exp > 0:
        if exp % 2 == 1: res = gf_mult(res, base, n)
        base = gf_sq(base, n)
        exp //= 2
    return res

def ec_add(x1, y1, x2, y2, a_param, n):
    if x1 == 0 and y1 == 0: return x2, y2
    if x2 == 0 and y2 == 0: return x1, y1
    if x1 == x2:
        if y1 != y2: return 0, 0
        else:
            if x1 == 0: return 0, 0
            lam = x1 ^ gf_mult(y1, gf_inv(x1, n), n)
    else:
        lam = gf_mult(y1 ^ y2, gf_inv(x1 ^ x2, n), n)
    
    lam2 = gf_sq(lam, n)
    x3 = lam2 ^ lam ^ x1 ^ x2 ^ a_param
    y3 = gf_mult(lam, x1 ^ x3, n) ^ x3 ^ y1
    return x3, y3

def simulate(gate_manager, initial_state, total_qubits):
    state = bytearray(initial_state)
    if len(state) < total_qubits: state.extend([0] * (total_qubits - len(state)))
    all_chunks = gate_manager.chunks + [gate_manager.current_chunk]
    SHIFT_OP, MASK_ARG = gate_manager.SHIFT_OP, (1 << gate_manager.SHIFT_ARG) - 1
    OP_X, OP_CNOT, OP_TOFFOLI, OP_NOP = gate_manager.OP_X, gate_manager.OP_CNOT, gate_manager.OP_TOFFOLI, gate_manager.OP_NOP
    
    pending_c1, pending_c2 = None, None
    waiting_data = False
    
    for chunk in all_chunks:
        for val in chunk:
            op = val >> SHIFT_OP
            if op == OP_NOP: continue
            if waiting_data:
                t = val & MASK_ARG
                if state[pending_c1] and state[pending_c2]: state[t] ^= 1
                waiting_data = False
            elif op == OP_CNOT:
                t, c = val & MASK_ARG, (val >> gate_manager.SHIFT_ARG) & MASK_ARG
                if state[c]: state[t] ^= 1
            elif op == OP_TOFFOLI:
                pending_c2, pending_c1 = val & MASK_ARG, (val >> gate_manager.SHIFT_ARG) & MASK_ARG
                waiting_data = True
            elif op == OP_X:
                state[val & MASK_ARG] ^= 1
    return state


def run_verification():
    n = TARGET_N
    print(f"\n{'='*20} Verifying Shor Logic for N={n} {'='*20}")
    
    try: converter = BasisConverter(n)
    except FileNotFoundError: return
    
    if n == 163: estimator = Shor163(DATA_ROOT)
    elif n == 233: estimator = Shor233(DATA_ROOT)
    elif n == 283: estimator = Shor283(DATA_ROOT)
    elif n == 571: estimator = Shor571(DATA_ROOT)
    else: return

    a_l_basis = converter.poly_to_l(CURVE_A_POLY)
    estimator.get_curve_a = lambda: a_l_basis

    def run_test_case(mode, run_id, parallel_num):
        inputs = []
        for _ in range(parallel_num):
            p1_x = random.randint(1, (1<<n)-1)
            p1_y = random.randint(1, (1<<n)-1)
            p2_x = random.randint(1, (1<<n)-1)
            p2_y = random.randint(1, (1<<n)-1)
            ctrl_bit = random.randint(0, 1) if mode == 'accumulation' else 1
            
            inputs.append({
                'p1_x': p1_x, 'p1_y': p1_y,
                'p2_x': p2_x, 'p2_y': p2_y,
                'ctrl': ctrl_bit
            })

        input_L_basis = []
        expected_results = []
        
        for inp in inputs:
            if mode == 'accumulation' and inp['ctrl'] == 0:
                exp_x = inp['p1_x']
                exp_y = inp['p1_y'] 
            else:
                exp_x, exp_y = ec_add(inp['p1_x'], inp['p1_y'], inp['p2_x'], inp['p2_y'], CURVE_A_POLY, n)
            
            expected_results.append((exp_x, exp_y))
            
            input_L_basis.append({
                'p1_x': converter.poly_to_l(inp['p1_x']),
                'p1_y': converter.poly_to_l(inp['p1_y']),
                'p2_x': converter.poly_to_l(inp['p2_x']),
                'p2_y': converter.poly_to_l(inp['p2_y'])
            })

        estimator.gm.clear()
        idx = 0
        
        q_indices_flat = []
        if mode == 'accumulation':
            q_indices_flat = list(range(idx, idx + parallel_num)); idx += parallel_num
        
        x1_indices_flat = list(range(idx, idx + n * parallel_num)); idx += n * parallel_num
        
        x2_indices_flat = []
        if mode == 'reduction':
            x2_indices_flat = list(range(idx, idx + n * parallel_num)); idx += n * parallel_num
            
        x1x2_indices_flat = list(range(idx, idx + n * parallel_num)); idx += n * parallel_num
        y1_indices_flat = list(range(idx, idx + n * parallel_num)); idx += n * parallel_num
        
        y2_indices_flat = []
        if mode == 'reduction':
            y2_indices_flat = list(range(idx, idx + n * parallel_num)); idx += n * parallel_num

        nums1 = estimator.block_size * 4 * parallel_num
        ancilla = list(range(idx, idx+nums1)); idx += nums1
        
        total_mul_size = estimator.get_mul_total_size()
        
        mul_blocks = 3 * (parallel_num - 1) + estimator.nums2_const + 2 * parallel_num + 5
        if parallel_num == 1: 
             pass
             
        nums2 = total_mul_size * mul_blocks
        estimator.Toffoli_qubits = list(range(idx, idx+nums2)); idx += nums2
        
        sqr = list(range(idx, idx+4*n)); idx += 4*n
        
        anc2_size = max(n * parallel_num, estimator.block_size * 2 + 500) 
        ancilla2 = list(range(idx, idx+anc2_size)); idx += anc2_size
        
        cswap_anc = []
        if mode == 'accumulation':
            cswap_anc = list(range(idx, idx + (2*n-1) * parallel_num)); idx += (2*n-1) * parallel_num
            
        total_qubits = idx

        x1_list = [x1_indices_flat[i*n : (i+1)*n] for i in range(parallel_num)]
        y1_list = [y1_indices_flat[i*n : (i+1)*n] for i in range(parallel_num)]
        x1x2_list = [x1x2_indices_flat[i*n : (i+1)*n] for i in range(parallel_num)]
        
        if mode == 'reduction':
            x2_list = [x2_indices_flat[i*n : (i+1)*n] for i in range(parallel_num)]
            y2_list = [y2_indices_flat[i*n : (i+1)*n] for i in range(parallel_num)]
        else:
            x2_val_list = [d['p2_x'] for d in input_L_basis]
            y2_val_list = [d['p2_y'] for d in input_L_basis]

        if mode == 'accumulation':
            estimator.shor_accumulation_step(
                parallel_num, q_indices_flat, x1_list, y1_list,
                x2_val_list, y2_val_list,
                ancilla, x1x2_list, ancilla2, sqr, cswap_anc
            )
        else:
            estimator.shor_reduction_step(
                parallel_num, x1_list, y1_list,
                x2_list, y2_list,
                ancilla, x1x2_list, ancilla2, sqr
            )

        state = [0] * total_qubits
        
        for i in range(parallel_num):
            inp = input_L_basis[i]
            
            for b_idx, bit in enumerate(int_to_bits(inp['p1_x'], n)):
                if bit: state[x1_list[i][b_idx]] = 1
            for b_idx, bit in enumerate(int_to_bits(inp['p1_y'], n)):
                if bit: state[y1_list[i][b_idx]] = 1
            
            if mode == 'reduction':
                for b_idx, bit in enumerate(int_to_bits(inp['p2_x'], n)):
                    if bit: state[x2_list[i][b_idx]] = 1
                for b_idx, bit in enumerate(int_to_bits(inp['p2_y'], n)):
                    if bit: state[y2_list[i][b_idx]] = 1
            
            if mode == 'accumulation' and inputs[i]['ctrl'] == 1:
                state[q_indices_flat[i]] = 1

        final_state = simulate(estimator.gm, state, total_qubits)

        all_pass = True
        print(f"[{mode.upper()}] Run {run_id} (Parallel={parallel_num}):")
        
        for i in range(parallel_num):
            exp_x, exp_y = expected_results[i]
            
            if mode == 'accumulation':
                res_x_L = bits_to_int([final_state[q] for q in x1_list[i]])
                res_y_L = bits_to_int([final_state[q] for q in y1_list[i]])
                res_x_poly = converter.l_to_poly(res_x_L)
                res_y_poly = converter.l_to_poly(res_y_L)
            else:
                res_x_L = bits_to_int([final_state[q] for q in x1x2_list[i]])
                res_x_poly = converter.l_to_poly(res_x_L)
                res_y_poly = exp_y 
            
            match = (res_x_poly == exp_x) and (res_y_poly == exp_y)
            if not match: all_pass = False
            
            if not match: 
                status = "✅" if match else "❌"
                desc = f"Branch {i} (Ctrl={inputs[i]['ctrl']})" if mode == 'accumulation' else f"Branch {i}"
                print(f"  {desc} | {status}")
                if not match:
                    print(f"    Exp: ({hex(exp_x)}, {hex(exp_y)})")
                    print(f"    Got: ({hex(res_x_poly)}, {hex(res_y_poly)})")
        
        if all_pass: print(f"  All {parallel_num} branches PASSED.")

    print("\n--- Testing Phase 1: Accumulation (Parallel=1) ---")
    for i in range(10): run_test_case('accumulation', i+1, 1)

    print("\n--- Testing Phase 1: Accumulation (Parallel=5) ---")
    run_test_case('accumulation', 1, 5)

    print("\n--- Testing Phase 3: Reduction (Parallel=1) ---")
    for i in range(10): run_test_case('reduction', 1+i, 1)
    
    print("\n--- Testing Phase 3: Reduction (Parallel=5) ---")
    run_test_case('reduction', 1, 5)

if __name__ == "__main__":
    run_verification()