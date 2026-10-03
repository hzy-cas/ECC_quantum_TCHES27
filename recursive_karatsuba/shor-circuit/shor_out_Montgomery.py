
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
import importlib
import math
import numpy as np
import gc
import random 
from init.stats_utils import get_exact_resources_optimized
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from init.basic_gates_163 import *
    from init.inversion_163 import *

TARGET_N =163
print(f"=== Initializing Backend for GF(2^{TARGET_N}) ===")
try:
    gates_module_name = f"init.basic_gates_{TARGET_N}"
    inversion_module_name = f"init.inversion_{TARGET_N}"
    gates_mod = importlib.import_module(gates_module_name)
    inv_mod = importlib.import_module(inversion_module_name)
    this_module = sys.modules[__name__]
    def inject_from_module(module):
        for attr in dir(module):
            if not attr.startswith('__'):
                setattr(this_module, attr, getattr(module, attr))
    inject_from_module(gates_mod)
    inject_from_module(inv_mod)
    globals().update({attr: getattr(gates_mod, attr) for attr in dir(gates_mod) if not attr.startswith('__')})
    globals().update({attr: getattr(inv_mod, attr) for attr in dir(inv_mod) if not attr.startswith('__')})
    GF_N = getattr(gates_mod, 'GF_N')
except ImportError as e:
    print(f"\n[Error] Failed to load backend for GF({TARGET_N}).")
    import traceback
    traceback.print_exc()
    sys.exit(1)

sys.setrecursionlimit(5000)


def CONST_ADD_n(target_reg, constant_val):
    """
    Add (XOR) the classical constant constant_val to the quantum register target_reg.
    """
    n = len(target_reg)
    for i in range(n):
        if (constant_val >> i) & 1:
            X(target_reg[i])

def copy_parallel(value, ancillas, n):
    """
    Binary tree copying (Fan-out)
    """
    divide = int(math.log2(n))
    last = n - 2 ** divide

    copy_list = []
    copy_list.append(value)
    for i in range(n - 1):
        copy_list.append(ancillas[i])

    for i in range(divide):
        for j in range(2 ** i):
            if (i == 0):
                CNOT(copy_list[0], copy_list[1])
            else:
                CNOT(copy_list[j], copy_list[2 ** i + j])

    for i in range(last):
        CNOT(copy_list[i], copy_list[2 ** divide + i])

    return copy_list


def cons_mul(xlist, count, ancilla):
    tree = []  
    tree.append(xlist[:])  
    num = KARATSUBA_ANC_BASE 
    current_level = xlist[:]
    n = GF_N 

    while len(tree[-1]) > 1:
        next_level = [] 
        i = 0 
        j = 0
        while i < len(current_level):
            if i + 1 < len(current_level):
                count = 0
                product, count, ancilla[num*j:num*(j+1)] = recursive_karatsuba(current_level[i], current_level[i+1], n, count, ancilla[num*j:num*(j+1)])
                product = Reduction(product)
                next_level.append(product)
                i += 2 
                j += 1
            else:
                next_level.append(current_level[i])
                i += 1 
        tree.append(next_level[:])
        current_level = next_level
    return tree[len(tree)-1][0], count, ancilla, tree

def Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, target_copy, inv_res):
    #Montgomery Trick
    n = GF_N
    num = KARATSUBA_ANC_BASE 
    parallel_num = len(x1x2_list)
    
    if parallel_num == 0:
        return [], ancilla

    start_ptr_inv = gm.current_pointer()
    
    if parallel_num == 1:
        count = 0
        x1x2_inv, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:] = Inverison_Itoh_Tsujii_based(
             x1x2_list[0], n, sqr[:n], sqr[n:2*n], sqr[2*n:3*n], sqr[3*n:4*n], count, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:]
        )
        mid_ptr_inv = gm.current_pointer()
        CNOT_n(x1x2_inv, inv_res[0:n])
    else:
        count = 0
        product_final, count, ancilla, tree = cons_mul(x1x2_list, count, ancilla)
        count = 0
        product_inv, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:] = Inverison_Itoh_Tsujii_based(
            product_final, n, sqr[:n], sqr[n:2*n], sqr[2*n:3*n], sqr[3*n:4*n], count, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:]
        )
        
        inv_tree = []
        inv_tree.append([product_inv])  
        
        for level in range(len(tree) - 2, -1, -1):
            current_level = tree[level]        
            current_inv_level = inv_tree[-1]   
            next_inv_level = []                
            inv_idx = 0  
            i = 0        
            idx=0
            while i < len(current_level):
                if i + 1 < len(current_level):
                    CNOT_n(current_inv_level[inv_idx], ancilla2[idx:idx+n])
                    copy = ancilla2[idx:idx+n]
                    count = 0
                    a_inv, count, ancilla[num*i:num*(i+1)] = recursive_karatsuba(
                        current_inv_level[inv_idx], current_level[i+1], n, count, ancilla[num*i:num*(i+1)]
                    )
                    a_inv = Reduction(a_inv)
                    count = 0
                    b_inv, count, ancilla[num*(i+1):num*(i+2)] = recursive_karatsuba(
                        copy, current_level[i], n, count, ancilla[num*(i+1):num*(i+2)]
                    )
                    b_inv = Reduction(b_inv)
                    CNOT_n(current_inv_level[inv_idx], copy)
                    next_inv_level.append(a_inv)
                    next_inv_level.append(b_inv)
                    idx += n
                    i += 2
                    inv_idx += 1
                else:
                    next_inv_level.append(current_inv_level[inv_idx])
                    i += 1
                    inv_idx += 1
            inv_tree.append(next_inv_level)
        
        inverses = inv_tree[-1]
        mid_ptr_inv = gm.current_pointer()
        for i in range(parallel_num):
            CNOT_n(inverses[i], inv_res[i*n : (i+1)*n])

    gm.replay_reverse(start_ptr_inv, mid_ptr_inv)
    gates_mod.cnt = 0
    
    start_ptr_mul = gm.current_pointer()
    lambda_list_dirty = []
    for i in range(parallel_num):
        count = 0
        current_inv = inv_res[i*n : (i+1)*n]
        lambda_i, count, ancilla[num*i:num*(i+1)] = recursive_karatsuba(y1_list[i], current_inv, n, count, ancilla[num*i:num*(i+1)])
        lambda_i = Reduction(lambda_i)
        lambda_list_dirty.append(lambda_i)

    mid_ptr_mul = gm.current_pointer()
    for i in range(parallel_num):
        target_slice = target_copy[i*n : (i+1)*n]
        CNOT_n(lambda_list_dirty[i], target_slice)

    gm.replay_reverse(start_ptr_mul, mid_ptr_mul)
    gates_mod.cnt = 0
    return target_copy, ancilla




def shor_level(parallel_num, q_list, x1_list, y1_list, x2_list_val, y2_list_val, ancilla, x1x2_list, ancilla2, sqr, level, lambda_res, mult_res, inv_res, cswap_ancilla):
    """
    Stream Accumulation
    """

    a = 3 
    num = KARATSUBA_ANC_BASE 
    n = GF_N 
    
    for i in range(parallel_num):
        CNOT_n(x1_list[i], x1x2_list[i])
        
        CONST_ADD_n(y1_list[i], y2_list_val[i])
        CONST_ADD_n(x1x2_list[i], x2_list_val[i])
    
    count = 0
    gates_mod.cnt = 0 
    lamba_flat, _ = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, lambda_res, inv_res)
    
    lamba_list = []
    for k in range(parallel_num):
        lamba_list.append(lamba_flat[k*n : (k+1)*n])
    gates_mod.cnt = 0 
    for i in range(parallel_num):
        
        CONST_ADD_n(y1_list[i], y2_list_val[i])
        CONST_ADD_n(x1x2_list[i], x2_list_val[i] ^ a)

        lamba_list[i], x1x2_list[i] = Squaring_plus_input(lamba_list[i] + x1x2_list[i], n) 

        start_mul = gm.current_pointer()
        count = 0
        result2, count, ancilla[num*i:num*(i+1)] = recursive_karatsuba(x1x2_list[i], lamba_list[i], n, count, ancilla[num*i:num*(i+1)])  
        result2 = Reduction(result2)
        end_mul = gm.current_pointer()

        target_y = mult_res[i*n : (i+1)*n] 
        CNOT_n(result2, target_y)
        gm.replay_reverse(start_mul, end_mul)

        CONST_ADD_n(x1x2_list[i], x2_list_val[i])
        CONST_ADD_n(target_y, y2_list_val[i])
        CNOT_n(x1x2_list[i], target_y)
        
        current_cswap_anc = cswap_ancilla[i * (2*n - 1) : (i + 1) * (2*n - 1)]
        copy_start_ptr = gm.current_pointer()
        
        control_bits = copy_parallel(q_list[i], current_cswap_anc, 2 * n)
        
        copy_end_ptr = gm.current_pointer()
        
        for bit_idx in range(n):
            CSWAP(control_bits[bit_idx], x1_list[i][bit_idx], x1x2_list[i][bit_idx])
            
        for bit_idx in range(n):
            CSWAP(control_bits[n + bit_idx], y1_list[i][bit_idx], target_y[bit_idx])
            
        gm.replay_reverse(copy_start_ptr, copy_end_ptr)

def shor(parallel_num, level):
    global q, nums1, n
    n = GF_N
    gm.clear()
    gates_mod.cnt = 0 
    
    nums1 = KARATSUBA_ANC_BASE * parallel_num 
    nums2 = TOFFOLI_BASE * (3 * (parallel_num - 1) + TOFFOLI_OFFSET) 
    
    q_list = [i for i in range(parallel_num)]
    offset = parallel_num

    x1 = [offset + i for i in range(n*parallel_num)]
    offset += n*parallel_num
    
    x1x2 = [offset + i for i in range(n*parallel_num)] 
    offset += n*parallel_num
    
    y1 = [offset + i for i in range(n*parallel_num)]
    offset += n*parallel_num
    
    x1_list, y1_list, x1x2_list = [], [], []
    for i in range(parallel_num):
        x1_list.append(x1[i*n:(i+1)*n])
        y1_list.append(y1[i*n:(i+1)*n])
        x1x2_list.append(x1x2[i*n:(i+1)*n])

    full_ones = (1 << n) - 1
    x2_list_val = [full_ones for _ in range(parallel_num)]
    y2_list_val = [full_ones for _ in range(parallel_num)]

    ancilla = [offset + i for i in range(int(nums1))]
    length = offset + len(ancilla)
    
    gates_mod.Toffoli_qubits = [length+i for i in range(nums2)]
    length += len(gates_mod.Toffoli_qubits)

    sqr = [length + i for i in range(4*n)]
    length += 4*n
    
    ancilla2 = [length + i for i in range(n*int(parallel_num/2))] 
    length += int(parallel_num/2)*n
    
    lambda_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    mult_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    inv_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    cswap_ancilla_size = (2 * n - 1) * parallel_num
    cswap_ancilla = [length + i for i in range(cswap_ancilla_size)]
    length += cswap_ancilla_size

    
    shor_level(parallel_num, q_list, x1_list, y1_list, x2_list_val, y2_list_val, ancilla, x1x2_list, ancilla2, sqr, level, lambda_res, mult_res, inv_res, cswap_ancilla)
    
    stats, full_depth, current_depth,toffoli_depth = get_exact_resources_optimized(gm, length)
    
    toffoli = stats['Toffoli_count']
    cnot = stats['CNOT_count']

    return length, toffoli, cnot, full_depth, current_depth,toffoli_depth

def shor_reduction_level(parallel_num, x1_list, y1_list, x2_list, y2_list, ancilla, x1x2_list, ancilla2, sqr, lambda_res, mult_res, inv_res):
    """
    Binary Tree Reduction
    """
    a = 3  
    num = KARATSUBA_ANC_BASE 
    n = GF_N 
    for i in range(parallel_num):
        
        CNOT_n(x1_list[i], x1x2_list[i])
        CNOT_n(x2_list[i], x1x2_list[i])
        CNOT_n(y2_list[i], y1_list[i])
    
    count = 0
    gates_mod.cnt = 0 
    lamba_flat, _ = Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr, lambda_res, inv_res)
    
    lamba_list = []
    for k in range(parallel_num):
        lamba_list.append(lamba_flat[k*n : (k+1)*n])
    gates_mod.cnt = 0 

    for i in range(parallel_num):
        CONST_ADD_n(x1x2_list[i], a)
        lamba_list[i], x1x2_list[i] = Squaring_plus_input(lamba_list[i] + x1x2_list[i], n) 
        CNOT_n(x2_list[i], x1x2_list[i])
        start_mul = gm.current_pointer()
        count = 0
        result2, count, ancilla[num*i:num*(i+1)] = recursive_karatsuba(x1x2_list[i], lamba_list[i], n, count, ancilla[num*i:num*(i+1)])  
        result2 = Reduction(result2)
        end_mul = gm.current_pointer()

        target_y = mult_res[i*n : (i+1)*n] 
        CNOT_n(result2, target_y) 
        gm.replay_reverse(start_mul, end_mul)

        CNOT_n(y2_list[i], target_y)
        CNOT_n(x2_list[i], x1x2_list[i]) 
        CNOT_n(x1x2_list[i], target_y)

def shor_reduction(parallel_num):
    n = GF_N
    gm.clear()
    gates_mod.cnt = 0 
    
    nums1 = KARATSUBA_ANC_BASE * parallel_num 
    nums2 = TOFFOLI_BASE * (3 * (parallel_num - 1) + TOFFOLI_OFFSET) 
    
    current_idx = 0
    
    x1 = [current_idx + i for i in range(n*parallel_num)]
    current_idx += n*parallel_num
    
    x2 = [current_idx + i for i in range(n*parallel_num)]
    current_idx += n*parallel_num
    
    x1x2 = [current_idx + i for i in range(n*parallel_num)]
    current_idx += n*parallel_num

    y1 = [current_idx + i for i in range(n*parallel_num)]
    current_idx += n*parallel_num
    
    y2 = [current_idx + i for i in range(n*parallel_num)]
    current_idx += n*parallel_num

    x1_list, y1_list, x1x2_list, x2_list, y2_list = [], [], [], [], []
    for i in range(parallel_num):
        x1_list.append(x1[i*n:(i+1)*n])
        y1_list.append(y1[i*n:(i+1)*n])
        x1x2_list.append(x1x2[i*n:(i+1)*n])
        x2_list.append(x2[i*n:(i+1)*n])
        y2_list.append(y2[i*n:(i+1)*n])

    ancilla = [current_idx + i for i in range(int(nums1))]
    length = current_idx + len(ancilla)
    
    gates_mod.Toffoli_qubits = [length+i for i in range(nums2)]
    length += len(gates_mod.Toffoli_qubits)

    sqr = [length + i for i in range(4*n)]
    length += 4*n
    
    ancilla2 = [length + i for i in range(n*int(parallel_num/2))] 
    length += int(parallel_num/2)*n
    
    lambda_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    mult_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    inv_res = [length + i for i in range(n * parallel_num)]
    length += n * parallel_num
    
    shor_reduction_level(parallel_num, x1_list, y1_list, x2_list, y2_list, ancilla, x1x2_list, ancilla2, sqr, lambda_res, mult_res, inv_res)
    
    stats, full_depth, current_depth,toffoli_depth =  get_exact_resources_optimized(gm, length)
    
    toffoli = stats['Toffoli_count']
    cnot = stats['CNOT_count']

    return length, toffoli, cnot, full_depth, current_depth,toffoli_depth

def accumulate_and_reduce(total_inputs, max_parallel_ops):
    """
    Phase 1: 
      - Depth += 2*ceil(log2(2n)) + 1 + d
      - CNOT += 2*ops*(2n-1) + 2n*ops + c
      - Qubit = length - ops 
      - Clear = length - 6n * ops
    Phase 2:
      - Demand = length - 2n*ops' - ops'
      - If Demand > Clear: Peak += (Demand - Clear), Clear = length - 6n*ops'
      - Else: Clear = length - 6n*ops' + Over
    Phase 3:
      - Demand = length - 4n*ops''
      - If Demand > Clear: Peak += (Demand - Clear), Clear = length - 8n*ops''
      - Else: Clear = length - 8n*ops'' + Over
    Final:
      - T_final = 2*T
      - C_final = 2*C + 2n
      - D_final = 2*D + 2
      - Q_final = Peak + 4n + 2
    """
    total_Toffoli = 0
    total_CNOT = 0
    total_fDepth = 0
    total_cDepth = 0
    total_tDepth = 0
    peak_qubits = 0
    clear = 0  
    
    n = GF_N 
    current_pool = total_inputs
    
    print(f"{'Accumulator Strategy (Strict Rules)':^60}")
    print(f"{'Input=' + str(total_inputs) + ', Batch Limit=' + str(max_parallel_ops):^60}")
    ops = max_parallel_ops
    print(f"\n[Phase 1] Initialization with ops={ops}")
    
    length, t, c, fd,cd,td = shor(ops, 0)

    fanout_depth = 2 * math.ceil(math.log2(2 * n))
    
    phase1_cnot = 2 * ops * (2 * n - 1) + 2 * n * ops + c
    
    total_Toffoli += t
    total_CNOT += phase1_cnot
    total_fDepth += fanout_depth + 1 + fd
    total_cDepth += fanout_depth + 1 + cd
    total_tDepth += td
    peak_qubits = length - ops
    clear = length - 6 * n * ops - ops

    print(f"  Phase 1 Stats: Depth={fanout_depth + 1 + fd}, CNOT={phase1_cnot}, Qubits={peak_qubits}, Clear={clear}")
    
    current_pool -= (2 * ops)

    print(f"\n[Phase 2] Accumulation Loop")
    batch_count = 0
    
    while current_pool > 0:
        ops_prime = min(max_parallel_ops, current_pool)
        batch_count += 1

        length, t, c, fd,cd,td = shor(ops_prime, 0)

        total_Toffoli += t
        total_CNOT += c
        total_fDepth += fd
        total_cDepth += cd
        total_tDepth += td
        demand = length - 2 * n * ops_prime - ops_prime
        
        needed_new = 0
        if demand > clear:
            needed_new = demand - clear
            peak_qubits += needed_new
            
            clear = length - 6 * n * ops_prime - ops_prime
        else:
            
            over = clear - demand
            
            clear = length - 6 * n * ops_prime - ops_prime + over
            
        current_pool -= ops_prime
        gc.collect()
        print(f"  Batch {batch_count} (ops={ops_prime}): New Qubits={needed_new}, Total Peak={peak_qubits}, Clear={clear}")
    print(f"\n[Phase 3] Tree Reduction")
    current_items = max_parallel_ops
    tree_level = 0
    
    while current_items > 1:
        ops_double_prime = current_items // 2
        leftover = current_items % 2

        length, t, c, fd,cd,td = shor_reduction(ops_double_prime)
        
        total_Toffoli += t
        total_CNOT += c
        total_fDepth += fd
        total_cDepth += cd
        total_tDepth += td

        demand = length - 4 * n * ops_double_prime
        
        needed_new = 0
        if demand > clear:
            needed_new = demand - clear
            peak_qubits += needed_new
            
            clear = length - 8 * n * ops_double_prime
        else:
            over = clear - demand
            
            clear = length - 8 * n * ops_double_prime + over
        
        current_items = ops_double_prime + leftover
        tree_level += 1
        gc.collect()
        print(f"  Level {tree_level} (ops={ops_double_prime}): New Qubits={needed_new}, Total Peak={peak_qubits}, Clear={clear}")
    final_Toffoli = 2 * total_Toffoli
    final_CNOT = 2 * total_CNOT + 2 * n
    final_fDepth = 2 * total_fDepth + 2
    final_cDepth = 2 * total_cDepth + 2
    final_tDepth = 2 * total_tDepth
    final_Qubits = peak_qubits + 4 * n + 2

    return final_Toffoli, final_CNOT, final_fDepth, final_cDepth,final_tDepth,final_Qubits

if __name__ == "__main__":
    n = GF_N 
    init_dir = os.path.join(os.path.dirname(__file__), 'init')
    if not os.path.exists(os.path.join(init_dir, '__init__.py')):
        os.makedirs(init_dir, exist_ok=True)
        with open(os.path.join(init_dir, '__init__.py'), 'w') as f: pass
    try:
        Toffoli_nums, CNOT_nums, full_depth, current_depth,toffoli_depth, num_qubit = accumulate_and_reduce(328,30)
        print("\n" + "="*40)
        print(f"{'Quantum Resource Estimates':^40}")
        print("="*40)
        print(f"{'Metric':<20} | {'Value':>15}")
        print("-" * 40)
        print(f"{'Toffoli Gates':<20} | {Toffoli_nums:>15,}")
        print(f"{'CNOT Gates':<20} | {CNOT_nums:>15,}")
        print(f"{'Full Depth':<20} | {full_depth:>15}")
        print(f"{'Current Depth':<20} | {current_depth:>15}")
        print(f"{'Toffoli Depth':<20} | {toffoli_depth:>15}")
        print(f"{'Qubits':<20} | {num_qubit:>15,}")
        print("-" * 40)
    except Exception as e:
        import traceback
        traceback.print_exc()