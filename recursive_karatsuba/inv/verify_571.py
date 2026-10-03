import random
import sys
import os

try:
    import init.basic_gates_571 as basic_gates_571
    from init.basic_gates_571 import Square, recursive_karatsuba, Reduction, CNOT_n, GF_N
    from init.qasm import FlatGateManager
    from init.inversion_571 import Inverison_Itoh_Tsujii_based
    
    print("The [System] module was successfully imported.")
except ImportError as e:
    print("[Error] Unable to import module. Please ensure that inversion_571.py, basic_gates_571.py, and qasm.py are in the same directory.")
    print(f"{e}")
    sys.exit(1)

# ==============================================================================
# 2. Classical Algorithm Implementation 
# Irreducible polynomial: P(x)= x^571 + x^10 + x^5 + x^2 + 1
# ==============================================================================

class GF2_571_Classic:
    """Classic GF(2^163) arithmetic operations"""
    
    POLY = (1 << 571) | (1 << 10) | (1 << 5) | (1 << 2) | 1

    @staticmethod
    def multiply(a, b):
        p = 0
        for i in range(571):
            if (b >> i) & 1:
                p ^= (a << i)
        
        deg_poly = 571
        for i in range(1142, deg_poly - 1, -1):
            if (p >> i) & 1:
                p ^= (GF2_571_Classic.POLY << (i - deg_poly))
        return p

    @staticmethod
    def inverse(a):
        if a == 0: return 0
        exponent = (1 << 571) - 2
        result = 1
        base = a
        while exponent > 0:
            if exponent % 2 == 1:
                result = GF2_571_Classic.multiply(result, base)
            base = GF2_571_Classic.multiply(base, base)
            exponent //= 2
            
        return result

# ==============================================================================
# 3. Quantum Circuit Simulator 
# ==============================================================================

def run_simulation(gate_manager, total_qubits, initial_state_int, n_bits=571):
    print(f"[Simulator] Start simulation, total number of sub-bits: {total_qubits} ...")
    
    state = [(initial_state_int >> i) & 1 if i < n_bits else 0 for i in range(total_qubits)]

    SHIFT_OP = gate_manager.SHIFT_OP
    SHIFT_ARG = gate_manager.SHIFT_ARG
    MASK_ARG = (1 << SHIFT_ARG) - 1
    
    OP_X = gate_manager.OP_X
    OP_CNOT = gate_manager.OP_CNOT
    OP_TOFFOLI = gate_manager.OP_TOFFOLI

    all_chunks = gate_manager.chunks + [gate_manager.current_chunk]

    pending_c1 = None
    pending_c2 = None
    waiting_toffoli_target = False
    
    op_count = 0
    
    for chunk in all_chunks:
        for val in chunk:
            op_code = val >> SHIFT_OP
            op_count += 1
            if waiting_toffoli_target:
                t = val & MASK_ARG
                c1 = pending_c1
                c2 = pending_c2
                state[t] ^= (state[c1] & state[c2])
                waiting_toffoli_target = False
                continue

            if op_code == OP_CNOT:
                c = (val >> SHIFT_ARG) & MASK_ARG
                t = val & MASK_ARG
                state[t] ^= state[c]
                
            elif op_code == OP_X:
                t = val & MASK_ARG
                state[t] ^= 1
                
            elif op_code == OP_TOFFOLI:
                c1 = (val >> SHIFT_ARG) & MASK_ARG
                c2 = val & MASK_ARG
                pending_c1 = c1
                pending_c2 = c2
                waiting_toffoli_target = True
            
            else:
                pass

    print(f"[Simulator] Simulation complete.")
    return state

# ==============================================================================
# 4. Main Authentication Process
# ==============================================================================

def main_verification():
    N = 571             
    NUM_TESTS = 20     

    print("=======================================================")
    print(f" GF(2^{N}) Itoh-Tsujii Inverse Circuit Batch Verification")
    print(f" Planned number of tests: {NUM_TESTS}")
    print("=======================================================")

    passed_count = 0
    
    for test_idx in range(1, NUM_TESTS + 1):
        print(f"\n[Testing Progress] Case {test_idx}/{NUM_TESTS}")
        
        input_val = random.getrandbits(N)
        while input_val == 0:
            input_val = random.getrandbits(N)
            
        expected_result = GF2_571_Classic.inverse(input_val) 


        basic_gates_571.gm.clear()  
        basic_gates_571.cnt = 0     

        qubit_idx_counter = 0
        def allocate_qubits(n):
            nonlocal qubit_idx_counter
            indices = list(range(qubit_idx_counter, qubit_idx_counter + n))
            qubit_idx_counter += n
            return indices
        
        # Assign variables
        a_qubits = allocate_qubits(N)
        sqr1 = allocate_qubits(N)
        sqr2 = allocate_qubits(N)
        sqr3 = allocate_qubits(N)
        sqr4 = allocate_qubits(N)
        
        # Auxiliary bits 
        ANCILLA_SIZE = 1500000
        ancilla1 = allocate_qubits(ANCILLA_SIZE)
        ancilla2 = allocate_qubits(ANCILLA_SIZE)
        
        global_toffoli_pool = allocate_qubits(3000000)
        basic_gates_571.Toffoli_qubits = global_toffoli_pool 
        
        total_qubits = qubit_idx_counter

        for i in range(N):
            if (input_val >> i) & 1:
                basic_gates_571.X(a_qubits[i]) 

        try:
            result_qubits_list, _, _ = Inverison_Itoh_Tsujii_based(
                a_qubits, N, sqr1, sqr2, sqr3, sqr4, 
                0, ancilla1, ancilla2
            )
        except Exception as e:
            print(f"❌[Error] Circuit construction failed: {e}")
            break

        final_state = run_simulation(basic_gates_571.gm, total_qubits, 0) 

        simulated_val = 0
        for i, q_idx in enumerate(result_qubits_list):
            if final_state[q_idx] == 1:
                simulated_val |= (1 << i)

        if simulated_val == expected_result:
            print(f"✅ Case {test_idx} Passed.")
            passed_count += 1
        else:
            print(f"❌ Case {test_idx} Failed!")
            print(f"   Input:    {hex(input_val)}")
            print(f"   Expected: {hex(expected_result)}")
            print(f"   Got:      {hex(simulated_val)}")

    print("\n=======================================================")
    print(f"Test complete.")
    print(f"Pass rate: {passed_count}/{NUM_TESTS}")
    if passed_count == NUM_TESTS:
        print("🎉All test cases passed!")
    else:
        print("⚠️ error")
    print("=======================================================")
if __name__ == "__main__":
    main_verification()