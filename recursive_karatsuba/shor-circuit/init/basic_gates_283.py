import os
import math
from .qasm import FlatGateManager, MatrixCache

gm = FlatGateManager()
matrix_loader = MatrixCache()
GF_N = 283                  
KARATSUBA_ANC_BASE = 19980 
TOFFOLI_BASE = 10273        
TOFFOLI_OFFSET = 11          
Toffoli_qubits = []

cnt = 0
def get_matrix_data_safe(file_path, n):
    if os.path.exists(file_path):
        return matrix_loader.get_matrix_data(file_path, n)
    else:
        print(f"CRITICAL WARNING: Matrix file {os.path.basename(file_path)} not found at {file_path}. Using dummy ops.")
        dummy_ops = []
        for i in range(n):
            dummy_ops.append((i, i, None))
        return dummy_ops, {}
    
def X(a):
    gm.add_X(a)
    return [a]

def Toffoli_gate(a, b, c):
    gm.add_Toffoli(a, b, c)
    return [a, b, c]

def Round_constant_XOR(k, rc, bit):
    for i in range(bit):
        if (rc >> i & 1):
            X(k[i])

def CNOT(a, b):
    gm.add_CNOT(a, b)
    return [a, b]

def CNOT_n(a, b):
    for i in range(GF_N):
        CNOT(a[i], b[i])

def CTRL_CNOT_n(copylist2, a, b):
    for i in range(GF_N):
        Toffoli_gate(copylist2[i], a[i], b[i])

def CSWAP(a, b, c):
    CNOT(c, b)
    Toffoli_gate(a, b, c)
    CNOT(c, b)

def SWAP(a, b):
    CNOT(a, b)
    CNOT(b, a)
    CNOT(a, b)

def copy_parallel(value, ancillas, n):
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

def Squaring_plus_input(x, n):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    div_dir = os.path.dirname(current_dir)
    project_root = os.path.dirname(div_dir)
    file_path = os.path.join(project_root, "square", "square_283", "square_Matrix_Squaring_plus.txt")
    ops, y_exists_map = get_matrix_data_safe(file_path, n)
    y = {} 
    for target, control, y_index in ops:
        if control < len(x) and target < len(x):
            CNOT(x[control], x[target])
            if y_index is not None:
                y[y_index] = x[target]
    if not y:
        return x[:n], x[:n]
    max_y_index = max(y.keys())
    y_list = [0] * (max_y_index + 1)
    for i, val in y.items():
        y_list[i] = val
    for i in range(len(y_list)):
        if y_list[i] == 0 and i not in y:
            if i < len(x):
                y_list[i] = x[i]
    if len(y_list) < 2*n:
        y_list.extend([0] * (2*n - len(y_list)))
    return y_list[:n], y_list[-n:]

def Square(x, n, power):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    div_dir = os.path.dirname(current_dir)
    project_root = os.path.dirname(div_dir)
    file_path = os.path.join(project_root, "square", "square_283", f"square_Matrix_2_{power}.txt")
    ops, y_exists_map = get_matrix_data_safe(file_path, n)
    y = {}
    for target, control, y_index in ops:
        if control < len(x) and target < len(x):
            CNOT(x[control], x[target])
            if y_index is not None:
                y[y_index] = x[target]
    if not y: 
        return x[:n], x[:n]
    max_y_index = max(y.keys())
    y_list = [0] * (max_y_index + 1)
    for i, val in y.items(): y_list[i] = val
    for i in range(len(y_list)):
        if y_list[i] == 0 and i not in y:
            if i < len(x):
                y_list[i] = x[i]
    if len(y_list) < 2*n:
        y_list.extend([0] * (2*n - len(y_list)))
    return y_list[:n], y_list[-n:]

def Reduction(result): 
    n = 283
    for i in range(n - 1):
        if (i < n): CNOT(result[i + n], result[i])
    for i in range(n - 1):
        if (i + 5 < n): CNOT(result[i + n], result[i + 5])
    for i in range(n - 1):
        if (i + 7 < n): CNOT(result[i + n], result[i + 7])
    for i in range(n - 1):
        if (i + 12 < n): CNOT(result[i + n], result[i + 12])
    Modular_small(result[561:], result, 4)
    Modular_small(result[559:], result, 6)
    Modular_small(result[554:], result, 11)
    return result[0:n]

def Modular_small(input, result, size):
    for i in range(size):
        CNOT(input[i], result[0 + i])   
    for i in range(size):
        CNOT(input[i], result[5 + i])   
    for i in range(size):
        CNOT(input[i], result[7 + i])   
    for i in range(size):
        CNOT(input[i], result[12 + i])  

def combine(a, b, r, n):
    if (n % 2 != 0):
        for i in range(n): CNOT(a[i], r[i])
        for i in range(n - 2): CNOT(b[i], r[i])
        for i in range(n // 2): CNOT(a[n // 2 + 1 + i], r[i])
        for i in range(n // 2): CNOT(b[i], r[n // 2 + 1 + i])
        out = []
        for i in range(n // 2 + 1): out.append(a[i])
        for i in range(n): out.append(r[i])
        for i in range((2 * n - 1) - n // 2 - 1 - n): out.append(b[n // 2 + i])
        return out
    half_n = int(n/2)
    for i in range(n-1):
        CNOT(a[i], r[i])
        CNOT(b[i], r[i])
    for i in range(half_n-1):
        CNOT(a[half_n+i], r[i])
        CNOT(b[i], r[half_n+i])
    result = []
    for i in range(half_n): result.append(a[i])
    for i in range(n-1): result.append(r[i])
    for i in range(half_n): result.append(b[half_n-1+i])
    return result

def recursive_karatsuba(a, b, n, count, ancilla):
    global cnt
    if(n == 1):
        c = Toffoli_qubits[cnt]
        Toffoli_gate(a[0], b[0], c)
        cnt += 1
        return [c], count, ancilla
    r_low = n//2
    if(n % 2 != 0): r_low = r_low +1
    r_a = ancilla[count:count + r_low]
    count = count + r_low
    r_b = ancilla[count:count + r_low]
    count = count + r_low
    start_ptr = gm.current_pointer()
    for i in range(r_low): CNOT(a[i], r_a[i])
    for i in range(n//2): CNOT(a[r_low + i], r_a[i])
    for i in range(r_low): CNOT(b[i], r_b[i])
    for i in range(n//2): CNOT(b[r_low + i], r_b[i])
    end_ptr = gm.current_pointer()
    if(r_low == 1):
        c = Toffoli_qubits[cnt:cnt+3]
        cnt += 3
        Toffoli_gate(a[0], b[0], c[0])
        Toffoli_gate(a[1], b[1], c[2])
        Toffoli_gate(r_a[0], r_b[0], c[1])
        CNOT(c[0], c[1])
        CNOT(c[2], c[1])
        gm.replay_reverse(start_ptr, end_ptr)
        return c, count, ancilla
    c_a, count, ancilla = recursive_karatsuba(a[0:r_low], b[0:r_low], r_low, count, ancilla)
    c_b, count, ancilla = recursive_karatsuba(a[r_low:n], b[r_low:n], n//2, count, ancilla)
    c_r, count, ancilla = recursive_karatsuba(r_a[0:r_low], r_b[0:r_low], r_low, count, ancilla)
    gm.replay_reverse(start_ptr, end_ptr)
    result = combine(c_a, c_b, c_r, n)
    return result, count, ancilla