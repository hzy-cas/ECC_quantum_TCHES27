import os

def CNOT_n(gm, a, b, n):
    for i in range(n):
        gm.add_CNOT(a[i], b[i])

def Toffoli_gate(gm, a, b, c):
    gm.add_Toffoli(a, b, c)

def mul_matrix(gm, matrix_loader, x, n, file_path):
    ops, matrix_rows = matrix_loader.get_matrix_data(file_path, n)
    y = {}
    for target, control, y_index in ops:
        gm.add_CNOT(x[control], x[target])
        if y_index is not None: y[y_index] = x[target]
    
    if not y: return []
    y_list = [0] * n
    for i, val in y.items(): y_list[i] = val
    for i in range(len(y_list)):
        if y_list[i] == 0 and i not in y:
            if i < len(matrix_rows) and matrix_rows[i]:
                j = matrix_rows[i][0]
                if j < len(x): y_list[i] = x[j]
    return y_list

def mul(gm, matrix_loader, a0, b0, c0, n, path_config):
    begin = gm.current_pointer()
    
    dim1 = path_config['mul_dim_1']
    dim2 = path_config['mul_dim_2']
    
    a1 = mul_matrix(gm, matrix_loader, a0[:dim1], dim1, path_config['TD'])
    b1 = mul_matrix(gm, matrix_loader, b0[:dim1], dim1, path_config['TD'])
    
    a2 = mul_matrix(gm, matrix_loader, a1 + a0[dim1:dim2], dim2, path_config['A'])
    b2 = mul_matrix(gm, matrix_loader, b1 + b0[dim1:dim2], dim2, path_config['A'])
    
    end = gm.current_pointer()
    
    for i in range(dim2):
        Toffoli_gate(gm, a2[i], b2[i], c0[i])
        
    c1 = mul_matrix(gm, matrix_loader, c0, dim2, path_config['C'])
    c2 = mul_matrix(gm, matrix_loader, c1[:dim1], dim1, path_config['Inv'])
    
    gm.replay_reverse(begin, end)
    return c2[:n]

def Square(gm, matrix_loader, x, n, power, base_path):
    file_path = os.path.join(base_path, f"square_{n}", f"square_Matrix_2_{power}.txt")
    ops, matrix_rows = matrix_loader.get_matrix_data(file_path, 2*n) 
    y = {}
    for target, control, y_index in ops:
        gm.add_CNOT(x[control], x[target])
        if y_index is not None: y[y_index] = x[target]
    
    if not y: return [], []
    y_list = [0] * (2*n)
    for i, val in y.items(): y_list[i] = val
    for i in range(len(y_list)):
        if y_list[i] == 0 and i not in y:
            if i < len(matrix_rows) and matrix_rows[i]:
                j = matrix_rows[i][0]
                if j < len(x): y_list[i] = x[j]
    return y_list[:n], y_list[n:]


def Inversion_163(gm, ml, a, n, sqr1, sqr2, sqr3, sqr4, count, ancilla2, ancilla3, toffoli_qubits, path_cfg, base_path):
    def S(x, p): return Square(gm, ml, x, n, p, base_path)
    def M(x, y, c_off): 
        res = mul(gm, ml, x, y, toffoli_qubits[c_off:c_off+906], n, path_cfg)
        return res, c_off + 906
    def C(x, y): CNOT_n(gm, x, y, n)

    # 1 -> 2
    a, A_1_1 = S(a+sqr1, 1)
    A_2_0, count = M(a+ancilla2[:743], A_1_1+ancilla2[743:], count)
    
    # Chain 1 (2->4)
    A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
    a, A_1_1 = S(a+A_1_1, 1) # clear sqr1
    sqr1 = A_1_1
    A_4_0, count = M(A_2_2+ancilla2[:743], A_2_0+ancilla2[743:], count)
    
    # 4 -> 8
    A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
    sqr2 = A_2_2
    A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
    A_8_0, count = M(A_4_0+ancilla2[:743], A_4_4+ancilla2[743:], count)

    # 8 -> 16
    A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
    sqr1 = A_4_4
    A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
    A_16_0, count = M(A_8_0+ancilla2[:743], A_8_8+ancilla2[743:], count)

    # 16 -> 32
    A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
    sqr2 = A_8_8
    A_16_0, A_16_16 = S(A_16_0+sqr1, 16)
    A_32_0, count = M(A_16_0+ancilla2[:743], A_16_16+ancilla2[743:], count)
    
    # S - Chain 2 (Prep 2->34 side branch)
    C(A_2_0, sqr3)
    sqr3, A_2_32 = S(sqr3+sqr4, 32)
    C(A_2_0, sqr3)

    # Chain 1 (32 -> 64)
    A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
    sqr1 = A_16_16
    A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
    
    # Merge (32 + 2 = 34)
    C(A_32_0, sqr3)
    # Main
    A_64_0, count = M(A_32_0+ancilla2[:743], A_32_32+ancilla2[743:], count)
    # Side (A_2_32 * sqr3) -> A_34_0
    A_34_0, count = M(A_2_32+ancilla3[:743], sqr3+ancilla3[743:], count)
    C(A_32_0, sqr3)
    
    # 64 -> 128
    A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
    sqr2 = A_32_32
    A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
    A_128_0, count = M(A_64_0+ancilla2[:743], A_64_64+ancilla2[743:], count)

    # Chain 1 Cleanup
    A_64_0, A_64_64 = S(A_64_0+A_64_64, 64)
    sqr1 = A_64_64 # clear sqr1
    
    # Final Merge Prep (128->1, 34->129)
    A_128_0, A_128_1 = S(A_128_0+sqr2, 1)
    A_34_0, A_34_129 = S(A_34_0+sqr3, 129)
    
    # Final Mul
    A_162_1, count = M(A_34_129+ancilla3[:743], A_128_1+ancilla3[743:], count)
    
    # Cleanup A_2_0 (optional based on logic, but good for completeness)
    A_2_0, sqr4 = S(A_2_0+A_2_32, 32)
    
    return A_162_1

def Inversion_233(gm, ml, a, n, sqr1, sqr2, sqr3, sqr4, count, ancilla2, ancilla3, toffoli_qubits, path_cfg, base_path):
    def S(x, p): return Square(gm, ml, x, n, p, base_path)
    def M(x, y, c_off): 
        res = mul(gm, ml, x, y, toffoli_qubits[c_off:c_off+1341], n, path_cfg)
        return res, c_off + 1341
    def C(x, y): CNOT_n(gm, x, y, n)
    
    # 1 -> 2
    a, A_1_1 = S(a+sqr1, 1)
    A_2_0, count = M(a+ancilla2[:1108], A_1_1+ancilla2[1108:], count)
    
    # 2 -> 4
    A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
    a, A_1_1 = S(a+A_1_1, 1) 
    sqr1 = A_1_1
    A_4_0, count = M(A_2_2+ancilla2[:1108], A_2_0+ancilla2[1108:], count)

    # 4 -> 8
    A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
    sqr2 = A_2_2
    A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
    A_8_0, count = M(A_4_0+ancilla2[:1108], A_4_4+ancilla2[1108:], count)

    # Chain 1 (8 -> 16)
    A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
    sqr1 = A_4_4
    A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
    A_16_0, count = M(A_8_0+ancilla2[:1108], A_8_8+ancilla2[1108:], count)

    # 16 -> 32
    A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
    sqr2 = A_8_8
    A_16_0, A_16_16 = S(A_16_0+sqr1, 16)
    A_32_0, count = M(A_16_0+ancilla2[:1108], A_16_16+ancilla2[1108:], count)

    # S - Chain 2 (Prep 8 -> 40 side branch)
    C(A_8_0, sqr3)
    sqr3, A_8_32 = S(sqr3+sqr4, 32)
    C(A_8_0, sqr3)

    # Chain 1 (32 -> 64)
    A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
    sqr1 = A_16_16
    A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
    
    # Merge (32 + 8 = 40)
    C(A_32_0, sqr3)
    # Main
    A_64_0, count = M(A_32_0+ancilla2[:1108], A_32_32+ancilla2[1108:], count)
    # Side (A_8_32 * sqr3) -> A_40_0
    A_40_0, count = M(A_8_32+ancilla3[:1108], sqr3+ancilla3[1108:], count)
    C(A_32_0, sqr3)

    # S - Chain 2 (Prep 40 -> 104 side branch)
    A_8_0, sqr4 = S(A_8_0+A_8_32, 32) # clear sqr4
    A_40_0, A_40_64 = S(A_40_0+sqr3, 64)
    
    # Chain 1 (64 -> 128)
    A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
    sqr2 = A_32_32
    A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
    
    # Merge (64 + 40 = 104)
    C(A_64_0, sqr4)
    # Main
    A_128_0, count = M(A_64_0+ancilla2[:1108], A_64_64+ancilla2[1108:], count)
    # Side (A_40_64 * sqr4) -> A_104_0
    A_104_0, count = M(A_40_64+ancilla3[:1108], sqr4+ancilla3[1108:], count)
    C(A_64_0, sqr4)
    
    # Chain 1 cleanup & Final Merge Prep
    A_64_0, sqr1 = S(A_64_0+A_64_64, 64) # clear sqr1
    A_128_0, A_128_1 = S(A_128_0+sqr2, 1)
    A_104_0, A_104_129 = S(A_104_0+sqr4, 129)

    # Final Mul
    A_232_1, count = M(A_104_129+ancilla3[:1108], A_128_1+ancilla3[1108:], count)
    
    A_40_0, sqr3 = S(A_40_0+A_40_64, 64) # clear sqr3
    
    return A_232_1

def Inversion_283(gm, ml, a, n, sqr1, sqr2, sqr3, sqr4, count, ancilla2, ancilla3, toffoli_qubits, path_cfg, base_path):
    ANC_LIMIT = 1385
    TOF_COUNT = 1668

    def S(x, p): return Square(gm, ml, x, n, p, base_path)
    def M(x, y, c_off): 
        res = mul(gm, ml, x, y, toffoli_qubits[c_off:c_off+TOF_COUNT], n, path_cfg)
        return res, c_off + TOF_COUNT
    def C(x, y): CNOT_n(gm, x, y, n)

    # 1 -> 2
    a, A_1_1 = S(a+sqr1, 1) # a^2
    A_2_0, count = M(a+ancilla2[:ANC_LIMIT], A_1_1+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (2 -> 4)
    A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
    a, A_1_1 = S(a+A_1_1, 1) # update sqr1
    sqr1 = A_1_1
    A_4_0, count = M(A_2_2+ancilla2[:ANC_LIMIT], A_2_0+ancilla2[ANC_LIMIT:], count)

    # 4 -> 8
    A_2_0, A_2_2 = S(A_2_0+A_2_2, 2)
    sqr2 = A_2_2 # update sqr2
    A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
    A_8_0, count = M(A_4_0+ancilla2[:ANC_LIMIT], A_4_4+ancilla2[ANC_LIMIT:], count)

    # S - Chain 2 (Prep Side Branch)
    C(A_2_0, sqr3)
    sqr3, A_2_8 = S(sqr3+sqr4, 8)
    C(A_2_0, sqr3)

    # Chain 1 (8 -> 16)
    A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
    A_8_0, A_8_8 = S(A_8_0+sqr2, 8)
    
    # Merge
    C(A_8_0, sqr3)
    # Main Branch
    A_16_0, count = M(A_8_0+ancilla2[:ANC_LIMIT], A_8_8+ancilla2[ANC_LIMIT:], count)
    # Side Branch
    A_10_0, count = M(A_2_8+ancilla3[:ANC_LIMIT], sqr3+ancilla3[ANC_LIMIT:], count)
    C(A_8_0, sqr3)

    # S - Chain 2
    A_2_0, sqr4 = S(A_2_0+A_2_8, 8) # clear sqr4
    A_10_0, A_10_16 = S(A_10_0+sqr3, 16)

    # Chain 1
    A_8_0, sqr2 = S(A_8_0+A_8_8, 8) # update sqr2? logic says: A_8_0, sqr2 = Square...
    A_16_0, A_16_16 = S(A_16_0+sqr1, 16)

    # Merge
    C(A_16_0, sqr4)
    # Main Branch (16->32)
    A_32_0, count = M(A_16_0+ancilla2[:ANC_LIMIT], A_16_16+ancilla2[ANC_LIMIT:], count)
    # Side Branch (10->26)
    A_26_0, count = M(A_10_16+ancilla3[:ANC_LIMIT], sqr4+ancilla3[ANC_LIMIT:], count)
    C(A_16_0, sqr4)

    # Chain 1 (32 -> 64)
    A_16_0, A_16_16 = S(A_16_0+A_16_16, 16)
    sqr1 = A_16_16
    A_32_0, A_32_32 = S(A_32_0+sqr2, 32)
    A_64_0, count = M(A_32_0+ancilla2[:ANC_LIMIT], A_32_32+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (64 -> 128)
    A_32_0, A_32_32 = S(A_32_0+A_32_32, 32)
    sqr2 = A_32_32
    A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
    A_128_0, count = M(A_64_0+ancilla2[:ANC_LIMIT], A_64_64+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (128 -> 256)
    A_64_0, A_64_64 = S(A_64_0+A_64_64, 64)
    sqr1 = A_64_64
    A_128_0, A_128_128 = S(A_128_0+sqr2, 128)
    A_256_0, count = M(A_128_0+ancilla2[:ANC_LIMIT], A_128_128+ancilla2[ANC_LIMIT:], count)

    # Final Prep
    A_128_0, sqr2 = S(A_128_0+A_128_128, 128) # clear sqr2
    A_256_0, A_256_1 = S(A_256_0+sqr1, 1)
    A_26_0, A_26_257 = S(A_26_0+sqr4, 257)

    # Final Mul (Combine 256 + 26 = 282)
    # Note: Logic file says recursive_karatsuba(A_26_257, A_256_1, ..., ancilla2)
    # Mapping ancilla2(logic) -> ancilla3(circuit)
    A_282_1, count = M(A_26_257+ancilla3[:ANC_LIMIT], A_256_1+ancilla3[ANC_LIMIT:], count)

    # Cleanup sqr3 (as per logic file comments)
    A_10_0, sqr3 = S(A_10_0+A_10_16, 16)

    return A_282_1

def Inversion_571(gm, ml, a, n, sqr1, sqr2, sqr3, sqr4, count, ancilla2, ancilla3, toffoli_qubits, path_cfg, base_path):
    ANC_LIMIT = 2998  
    TOF_COUNT = 3569  

    def S(x, p): return Square(gm, ml, x, n, p, base_path)
    def M(x, y, c_off): 
        # ancilla2 is used for main branch
        # ancilla3 is used for side branch
        res = mul(gm, ml, x, y, toffoli_qubits[c_off:c_off+TOF_COUNT], n, path_cfg)
        return res, c_off + TOF_COUNT
    def C(x, y): CNOT_n(gm, x, y, n)

    # 1 -> 2
    a, A_1_1 = S(a+sqr1, 1) 
    A_2_0, count = M(a+ancilla2[:ANC_LIMIT], A_1_1+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (2 -> 4)
    A_2_0, A_2_2 = S(A_2_0+sqr2, 2)
    a, sqr1 = S(a+A_1_1, 1) 
    A_4_0, count = M(A_2_2+ancilla2[:ANC_LIMIT], A_2_0+ancilla2[ANC_LIMIT:], count)

    # 4 -> 8
    A_2_0, sqr2 = S(A_2_0+A_2_2, 2)
    A_4_0, A_4_4 = S(A_4_0+sqr1, 4)
    A_8_0, count = M(A_4_0+ancilla2[:ANC_LIMIT], A_4_4+ancilla2[ANC_LIMIT:], count)

    # S - Chain 2 (Prep Side)
    C(A_2_0, sqr3)
    sqr3, A_2_8 = S(sqr3+sqr4, 8)
    C(A_2_0, sqr3)

    # Chain 1 (8 -> 16)
    A_4_0, A_4_4 = S(A_4_0+A_4_4, 4)
    A_8_0, A_8_8 = S(A_8_0+sqr2, 8)

    # Merge
    C(A_8_0, sqr3)
    # Main
    A_16_0, count = M(A_8_0+ancilla2[:ANC_LIMIT], A_8_8+ancilla2[ANC_LIMIT:], count)
    # Side
    A_10_0, count = M(A_2_8+ancilla3[:ANC_LIMIT], sqr3+ancilla3[ANC_LIMIT:], count)
    C(A_8_0, sqr3)

    # S - Chain 2
    A_2_0, sqr4 = S(A_2_0+A_2_8, 8) # clear sqr4
    A_10_0, A_10_16 = S(A_10_0+sqr3, 16)

    # Chain 1
    A_8_0, A_8_8 = S(A_8_0+A_8_8, 8)
    sqr2 = A_8_8
    A_16_0, A_16_16 = S(A_16_0+sqr1, 16)

    # Merge
    C(A_16_0, sqr4)
    # Main
    A_32_0, count = M(A_16_0+ancilla2[:ANC_LIMIT], A_16_16+ancilla2[ANC_LIMIT:], count)
    # Side
    A_26_0, count = M(A_10_16+ancilla3[:ANC_LIMIT], sqr4+ancilla3[ANC_LIMIT:], count)
    C(A_16_0, sqr4)

    # S - Chain 2
    A_26_0, A_26_32 = S(A_26_0+sqr4, 32)
    A_10_0, sqr3 = S(A_10_0+A_10_16, 16) # clear sqr3

    # Chain 1 (32 -> 64)
    A_16_0, sqr1 = S(A_16_0+A_16_16, 16)
    A_32_0, A_32_32 = S(A_32_0+sqr2, 32)

    # Merge
    C(A_32_0, sqr3)
    # Main
    A_64_0, count = M(A_32_0+ancilla2[:ANC_LIMIT], A_32_32+ancilla2[ANC_LIMIT:], count)
    # Side (Note: A_26_32 + sqr3)
    A_58_0, count = M(A_26_32+ancilla3[:ANC_LIMIT], sqr3+ancilla3[ANC_LIMIT:], count)
    C(A_32_0, sqr3)

    # Chain 1 (64 -> 128)
    A_32_0, sqr2 = S(A_32_0+A_32_32, 32)
    A_64_0, A_64_64 = S(A_64_0+sqr1, 64)
    A_128_0, count = M(A_64_0+ancilla2[:ANC_LIMIT], A_64_64+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (128 -> 256)
    A_64_0, sqr1 = S(A_64_0+A_64_64, 64)
    A_128_0, A_128_128 = S(A_128_0+sqr2, 128)
    A_256_0, count = M(A_128_0+ancilla2[:ANC_LIMIT], A_128_128+ancilla2[ANC_LIMIT:], count)

    # Chain 1 (256 -> 512)
    A_128_0, sqr2 = S(A_128_0+A_128_128, 128)
    A_256_0, A_256_256 = S(A_256_0+sqr1, 256)
    A_512_0, count = M(A_256_0+ancilla2[:ANC_LIMIT], A_256_256+ancilla2[ANC_LIMIT:], count)

    # Final Prep
    A_256_0, sqr1 = S(A_256_0+A_256_256, 256) # clear sqr1
    A_512_0, A_512_1 = S(A_512_0+sqr2, 1)
    A_58_0, A_58_513 = S(A_58_0+sqr3, 513)

    # Final Mul (Combine 58 + 512 = 570)
    # Logic file uses ancilla2 for this step
    A_570_1, count = M(A_58_513+ancilla3[:ANC_LIMIT], A_512_1+ancilla3[ANC_LIMIT:], count)

    # Cleanup sqr4 (as per logic file)
    A_26_0, sqr4 = S(A_26_0+A_26_32, 32)

    return A_570_1