import os

# ----------------- Common gate and arithmetic wrappers -----------------
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
    # path_config must contain 'TD', 'A', 'C', 'Inv', 'mul_dim_1', and 'mul_dim_2'.
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
    ops, matrix_rows = matrix_loader.get_matrix_data(file_path, 2*n) # Square normally uses width 2n.
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
