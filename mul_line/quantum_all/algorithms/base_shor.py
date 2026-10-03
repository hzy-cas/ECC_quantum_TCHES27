import math
import sys
import gc
import os
import random
from init.qasm import FlatGateManager, MatrixCache
from init.stats_utils import get_exact_resources_optimized

class BaseShor:

    def __init__(self, n, block_size, nums2_const, data_root):
        self.n = n
        self.block_size = block_size
        self.nums2_const = nums2_const  
        self.data_root = data_root     
        
        self.gm = FlatGateManager()     
        self.matrix_loader = MatrixCache()
        self.Toffoli_qubits = []      

    # ================= Gate Wrappers=================
    def X(self, a):
        self.gm.add_X(a)
    
    def Toffoli_gate(self, a, b, c):
        self.gm.add_Toffoli(a, b, c)

    def CNOT(self, a, b):
        self.gm.add_CNOT(a, b)

    def CNOT_N(self, a, b):
        for i in range(self.n): self.gm.add_CNOT(a[i], b[i])

    def CONST_ADD(self, target_reg, constant_val):
        for i in range(self.n):
            if (constant_val >> i) & 1:
                self.X(target_reg[i])

    def CSWAP(self, c, a, b):
        self.CNOT(b, a)
        self.Toffoli_gate(c, a, b)
        self.CNOT(b, a)

    def CSWAP_N(self, ctrl, reg_a, reg_b):
        for i in range(self.n):
            self.CSWAP(ctrl, reg_a[i], reg_b[i])

    def Round_constant_XOR(self, k, rc, bit):
        for i in range(bit):
            if (rc >> i & 1): self.X(k[i])
    def copy_parallel(self, ctrl, ancilla, total_count):
        controls = [ctrl]
        anc_idx = 0
        while len(controls) < total_count:
            current_sources = controls[:]
            for src in current_sources:
                if len(controls) >= total_count: break
                if anc_idx >= len(ancilla):
                    raise ValueError(f"Not enough ancilla: need {total_count-1}, has {len(ancilla)}")
                target = ancilla[anc_idx]
                anc_idx += 1
                self.CNOT(src, target)
                controls.append(target)
        return controls
    # =================Basic Ops=================
    def get_matrix_path(self, folder, filename):
        return os.path.join(self.data_root, f"quantum_{self.n}", folder, filename)

    def mul_matrix(self, x, n_rows, file_path):
        ops, matrix_rows = self.matrix_loader.get_matrix_data(file_path, n_rows)
        y = {}

        for target, control, y_index in ops:
            self.CNOT(x[control], x[target])
            if y_index is not None: y[y_index] = x[target]
        
        if not y: return []
        y_list = [0] * n_rows
        for i, val in y.items(): y_list[i] = val
        for i in range(len(y_list)):
            if y_list[i] == 0 and i not in y:
                if i < len(matrix_rows):
                    if matrix_rows[i]:
                        j = matrix_rows[i][0]
                        if j < len(x): y_list[i] = x[j]
        return y_list

    def Square(self, x, n_rows, power):
        file_path = self.get_matrix_path(f"square_{self.n}", f"square_Matrix_2_{power}.txt")
        y_list = self.mul_matrix(x, n_rows, file_path)
        return y_list[:int(n_rows/2)], y_list[-int(n_rows/2):]

    def Squaring_plus_input(self, x, n_rows):
        file_path = self.get_matrix_path(f"square_{self.n}", "square_Matrix_Squaring_plus.txt")
        y_list = self.mul_matrix(x, n_rows, file_path)
        return y_list[:int(n_rows/2)], y_list[-int(n_rows/2):]

    def mul(self, a0, b0, c0, toffoli_offset):
        path_td = self.get_matrix_path("CNOT_mul", f"result_{self.n}_TD_seq.txt")
        path_a = self.get_matrix_path("CNOT_mul", f"result_{self.n}_A.txt")
        path_c = self.get_matrix_path("CNOT_mul", f"result_{self.n}_C.txt")
        path_inv = self.get_matrix_path("CNOT_mul", f"result_{self.n}_CT2D_inv_seq.txt")
        
        slice1 = 2 * self.n - 1
        total_mul_size = self.get_mul_total_size()

        begin = self.gm.current_pointer()

        a1 = self.mul_matrix(a0[:slice1], slice1, path_td)
        b1 = self.mul_matrix(b0[:slice1], slice1, path_td)
        
        a2 = self.mul_matrix(a1 + a0[slice1:total_mul_size], total_mul_size, path_a)
        b2 = self.mul_matrix(b1 + b0[slice1:total_mul_size], total_mul_size, path_a)
        
        end = self.gm.current_pointer()

        for i in range(total_mul_size):
            self.Toffoli_gate(a2[i], b2[i], c0[i])

        c1 = self.mul_matrix(c0[:total_mul_size], total_mul_size, path_c)
        c2 = self.mul_matrix(c1[:slice1], slice1, path_inv)

        self.gm.replay_reverse(begin, end)
        return c2[:self.n]

    # =================Inversion=================
    def cons_mul(self, xlist, count, ancilla):
        tree = []
        tree.append(xlist[:])
        num = self.block_size * 2
        total_mul_size = self.get_mul_total_size()
        
        current_level = xlist[:]
        while len(tree[-1]) > 1:
            next_level = []
            i = 0
            j = 0
            while i < len(current_level):
                if i + 1 < len(current_level):
                    product = self.mul(
                        current_level[i] + ancilla[num*j : num*j+self.block_size], 
                        current_level[i+1] + ancilla[num*j+self.block_size : num*j+self.block_size*2], 
                        self.Toffoli_qubits[count : count+total_mul_size], 
                        0
                    )
                    count += total_mul_size
                    next_level.append(product)
                    i += 2
                    j += 1
                else:

                    next_level.append(current_level[i])
                    i += 1
            tree.append(next_level[:])
            current_level = next_level
        return tree[-1][0], count, ancilla, tree

    def Montgomerytrick(self, x1x2_list, y1_list, count, ancilla, ancilla2, sqr):
        num = self.block_size * 2
        parallel_num = len(x1x2_list)
        total_mul_size = self.get_mul_total_size()
        
        if parallel_num == 0: return [], ancilla, count

        if parallel_num == 1:
            x1x2_inv, _, count = self.inversion_logic(
                x1x2_list[0], 
                sqr, 
                ancilla[:num],    
                ancilla[-num:],  
                count
            )
            lambda_result = self.mul(
                y1_list[0] + ancilla[:self.block_size], 
                x1x2_inv + ancilla[self.block_size : self.block_size*2], 
                self.Toffoli_qubits[count : count+total_mul_size],
                0
            )
            count += total_mul_size
            return [lambda_result], ancilla, count
        product_final, count, ancilla, tree = self.cons_mul(x1x2_list, count, ancilla)
        product_inv, _, count = self.inversion_logic(
            product_final, 
            sqr, 
            ancilla[:num],  
            ancilla[-num:], 
            count
        )
        
        inv_tree = [[product_inv]]

        for level in range(len(tree) - 2, -1, -1):
            current_level = tree[level]
            current_inv_level = inv_tree[-1]
            next_inv_level = []
            inv_idx = 0
            i = 0
            idx = 0
            while i < len(current_level):
                if i + 1 < len(current_level):
                    begin = self.gm.current_pointer()
                    self.CNOT_N(current_inv_level[inv_idx], ancilla2[idx:idx+self.n])
                    copy = ancilla2[idx:idx+self.n]
                    idx += self.n
                    end = self.gm.current_pointer()
                    a_inv = self.mul(
                        current_inv_level[inv_idx] + ancilla[num*i : num*i+self.block_size],
                        current_level[i+1] + ancilla[num*i+self.block_size : num*i+self.block_size*2],
                        self.Toffoli_qubits[count : count+total_mul_size],
                        0
                    )
                    count += total_mul_size
                    b_inv = self.mul(
                        copy + ancilla[num*(i+1) : num*(i+1)+self.block_size],
                        current_level[i] + ancilla[num*(i+1)+self.block_size : num*(i+1)+self.block_size*2],
                        self.Toffoli_qubits[count : count+total_mul_size],
                        0
                    )
                    count += total_mul_size

                    self.gm.replay_reverse(begin, end)
                    next_inv_level.append(a_inv)
                    next_inv_level.append(b_inv)
                    i += 2
                    inv_idx += 1
                else:
                    next_inv_level.append(current_inv_level[inv_idx])
                    i += 1
                    inv_idx += 1
            inv_tree.append(next_inv_level)
        
        inverses = inv_tree[-1]
        lambda_list = []
        for i in range(parallel_num):
            lambda_i = self.mul(
                y1_list[i] + ancilla[num*i : num*i+self.block_size],
                inverses[i] + ancilla[num*i+self.block_size : num*i+self.block_size*2],
                self.Toffoli_qubits[count : count+total_mul_size],
                0
            )
            count += total_mul_size
            lambda_list.append(lambda_i)
        return lambda_list, ancilla, count
    def shor_accumulation_step(self, parallel_num, q_list, x1_list, y1_list, x2_val_list, y2_val_list, ancilla, x1x2_list, ancilla2, sqr, cswap_ancilla):
        curve_a = self.get_curve_a()
        total_mul_size = self.get_mul_total_size()
        num = self.block_size * 2
        n = self.n

        for i in range(parallel_num):
            self.CNOT_N(x1_list[i], x1x2_list[i])       
            self.CONST_ADD(y1_list[i], y2_val_list[i])  
            self.CONST_ADD(x1x2_list[i], x2_val_list[i])
        count = 0
        lamba, ancilla, count = self.Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr)

        for i in range(parallel_num):

            self.CONST_ADD(x1x2_list[i], x2_val_list[i] ^ curve_a) 
            self.CONST_ADD(y1_list[i], y2_val_list[i]) 

            lamba[i], x1x2_list[i] = self.Squaring_plus_input(lamba[i] + x1x2_list[i], 2*self.n)
            y_target = self.mul(
                x1x2_list[i] + ancilla[num*i : num*i+self.block_size],
                lamba[i] + ancilla[num*i+self.block_size : num*i+self.block_size*2],
                self.Toffoli_qubits[count : count+total_mul_size],
                0
            )
            count += total_mul_size

            self.CONST_ADD(x1x2_list[i], x2_val_list[i])
            self.CONST_ADD(y_target, y2_val_list[i])     
            self.CNOT_N(x1x2_list[i], y_target)       

            current_cswap_anc = cswap_ancilla[i * (2 * n - 1) : (i + 1) * (2 * n - 1)]
            
            copy_start_ptr = self.gm.current_pointer()

            control_bits = self.copy_parallel(q_list[i], current_cswap_anc, 2 * n)
            
            copy_end_ptr = self.gm.current_pointer()
            for bit_idx in range(n):
                self.CSWAP(control_bits[bit_idx], x1_list[i][bit_idx], x1x2_list[i][bit_idx])
            
            for bit_idx in range(n):
                self.CSWAP(control_bits[n + bit_idx], y1_list[i][bit_idx], y_target[bit_idx])
                
            # Uncompute Fan-out
            self.gm.replay_reverse(copy_start_ptr, copy_end_ptr)

    # ================= Phase 3 =================
    def shor_reduction_step(self, parallel_num, x1_list, y1_list, x2_list, y2_list, ancilla, x1x2_list, ancilla2, sqr):
        curve_a = self.get_curve_a()
        total_mul_size = self.get_mul_total_size()
        num = self.block_size * 2
        for i in range(parallel_num):
            self.CNOT_N(x1_list[i], x1x2_list[i])
            self.CNOT_N(x2_list[i], x1x2_list[i]) 
            self.CNOT_N(y2_list[i], y1_list[i])  
        count = 0
        lamba, ancilla, count = self.Montgomerytrick(x1x2_list, y1_list, count, ancilla, ancilla2, sqr)

        for i in range(parallel_num):
            self.CONST_ADD(x1x2_list[i], curve_a)

            lamba[i], x1x2_list[i] = self.Squaring_plus_input(lamba[i] + x1x2_list[i], 2*self.n)

            self.CNOT_N(x2_list[i], x1x2_list[i]) 
            
            y_target = self.mul(
                x1x2_list[i] + ancilla[num*i : num*i+self.block_size],
                lamba[i] + ancilla[num*i+self.block_size : num*i+self.block_size*2],
                self.Toffoli_qubits[count : count+total_mul_size],
                0
            )
            count += total_mul_size

            self.CNOT_N(y2_list[i], y_target) 
            self.CNOT_N(x2_list[i], x1x2_list[i]) 
            self.CNOT_N(x1x2_list[i], y_target) 
    def run_shor_logic(self, parallel_num, mode='accumulation'):
        self.gm.clear()
        gc.collect()
        
        n = self.n
        nums1 = self.block_size * 4 * parallel_num
        total_mul_size = self.get_mul_total_size()
        nums2 = total_mul_size * (3 * (parallel_num - 1) + self.nums2_const + 2 * parallel_num)

        current_idx = 0

        if mode == 'accumulation':
            q_list = [current_idx + i for i in range(parallel_num)]
            current_idx += parallel_num
        else:
            q_list = [] 

        x1 = [current_idx + i for i in range(n*parallel_num)]
        current_idx += n*parallel_num

        if mode == 'reduction':
            x2 = [current_idx + i for i in range(n*parallel_num)]
            current_idx += n*parallel_num
        else:
            x2 = [] 

        x1x2 = [current_idx + i for i in range(n*parallel_num)]
        current_idx += n*parallel_num

        y1 = [current_idx + i for i in range(n*parallel_num)]
        current_idx += n*parallel_num

        if mode == 'reduction':
            y2 = [current_idx + i for i in range(n*parallel_num)]
            current_idx += n*parallel_num
        else:
            y2 = []

        x1_list = [x1[i*n:(i+1)*n] for i in range(parallel_num)]
        y1_list = [y1[i*n:(i+1)*n] for i in range(parallel_num)]
        x1x2_list = [x1x2[i*n:(i+1)*n] for i in range(parallel_num)]
        
        x2_list_q = []
        y2_list_q = []
        x2_list_val = []
        y2_list_val = []

        if mode == 'reduction':
            x2_list_q = [x2[i*n:(i+1)*n] for i in range(parallel_num)]
            y2_list_q = [y2[i*n:(i+1)*n] for i in range(parallel_num)]
        else:
            full_ones = (1 << n) - 1
            x2_list_val = [full_ones for _ in range(parallel_num)]
            y2_list_val = [full_ones for _ in range(parallel_num)]

        ancilla = [current_idx + i for i in range(int(nums1))]
        length = current_idx + len(ancilla)
        
        self.Toffoli_qubits = [length+i for i in range(int(nums2))]
        length += len(self.Toffoli_qubits)
        
        sqr = [length + i for i in range(4*n)]
        length += 4*n
        
        ancilla2 = [length + i for i in range(n*int(parallel_num/2))]
        length += int(parallel_num/2)*n

        if mode == 'accumulation':
            cswap_anc_size = (2 * n - 1) * parallel_num
            cswap_ancilla = [length + i for i in range(cswap_anc_size)]
            length += cswap_anc_size

        if mode == 'accumulation':
            self.shor_accumulation_step(
                parallel_num, q_list, x1_list, y1_list, 
                x2_list_val, y2_list_val, 
                ancilla, x1x2_list, ancilla2, sqr, cswap_ancilla
            )
        else:
            self.shor_reduction_step(
                parallel_num, x1_list, y1_list, 
                x2_list_q, y2_list_q, 
                ancilla, x1x2_list, ancilla2, sqr
            )
        
        stats, full_depth, current_depth,toffoli_depth = get_exact_resources_optimized(self.gm, length)
        t_count = stats['Toffoli_count']
        c_count = stats['CNOT_count']
        
        self.gm.clear()
        return length, t_count, c_count, full_depth, current_depth,toffoli_depth
    def estimate_total_resources(self, total_inputs):
        total_Toffoli = 0
        total_CNOT = 0
        total_fDepth = 0
        total_cDepth = 0
        total_tDepth = 0
        peak_qubits = 0
        clear = 0 
        
        n = self.n
        current = total_inputs
        level = 0
        
        print(f"{'Estimation Strategy':^60}")
        print(f"{'Strict Rules Applied':^60}")
        
        while current > 1:
            ops = current // 2
            if level == 0:
                length, t, c, fd,cd,td = self.run_shor_logic(ops, mode='accumulation')
                
                # 1. Depth: 2*ceil(log(2n)) + 1 + d
                log_2n = math.ceil(math.log2(2 * n))
                depth_overhead = 2 * log_2n + 1
                total_fDepth += (depth_overhead + fd)
                total_cDepth += (depth_overhead + cd)
                total_tDepth +=  td
                
                # 2. CNOT: 2*ops*(2n-1) + 2n*ops + c
                cnot_overhead = 2 * ops * (2 * n - 1) + 2 * n * ops
                total_CNOT += (cnot_overhead + c)
                
                # 3. Toffoli:
                total_Toffoli += t
                
                # 4. Qubit: length - ops 
                peak_qubits = length - ops
                
                # 5. Clear: 
                #  2n + n*int(ops/2) + (2n - 1)*ops
                clear = 2 * n + n * int(ops / 2) + (2 * n - 1) * ops +self.block_size * 4 * ops
                
                print(f"Level {level} (ops={ops}): Depth={depth_overhead+fd}, Qubits={peak_qubits}, Clear={clear}")

            else:
                length, t, c, fd,cd,td = self.run_shor_logic(ops, mode='reduction')
                
                total_fDepth += fd
                total_cDepth += cd
                total_tDepth +=  td
                total_CNOT += c
                total_Toffoli += t
                
                # 1. Demand : length - 4n*ops
                demand = length - 4 * n * ops
                
                # 2. Qubit & Clear 
                if demand > clear:
                    peak_qubits += (demand - clear)
                    # new Clear: 2n + n*int(ops/2)
                    clear = 2 * n + n * int(ops / 2)+self.block_size * 4 * ops
                else:
                    over = clear - demand
                    # new Clear: 2n + n*int(ops/2) + Over
                    clear = 2 * n + n * int(ops / 2) + over +self.block_size * 4 * ops
            
                print(f"Level {level} (ops={ops}): cDepth={cd},fDepth={fd}, Demand={demand}, Qubits={peak_qubits}, Clear={clear}")
            
            current = ops + (current % 2)
            level += 1
            
        # === Final Calculation (Uncompute & Result Copy) ===
        # 1. Toffoli: 2 * Total
        final_Toffoli = 2 * total_Toffoli
        
        # 2. CNOT: 2 * Total + 2n
        final_CNOT = 2 * total_CNOT + 2 * n
        
        # 3. Depth: 2 * Total + 2 
        final_cDepth = 2 * total_cDepth + 2
        final_fDepth = 2 * total_fDepth + 2
        final_tDepth = 2 * total_tDepth
        
        # 4. Qubits: Peak + 4n + 2 (Result Register + q)
        final_Qubits = peak_qubits + 4 * n + 2
        return final_Toffoli, final_CNOT, final_fDepth,final_cDepth,final_tDepth,  final_Qubits

    def get_curve_a(self):
        return 3

    def get_mul_total_size(self):
        raise NotImplementedError