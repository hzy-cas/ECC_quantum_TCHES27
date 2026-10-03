import collections
import numpy as np
import array
import re

class FlatGateManager:
    def __init__(self):
        self.chunks = []
        self.current_chunk = array.array('Q')
        self.CHUNK_LIMIT = 10_000_000 
    
        self.SHIFT_OP = 60
        self.SHIFT_ARG = 30
        self.OP_DATA = 0
        self.OP_X = 1
        self.OP_CNOT = 2
        self.OP_TOFFOLI = 3
        
        self.OP_NOP = 15

    def _append(self, val):
        if len(self.current_chunk) >= self.CHUNK_LIMIT:
            self.chunks.append(self.current_chunk)
            self.current_chunk = array.array('Q')
        self.current_chunk.append(val)

    def add_X(self, t):
        val = (self.OP_X << self.SHIFT_OP) | t
        self._append(val)

    def add_CNOT(self, c, t):
        val = (self.OP_CNOT << self.SHIFT_OP) | (c << self.SHIFT_ARG) | t
        self._append(val)

    def add_Toffoli(self, c1, c2, t):
        val1 = (self.OP_TOFFOLI << self.SHIFT_OP) | (c1 << self.SHIFT_ARG) | c2
        self._append(val1)
        val2 = t 
        self._append(val2)

    def get_stats(self):
        cnot_count = 0
        toffoli_count = 0
        
        SHIFT_OP = self.SHIFT_OP
        OP_CNOT = self.OP_CNOT
        OP_TOFFOLI = self.OP_TOFFOLI
        
        all_data_blocks = self.chunks + [self.current_chunk]
        
        for chunk in all_data_blocks:
            for val in chunk:
                op = val >> SHIFT_OP
                if op == OP_CNOT:
                    cnot_count += 1
                elif op == OP_TOFFOLI:
                    toffoli_count += 1
                    
        return toffoli_count, cnot_count

    def current_pointer(self):
        total = sum(len(c) for c in self.chunks) + len(self.current_chunk)
        return total

    def calculate_depth(self, width):
        if width == 0: return 0
        wire_depths = array.array('I', [0] * width)
        
        SHIFT_OP = self.SHIFT_OP
        SHIFT_ARG = self.SHIFT_ARG
        MASK_ARG = (1 << SHIFT_ARG) - 1
        OP_X, OP_CNOT, OP_TOFFOLI, OP_NOP = self.OP_X, self.OP_CNOT, self.OP_TOFFOLI, self.OP_NOP
        
        all_chunks = self.chunks + [self.current_chunk]
        
        pending_toffoli_c1 = None
        pending_toffoli_c2 = None
        waiting_for_toffoli_data = False

        for chunk in all_chunks:
            for val in chunk:
                op = val >> SHIFT_OP
                
                if op == OP_NOP:
                    continue

                if waiting_for_toffoli_data:
                    t = val & MASK_ARG 
                    
                    c1 = pending_toffoli_c1
                    c2 = pending_toffoli_c2
                    
                    d_c1 = wire_depths[c1]
                    d_c2 = wire_depths[c2]
                    d_t = wire_depths[t]
                    
                    mx = d_c1 if d_c1 > d_c2 else d_c2
                    if d_t > mx: mx = d_t
                    new_depth = mx + 1
                    
                    wire_depths[c1] = new_depth
                    wire_depths[c2] = new_depth
                    wire_depths[t] = new_depth
                    
                    waiting_for_toffoli_data = False
                    
                elif op == OP_CNOT:
                    t = val & MASK_ARG
                    c = (val >> SHIFT_ARG) & MASK_ARG
                    d_c = wire_depths[c]
                    d_t = wire_depths[t]
                    new_depth = (d_c if d_c > d_t else d_t) + 1
                    wire_depths[c] = new_depth
                    wire_depths[t] = new_depth
                    
                elif op == OP_TOFFOLI:
                    c2 = val & MASK_ARG
                    c1 = (val >> SHIFT_ARG) & MASK_ARG
                    pending_toffoli_c1 = c1
                    pending_toffoli_c2 = c2
                    waiting_for_toffoli_data = True
                    
                elif op == OP_X:
                    t = val & MASK_ARG
                    wire_depths[t] += 1
        
        return max(wire_depths) if wire_depths else 0

    def replay_reverse(self, start_idx, end_idx):
        curr_ptr = self.current_pointer()
        if start_idx >= end_idx: return

        all_chunks = self.chunks + [self.current_chunk]
        
        SHIFT_OP = self.SHIFT_OP
        SHIFT_ARG = self.SHIFT_ARG
        MASK_ARG = (1 << SHIFT_ARG) - 1
        OP_X, OP_CNOT, OP_TOFFOLI = self.OP_X, self.OP_CNOT, self.OP_TOFFOLI
        
        idx = curr_ptr
        
        for chunk_idx, chunk in enumerate(reversed(all_chunks)):
            real_chunk_idx = len(all_chunks) - 1 - chunk_idx
            
            length = len(chunk)
            chunk_start_abs = idx - length
            
            if chunk_start_abs >= end_idx:
                idx -= length
                continue
            
            if idx <= start_idx:
                break
                
            for i in range(length - 1, -1, -1):
                abs_pos = chunk_start_abs + i
                
                if abs_pos >= end_idx: continue
                if abs_pos < start_idx: break 
                
                val = chunk[i]
                op = val >> SHIFT_OP
                
                if op == 0: # OP_DATA (Toffoli t)
                    t = val
                    prev_i = i - 1
                    
                    val_prev = None
                    if prev_i >= 0:
                        val_prev = chunk[prev_i]
                    elif real_chunk_idx > 0: 
                        prev_chunk = all_chunks[real_chunk_idx - 1]
                        if prev_chunk:
                            val_prev = prev_chunk[-1]
                    
                    if val_prev is not None:
                        c2 = val_prev & MASK_ARG
                        c1 = (val_prev >> SHIFT_ARG) & MASK_ARG
                        self.add_Toffoli(c1, c2, t)
                
                elif op == OP_TOFFOLI:
                    pass 

                elif op == OP_CNOT:
                    t = val & MASK_ARG
                    c = (val >> SHIFT_ARG) & MASK_ARG
                    self.add_CNOT(c, t)
                    
                elif op == OP_X:
                    t = val & MASK_ARG
                    self.add_X(t)
            
            idx -= length

    def clear(self):
        self.chunks = []
        self.current_chunk = array.array('Q')


class MatrixCache:
    def __init__(self):
        self.cache = {}

    def get_matrix_data(self, file_path, n):
        if file_path in self.cache:
            return self.cache[file_path]
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            return [], []

        ops = []
        start_index = n
        
        for i, line in enumerate(lines):
            if "CNOT Depth" in line:
                start_index = i + 1
                break
        for line in lines[start_index:]:
            line = line.strip()
            if not line or '=' not in line: continue
            try:
                left, right = line.split('=')
                target = int(left.split('[')[1].split(']')[0])
                parts = right.split('^')
                if len(parts) < 2: continue
                control = int(parts[1].strip().split('[')[1].split(']')[0])
                y_index = None
                if 'y[' in line:
                    y_index = int(line.split('y[')[1].split(']')[0])
                ops.append((target, control, y_index))
            except:
                continue

        matrix_rows = []
        line_count = 0
        for line in lines:
            if line_count >= n: break
            stripped = line.strip()
            if stripped and not stripped.startswith('Original') and not stripped.startswith('Reduced'):
                ones_indices = [j for j, char in enumerate(stripped) if char == '1']
                matrix_rows.append(ones_indices)
                line_count += 1
        
        data = (ops, matrix_rows)
        self.cache[file_path] = data
        return data