"""Read AC linear gate streams and cache their matrix metadata."""

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
            # Preserve legacy behavior: a missing file yields empty data.
            return [], []

        ops = []
        start_index = n
        
        # 1. Find the CNOT Depth marker and locate the first operation.
        for i, line in enumerate(lines):
            if "CNOT Depth" in line:
                start_index = i + 1
                break
        
        # 2. Parse CNOT operations.
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

        # 3. Parse the leading binary matrix (required by the Shor circuit).
        matrix_rows = []
        line_count = 0
        for line in lines:
            if line_count >= n: break
            stripped = line.strip()
            # Ignore non-matrix lines.
            if stripped and not stripped.startswith('Original') and not stripped.startswith('Reduced'):
                ones_indices = [j for j, char in enumerate(stripped) if char == '1']
                matrix_rows.append(ones_indices)
                line_count += 1
        
        data = (ops, matrix_rows)
        self.cache[file_path] = data
        return data
