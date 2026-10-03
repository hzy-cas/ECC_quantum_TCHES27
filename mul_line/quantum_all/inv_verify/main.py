import os
import time
import random
import math


from qasm import FlatGateManager,MatrixCache
from poly_utils import GF2n
import circuit_library as lib

# ==========================================

# ==========================================
class BasisTransformer:
    def __init__(self, n, l2x_path, x2l_path):
        self.n = n
        self.matrix_l2x = self._load_matrix(l2x_path)
        self.matrix_x2l = self._load_matrix(x2l_path)

    def _load_matrix(self, path):
        """Read a binary matrix stored as rows."""
        mat = []
        if not os.path.exists(path):
            print(f"[Warning] Basis matrix not found: {path}")
            return None
        try:
            with open(path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line: continue
                    row = [int(c) for c in line]
                    if len(row) != self.n:
                        print(f"[Error] Matrix row length mismatch in {path}")
                        return None
                    mat.append(row)
            if len(mat) != self.n:
                print(f"[Error] Matrix height mismatch in {path}")
                return None
            return mat
        except Exception as e:
            print(f"[Error] Failed to load matrix {path}: {e}")
            return None

    def apply(self, vec_bits, direction='x2l'):
        """Compute M times a column vector over GF(2), with coordinates in file order."""
        matrix = self.matrix_x2l if direction == 'x2l' else self.matrix_l2x
        if matrix is None:
            return vec_bits # Fallback: return as is if matrix missing



        res = [0] * self.n
        for i in range(self.n):
            row = matrix[i]
            val = 0
            for j in range(self.n):
                if row[j] == 1 and vec_bits[j] == 1:
                    val ^= 1
            res[i] = val
        return res

# ==========================================

# ==========================================
def run_simulation(gm, width, input_bits, output_indices):
    print(f"Running simulation on {width} qubits...")
    start_t = time.time()
    Q = [0] * width


    for i, bit in enumerate(input_bits):
        Q[i] = bit


    SHIFT_OP = 60
    SHIFT_ARG = 30
    MASK_ARG = (1 << 30) - 1
    OP_X, OP_CNOT, OP_TOFFOLI, OP_NOP = 1, 2, 3, 15

    all_chunks = gm.chunks + [gm.current_chunk]

    for chunk in all_chunks:
        i = 0
        length = len(chunk)
        while i < length:
            val = chunk[i]
            op = val >> SHIFT_OP

            if op == OP_NOP:
                i += 1
            elif op == OP_CNOT:
                t = val & MASK_ARG
                c = (val >> SHIFT_ARG) & MASK_ARG
                Q[t] ^= Q[c]
                i += 1
            elif op == OP_X:
                t = val & MASK_ARG
                Q[t] ^= 1
                i += 1
            elif op == OP_TOFFOLI:
                c2 = val & MASK_ARG
                c1 = (val >> SHIFT_ARG) & MASK_ARG
                if i + 1 < length:
                    t = chunk[i+1]
                    Q[t] ^= (Q[c1] & Q[c2])
                    i += 2
                else:
                    break
            else:
                i += 1

    print(f"Simulation finished in {time.time() - start_t:.2f}s")
    return [Q[idx] for idx in output_indices]

def int_to_bits(val, n):
    return [(val >> i) & 1 for i in range(n)]

def bits_to_int(bits):
    val = 0
    for i, b in enumerate(bits):
        if b: val |= (1 << i)
    return val
# ==========================================

# ==========================================
def main():


    TARGET_N = 283


    BASE_PATH = "/home/kekai/home/mul_line/quantum_all"



    CONFIGS = {
        163: {
            'poly': [163, 7, 6, 3, 0],
            'ancilla_mul': 743,
            'toffoli_block': 906,
            'mul_dims': (325, 906),
            'func': lib.Inversion_163
        },
        233: {
            'poly': [233, 74, 0],
            'ancilla_mul': 1108,
            'toffoli_block': 1341,
            'mul_dims': (465, 1341),
            'func': lib.Inversion_233
        },
        283: {
            'poly': [283, 12, 7, 5, 0],
            'ancilla_mul': 1385,
            'toffoli_block': 1668,
            'mul_dims': (565, 1668),
            'func': lib.Inversion_283
        },
        571: {
            'poly': [571, 10, 5, 2, 0],
            'ancilla_mul': 2998,
            'toffoli_block': 3569,
            'mul_dims': (1141, 3569),
            'func': lib.Inversion_571
        }
    }

    if TARGET_N not in CONFIGS:
        print(f"Error: Unsupported N={TARGET_N}")
        return

    cfg = CONFIGS[TARGET_N]
    n = TARGET_N
    print(f"=== Target: GF(2^{n}) [Basis Transformed] ===")

    gm = FlatGateManager()
    ml = MatrixCache()



    base_dir = os.path.join(BASE_PATH,  f"quantum_{n}")
    if not os.path.exists(base_dir):
        base_dir = os.path.join(BASE_PATH, f"quantum_{n}")

    l2x_path = os.path.join(base_dir, f"{n}_L2X.txt")
    x2l_path = os.path.join(base_dir, f"{n}_X2L.txt")

    transformer = BasisTransformer(n, l2x_path, x2l_path)


    print("Building circuit...")

    path_config = {
        'TD': os.path.join(base_dir, "CNOT_mul", f"result_{n}_TD_seq.txt"),
        'A': os.path.join(base_dir, "CNOT_mul", f"result_{n}_A.txt"),
        'C': os.path.join(base_dir, "CNOT_mul", f"result_{n}_C.txt"),
        'Inv': os.path.join(base_dir, "CNOT_mul", f"result_{n}_CT2D_inv_seq.txt"),
        'mul_dim_1': cfg['mul_dims'][0],
        'mul_dim_2': cfg['mul_dims'][1]
    }

    a = [i for i in range(n)]
    sqr1 = [n + i for i in range(n)]
    sqr2 = [2*n + i for i in range(n)]
    sqr3 = [3*n + i for i in range(n)]
    sqr4 = [4*n + i for i in range(n)]

    offset = 5 * n
    toffoli_size = cfg['toffoli_block'] * 20
    toffoli_qubits = [offset + i for i in range(toffoli_size)]
    offset += len(toffoli_qubits)

    anc_mul_size = cfg['ancilla_mul']
    ancilla2 = [offset + i for i in range(anc_mul_size * 2)]
    offset += len(ancilla2)
    ancilla3 = [offset + i for i in range(anc_mul_size * 2)]
    offset += len(ancilla3)

    count = 0


    result_qubits = cfg['func'](
        gm, ml, a, n, sqr1, sqr2, sqr3, sqr4, count,
        ancilla2, ancilla3, toffoli_qubits, path_config, base_dir
    )

    t_count, c_count = gm.get_stats()
    print("-" * 30)
    print(f"Toffoli: {t_count}, CNOT: {c_count}, Qubits: {offset}")


    print("\nVerifying correctness with Basis Transformation...")


    input_val = random.randint(1, (1 << n) - 1)
    print(f"Input (Std Basis): {input_val}")


    gf = GF2n(n, cfg['poly'])
    # exponent = (1 << 128) - 1  # 2^n - 2
    # expected_val = gf.power(input_val, exponent)
    expected_val = gf.inverse(input_val)
    print(f"Expected (Std)   : {expected_val}")


    input_bits_std = int_to_bits(input_val, n)
    if transformer.matrix_x2l:
        print("Applying X -> L transformation to input...")
        input_bits_circuit = transformer.apply(input_bits_std, 'x2l')
    else:
        print("[Warn] No X->L matrix found, using standard basis for input.")
        input_bits_circuit = input_bits_std


    output_bits_circuit = run_simulation(gm, offset, input_bits_circuit, result_qubits)


    if transformer.matrix_l2x:
        print("Applying L -> X transformation to output...")
        output_bits_std = transformer.apply(output_bits_circuit, 'l2x')
    else:
        print("[Warn] No L->X matrix found, assuming standard basis output.")
        output_bits_std = output_bits_circuit

    output_val = bits_to_int(output_bits_std)
    print(f"Output (Std)     : {output_val}")


    if output_val == expected_val:
        print("\n[SUCCESS] Result matches!")
    else:
        print("\n[FAILURE] Result mismatch.")
        # Debug info
        prod = gf.multiply(input_val, output_val)
        print(f"Check a * a^-1 (Std) = {prod}")

    gm.clear()

if __name__ == "__main__":
    main()
