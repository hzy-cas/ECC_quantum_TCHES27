import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from algorithms.implementations import Shor163, Shor233, Shor283, Shor571
from init.stats_utils import get_exact_resources_optimized

def measure_single_multiplication(target_n=163):
    DATA_ROOT = os.path.dirname(os.path.abspath(__file__))
    
    print(f"Starting Single Multiplication Resource Estimation for N={target_n}...")
    if target_n == 163:
        estimator = Shor163(DATA_ROOT)
    elif target_n == 233:
        estimator = Shor233(DATA_ROOT)
    elif target_n == 283:
        estimator = Shor283(DATA_ROOT)
    elif target_n == 571:
        estimator = Shor571(DATA_ROOT)
    else:
        print(f"Unsupported N={target_n}")
        return

    mul_size = estimator.get_mul_total_size()
    print(f"Required input register size: {mul_size} qubits")

    current_idx = 0

    reg_a = [current_idx + i for i in range(mul_size)]
    current_idx += mul_size

    reg_b = [current_idx + i for i in range(mul_size)]
    current_idx += mul_size

    reg_c = [current_idx + i for i in range(mul_size)]
    current_idx += mul_size
    
    total_qubits_allocated = current_idx

    try:
        estimator.gm.clear()
        
        result = estimator.mul(reg_a, reg_b, reg_c, 0)
        
    except FileNotFoundError:
        print("\n[WARNING] Matrix data file not found!")
        print("Please ensure the 'quantum_xxx' folder is in the current directory.")
        print("Missing files will cause inaccurate CNOT quantity and depth statistics (possibly 0).")
    except Exception as e:
        print(f"\n[Error] An exception occurred while performing multiplication: {e}")
        import traceback
        traceback.print_exc()
        return

    stats, full_depth, current_depth,toffoli_depth = get_exact_resources_optimized(estimator.gm, total_qubits_allocated)
    t_gates = stats['Toffoli_count']
    c_gates = stats['CNOT_count']
        
    print("\n" + "="*40)
    print("="*40)
    print(f"{'Metric':<20} | {'Value':>15}")
    print("-" * 40)
    print(f"{'Toffoli Gates':<20} | {t_gates:>15}")
    print(f"{'CNOT Gates':<20} | {c_gates:>15}")
    print(f"{'Full Depth':<20} | {full_depth:>15}")
    print(f"{'Current Depth':<20} | {current_depth:>15}")
    print(f"{'Toffoli Depth':<20} | {toffoli_depth:>15}")
    print(f"{'Qubits':<20} | {total_qubits_allocated:>15}")
    print("-" * 40)
        

if __name__ == "__main__":
    measure_single_multiplication(233)