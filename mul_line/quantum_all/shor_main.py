import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from algorithms.implementations import Shor163, Shor233, Shor283, Shor571

def main():
    TARGET_N = 163
    DATA_ROOT = os.path.dirname(os.path.abspath(__file__))
    OPS_COUNT = 2  
    # ===========================================

    print(f"Starting Resource Estimation for N={TARGET_N}...")
    
    estimator = None
    if TARGET_N == 163:
        estimator = Shor163(DATA_ROOT)
        OPS_COUNT = 328 
    elif TARGET_N == 233:
        estimator = Shor233(DATA_ROOT)
        OPS_COUNT = 468
    elif TARGET_N == 283:
        estimator = Shor283(DATA_ROOT)
        OPS_COUNT = 568
    elif TARGET_N == 571:
        estimator = Shor571(DATA_ROOT)
        OPS_COUNT = 1144
    else:
        print(f"Error: Unsupported N={TARGET_N}")
        return

    try:
        t_gates, c_gates, full_depth, current_depth,toffoli_depth, qubits = estimator.estimate_total_resources(OPS_COUNT)
        
        print("\n" + "="*40)
        print(f"{'Final Results (N=' + str(TARGET_N) + ')':^40}")
        print("="*40)
        print(f"{'Metric':<20} | {'Value':>15}")
        print("-" * 40)
        print(f"{'Toffoli Gates':<20} | {t_gates:>15}")
        print(f"{'CNOT Gates':<20} | {c_gates:>15}")
        print(f"{'Full Depth':<20} | {full_depth:>15}")
        print(f"{'Current Depth':<20} | {current_depth:>15}")
        print(f"{'Toffoli Depth':<20} | {toffoli_depth:>15}")
        print(f"{'Qubits':<20} | {qubits:>15}")
        print("-" * 40)
        
    except FileNotFoundError as e:
        print(f"\nError: Data files not found. Please check DATA_ROOT path.")
        print(f"Details: {e}")
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()