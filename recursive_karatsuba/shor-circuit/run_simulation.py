import multiprocessing
import pandas as pd
import time
import os
import sys

try:
    sys.path.append(os.path.join(os.path.dirname(__file__), "gf163"))
    import shor_out_Montgomery as shor_core
    
except ImportError:
    print("erro")
    sys.exit(1)

def worker_task(x):
    input_total = 33
    try:
        start_time = time.time()
        toffoli, cnot, depth, qubits = shor_core.accumulate_and_reduce(input_total, x)
        end_time = time.time()
        duration = end_time - start_time
        dw_cost = depth * qubits
        return {
            "Total_Inputs": input_total,
            "Parallel_Ops (x)": x,
            "Toffoli": toffoli,
            "CNOT": cnot,
            "Depth": depth,
            "Qubits": qubits,
            "DW_Cost": dw_cost,
            "Time_Seconds": round(duration, 4),
            "Status": "Success"
        }
        
    except Exception as e:
        
        return {
            "Parallel_Ops (x)": x,
            "Status": "Failed",
            "Error_Message": str(e)
        }

def main():
    x_values = list(range(2, 30, 1))
    num_cores = 20
    output_file = "163_Results.xlsx"
    
    print(f"Start simulation...")
    print(f"Parameter settings: Total Inputs=328, x range=[2, 4, ..., 128] (total {len(x_values)} groups)")
    print(f"Number of processes enabled: {num_cores}")
    
    global_start = time.time()
    
    results = []
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.map(worker_task, x_values)

    global_end = time.time()
    print(f"\ntime: {global_end - global_start:.2f} ")
    df = pd.DataFrame(results)

    if "Parallel_Ops (x)" in df.columns:
        df = df.sort_values("Parallel_Ops (x)")

    print("\nresult:")
    print(df.head())
    try:
        df.to_excel(output_file, index=False)
        print(f"\nSuccess: Results saved to {output_file}")
    except Exception as e:
        print(f"\nError: Unable to save Excel file ({e}). Attempting to save as CSV...")
        df.to_csv("Shor_Simulation_Results.csv", index=False)

if __name__ == "__main__":
    
    multiprocessing.freeze_support()
    main()