# Inverse resource estimation
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
import importlib
import math
import numpy as np
import gc
import random 
from init.stats_utils import get_exact_resources_optimized
# ==========================================
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from init.basic_gates_163 import *
    from init.inversion_163 import *

TARGET_N =163
print(f"=== Initializing Backend for GF(2^{TARGET_N}) ===")
try:
    gates_module_name = f"init.basic_gates_{TARGET_N}"
    inversion_module_name = f"init.inversion_{TARGET_N}"
    gates_mod = importlib.import_module(gates_module_name)
    inv_mod = importlib.import_module(inversion_module_name)
    this_module = sys.modules[__name__]
    def inject_from_module(module):
        for attr in dir(module):
            if not attr.startswith('__'):
                setattr(this_module, attr, getattr(module, attr))
    inject_from_module(gates_mod)
    inject_from_module(inv_mod)
    globals().update({attr: getattr(gates_mod, attr) for attr in dir(gates_mod) if not attr.startswith('__')})
    globals().update({attr: getattr(inv_mod, attr) for attr in dir(inv_mod) if not attr.startswith('__')})
    GF_N = getattr(gates_mod, 'GF_N')
except ImportError as e:
    print(f"\n[Error] Failed to load backend for GF({TARGET_N}).")
    import traceback
    traceback.print_exc()
    sys.exit(1)

sys.setrecursionlimit(5000)




n = GF_N
gm.clear()
gates_mod.cnt = 0 
    
nums1 = KARATSUBA_ANC_BASE*2
nums2 = TOFFOLI_BASE * ( TOFFOLI_OFFSET) 
    
a = [i for i in range(n)]
sqr=[n+i for i in range(4*n)]
ancilla = [5*n+i for i in range(nums1)]
gates_mod.Toffoli_qubits = [5*n+nums1+i for i in range(nums2)]
count=0
result, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:] = Inverison_Itoh_Tsujii_based(
            a, n, sqr[:n], sqr[n:2*n], sqr[2*n:3*n], sqr[3*n:4*n], count, ancilla[:KARATSUBA_ANC_BASE], ancilla[-KARATSUBA_ANC_BASE:]
        )

stats, full_depth, current_depth,toffoli_depth = get_exact_resources_optimized(gm, 5*n+nums1+nums2)
    
toffoli = stats['Toffoli_count']
cnot = stats['CNOT_count']

print("\n" + "="*40)
print(f"{'Quantum Resource Estimates':^40}")
print("="*40)
print(f"{'Metric':<20} | {'Value':>15}")
print("-" * 40)
print(f"{'Toffoli Gates':<20} | {toffoli:>15,}")
print(f"{'CNOT Gates':<20} | {cnot:>15,}")
print(f"{'Full Depth':<20} | {full_depth:>15}")
print(f"{'Current Depth':<20} | {current_depth:>15}")
print(f"{'Toffoli Depth':<20} | {toffoli_depth:>15}")
print(f"{'Qubits':<20} | {5*n+nums1+nums2:>15,}")
print("-" * 40)