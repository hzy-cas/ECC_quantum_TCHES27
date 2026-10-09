# Shared NCT circuit backend

This directory is the single implementation location for the main circuits'
gate-stream construction, inverse replay, cancellation, and depth/count analysis.
Field arithmetic, point-addition schedules, register allocation, and whole-stage
resource composition remain in their respective circuit modules.

## Implementations and consumers

| Implementation | Consumers |
|---|---|
| `cpp/GateManager.h`, `cpp/GateManager.cpp` | AC Balanced, AC T-optimal, Karatsuba Balanced, and Window |
| `python/qasm.py`, `python/stats_utils.py` | AC in-place and JSB25 in-place/out-of-place |
| `python/matrix_cache.py` | AC linear-matrix and CNOT-stream loading |

The old `GateManager.*`, `qasm.py`, and `stats_utils.py` paths are compatibility
entry points. They contain includes/imports rather than additional copies of the
implementations. Existing module commands continue to work from their documented
directories. Keep `circuit_backend/` beside the main circuit directories when
copying the repository.

Python usage from the repository root:

```python
from circuit_backend.python.qasm import FlatGateManager
from circuit_backend.python.stats_utils import get_exact_resources_optimized

gm = FlatGateManager()
gm.add_Toffoli(0, 1, 2)
gm.add_Toffoli(1, 0, 2)
counts, full_depth, current_depth, toffoli_depth = get_exact_resources_optimized(gm, 3)
assert counts["Toffoli_count"] == 0
```

`get_stats()` returns `(Toffoli, CNOT, X)` without optimization.
`count_stream(width)` reports raw gate-stream counts and dependency depths.
`get_exact_resources_optimized()` cancels gates in place and reports optimized
counts and depths. `calculate_depth(width)` remains available as a full-depth
compatibility method. Both Python depth arrays and C++ counters use 64-bit storage.

For C++, include `cpp/GateManager.h` and compile `cpp/GateManager.cpp` once.
The existing module build files compile their forwarding source, which includes
this implementation once per executable. CMake tracks the included files, and
the Makefiles explicitly depend on the shared sources.

## Gate and depth semantics

- X and CNOT occupy one packed word; Toffoli occupies two.
- Toffoli controls are normalized because their order does not affect the gate.
  Python cancellation also recognizes exchanged controls in imported old streams.
- C++ keeps both Toffoli words in one chunk. Python uses absolute word addresses
  and supports Toffoli gates split across chunks.
- Cancellation only removes matching self-inverse gates when the participating
  wire stacks permit it. Intervening dependent operations block cancellation.
- `full_depth` is the dependency depth of the remaining NCT stream, with unit
  cost for each X, CNOT, or Toffoli gate.
- `toffoli_depth` gives Toffoli gates weight one and Clifford gates weight zero,
  while retaining the dependencies induced by the emitted gate order.
- `current_depth` reconstructs layers grouped by Toffoli depth. It is a separate
  metric from `full_depth`.
- Register width is supplied by the circuit builder. QROM formulas and
  whole-stage composition are outside the backend.

## Regression checks

```bash
python3 -m unittest discover -s circuit_backend/tests -v
```

The tests include exchanged-control cancellation, both cross-chunk cancellation
cases, inverse replay, exhaustive small truth tables, compatibility-import
identity, and 303 common C++/Python gate sequences at three C++ chunk capacities.
A C++17 compiler is needed for the cross-language checks.

AC Balanced generated caches now use `ac-balanced-clean-nct-v2`. Its estimator
rejects incompatible caches. The export-only script can read explicitly selected
historical v1 records or new records, but rejects files mixing both versions.
Window fingerprints include the shared C++ source and header; changes invalidate
old build stamps and layer identities.

## n=163 reevaluation

The saved results in `results/n163/` cover all five paper configurations and both
JSB25 modes. See `results/n163/README.md` for the evaluation scope and resource records.

To reproduce the evaluation into a new directory:

```bash
python3 circuit_backend/evaluate_n163.py --output /tmp/ecc-n163-new
```

This builds the C++ estimators, recomputes all 192 AC Balanced layers needed for
`w=1..128`, evaluates AC in-place and JSB25, evaluates Karatsuba at `w=30`, runs
the AC T-optimal stage, and verifies the n=163 Window TDW candidates. It writes
`resources.csv` and `summary.json`. Existing output directories containing files
are rejected to avoid mixing evaluations.

To aggregate an existing complete evaluation without regenerating circuits:

```bash
python3 circuit_backend/evaluate_n163.py --output /tmp/ecc-n163-new --summarize-only
```
