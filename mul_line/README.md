# AC-based quantum arithmetic and point addition

This directory contains the AC-based multiplication construction, its quantum
arithmetic circuits, and the three non-windowed AC point-addition configurations
in Table 6. The standalone Window implementation remains at
`../windowed_qrom_toptimal`.

| Configuration | Code | Workspace strategy |
|---|---|---|
| In-place+AC | `algorithm1_inplace_point_add/` | In-place point update with clean arithmetic |
| Balanced+AC | `ac_balanced_dw_search/` | Copy arithmetic outputs, uncompute, and reuse workspace |
| T-optimal+AC | `C++/`, with Python reference in `quantum_all/` | Retain nonlinear arithmetic intermediates |

## T-optimal retained implementation

From the repository root:

```bash
cmake -S mul_line/C++ -B mul_line/C++/build -DCMAKE_BUILD_TYPE=Release
cmake --build mul_line/C++/build -j 2
mul_line/C++/build/ShorEstimator --n 163
```

`--n` supports 163, 233, 283, and 571. The default number of scalar-controlled
inputs is `2n+2`; `--inputs 4` is a small schedule smoke check, not a Table 6 run.
The executable finds the data directory configured at build time. Override it
with `--data-root /path/to/mul_line/C++/data` if needed. Full runs can require
substantial memory. The existing schedule pairs the `2n+2` inputs in its first
layer, hence starts with `n+1` active additions.

## Balanced implementation

```bash
cmake -S mul_line/ac_balanced_dw_search -B mul_line/ac_balanced_dw_search/build -DCMAKE_BUILD_TYPE=Release
cmake --build mul_line/ac_balanced_dw_search/build -j 2
ctest --test-dir mul_line/ac_balanced_dw_search/build --output-on-failure
cd mul_line/ac_balanced_dw_search
python3 -m unittest discover -s tests -v
python3 scripts/export_full_shor.py --max-w-571 128
```

The bundled caches reproduce the four Balanced+AC rows at `w=56,59,64,128`.
The last command aggregates existing gate-level resource records and does not generate gates.
`--max-w-571 128` is required for the cache range that produced the paper row.

## In-place implementation

```bash
cd mul_line/algorithm1_inplace_point_add
python3 resources/point_addition_resources.py --sizes 163 --inversion optimal_depth
python3 tests/test_point_addition.py
python3 verification/verify_point_addition.py --sizes 163 233 283 571 --algebra-only
```

The deliverable is one controlled point-addition circuit, with its resource
estimator and verification code. The circuit implementation is
`algorithm1_inplace_point_add/circuit/point_addition.py` (`Point_addition`).
No complete Shor implementation or whole-stage export wrapper is required.
For comparison with Table 6, gate counts and depths are multiplied by `2n+2`,
while width remains unchanged.

## Construction and verification support

- `classcial_mul/`: classical AC multiplication for GF(2^163), GF(2^233),
  GF(2^283), and GF(2^571). The historical directory spelling is preserved.
- `Construction_Strategy.ipynb`, `choose_basis.ipynb`, `local_expansion.ipynb`,
  `matrix_full.ipynb`, and `Matrix_txt/`: multiplication/basis construction.
- `quantum_all/inv_verify/`, `quantum_all/mul.py`, and
  `quantum_all/verify_shor_logic.py`: Python arithmetic verification and resources.
- `C++/data/`: matrices and optimized linear gate streams for the C++ estimators.

Classical `mul_1.txt`/`mul_2.txt` expressions have representation-specific
coordinates; they are not automatically interchangeable with polynomial-basis
inputs. The quantum estimators consume their own matrix/gate-stream files.

## Shared resource backend

Gate-stream construction and resource analysis are maintained in the repository
root `circuit_backend/` directory. Local backend files are compatibility entry
points. Keep the repository directory structure when copying these modules.
See `circuit_backend/results/n163/` for the shared-backend reevaluation.
