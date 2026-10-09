# Quantum Arithmetic over Binary Fields

Code and data for the manuscript *Quantum Arithmetic over Binary Fields with Its
application*. 

## Table 6 implementations

| Paper configuration | Implementation | Entry point |
|---|---|---|
| In-place+AC | [mul_line/algorithm1_inplace_point_add](mul_line/algorithm1_inplace_point_add) | `resources/point_addition_resources.py --inversion optimal_depth` |
| Balanced+Kara. | [recursive_karatsuba/C++](recursive_karatsuba/C++) | `karatsuba_balanced_estimator` |
| Balanced+AC | [mul_line/ac_balanced_dw_search](mul_line/ac_balanced_dw_search) | `ac_balanced_estimator` and `scripts/export_full_shor.py` |
| T-optimal+AC | [mul_line/C++](mul_line/C++) | `ShorEstimator` |
| Window+AC | [windowed_qrom_toptimal](windowed_qrom_toptimal) | `run.py --objective tdw` on the `plan`, `run`, or `scan` subcommand |
| JSB25 in-place/out-of-place baselines | [jsb25](jsb25) | `shor_table_resources.py` |


Balanced arithmetic uses compute-copy-uncompute to release and reuse arithmetic
workspace. T-optimal arithmetic retains nonlinear intermediate results; linear
cleanup remains. Both complete-stage estimators account for the outer
compute-copy-uncompute composition. QFT, measurement, and classical post-processing
are outside the reported resource scope.

## Shared gate-stream backend

The five main configurations and the JSB25 comparison circuits use the shared
C++ or Python implementations in [circuit_backend](circuit_backend). Their
original backend paths are compatibility entry points. See that directory for
regression tests and the n=163 resource reevaluation.

## Arithmetic and synthesis support

- `mul_line/classcial_mul`, the construction notebooks, and `mul_line/Matrix_txt`
  contain the AC multiplication construction and classical supporting data.
- `mul_line/quantum_all` contains the Python arithmetic/retained-circuit reference
  and verification code. `mul_line/C++/data` is the canonical data location used
  by the new AC Balanced and in-place entry points.
- `recursive_karatsuba/inv`, `square`, and `shor-circuit` contain the arithmetic
  reference and verification code used with Karatsuba multiplication.
- `DaCSynth` and `move_CNOT_depth` implement the linear-layer synthesis/optimization
  used in the arithmetic circuits and remain part of the paper artifact.

The Python estimators require Python 3.10 or newer. C++ builds require CMake 3.16
or newer and a C++17 compiler. AC Balanced workbook export additionally requires
`openpyxl` (`mul_line/ac_balanced_dw_search/requirements.txt`). The classical
construction notebooks use SageMath; it is not needed for the C++ or pure-Python
point-addition entry points.