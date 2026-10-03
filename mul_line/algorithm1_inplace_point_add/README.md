# In-Place Point Addition

This pure-Python implementation reproduces the [JSB25] in-place point-addition structure
with AC-based multiplication. It provides two inversion backends:

- `binary_ecc`: the low-width [JSB25] inversion chain;
- `optimal_depth`: the parallel minimum-Toffoli-depth inversion used in our work.

It supports `n=163,233,283,571` and does not require ProjectQ. The code and optimized
squaring streams are contained here; multiplication matrices and basis transforms are
read from the repository's `mul_line/C++/data` directory, which is located automatically.

## Scope and circuit entry point

This package provides one controlled point-addition circuit, its resource
estimator, and verification code. The circuit entry point is `Point_addition`
in `circuit/point_addition.py`, implementing
`|q>|P>|0> -> |q>|P+qP2>|0>` on the supported affine-addition branch.
The paper configuration selects `--inversion optimal_depth`. A complete Shor
implementation or whole-stage export wrapper is outside this package's scope.

## Run

From `mul_line/algorithm1_inplace_point_add`:

```bash
# Fast algebraic verification for every supported field size
python3 verification/verify_point_addition.py \
  --sizes 163 233 283 571 --algebra-only

# Complete n=163 gate-stream verification with both inversion backends
python3 verification/verify_point_addition.py --sizes 163 --inversion both

# Regression tests
python3 tests/test_point_addition.py

# Single-point-addition resources for the paper configuration
python3 resources/point_addition_resources.py --sizes 163 --inversion optimal_depth
```
