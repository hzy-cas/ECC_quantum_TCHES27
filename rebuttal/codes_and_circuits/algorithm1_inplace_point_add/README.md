# In-Place Point Addition

This pure-Python implementation reproduces the [JSB25] in-place point-addition structure
with AC-based multiplication. It provides two inversion backends:

- `binary_ecc`: the low-width [JSB25] inversion chain;
- `optimal_depth`: the parallel minimum-Toffoli-depth inversion used in our work.

It supports `n=163,233,283,571` and does not require ProjectQ. The code and optimized
squaring streams are contained here; multiplication matrices and basis transforms are
read from the repository's `mul_line/C++/data` directory, which is located automatically.

## Run

From `rebuttal/codes_and_circuits/algorithm1_inplace_point_add`:

```bash
# Fast algebraic verification for every supported field size
python3 verification/verify_point_addition.py \
  --sizes 163 233 283 571 --algebra-only

# Complete n=163 gate-stream verification with both inversion backends
python3 verification/verify_point_addition.py --sizes 163 --inversion both

# Regression tests
python3 tests/test_point_addition.py

# Single-point-addition resource estimates
python3 resources/point_addition_resources.py --sizes 163 --inversion both
```
