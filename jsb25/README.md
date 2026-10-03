# Binary ECC Full-Stream Point Addition

This pure-Python implementation emits complete NCT gate streams for the [JSB25] in-place and out-of-place point-addition circuits. It supports
`n=163,233,283,571` and requires neither ProjectQ nor external matrix files.

## Run

From `jsb25`:

```bash
# Functional, ancilla-clearing, and compute-copy-uncompute verification
python3 verify_point_addition.py --n 163 --mode both

# Single-point-addition resources
python3 point_addition_resources.py --n 163 --mode both

# Complete controlled-point-addition stage in Shor's algorithm
python3 shor_table_resources.py --n 163 --mode both

# Default regression suite
python3 -m unittest discover -s tests -v
```

Valid modes are `inplace`, `outofplace`, and `both`. 

## Reported depths

The estimator reports gate counts, width, strict NCT ASAP depth (`full_depth`), depth after
reconstruction by Toffoli layer (`current_depth`), and ordered Toffoli depth
(`toffoli_depth`). Gates assigned to the same Toffoli layer are processed together for
automatic depth estimation.

For complete-Shor estimates, the number of controlled point additions is `2n+2`. FLT-out
uses compute-copy-uncompute to clear retained garbage; the script applies this composition
and the final output-copy CNOT layer automatically.
