# Balanced point addition with Karatsuba arithmetic

`C++/QuantumShor.cpp` is the **Balanced** implementation used for the
Balanced+Kara. configuration in Table 6. It already implements arithmetic
compute-copy-uncompute: Montgomery inversion, the slope product, and the coordinate
product are copied into clean output registers before their workspaces are
uncomputed and reused. A separate retained Karatsuba implementation is not needed
for the five paper configurations.


## Build and run

From the repository root:

```bash
cmake -S recursive_karatsuba/C++ -B recursive_karatsuba/C++/build -DCMAKE_BUILD_TYPE=Release
cmake --build recursive_karatsuba/C++/build -j 2
recursive_karatsuba/C++/build/karatsuba_balanced_estimator --n 163 --w 30 --output /tmp/karatsuba_n163.csv
```

The driver is sequential to bound memory use and needs no OpenMP runtime.
For a scan, replace `--w 30` with `--w-start 1 --w-end 64`. Bounds are inclusive,
and `1 <= w <= n+1` is enforced. Use `--data-root PATH` to point to a directory
containing `square_163/`, `square_233/`, `square_283/`, and `square_571/`.

The saved paper candidate widths are:

| n | w |
|---:|---:|
| 163 | 30 |
| 233 | 32 |
| 283 | 30 |
| 571 | 64 |

## Arithmetic support

- `inv/`: Karatsuba-based inversion estimation and verification.
- `square/`: optimized squaring data and construction support.
- `shor-circuit/`: Python point-addition reference and verification programs.
- `C++/square_*/`: squaring streams used by the C++ Balanced circuit.
