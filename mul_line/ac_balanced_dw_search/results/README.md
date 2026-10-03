# Result-File Status

- `layers_n163.csv`, `layers_n233.csv`, and `layers_n283.csv` contain gate-level resource records covering controlled accumulation lanes 1–128 and uncontrolled reduction lanes
  1–64.
- `layers_n571.csv` currently covers accumulation lanes 1–128 and reduction lanes 1–64.
  It is sufficient for exact `w=1,...,128` estimates but is not the complete cache required
  for the extended `w=1,...,256` scan.
- `ac_balanced_full_shor.xlsx` aggregates those gate-level resource records into complete controlled
  point-addition maps for Shor's algorithm. All four sheets currently cover
  `w=1,...,128`. The `n=571` optimum lies on the old scan boundary and must not be treated
  as the final optimum over 1–256.

Rebuild the current workbook with:

```bash
python3 scripts/export_full_shor.py --max-w-571 128
```

After generating accumulation layers 129–256 and reduction layers 65–128, omit that option
to extend the `571` sheet to 256 rows. See Section 4.1 of the project README for the scan
commands.

The resource scope matches Table 6 of `main.pdf`: it includes the outer
compute/copy/uncompute but excludes QFT, measurement, and classical post-processing. The
project README and `BalancedResource.cpp` document the aggregation formulas.
