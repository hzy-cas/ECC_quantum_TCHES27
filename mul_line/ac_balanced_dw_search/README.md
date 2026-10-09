# AC-Based Balanced Resource Estimator

This project combines AC-based multiplication and Itoh--Tsujii inversion with the
Balanced point-addition schedule. It scans the parallelism parameter `w` for the minimum
`DW` and `TDW`, where:

- `DW = full NCT depth × width`;
- `TDW = Toffoli depth × width`.

Supported field sizes are `n=163,233,283,571`.

## Data and build

The source code is contained here. Arithmetic matrices are read from the repository's
`mul_line/C++/data` directory, which is located automatically. Use
`--data-root PATH` only to override that location.

```bash
cd mul_line/ac_balanced_dw_search
make
make test
```

Inspect individual circuit components:

```bash
build/ac_balanced_estimator primitive --n 163 --kind multiplication
build/ac_balanced_estimator layer --n 163 --lanes 1 --mode accumulation
```

## Search and export

Scan an interval of `w` values:

```bash
python3 scripts/run_experiment.py --n 163 --w-start 1 --w-end 32
python3 scripts/run_experiment.py --n 233 --w-start 1 --w-end 16 --jobs 2
```

Useful options are `--g` for the stride, `--jobs` for independent workers, `--force` for
recomputation, and `--no-build` when the executable is already built. 

Outputs are written to:

```text
results/ac_balanced_n<N>.csv
results/ac_balanced_n<N>.xlsx
cache/layers_n<N>.csv
```

Rebuild an Excel file without running circuits:

```bash
python3 scripts/run_experiment.py --n 163 --export-only
```

Export complete-Shor estimates from existing layer caches:

```bash
python3 scripts/export_full_shor.py
```

The included `n=571` cache currently covers `w<=128`; use
`python3 scripts/export_full_shor.py --max-w-571 128` until layers for `w=129,...,256`
have been generated. The combined workbook is `results/ac_balanced_full_shor.xlsx`.

## Shared resource backend

Gate-stream construction and resource analysis are maintained in the repository
root `circuit_backend/` directory. Local backend files are compatibility entry
points. Keep the repository directory structure when copying these modules.
See `circuit_backend/results/n163/` for the shared-backend reevaluation.
