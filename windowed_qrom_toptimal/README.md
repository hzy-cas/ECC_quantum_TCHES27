# Window-QROM resource estimation

This package supports GF(2^163), GF(2^233), GF(2^283), and GF(2^571).
It includes its arithmetic C++17 code, 84 data files, basis transforms, and
Python tests. Gate storage and resource analysis use `../circuit_backend/cpp`. Runtime dependencies are Python 3.10+, CMake 3.16+, and a C++17
compiler on Linux or macOS. No external repository, SageMath, OpenSSL, network
access, or third-party Python package is needed.

**The updated paper's Table 6 Window+AC configuration uses `--objective tdw`.**
The historical directory/executable names do not select that objective.

## Objectives and resource scope

- `--objective td` (the default): retained/T-optimal arithmetic keeps nonlinear
  intermediate results and minimizes Toffoli depth, breaking ties by width.
- `--objective tdw`: Balanced arithmetic uses compute-copy-uncompute to reclaim
  arithmetic workspace and minimizes Toffoli depth times width over the full
  specified search range. It does not restrict the search to minimum-depth points.

Both modes use the same lane-first schedule, unary-iteration QROM, early unlookup,
AC data, and outer compute-copy-uncompute. The Balanced and retained families use
different workspace lifetimes; optimality is restricted to the selected family,
allocation policy, module boundaries, and search range.

For N=2n+2 and p=w, lane j receives `S_j, S_(p+j), S_(2p+j), ...`, then groups up
to s entries into each window. Entry t of window q has index `(q*s+t)*p+j`.
There is no P/Q boundary in this lane-first grouping. The obsolete scheme that
first windowed the two scalars separately is not an option in this package.

Arithmetic layers use real optimized NCT streams with per-wire ASAP scheduling.
QROM costs use the conservative table-independent `h_u=2n` model. Network depth
is composed across explicit module barriers, not optimized as one monolithic
stream. Width includes the address, copied final output, retained intermediate
states, and the peak reusable workspace. TD mode retains squaring workspace;
TDW mode also reuses cleaned nonlinear and squaring workspace but retains copied
inverse/slope/coordinate outputs and old point states. Both include the outer
compute, output copy, and uncompute, excluding QFT and classical post-processing.
The circuits implement the ordinary affine-addition branch, with the existing
shifted-table/exceptional-input assumptions; they do not add a universal circuit
for infinity or equal-point cases. See [METHOD.md](METHOD.md) for exact lifetimes.

## Build and verify data

Copy the repository, or this directory together with the sibling
`circuit_backend/` directory, excluding `build/`, and rebuild:

```bash
cd windowed_qrom_toptimal
python3 run.py verify-data
python3 run.py build --jobs 2
python3 -m unittest discover -s tests -v
```

`build --jobs` controls compilation. `run --jobs` controls independent estimator
processes. An incompatible build is preserved and rebuilt automatically when a
run needs new layers. The program does not install system packages.

## Reproduce the paper objective

```bash
python3 run.py plan --n 163 233 283 571 --objective tdw --scope winners
python3 run.py run --n 163 233 283 571 --objective tdw --scope winners --jobs 4 --memory-budget-gb 200 --output results_tdw
```

The default grid is `s=2..18`, `w=1..2n+2`. It has 53 tied TDW candidates requiring
25 distinct Balanced reduction layers. The largest job's scheduling estimate is
191 GiB, not a measured upper bound. Each generated layer must pass Toffoli-count,
Toffoli-depth, and register-layout checks before complete resource rows are emitted.
Add `--verify` for deterministic basis-state functional and cleanup checks.

The representative paper parameters are:

| n | s | w | Toffoli | Width | Toffoli depth |
|---:|---:|---:|---:|---:|---:|
| 163 | 6 | 55 | 1,122,328 | 198,025 | 608 |
| 233 | 5 | 94 | 2,768,248 | 501,605 | 560 |
| 283 | 6 | 95 | 3,532,776 | 624,661 | 716 |
| 571 | 6 | 191 | 14,751,392 | 2,658,604 | 848 |

For a small local check:

```bash
python3 run.py run --n 163 --objective tdw --scope smoke --jobs 1 --memory-budget-gb 6 --verify --output results_tdw
python3 run.py run --n 163 --objective tdw --scope winners --jobs 1 --memory-budget-gb 6 --verify --output results_tdw
```

`smoke` uses k=1,2,3 and does not establish a grid optimum. For paper candidates,
check `formula_minimizers_all_available=true`. That flag covers all formula
minimizers of the selected objective, not every grid point.

## Retained mode and full scans

The retained comparison uses `--objective td` and a separate `results/` directory:

```bash
python3 run.py plan --n 163 233 283 571 --objective td --scope winners
python3 run.py run --n 163 233 283 571 --objective td --scope winners --jobs 4 --memory-budget-gb 200 --output results
```

Its default winners need 74 distinct layers (971 tied candidates), with a largest
scheduling estimate of 193 GiB. These are different results from the paper's
Balanced Window row. Larger job counts do not change the circuit parallelism w.

`--scope full` requires 2,089 layers and evaluates 42,636 configurations. The
largest n=571 estimates are 573 GiB for TD and 1,145 GiB for TDW. A 200 GiB or
448 GiB budget therefore rejects those full scans; winners do not require them.
Use `plan` before allocating a full scan. For n=163 the retained full-scan estimate
is 14.12 GiB. All memory estimates are scheduling heuristics; observe actual RSS.
Linux additionally supports `--job-memory-gb` (RLIMIT_AS, virtual address space,
not RSS). `--timeout` limits each C++ job in seconds.

## Resume, outputs, and provenance

Repeat the same run command to resume. Successful layer JSON records are written
atomically and retained after interruptions. Failed jobs return a nonzero status.
Use separate directories for the two objectives; mixed caches are rejected.
One coordinator may run per output directory, with internal `--jobs` concurrency.
A stale `.lock` may be removed only after confirming no task owns it.

- `results_tdw/layers/*.json`: gate-level resource records, register layout, RSS, and elapsed time.
- `results_tdw/scan.csv`: grid rows; missing-layer rows have `available=False` and
  blank gate-level resource fields. `predicted_*` columns are formulas.
- `results_tdw/summary.json`: selected optimum, objective, and resource coverage.
- `results_tdw/run_progress.json` and `logs/`: job progress and failures.

To compose existing records without generating gates:

```bash
python3 run.py scan --n 163 --objective tdw --output results_tdw
python3 run.py scan --n 163 --objective td --output results
```

TD uses `window-retained-lane-first-v2`; TDW uses
`window-balanced-lane-first-tdw-v1`. Read `selected_optimum` for the selected
objective. `minimum_toffoli_depth_over_available` always means the lowest TD
among available rows, including in TDW mode; it is not the TDW winner.
Network and arithmetic-layer fingerprints are separate. Compatible retained
layers can be reused without rewriting their JSON, while network summaries are
recomputed. Balanced records cannot be fabricated by doubling retained numbers
or editing fingerprints. Rebuild after arithmetic/data changes and use a new
output directory.

`backend/` contains the C++ implementations; `model.py` composes schedules and
resources; `run.py` manages builds and jobs; `support/` provides local field,
basis, QROM, and fixture utilities. `data_manifest.json` checks all 84 data files.

The shared C++ backend is included in the arithmetic-layer fingerprint. Use a
new output directory when changing backend versions. The n=163 shared-backend
reevaluation is stored in `../circuit_backend/results/n163/`.
