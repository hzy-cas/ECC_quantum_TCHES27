# Window-QROM resource-estimation artifact

This directory contains the code needed to reproduce and audit the coherent
NCT Window-QROM architecture described in
`../../results_and_tables/window.pdf`.

The artifact provides two complementary paths:

1. `analytic-scan` reproduces the reported module-composed resource table.  It
   combines the bundled gate-level AC reduction-layer records with the closed
   QROM formulas and sets every QROM word weight to `h_u=2n`.
2. The table generator, Python QROM simulator, and bundled C++ AC backend
   implement coherent `U_T`, its literal inverse on the same wires, clean
   `z=x_R+x_T` and `v=y_R+y_T`, early unlookup, Montgomery inversion, and the
   quantum--quantum point-addition tail.

Measurement-assisted uncomputation is outside the scope of this artifact.

## Reproduce the resource scan

Run the following commands from the repository's
`rebuttal/codes_and_circuits/` directory:

```bash
python3 -m windowed_qrom_nct.cli analytic-scan \
  --output-dir /tmp/windowed_qrom_scan
```

The command scans every `s=2,...,18` and every `p=1,...,L_s` point for
`n=163,233,283,571`.  It emits:

- `windowed_theory_upper.csv`, containing all 6,304 candidate rows and an
  explicit availability flag;
- `windowed_theory_upper_minima.csv`, containing the minimum Toffoli count,
  TDW, and DW among available points;
- `summary.json`, containing coverage and missing-lane records.

The bundled caches cover all 832, 1,178, and 1,428 configurations for
`n=163,233,283`, respectively.  For `n=571`, they cover 1,734 of the 2,866
configurations.  Accordingly, the `n=571` rows are the best evaluated points,
not certified global optima.

The scan emits no quantum gates.  QROM CNOT count and per-line ASAP depth use
the table-independent upper-weight model `h_u=2n`.  The AC point-addition
modules are gate-level records from the uncontrolled Phase-II reduction
implementation.  Network depth is aggregated with explicit module barriers
and width with the documented liveness model.  These results are analytical
module-composed estimates, not monolithic whole-network gate-stream results.

## Build the coherent AC backend

```bash
cmake -S windowed_qrom_nct/ac_backend \
      -B windowed_qrom_nct/ac_backend/build \
      -DCMAKE_BUILD_TYPE=Release
cmake --build windowed_qrom_nct/ac_backend/build -j
ctest --test-dir windowed_qrom_nct/ac_backend/build --output-on-failure
```

The executable is
`windowed_qrom_nct/ac_backend/build/ac_balanced_estimator`.  Its Window-QROM
commands are `window-init`, `qrom-pair`, `early-tail`, and `window-layer`.

## Frozen shifted tables

The standardized generator is `P`, and the default reproducible test point is
`Q=[7]P`.  This is a gate-generation instance, not an ECDLP challenge.  The
public seed is derived from
`Binary_ECC/windowed-QROM/shifted-tables/v1` and has SHA-256 digest

```text
de0faddf37eaf09c24fe74b8dc130a5f75e73234f62f1ca55dc66b7607a67340
```

Generate a descriptor-only manifest quickly:

```bash
python3 -m windowed_qrom_nct.cli manifest \
  --n 163 --s 10 --descriptors-only --output /tmp/n163_s10.json
python3 -m windowed_qrom_nct.cli verify-manifest \
  --manifest /tmp/n163_s10.json
```

Omitting `--descriptors-only` also computes every concrete table hash and
Hamming weight and therefore takes substantially longer for large windows.

An arbitrary public subgroup point can be supplied by affine coordinates;
the generator validates curve and subgroup membership and never requires or
records its discrete logarithm relative to `P`:

```bash
python3 -m windowed_qrom_nct.cli manifest \
  --n 163 --s 10 --q-x 0xPUBLIC_X --q-y 0xPUBLIC_Y \
  --descriptors-only --output /tmp/n163_s10_custom_q.json
```

The same `--q-x` and `--q-y` options are accepted by the `table` command.
For a user-supplied `Q`, each Q-stream offset and window step are represented
as known scalar multiples of the public point itself.  This permits exact
infinity rejection without knowing a scalar relating `Q` to `P`.

Generate one AC-L-basis table and synthesize its coherent lookup pair:

```bash
python3 -m windowed_qrom_nct.cli table \
  --n 163 --s 2 --stream P --j 0 --basis l --output /tmp/qrom.bin

windowed_qrom_nct/ac_backend/build/ac_balanced_estimator qrom-pair \
  --n 163 --tables 2:/tmp/qrom.bin --data-root ../../mul_line/C++/data
```

## Validation

```bash
python3 -m unittest discover -s windowed_qrom_nct/tests -v
```

The tests cover the standardized curves, deterministic offsets, short top
windows, basis conversion, collision-preserving translation metadata, literal
QROM inversion, balanced fanout counts, clean cache removal after forming
`z,v`, fail-closed search coverage, and exact reproduction of the reported
upper-model minima.
