# n=163 shared-backend reevaluation

Evaluation date: 2026-10-05. Backend: `nct-v2`.
Reference: Table 6, page 21 of the updated 24-page manuscript.

All five paper configurations were evaluated with the shared C++ or Python
backend. Every reported gate count, width, Toffoli depth, NCT depth, DW, and TDW
for these five n=163 rows is unchanged from Table 6.

| Configuration | Toffoli | CNOT | Width | Toffoli depth | NCT depth | DW |
|---|---:|---:|---:|---:|---:|---:|
| In-place + AC | 13,182,320 | 2,317,467,600 | 13,174 | 13,776 | 9,798,016 | 129,079,062,784 |
| Balanced + AC | 6,341,768 | 857,859,900 | 512,450 | 810 | 577,596 | 295,989,070,200 |
| Balanced + Karatsuba | 30,464,596 | 265,802,254 | 873,913 | 1,104 | 73,788 | 64,484,292,444 |
| T-optimal + AC | 3,167,396 | 571,478,060 | 1,643,885 | 322 | 233,752 | 384,261,406,520 |
| Window + AC | 1,122,328 | 167,316,138 | 198,025 | 608 | 264,098 | 52,298,006,450 |

## Evaluation scope

- **AC Balanced:** regenerated all 192 layer records needed for `w=1..128`,
  including 128 accumulation sizes and 64 reduction sizes. All six layer metrics
  match the corresponding saved pre-change records. The minimum DW over this
  interval remains at `w=56`. This evaluation does not extend the search interval.
- **AC in-place:** generated the complete single-point stream using
  `optimal_depth`. The point resources remain width 13,174, Toffoli count 40,190,
  CNOT count 7,065,450, Toffoli depth 42, and full depth 29,872. Gate counts and
  depths are multiplied by 328 for the table; width is unchanged.
- **Karatsuba Balanced:** reran the complete stage at `w=30`.
- **AC T-optimal:** reran the complete stage with the default 328 inputs.
- **Window:** regenerated and basis-verified six Balanced reduction layers,
  covering all 13 n=163 formula minimizers of TDW. The selected configuration
  remains `s=6, w=55`. This evaluates every formula minimizer, not the full gate-level
  resource grid. QROM costs retain the existing table-independent model.
- **JSB25:** reran both in-place and out-of-place complete-stage estimators. Both
  match their saved pre-change repository results.

## Files

- `resources.csv`: resource metrics for each evaluated configuration.
- `summary.json`: resource metrics, selected parameters, and shared-source SHA-256 hashes.
- `ac_balanced_layers.csv`, `ac_balanced_scan.csv`: fresh layers and all 128 scan rows.
- `ac_balanced_w56.json`: direct estimate of the paper's AC Balanced configuration.
- `ac_inplace_point.json`: single-point AC in-place resources.
- `karatsuba_w30.csv`, `retained_full.txt`: complete-stage C++ outputs.
- `jsb25_shor.json`: both JSB25 complete-stage outputs.
- `window_tdw/`: verified Window layers, candidate scan, and summary.

These are gate-level resource records and their configured composition results.
See `../../evaluate_n163.py` and `../../README.md` for reproduction commands.
