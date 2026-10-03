# Validation record

This document translates the validation history shipped with the package. The
historical gate-level checks below were run at **n=163 only**. Predictions and
planning checks at n=233,283,571 are not complete gate-level resource results.

## Dual-objective version, 2026-09-30

- Release CMake build and CTest passed; all 31 Python tests passed.
- Balanced tail/reduction layers were checked at k=1,2,3,4,5,7,8,9, including odd
  sizes and powers-of-two boundaries. Checks covered ordinary point-add outputs,
  zeroed reusable workspace, and outer compute-copy-uncompute.
- All n=163 TDW winners' reduction layers k=1,2,3,7,14,27 were generated and
  verified. The 13 tied candidates have s=6..18 and w=55.
- With smoke caches, complete resources cover 71/5,576 rows; all formula minima
  are available, but the full resource grid is incomplete.
- Maximum observed RSS for these new winner jobs was about 1.68 GiB at k=27.
- All 25 retained winner layer caches were reused without changing their JSON;
  retained smoke verification also passed. Cross-objective output mixing was
  rejected before writing, and resumed runs did not rewrite successful layers.
- Across all four fields, the TDW plan has 53 candidates and 25 distinct layers
  (largest scheduling estimate 191 GiB); TD has 971 candidates and 74 layers
  (largest estimate 193 GiB). These are not measured memory upper bounds.
- A synthetic-record test checks the full n=163 grid's composer/minimizer logic.
  Synthetic records are confined to tests and are never exported as
  gate-level resource records.
- QROM and lane-first scheduling regressions passed. A read-only comparison of
  256 historical Balanced reduction records checked Toffoli counts/depths and
  allocated widths; those records were not imported as gate-level resource records
  from the new backend.

| n=163 metric | Retained, minimum TD | Balanced, minimum TDW |
|---|---:|---:|
| (s,w) | (3,110) | (6,55) |
| Toffoli | 1,066,252 | 1,122,328 |
| CNOT | 201,599,530 | 167,316,138 |
| Width | 611,406 | 198,025 |
| Toffoli depth | 248 | 608 |
| NCT depth with module barriers | 164,598 | 264,098 |
| DW | 100,636,204,788 | 52,298,006,450 |
| TDW | 151,628,688 | 120,399,200 |

Balanced exchanges more Toffoli depth for less width, reducing TDW by about 20.6%
in this comparison. The estimates use `h_u=2n` QROM costs and module barriers;
they are not a monolithic whole-network optimum.

## Other fields: formula candidates

| n | (s,w) | Predicted Toffoli | Predicted TD | Predicted width | Predicted TDW |
|---:|---:|---:|---:|---:|---:|
| 233 | (5,94) | 2,768,248 | 560 | 501,605 | 280,898,800 |
| 283 | (6,95) | 3,532,776 | 716 | 624,661 | 447,257,276 |
| 571 | (6,191) | 14,751,392 | 848 | 2,658,604 | 2,254,496,192 |

## Historical lane-first v2 checks, 2026-09-27

- Release build, CTest, and 21 Python tests passed, including real C++ simulation,
  two-process execution, resume, compatible layer reuse, new network aggregation,
  and unknown-cache rejection.
- Integer schedule tests at all four sizes checked that every original addend is
  used once, at index `(q*s+t)*p+j`, with `ceil((2n+2)/(p*s))-1` accumulation
  batches and a search including p=2n+2.
- At n=233,s=4,p=117 the schedule has 117 four-bit tables; p=118 has 114 four-bit
  and four three-bit tables. The old separately windowed scalars gave 116
  four-bit and two two-bit tables and are not the current schedule.
- The four-field winners plan required 74 layer types and 971 candidates, with a
  largest scheduling estimate of 193 GiB. CMake 3.16-compatible build and foreign
  build-preservation tests passed.
- All 84 bundled data hashes passed, with no data symlinks.
- n=163 reduction and tail layers at k=1,2,3 passed basis-state tests. All 25
  retained winner reduction layers passed (largest k=64), including 18 previously
  unverified cached layers.
- Tests checked outputs, reusable workspace cleanup, input restoration, and
  copied-output preservation. They are deterministic sample checks, not a formal
  proof for all inputs.
- The existing 164 n=163 layer records remained usable and covered 2,765/5,576
  lane-first grid points, including all minimum-depth candidates. Full-grid
  completion was not claimed.
- Peak RSS in the supplemental candidate verification was about 2.33 GiB at
  k=64; historical complete n=163 layer records reached about 2.98 GiB.

The retained n=163 optimum over s=2..18,w=1..328 has 304 tied parameter pairs
(s=3..18,w=110..128). Tie-breaking selects (3,110), with the resources shown
above. Its agreement with an earlier representative row does not make other
field sizes exempt from recomputation. Width is the documented allocation peak,
not a globally width-optimal result.
