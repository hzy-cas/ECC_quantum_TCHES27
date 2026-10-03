# Third-party notices

## Binary_ECC

- Project: `starj1023/Binary_ECC`
- Source: <https://github.com/starj1023/Binary_ECC>
- Audited point-addition revision:
  `e844e669c1f0cb8ccfb54144b9c911ba964f405d`
- License: Eclipse Public License 2.0; see `LICENSE-EPL-2.0.txt`.

The point-addition schedules, Itoh--Tsujii chains, reduction order, and
literal recursive-Karatsuba gate order are Python reimplementations of the
corresponding ProjectQ source.  ProjectQ itself is neither included nor used.

## Local Python backend lineage

The packed `FlatGateManager`, `replay_reverse()` construction pattern, global
arithmetic backend, and three-stage resource interface follow the local
`recursive_karatsuba` Python framework.  The implementation in this directory
is self-contained and does not import that sibling project at runtime.

