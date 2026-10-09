# Resource model implemented by this artifact

Let `m=n+1`, `ell=ceil(m/s)`, and `L_s=2 ell`.  Separate `P` and `Q` windows
are ordered as

```text
P_0, Q_0, P_1, Q_1, ..., P_(ell-1), Q_(ell-1).
```

The top window has its true width

```text
b_j = min(s, m-sj),    K_j = 2^b_j.
```

Each shifted table contains

```text
T_tilde^S_j[u] = [u 2^(sj)]S + Delta^S_j.
```

The fixed offsets are derived by domain-separated SHAKE-256 and rejected if
any represented table address would yield the point at infinity.  The circuit
therefore computes

```text
F'(k,l) = [k]P + [l]Q + Delta_total.
```

No final subtraction is performed.  Translation by a fixed group element is
a bijection, so `F'` has exactly the same collision relation and hidden
subgroup as the unshifted oracle.

The default reproducibility instance uses `Q=[7]P`.  Alternatively, the
manifest and table commands accept the public affine coordinates of any
finite point in the declared prime-order subgroup.  In that mode the Q-stream
offset is generated as a known multiple of public `Q`, so neither table
generation nor infinity rejection uses the discrete logarithm of `Q` relative
to `P`.  The manifest records the coordinates and explicitly records the
relative scalar as unknown.

## Coherent QROM formulas

For one `K=2^b` table, let `h_u` be the Hamming weight of the `u`-th `2n`-bit
affine word.  The coherent unary-iteration lookup uses

```text
T_QROM = 2(K-2),
C_QROM = (K-2) + sum_u (3h_u-2),
D_QROM = 3K-4 + sum_u (2 ceil(log2 h_u)+1).
```

The paper scan sets `h_u=2n` for every row.  This upper-bounds the CNOT count
and per-line ASAP depth of the implemented unary-iteration and balanced-
fanout QROM module.  It does not change selector Toffoli count or allocated
QROM width.

After internal self-inverse cancellation, the lookup, clean `z,v`
precomputation, and literal unlookup pair uses

```text
T_pair = 4K-2b-6,
C_pair = 2 C_QROM + 4n - 2(h_(K-1)-1),
D_pair = 2 D_QROM - 2b - 2 ceil(log2 h_(K-1)) + 3.
```

## Window schedule

The first `p` QROM outputs initialize the accumulators.  Later Phase-I batches
use coherent lookup followed by uncontrolled quantum--quantum point addition.
Phase II contains only uncontrolled quantum--quantum reduction layers.  The
legacy controlled constant-point layer is not used by the windowed scan.

The forward network contains

```text
accumulation layers = ceil(L_s/p)-1,
tree layers         = ceil(log2 p),
point additions     = L_s-1.
```

## Early unlookup and width

For every post-initialization window addition, the gate order is

```text
U_T -> clean z,v -> literal U_T^dagger -> Montgomery inversion -> PA tail.
```

The tail evaluates

```text
lambda = v z^(-1),
x3 = lambda^2 + lambda + z + a,
y3 = lambda (x_R+x3) + x3 + y_R.
```

It never reads `x_T,y_T` after unlookup.  The `z,v` registers coexist with the
QROM cache and are counted explicitly.  After the cache is cleared, its
scratch lines may be reused by the Montgomery and PA workspaces.  The
module-level peak therefore uses the maximum of those two scratch demands,
not their sum.

## Scope of reported results

The analytical scan composes concrete AC reduction-layer records with the
closed QROM formulas.  It adds module depths layer by layer and applies the
documented liveness recurrence for width.  It cannot capture cancellation or
parallel execution across module boundaries.  The outer accounting follows
the existing Table-6 convention: forward compute, output copy, and reverse
compute.  QFT, final measurement, and classical post-processing are excluded.

The bundled C++ backend provides the separate gate-stream implementation used
to audit coherent lookup, literal inverse, early unlookup, and the AC point-
addition tail.  A complete monolithic whole-network synthesis remains a
different and substantially more expensive calculation than the analytical
scan reported in `../../results_and_tables/window.pdf`.
