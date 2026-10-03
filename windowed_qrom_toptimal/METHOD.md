# Window-QROM resource models: retained TD and Balanced TDW

The default objective `td` preserves the retained model below. Objective
`tdw` switches ONLY the point-add arithmetic/liveness policy to Balanced,
and minimizes the full-network Toffoli-depth times allocated peak width.
Both objectives share the exact same lane-first schedule and QROM formulas.

## Lane-first partition (matches the manuscript)

Let N=2n+2 and p=w. Enumerate the original constant addends as
S_0,...,S_(N-1), with the P powers followed by the Q powers. FIRST assign
S_j,S_(p+j),S_(2p+j),... to lane j, for j=0,...,p-1. Its length is
M_j=floor((N-1-j)/p)+1. THEN divide that lane into L_j=ceil(M_j/s) windows.
Window q has b_(j,q)=min(s,M_j-q*s) bits and table

T_(j,q)[u] = sum_(t=0)^(b_(j,q)-1) u_t S_((q*s+t)*p+j).

Thus one table can mix P and Q addends. There is no special split at n+1.
As in the manuscript, shift each table by a constant offset to avoid infinity;
the table-independent h_u=2n bound does not depend on the selected offsets.
The fixed total output translation does not change the collision relation.

Initialize each lane with its q=0 table. At round q>=1 process only lanes
with q<L_j, preserving their identities. For this strided partition they
form a prefix. Each round uses QROM, clean z,v, literal unlookup and the
retained point-add tail. Finally reduce p accumulators pairwise, carrying
an odd item forward. Phase I has ceil(N/(p*s))-1 point-add rounds; Phase II
has ceil(log2 p) rounds. The total number of tables is sum_j L_j, NOT
2 ceil((n+1)/s). The forward point-add count is sum_j L_j-1.

Search every s=2,...,18 and p=1,...,N. Here s is a maximum window size, so
many different caps can describe the same effective circuit. All tied
depth candidates are retained; their distinct arithmetic layers are deduplicated.
Only the network model was changed to window-retained-lane-first-v2.
The arithmetic module schema remains window-retained-clean-zv-v1.

## Retained arithmetic

For k simultaneous point additions, M_n is the AC multiplication count and
r_n is the number of forward inversion multiplications:

| n | M_n | r_n | inversion Toffoli depth d_n |
|---:|---:|---:|---:|
| 163 | 906 | 9 | 8 |
| 233 | 1341 | 10 | 8 |
| 283 | 1668 | 11 | 9 |
| 571 | 3569 | 13 | 10 |

Montgomery uses k-1 product-tree multiplications, one inversion, and 2(k-1)
recovery multiplications. Each lane adds one multiplication for lambda and
one for y3. Retain their output blocks; restore multiplication input linear
transforms and uncopy recovery fanout, as in the original T-optimal code.

T_PA(n,k) = (r_n + 5k - 3) M_n.
D_PA(n,k) = d_n + 2 ceil(log2 k) + 2.

These are forward module costs. The layer generator verifies both formulas
against its actual optimized gate stream and fails if either differs.
The formula-only grid is a screening step; it does not fabricate unavailable
gate-level counts, full depths, or widths.

## Balanced arithmetic for objective TDW

Use the existing `montgomery_clean` compute-copy-uncompute routine to produce
copied inverses and slopes. Each coordinate multiplication similarly copies
its result and reverses its multiplication gates. The reusable Toffoli target
pool has R_B=M_n(r_n+3(k-1)) wires, sufficient also for the later k independent
lambda/coordinate multiplications. This yields forward module formulas

T_PA_balanced(n,k) = 2 (r_n+5k-3) M_n,
D_PA_balanced(n,k) = 2 (d_n+2 ceil(log2 k)+2).

These are not substituted for the CNOT/full-depth costs in gate-level resource
records: the backend generates actual gates and validates the Toffoli formulas for every emitted
module. `--verify` additionally runs ordinary affine-point fixtures, checks
all declared reusable blocks are zero, and checks the external
compute-copy-uncompute restoration with the correct copied output.

Both Balanced interfaces have 4nk input wires. Tail inputs are x_R,y_R,z,v
from the unchanged QROM pair. Reduction inputs are x_R,y_R,x_T,y_T; as in
the original Balanced reduction, the numerator uses y_R ^= y_T, with a
fresh nk-wire z block. Its y-coordinate formula then uses the unchanged T
point. Tail instead uses R. Both evaluate the same affine sum and use the
actual curve coefficient in the bundled AC basis.

The Balanced clean/reusable block is

S_B = 4(M_n-n)k + R_B + 4n + n floor(k/2).

Only linear ancillas, the uncomputed multiplication-target pool, the
uncomputed square workspace, and recovery copies enter S_B. Copied
inverses, copied slopes, and coordinate outputs are NOT zero scratch:

G_tail_balanced = 3nk,
G_reduction_balanced = 4nk.

Thus module width is 4nk+G+S_B. Tail x3 occupies its existing z input;
reduction x3 occupies its new z block. Old inputs remain live, even after
leaving the active frontier. Tail QROM pairs still introduce 2nk fresh z,v
wires before arithmetic. The network uses the SAME live/peak recurrence
described below, replacing only G and S with these Balanced values.

For TDW search, compute both D_T(s,w) and this exact allocation-policy peak
W(s,w) for EVERY point in the specified grid. Select all ties minimizing
D_T*W, then generate their required modules and check their resources.
Rows without gate-level resource records expose only explicitly named
`predicted_*` fields, not complete resources. The result is a minimum within the Balanced family,
not necessarily a minimum across both families or all possible circuits.
There is still an outer inverse of the whole network in addition to the
internal Balanced uncomputation; these are distinct operations, both counted.

Each lane uses z=x_R+x_T and v=y_R+y_T, followed by

lambda = v/z,
x3 = lambda^2 + lambda + z + a,
y3 = lambda (x_R+x3) + x3 + y_R.

`reduction` allocates clean z,v and computes them from both input points.
`tail` begins with existing clean z,v and does not read QROM coordinates.
The x3 register aliases z; y3 lives inside a retained multiplication output
block. The curve coefficient a is converted from polynomial to AC L basis.
All square workspace remains allocated until the outer inverse, whether or
not individual square wires could be proven reusable by a stronger analysis.

## Exact allocation policy and module interface

Write B=M_n-n and R=M_n(r_n+5k-3). Both module types use

- 4nk input wires;
- R retained multiplication target wires;
- 4n retained square-workspace wires;
- S=4Bk+n floor(k/2) clean/reusable linear and recovery-copy wires.

A reduction additionally allocates 2nk fresh z,v wires. Its new retained
allocation is G_red=R+4n+2nk. A tail's new retained allocation is G_tail=R+4n;
its z,v are already in the input interface. Total module width is
4nk+G+S. Output width is 2nk; n k of the tail output wires alias its input z.
The input count is an interface size, not an assertion that all inputs are
unchanged: tail z becomes x3. Existing state and retained garbage must not
be treated as zero scratch just because they are no longer active points.

The generator checks S wires are zero on the deterministic basis fixtures.
The multiplication implementation restores its linear input workspaces;
recovery copies are explicitly undone. No inversion/product/output blocks
are recycled between arithmetic layers.

## QROM model

For K=2^b entries, all words have upper-model weight h=2n:

T_lookup=D_T_lookup=2(K-2),
C_lookup=K-2+K(3h-2),
D_lookup=3K-4+K(2 ceil(log2 h)+1).

After internal cancellations in lookup -> clean z,v -> literal unlookup:

T_pair=D_T_pair=4K-2b-6,
C_pair=2 C_lookup+4n-2(h-1),
D_pair=2 D_lookup-2b-2 ceil(log2 h)+3.

Initialization outputs occupy 2n persistent wires per table; scratch occupies
2n+b-2 wires excluding the existing address and the persistent output.
A lookup pair introduces 2n persistent z,v wires. Its temporary table,
selector and data-fanout workspace is 4n+b-2 wires, excluding the address,
the existing accumulator and the new z,v.

Parallel counts and allocations sum; depths take the maximum. After early
unlookup, all QROM-only scratch can be reused by the arithmetic tail.

## Network liveness

Reserve 2m address wires and 2n final result-copy wires throughout.
Maintain live retained width W and allocated peak P. Initialization adds
2np live wires, with its disjoint QROM scratch counted at the peak.
For each QROM pair batch first add 2nk to W and evaluate its scratch peak.
For each arithmetic layer add G to W, and update P=max(P,W+S).
Only zero scratch is recycled. Older input points and nonlinear garbage
remain part of W even after leaving the active accumulation/tree frontier.
This is equivalent to an allocator that obtains new live wires from its
free-zero pool and grows only when a stage's simultaneous demand exceeds P.

Copy the final 2n output bits and reverse the entire forward construction.
The reversal uses the same wires, restores z before replaying earlier
QROM pairs, and returns all work to its initial state. Count and Toffoli
depth are twice the forward totals. Following the existing table convention,
full NCT depth is 2 D_forward+2, including a conservative output-copy boundary.
DW and TDW are computed with P, not a sum of all layer widths.

## Limits of the estimates

Only individual arithmetic modules are gate-exact. The overall network has
explicit module barriers. There is no cancellation, overlapping schedule,
or optimizer across module boundaries. QROM CNOT count/full depth are
table-independent upper estimates. A minimum reported here is a minimum
within this specified schedule/model and search range, not a lower bound
over all possible quantum circuits or QROM implementations.

No measurement-assisted uncomputation is used. The ordinary affine-add
validity assumptions are inherited; the test fixtures choose valid ordinary
additions and do not establish a new exceptional-input probability bound.
