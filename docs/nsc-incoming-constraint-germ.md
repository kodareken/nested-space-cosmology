# Two independent normal directions in the incoming constraints

On the restricted incoming plane
`u=delta(r_TTT)`, `v=delta(r_Tz)`, with all other jets fixed to their baseline,
the retained bulk functional has the structural form

$$
\mathcal E_N=S_N+c_u u+c_{v^2}v^2,
\qquad \mathcal E_\beta=S_\beta+c_v v.
$$

The numerical response gives approximately `c_u=0.05739928` and
`c_v=16.57598`, so the two directions are independently resolved. These
coefficients are determined responses of the existing action; they are
not adjustable gravitational, source or inheritance parameters.

## Reused functional and exact structure

The local Euler owner varies the locked compact action and charged light
restoration. Their terms have coordinate-derivative grade at most four.
Metric/radius denominators, logarithms and coefficients depend only on
undifferentiated fields. Euler variation lowers derivative grade by the
order of the differentiated field slot, then restores it through integration
by parts. The order-four reference projector likewise has formal order n
at derivative grade n. Its gap denominators depend on unchanged zeroth-order
fields. These are properties of the current owners, not a fitted polynomial.

Here u has grade three and v grade two. Thus `u*v`, `u^2`, and `v^3` cannot
occur. The scalar local actions are invariant under `z -> -z`, `beta -> -beta`.
For the reference Hamiltonian the simultaneous change
`(z,k,ell,beta) -> (-z,-k,-ell,-beta)` is implemented by conjugation with
sigma1. The reference recurrence and spatial Weyl product preserve it.
Both momentum signs and the inherited equal angular-sign weights are already
integrated. The lapse response is consequently even in v and the shift
response odd. The unchanged physical state need not have this reflection
symmetry; its contribution remains in the constants `S_N,S_beta`.

For this plane, zeroth and first metric jets on the surface are unchanged.
Hence `Delta P0=Delta P1=0` at the evaluation point. The
[existing reference power bounds](nsc-reference-band-bulk.md) give
`Pj=O(|k|^(-j-1))`, so the difference starts at `O(|k|^-3)` and the raw
vertices are at most `O(|k|)`. The resulting reference insertion has an
absolutely integrable `O(|k|^-2)` tail in each retained gapped channel.
The finite-k gap is nonzero; the LLL has its separately owned constant
chiral reference. This is an integrability statement for this restricted
finite-inventory plane, not a general twenty-jet or full4D state theorem.

## Bounded response experiment

Reuse the committed `v=0.01` response, the same local/reference owners and
the fixed C0 identification. Evaluate one labeled `u=0.01` response with
16/24 reference quadrature and the existing local contour refinement. The
physical matter insertion cancels exactly in these response differences,
so its unfinished numerical error budget does not prevent this rank check.
One independent mixed control `u=0.02,v=-0.005` audits the predicted form.
Stop after this structural/rank check; no optimizer or physical IV is used.

Risks are finite-contour aliasing, a numerically unresolved coefficient,
and confusing a point germ with a surface solution. The artifact retains
raw shift roundoff in the homogeneous probe and all new reference nodes.
Coefficient indicators remain numerical indicators. The exact grade/parity
claim is separate from their floating evaluation.

If the exact coefficients satisfy `c_u*c_v != 0`, the conditional relation is

$$
v_*=-S_\beta/c_v,\qquad
u_*=-\frac{S_N+c_{v^2}v_*^2}{c_u}.
$$

This defines how a local germ would depend on the actual fixed source.
The [record](../results/development/nsc-incoming-constraint-germ.json) does
not select these as physical initial data or claim a certified root. Numerical
source uncertainty must still be propagated for a numerical solution.
A surface solution additionally requires spatial normal-data functions whose
derivatives reproduce the mixed jets, constraints between sample points,
and a stated spatial domain with boundary/asymptotic conditions. Global
parent preparation, endpoint variation and extended stationarity remain OPEN.
No new full4D Hadamard/raw-heat gate, action term, duration or metric step is
introduced by this local relation.

```sh
python3 scripts/derive_nsc_incoming_constraint_germ.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_constraint_germ.py
```
