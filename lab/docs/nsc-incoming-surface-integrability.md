# Reference integrability on compatible incoming normal functions

For the local compatible family

\[
\delta r=T\,w(z)+T^3U(z)/6,\qquad\delta a=0,
\]

the included lapse and shift reference changes are absolutely integrable in
canonical momentum through the owned formal order four. The intrinsic
surface, normal/frame, lapse/shift and abstract C0 remain unchanged. This is
a local statement for finite compatible jets; it selects no functions,
spatial boundary conditions, history, duration or physical initial data.

Decision: reuse the existing projector recursion, raw metric vertices and
stored UV powers. The missing connection is the newly nonzero first-order
projector when `w=delta(r_T)` is nonzero. Check its two constraint contractions
and the spatial Weyl terms by exact2x2 algebra; then import the existing
higher-order decay result. Stop after these identities, provenance replay and
focused tests. No reference quadrature, numerical probe or radial generator
is required.

## First order changes, but its constraint insertions vanish

On the fixed incoming surface write

\[
H=-m\sigma_1+(\ell/r)\sigma_2+(k/a)\sigma_3,\qquad
\omega^2=m^2+\ell^2/r^2+k^2/a^2,
\quad P_0=(I-H/\omega)/2.
\]

Every pure spatial intrinsic derivative is zero at T=0. In particular,
`H_z=(P0)_z=0`, so the actual first Weyl idempotence term `G1` and the
first Weyl commutator in the transport equation vanish. This keeps the
spatial terms before imposing the surface identities. Mixed derivatives
such as `r_Tz=w′` do not alter those pure intrinsic derivatives.

The owned recursion therefore gives

\[
\Delta P_1=\frac{i[H,w\partial_r P_0]}{4\omega^2}
=\frac{\ell w}{4r^2\omega^3}
  \left(\frac{k}{a}\sigma_1+m\sigma_3\right).
\]

This matrix is generally nonzero. Nevertheless, with the existing raw
vertices `V_N=H` and `V_beta=-k I` at N=1,beta=0,

\[
\operatorname{tr}(\Delta P_1V_N)
=\operatorname{tr}(\Delta P_1V_\beta)=0
\]

exactly for each momentum and angular sign. No occupation symmetry is used.
For the bare intrinsic vertices, all spatial derivatives vanish and their
momentum dependence is affine. Their only possible Moyal correction in the
symmetric traced insertion is `i*tr[Pj_z,(V_B)_k]/4=0`; higher ones vanish
because second momentum derivatives of the vertex are zero. This statement
does not discard terms of a nonconstant test envelope pointwise: the bulk
Euler coefficient still uses the existing complete unweighted trace and
compact-test integration-by-parts domain.

## Reuse the higher-order tails

The [existing UV certificate](nsc-reference-band-bulk.md) applies to all
smooth finite metric jets with positive nondegenerate N,a,r. It gives
`Pj=O(|k|^(-j-1))`. For j=2,3,4 the at-most-linear raw vertices yield powers
`-2,-3,-4`, hence absolute momentum integrability. Zeroth order is unchanged,
and the first-order constraint contractions vanish exactly as above.
At finite momentum each retained non-LLL channel has nonzero normal gap.
Finite inventory and multiplicities preserve finiteness. Bounds are local
in z, or uniform on compact z sets with bounded jets; no global falloff is
introduced. The LLL keeps its owned fixed chiral reference charts and its
separate geometric allocation once.

The [record](../results/development/nsc-incoming-surface-integrability.json)
contains exact zero residuals and the imported order bounds. This does not
prove pressure-direction integrability, a full4D Hadamard statement, a
surface solution, transmitting endpoint completion, or extended stationarity.
The fixed rho=0 seam and all other declared action terms remain unchanged.

```sh
python3 scripts/derive_nsc_incoming_surface_integrability.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_integrability.py
```
