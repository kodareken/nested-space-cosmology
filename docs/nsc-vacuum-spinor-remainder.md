# Upstream vacuum spinor from its covariance remainder

The owned homogeneous vacuum calculation bounds the exact unit Bloch vector
by a finite polynomial: |n-n_M| <= C_n/E^M. At the one upstream slice choose
u0 = sqrt((1+nz)/2) > 0 and u1 = (nx+i ny)/(2u0).

This fixes one common phase per energy, held constant during subsequent
spacetime evolution. Imposing a new real-major condition at every radius would
introduce an additional connection; this construction does not do that.
The rank-one vacuum covariance is unchanged. Finite thermal/coherent terms
remain separate, as does matching to the actual three-column preparation.

Set x=1/E and evaluate the chart on n_M(x). On the whole interval [0,1/E0],
let a_min>0 bound both exact and approximate major components from below,
and let W bound the approximate transverse magnitude. Exact differences of
square roots and reciprocals give the conservative Euclidean bound

    |u-u(n_M)| <= [3/(4 a_min) + W/(8 a_min^3)] C_n x^M.

The truncated Bloch vector is not assumed to have unit norm. The code checks
chart separation after including the Bloch remainder. Expand the chart of
n_M through degree M-1. Arb Taylor coefficients of order M over the entire
inverse-energy interval bound its remainder by C_T x^M. The resulting constant
is C_u = C_T + [3/(4 a_min) + W/(8 a_min^3)] C_n.

This is a uniform upstream-column bound for every E>=E0, not a refinement
estimate. The upstream envelope is spatially homogeneous, so its axial
derivatives vanish there; the physical carrier derivatives do not vanish.
It does not yet give the propagated changed-history remainder or its full
Sobolev norm, and cannot fill the UV-tail budget by itself.

Leading coefficients reproduce the inherited minor coefficient and real-major
normalization term Re(A2)=-|A1|^2/2 in this fixed phase convention. Tests check
these coefficients, normalization through the retained order, an independent
full Bloch ODE comparison, and rejection of lost chart separation.
Independent Grok review of the core conversion found no under-enclosure and
ran its original three tests successfully. A supplementary review also accepted the truncation adapter, signed control
and record builder, and independently passed all seven tests.

Owners: `src/recursive_horizons/nsc_vacuum_spinor_remainder.py` and
`tests/test_nsc_vacuum_spinor_remainder.py`. Reuses the
[vacuum covariance remainder](nsc-vacuum-source-remainder.md).

Run from the repository root:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_vacuum_spinor_remainder.py
```

For the paired negative-energy vacuum, the occupied covariance is the
complement of the sigma3-conjugated positive covariance. If the positive
representative is (a,v) with a real, the negative representative is (v,a).
This swap preserves the error norm. Merely sigma3-conjugating the positive
occupied column selects the wrong vacuum projector; the signed-partner test
explicitly detects that mistake. This statement is about the vacuum term,
not a substitute for the actual finite-occupation signed source law.

## Connection to the existing fourth-order envelope recurrence

`truncated_upstream_envelope` returns this same phase-fixed polynomial through
A4, not a new set of zero boundary constants. It bounds the omitted degrees
5 through M-1 explicitly and adds the proved degree-M remainder. Multiplying
by E^4 and maximizing inverse powers at the cutoff gives a uniform constant.
For a homogeneous envelope on a period of length L, the initial Sobolev triple
is (sqrt(L) times that constant, 0, 0). These are envelope derivative norms;
physical carrier derivatives must still be included by the existing contraction.
The subsequent recurrence must use exactly the returned upstream coefficients.
No existing UV record is changed by this helper; propagated L0 A4 bounds remain
missing. The supplementary independent review accepted this initial-envelope connection.

The record builder `scripts/derive_nsc_vacuum_spinor_remainder.py` binds the
original group-1 and group-14 masses, both angular signs, the actual archived
upstream radius, and cutoffs 160/320. It stores coefficient intervals and
pointwise remainder constants, with the spatial norm domain and changed-history
remainder still unset. The v1 result record is generated only after the independent review. Its test
recomputes all eight cases and checks those scope markers.

An additional independent control evaluates the original metric directly at
the upstream radius. With l=a(m-i ell/r)/2, the first two minor coefficients
agree with l and -i a^2 (d l/d rho)/2 from the Dirac recurrence. Reversing the
operator i-sign fails the control. This verifies the initial coefficient
identification; it does not replace their subsequent retarded transport.
