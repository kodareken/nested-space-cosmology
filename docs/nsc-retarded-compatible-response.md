# Retarded compatible-radius restriction pilot

## Reuse, decision and stopping condition

Reuse the causal construction in [retarded compatible preparation](nsc-retarded-compatible-preparation.md),
its completed geometry/prepared-tangent adapters, and the authenticated
`14_1` global reference columns. No radial mode equation, scattering, old
source/control producer or metric equation is rerun. The named gap is the
actual rho1 response of a causally prepared off-shell radius variation,
rather than an arbitrary initial/preparation tangent.

This one-family pilot uses amplitude zero, normal window inner0.007 and
outer0.03, U=0, and the existing normalized smooth bump
`w(z)=profile((z-S(1)-0.15)/0.06)[0]`. Its first three derivatives are those
of this same function. Caller PG times run from0 to0.3. These are response
control coordinates and widths, not physical initial data or a duration
selected by the source. The radius perturbation is zero at rho1 but its
normal derivative there is w.

The radial support and axial/time support inequalities are checked before
propagation. The initial time precedes the entire perturbation, and the
right edge is upstream of its radial support with both characteristic
speeds negative. Thus initial and right incident parameter tangents are
explicitly zero by the unchanged-past/source construction. The original
coherent source covariance is fixed, so its derivative is explicitly zero.

The only solves are `(401 grid points,64 midpoint steps)`, `(801,64)` and
`(801,128)`. Total producer CPU budget is120 seconds; a signal timer stops
and saves completed cases when exhausted. There is no automatic refinement.
For each of final field tangent, rho1 trace tangent and full sampled
covariance tangent, use maximum-absolute entry norms. Predeclared diagnostic
acceptance requires spatial difference / fine response norm <0.25 and
temporal difference / fine response norm <0.10. Failure gives OPEN and
stops this pilot. These differences are indicators, not continuum error
bounds or stationarity tolerances.

## Authenticated fields and preparation

The stored history artifact `9b4f7b14572697ffe4afa14653ced9eef57d02b0df561b961476322f8202786a`
contains `14_1/x` and `14_1/reference_field`:801 nodes on[-2,1.2] and
1602 characteristic rows by12 columns. Stride-two sampling within each
spin block gives the401 grid without a radial solve. The exact rho1 rows
are(750,1551) and(375,776), respectively. This group was not a failed family
in the later selective-refinement artifact; that artifact has no14_1 grid.

The matching source-mode artifact `f6d34987d41fd9ea72abdc3ac48b16d403a5b74895f25ef8aa52ca4e65e95232`
supplies four signed real energies, inherited positive weights, the complete
three-source covariance and projectors, m=pi/2 and angular label+sqrt(5).
All four energies have absolute value below m, so the infinity column is
closed. The original coherent horizon blocks and closed column are retained
without occupation or phase replacement. Weights enter exactly once through
`sqrt(weight/(2*pi))` on each source column. This is a reduced positive-angular
family response, without a new multiplicity factor or a claim to the full
angular/spectral sum.

At tau0=0 the saved columns are the reference initial data. Incident fields
carry `exp(-i E tau)` repeated for their three source columns and the
existing right SAT coefficients. At rho1 the solver stores both actual
spin rows at each time node. The frozen characteristic/current frame and
normal/half-density restriction give F and dF, with z=tau+S(1). There is
no interpolation and no extra stationary exp(i E S(1)) factor.

The full finite sampled two-position kernel derivative is
`dF C_source F† + F C_source dF†`, with all coherences and cross-position
entries retained. It is not a diagonal stress integral. An independent
stationary reference trace is formed from the archived columns and their
known time phase. Drift of the unperturbed discrete field and trace from
that reference is reported separately and never subtracted as a physical
source or fitted correction.

## Result and interpretation

The [record](../results/development/nsc-retarded-compatible-response.json)
contains the actual response norms, flux/support/Hermiticity checks,
stationary-reference drift, resolution differences, elapsed CPU and wall
time, and PASS/OPEN according to the predeclared indicators. Its small
artifact preserves final fields/tangents, actual rho1 traces, and full
sampled two-position covariance tangents. Receipt replay computes metrics
from those saved arrays and does not propagate another field.

The complete canonical CAR kernel on the unchanged intrinsic Sigma has
zero variation after full source resolution. A four-energy projector kernel
need not: the same traces also give its sampled tangent `d(F P F†)`, recorded
as a diagnostic without enforcing zero or subtracting a source/stress
correction. A nonzero finite covariance response must not be mistaken for
a completed change of canonical state.

A nonzero numerically resolved result concerns this sampled prepared
component only. Four inherited control energies are not a complete
quadrature; omitted frequencies and the existing analytic subgap/tail
representations prevent a conclusion about the full incoming C0 equation.
The compact variation is not a global constraint solution. Full parent
matching, extended stationarity, stress and metric evolution remain OPEN.
The fixed rho0 seam, all Gamma allocations and the93.54264532195464 diagnostic
are unchanged; no extra action or counterforce is introduced.

```sh
python3 scripts/derive_nsc_retarded_compatible_response.py --run
python3 scripts/derive_nsc_retarded_compatible_response.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_retarded_compatible_response.py
```

## Completed bounded result

The three solves used5.544 CPU seconds and returned OPEN. Spatial relative
indicators were0.702 for the field,0.612 for the trace and0.630 for the
covariance, above the0.25 limit. Temporal indicators were0.00184,0.00311,
and0.00316, below0.10. No fourth solve was run. The fine covariance tangent
maximum was7.75998e-6; its sampled projector counterpart was1.18031e-5.
Fine stationary-reference field/trace drifts were2.15385e-4 and3.29398e-5.
The response is spatially unresolved under the declared check; the algebraic
flux and Hermiticity residuals do not override that conclusion.

A post-run audit found that the fixed-surface constructor received
`(N,beta,q,r)` instead of `(N,q,beta,r)`. The producer now uses explicit
arguments in the correct order. For this fixed surface, the used trace,
inverse and Gram arrays—and the final restriction matrix—are bitwise
identical under the correction: they depend only on N,r,q*beta, while the
external r*sqrt(q) divisor was already correct. The coordinate normal and
orbit metric are not claimed invariant. The artifact preserves the exact
executed producer/test/document bytes and execution signature separately
from current verification hashes, and authenticates every unchanged output
array. Only receipt/provenance was regenerated; no solve was repeated.
