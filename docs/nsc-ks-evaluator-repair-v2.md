# Evaluator repair and original-source controls, September 25

The source-fixed local gate remains OPEN. This delivery repairs the numerical
operator and cache identity, then measures the actual remaining value errors.
The physical action, source, history class and interval are unchanged.

## Propagation and replay

The augmented variables are the homogeneous reference A, the difference
D=X-A, and the retarded directions Y. The CF4 linear operator now supplies
the true Euclidean adjoint of all these blocks. Fourier differentiation
changes sign under adjunction, local matrices are conjugate-transposed,
the adjoint of spatial replication sums over z, and the tangent coupling
has both Y-to-D and Y-to-A contributions. Its analytic trace is supplied
to the exponential routine. The former all-ones norm substitution is gone.

The chronological two-exponential product is tested independently against
noncommuting time-dependent matrices, in both time directions. The observed
halving ratios are about 16 for CF4 and 4 for reversed ordering. These are
implementation controls, not a continuum error enclosure.

Production operator evaluations use primal-block DOP853 control. Optional
tangent tolerances act on each direction separately; adding zero directions
cannot dilute a nonzero direction's tolerance. The existing streaming
DOP853 trajectory path remains available for interval validation.

Version-2 value caches bind the actual source arrays, sampled history,
channel inventory, numerical settings, nodes, coefficients and first-party
implementation dependencies. Workers compare their current binding with the
parent checkpoint before evolving. A complete value record names its family
JSON inputs, numerical payloads and binding file. Old version-1 data are
historical replay only and are not reused for new production evaluations.
The old scalar residuals are not automatically re-ranked after these repairs.

## Fixed-history controls

The current history is `0b0e4ced…`. Four existing energy rows are selected
from each of positive families `14_1` and `32_-1`. Their original
columns, weights and complete three-fibre coherence blocks are retained;
the original negative-energy partners are also included. The rows span
approximately `0.0013–160` and `0.00058–320`. No source preparation
or source interpolation is performed in these controls.

Tight DOP853 uses `rtol=5e-14`, `atol=5e-19` and
`max_step=1/8192`. Every result below is a refinement indicator for
these selected rows, not a bound or a full-source residual.

| Comparison | N movement indicator | beta movement indicator |
|---|---:|---:|
| Refined CF4 against tight DOP853, 256 nodes, period 1.6 | 7.5454e-12 | 9.8855e-13 |
| Spatial grid 256 to 512, tight DOP853, period 1.6 | 1.4866e-7 | 6.3687e-8 |
| Spatial grid 512 to 1024, same period | 1.4317e-9 | 1.5430e-9 |
| Spatial grid 1024 to 2048, same period | 4.0679e-12 | 2.5107e-11 |
| Spatial grid 2048 to 4096, same period | 3.4970e-15 | 2.6370e-14 |
| Equal spacing: period 1.6/grid4096 versus period 0.4/grid1024 | 1.7737e-17 | 1.6696e-17 |
| Period 0.4, grid1024 to grid2048 | 1.0084e-16 | 9.7135e-17 |
| DOP853 time-step halving, grid512/period1.6 | 1.1096e-18 | 9.7606e-19 |

The comparison owner adds the two families' componentwise sup differences,
retaining both signs within each family before taking the maximum.
All arrays and comparisons are stored in versioned development records.
For example:

```sh
PYTHONPATH=src python3 scripts/derive_nsc_ks_source_control_v2.py --check --name dop-tight-g1024-l04
PYTHONPATH=src python3 scripts/derive_nsc_ks_source_control_comparison_v2.py --check --name boundary-matched-spacing
```

These observations motivate the shorter resolved domain for subsequent
complete-family tests. They do not declare search or certificate numerics
for the whole source.

## Complete families expose a reference-reconstruction floor

The subsequent controls use all retained source rows of families `14_1`
and `32_-1`, their separate negative-energy partners, and the
energy-interpolated operator. Increasing degree 32 to 48 changes N by
`2.1696e-11` and `7.2874e-11` respectively on grid1024/period0.4.
This is much larger than the ideal interpolation remainder alone.

A separate zero-history control has exactly the reference geometry. Its
matter change should be zero, but the old reconstruction gives
`(8.476267935e-11, 1.513078320e-11)` in family `32_-1`.
This identifies a numerical mismatch between the interpolated homogeneous
reference and the independently propagated reference used by the contraction.
It is not a physical compensating stress and is not subtracted from another
history's residual.

The subsequent numerical refinement uses the exact decomposition
`U_g = U_ref + (U_g-U_ref)`: propagate the unchanged reference at the
original source energies and interpolate only the history-dependent difference.
Fresh operators carry `reference_mode=direct-original-energies`;
historical operators without this marker retain the earlier interpolant.
The mode participates in the operator digest. The four-field and retarded
tangent controls against direct evolution pass with the unchanged source.

The repeated complete-family zero-history control now has matter change
`(0,0)` exactly. This is the expected reference identity, not a root:
the baseline constraints remain. On the current nonzero history, degree
32 to 48 changes family `32_-1` by only
`(2.52055321e-15, 2.42102345e-15)` and family `14_1` by
`(5.59956541e-17, 4.15486601e-17)`. The source matrices, columns,
phases and weights remain unchanged. No measured offset enters reconstruction.

These controls justify evaluating the full retained source at the candidate
settings. They do not establish a full source/field enclosure.

With the same direct-reference reconstruction, doubling the spatial grid
from 1024 to 2048 at period 0.4 and degree 48 changes the complete
family `14_1` contraction by `(2.712e-16, 5.260e-16)` and `32_-1`
by `(3.006e-14, 3.612e-14)`. The comparison explicitly verifies the
serializer-only code change against its pinned producer before comparing
the saved arrays. It does not waive a changed numerical/source dependency.

## Error components and interpretation

The successor budget preserves the V5 covered-baseline bound
`(8.916253820835393e-12, 3.882373567206641e-18)`. The version-3 budget also
includes the current difference-interpolation enclosure; seven components
remain missing. A field refinement difference is not assigned to unrelated
baseline components. The N baseline bound exceeds its provisional `1e-12`
allocation; no allocation has silently changed.

The new exact Chebyshev profile bound treats binary coefficients and the
actual domain map as rationals. It uses endpoint bounds for the polynomial
and its first two derivatives before applying the existing cutoff-product
bounds. This avoids the enormous overestimate caused by summing absolute
power coefficients of a high-degree Chebyshev polynomial. The endpoint
inequality follows from the standard Chebyshev/ultraspherical bounds and
derivative identities: [DLMF 18.14](https://dlmf.nist.gov/18.14),
[DLMF 18.9](https://dlmf.nist.gov/18.9).

For the current w/U history, the ideal degree-48 interpolation component
over all 60 positive families and 120 signed families is enclosed at about
`(1.208e-26, 9.460e-27)` for the full-operator interpolant. This does
not include numerical node solves, reconstruction arithmetic, source
preparation, quadrature, the changed-history UV tail or between-node
residual error. A difference-only interpolant needs its corresponding
majorant; the full-operator bound is not transferred by renaming it.

The version-3 difference-interpolant bound explicitly uses the triangle
inequality for `U_g-U_ref`. Its field and axial-field interpolation
errors are at most twice the earlier uniform majorants. The owned stress
error polynomial has nonnegative coefficients and degree at most two in
these errors, so four times its previous bound is conservative. The
degree-48 component for the current history is therefore enclosed at
`(4.831544944e-26, 3.784029101e-26)`. The direct reference solve's
numerical error is separate and remains to be enclosed. No tangent
interpolation, UV or source-preparation bound is inferred.

At degree 32 the corresponding conservative component is
`(2.213119788e-9, 1.733298516e-9)`, which exceeds its allocation.
The small observed degree-32 movement therefore does not justify using
degree 32 for certification. The degree-48 record replays with the pinned
`python-flint==0.9.0` validation environment:

```sh
PYTHONPATH=src .venv/validation/bin/python scripts/derive_nsc_ks_difference_interpolation_accuracy_v3.py --check --degree 48
PYTHONPATH=src .venv/validation/bin/python scripts/derive_nsc_ks_gate_budget_v3.py --check
```

Grok independently reviewed the exact Chebyshev/radius and K0/D/K1
mapping. Its review confirmed the analytic bounds and emphasized that
floating profile evaluation remains a separate error. The follow-up
request about the reference reconstruction did not return a review receipt;
the reconstruction claims here use the stated identity and the completed
direct-evolution controls.

The reusable evidence graph checks all declared descendants, including
pinned historical Git sources. Its integrity PASS is deliberately separate
from a physical claim. The componentwise gate arithmetic requires every
error term, includes the between-node remainder, and compares outward bounds
with the exact decimal tolerance `3/10^11`.

No physical EXISTENCE or NON-EXISTENCE result, manuscript update, or public
release follows from this repair/control delivery.

### Local-Fourier successor

The v4 field-method successor replaces the global derivative-L1 estimate of
the omitted compact-profile band. It directly integrates 256 positive Fourier
coefficients with Arb, encloses the exponentially flat endpoint strips, and
bounds the analytic profile minus the retained Fourier polynomial on local
Taylor panels. On the same saved family-`14_1` cell, the continuous normalized
order-zero/order-one residual enclosure is approximately
`(3.00849e-16, 1.36374e-12)`, rather than v3's
`(5.481e-4, 6.429e-1)`. See
[the v4 record](../results/development/nsc-ks-current-field-pilot-v4.json).
The result validates one method/cell only; it is not transferred to all time
cells or source families.

## Full-family serialization control and recovery

The first complete-inventory campaign stopped after writing 34 families.
The reader correctly rejected the writer's omission of the redundant
`evaluation_identity` inside each NPZ metadata block. Earlier cache tests
constructed fixtures containing that field, so they did not exercise this
writer defect. The new test calls the actual writer and strict reader.

The corrected writer and explicit successor converter are in `fd347e3d`.
The converter requires the original pinned producer, proves that code
outside the serializer is unchanged, verifies every other source dependency,
restores the original operator digest, and re-contracts each stored operator
against the unchanged source. It then verifies bit-identical numerical arrays
in a new payload. Original JSON and NPZ bytes remain in place. Converted
records declare zero new operator solves and retain their original record,
payload and producer identities. This conversion does not make unbound
historical caches eligible for a new production request.

## Continuous-proof pilots and dependency closure

The [phase pilot](nsc-ks-current-phase-accuracy-v2.md) compares the actual
192/384-point production contractions with directed true-gradient balls.
Its three nodes leave a maximum N enclosure of about `1.20e-12`, slightly
above the provisional allocation. That miss is dominated by the proof's
remainder, not the observed change when the quadrature is doubled.

The [field pilot](nsc-ks-current-field-pilot-v2.md) captures and validates
one current-history DOP853 cell. It forms the difference residual before
taking norms and includes all fourteen `w^p U^q` monomials. The first
Fourier configuration retains only index 64 and returns an unusably loose
profile-tail bound. This is a limitation of that enclosure configuration;
it excludes neither a tighter proof method nor any physical history.
The actual Fourier period comes from the stored grid spacing, including
its difference from the nominal period 0.4.

`nsc-local-gate-evidence-components-v3.json` binds the existing covered
baseline and the current difference-interpolation component through 1020
nodes, including four pinned historical versions. Its integrity replay
passes; it deliberately makes no physical claim. The older aggregate
bookkeeping chain separately names an unavailable historical hash of
`docs/nsc-ks-gate-value.md`. That gap is not silently accepted or repaired
by changing the old bytes.

The combined operator, cache, mutation, reconstruction, phase and bound-helper
suite passes 93 tests. Seven further focused tests cover the new difference
residual polynomial. These are implementation and conditional-bound checks;
the whole-source continuous field and remaining source errors are still open.

## Prepared-source error that remains to be bounded

The source inventory and retained upstream archive explicitly leave
`continuation_error_bound` and `physical_source_error_bound` unset.
Grok's read-only audit verified that the stored covariance reproduces the
owned covariance recipe bit for bit on families `14_1` and `32_-1`.
That checks the discrete recipe; it does not enclose its arithmetic error
relative to the exact formula or the error of the prepared field columns.

The missing input to the existing transport majorant is a directed bound
on those columns at the original
`rho_up=1159676904047903/1125899906842624`. Initial mode errors at
`rho=1` and the subsequent homogeneous continuation must be accounted for
with their actual domains. The V5 covered-background action bound and the
separate group32 vacuum-column certificate on `[24,32]` cannot be used as
this matrix error without an explicit identification and transport relation.

The relevant owners are `nsc_ks_source_inventory.py`,
`nsc_ks_retained_upstream_archive.py`, and the conditional transport in
`nsc_ks_difference_error.py`. This is a specific missing numerical proof
input for the unchanged source, not a new physical source or a global
cosmology requirement.

Review also found a domain issue in the new polynomial helper: a derivative
majorant valid on the cutoff support was being applied to arbitrary input
balls. The `T_2` control at `2 +/- 0.1` reproduced the failure. Commit
`316fca2d` restricts those intersections to their proven domain and uses
the complete finite Taylor enclosure outside it, with two regression cases.
The earlier field-v3 records retain their original bytes and source context
at `189b703b`; subsequent bounds must use the corrected helper and successor
records.

## Changed-history UV term

The current residual contains the baseline's existing background tail
values, the retained finite-source change, and the source-cutoff edge.
The edge is an integrated asymptotic **value**, as specified in
[the cutoff bridge](nsc-source-cutoff-bridge.md). Its nodal evaluation
error is the phase-value component already enclosed above.

The separate remainder `T_{B,Lambda}` after the retained cutoffs
160/320 is still unbounded. The bridge proves a finite inverse-energy
expansion and integrability, but its record has `numerical_C_M=None`.
The linear reference contractions and the profile Fourier remainder do
not supply that finite-history constant. Grok's audit traced this gap
to the missing executable envelope coefficients and their residual
constant, including the changed thermal/coherent remainder.

The next bounded construction uses the existing coefficient recurrence,
the unchanged upstream expansion and the actual N/beta vertices.
Restoring the source carrier lowers inverse-energy order, so a term
such as `h_s,z/E^2` is not treated as the complete leading current
without the accompanying higher density terms. No discarded tail is
set to zero from a formal order or from a finite-energy refinement.

That construction now supplies an additional analytic result. The
coefficient recurrence and density continuity give `n3=2s h_z` with the
normalized, spatially homogeneous upstream condition. Including the
carrier term and the actual N/beta vertices gives coefficients
`mu*s*delta(h_z)/(2*pi*a)` and `-mu*delta(h_z)/(2*pi)` at Sigma.
The h transport is odd in the signed angular eigenvalue. The original
30 angular pairs have equal multiplicities and matching source grids,
so their leading algebraic cutoff coefficient cancels. This is a scoped
vacuum-coefficient statement; `C_M`, the higher remainder and the thermal
contribution are still unbounded.

The current proof register is `nsc-ks-current-uv-coefficients-v4.json`,
with 34 coefficient identities and 24 Sigma identities checked exactly.
Its owner is `derive_nsc_ks_current_uv_coefficient_control_v4.py`.
The v2/v3 development records retain their bytes, but their working-tree
HEAD metadata was not a source snapshot. They are not proof inputs for
v4, which reconstructs the identities directly and refuses `--check`
when the actual record is missing.

## Consistent candidate measurements

All rows below use the corrected DOP853 value path, the same original
retained source, 129 target nodes, spatial grid 1024, period 0.4 and
energy degree 48. These are nodal measurements, not error-enclosed
continuum residuals.

| History | N maximum | beta maximum |
|---|---:|---:|
| Current `0b0e4ced…` | 9.640931165e-4 | 3.809601596e-4 |
| Older n64 `ad759424…` | 1.027056385e-3 | 6.836887627e-4 |
| Previous value leader `a9572341…` | 1.044453287e-3 | 3.565823073e-4 |

The current history minimizes the declared joint nodal merit among these
three. `nsc-ks-gate-reanchor-v2.json` replays all three source assemblies
before ranking them. No optimizer step or physical certificate follows
from that ranking. Certificate numerics remain undeclared while the six
remaining error components lack enclosures.
