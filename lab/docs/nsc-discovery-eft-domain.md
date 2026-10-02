# Actual-unit local-curvature EFT diagnostics

This finite consumer reads the uniform, pattern and coherent half-step
`nsc-discovery-dynamic-episode-v1` checkpoints and the balanced/lower/upper
half-step `nsc-discovery-initial-sector-balanced-v2` checkpoints. It computes
actual metric jets, four-dimensional curvature, and normal-clock field scales.
It makes no new Cauchy data, trajectory, source choice, error enclosure or gate.
The owner is [assess_nsc_discovery_eft_domain.py](../scripts/assess_nsc_discovery_eft_domain.py),
with [tests](../tests/test_nsc_discovery_eft_domain.py).

## Units and the existing validity requirement

The owning [curvature EFT note](nsc-curvature-eft.md#quantum-and-validity-requirements-remain-explicit)
requires all relevant `|c_i curvature|/A` and `|c_i T|/A^2` to be small, with
external scales below the matched local expansion scales. It also retains the
finite-regulator/state/boundary Jacobian and nonlocal remainder. The diagnostics
here do not establish all those requirements or a full nonlocal-theory domain.

For the actual conformal spherical metric

\[
g=r^2[Q^2(dt^2-dx^2)-d\Omega^2],\qquad N=q=rQ,
\]

the code's t/x/period are dimensionless chart coordinates in the locked reference
units. Areal radius r, N and q have physical length; Q and chi are dimensionless.
The raw two-dimensional `R_h` belongs to the dimensionless h metric. It must not
be multiplied directly by `C_W/A` as a physical four-dimensional curvature test.

| Quantity | Physical units |
|---|---|
| Einstein coefficient A | length^-2 |
| Weyl coefficient C_W | dimensionless |
| abs(C_W)/A | length^2 |
| R4, tetrad Riemann/Ricci components | length^-2 |
| Ricci square, Kretschmann K, Weyl square | length^-4 |
| Coordinate representative H level | dimensionless |
| Proper frequency, pole comparison mass | length^-1 |
| Normal-observer stress | length^-4 |

The locked coefficient ratio is `abs(C_W)/A = 0.01875`, in length-squared units.
The extra flat-mode comparison mass is

\[
m_2=\sqrt{A/(2|C_W|)}=5.163977794943227,\qquad
\epsilon_{\cal R}=|C_W|\,\mathcal R/A=\mathcal R/(2m_2^2).
\]

This flat scale is not a separately established matched cutoff and does not
determine the extra mode's physical health. Nothing changes the action
coefficients, magnetic flux, Dirac kappa, copy count or source covariance.

## Physical curvature and cancellation

The consumer calls the actual full-rate analytic acceleration and metric-jet
APIs in [nsc_discovery_tidal.py](../src/recursive_horizons/nsc_discovery_tidal.py).
It retains both

\[
R_4=(R_h-2)/r^2-6(r_{tt}-r_{xx})/(r^3Q^2),\qquad
C^2=(R_h-2)^2/(3r^4).
\]

Chi is not substituted for either curvature. The physical scale table includes
absolute R4; maximum absolute normal-tetrad Riemann, Ricci and tidal components;
and `sqrt(abs(Ricci2))`, `sqrt(abs(K))`, and `sqrt(Weyl2)`. These square roots
restore length^-2 units. Lorentzian scalar contractions can cancel and are not
positive norms of all tensor components. A small R4 or finite Weyl scalar does
not bound the normal observer's physical tides. The report retains the signed
invariants, individual components, scalar-route discrepancy, and discrepancy
from subtracting large invariant terms to obtain a small Weyl scalar.

The tetrad is the actual normal frame at each point, with `e0=N^-1 ∂t` and
`e1=N^-1 ∂x` as coordinate tangent vectors, and angular vectors
scaled by r. The default recorded worldline is x2. Its component maxima are
observer-specific, while the scalar contractions are invariants. Condition
indicators remain indicators, not arithmetic or continuum error certificates.

## Spectral versus source-state scales

The instantaneous representative H is diagonalized on the actual full fermion
band. Its six lowest positive levels are converted to normal-observer frequencies
by `omega=E/(rQ)`, with global ranges and the values at x2. These are instantaneous
frozen-generator measurements. The evolved covariance is not asserted to occupy
those six eigenvectors. Initially the uniform source uses six positive standing
modes; pattern uses the declared Q profile, and coherent rotates the fixed
unequal-weight reference columns. Their later source evolution is retained.

The exact angular Dirac shell term is `kappa/r`. For spatially uniform geometry,
the instantaneous dispersion is

\[
\omega_k=\sqrt{[k/(rQ)]^2+(\kappa/r)^2},\qquad
k=2\pi(m+1/2)/L.
\]

A decreasing angular term does not bound radial proper frequencies when N shrinks.
Multiplicity is counted in source forces, never as a multiplier of H frequencies.
The consumer separately reports the normal-clock temporal RMS of the actual
canonical columns:

\[
\omega_{\rm state}=\frac1{rQ}
\sqrt{\frac{\sum_j c_j|U_f\dot\phi_j|^2}
{\sum_j c_j|U_f\phi_j|^2}}.
\]

This is a source-state variation scale in the owned half-density representation,
not a particle-occupation or adiabatic-state certificate. No shell is reselected
and no state is reset. External background and higher-derivative scales, the
stress budget and the nonlocal remainder remain additional validity requirements.

## What the saved stations show

The original preparations have physical component ratios around `0.0024` at T0,
`0.0125`–`0.0156` at T0.3, and `0.43`–`0.50` at T1. At T1 the square-root
K ratios are already `1.66`–`1.95`; this is not a strong local curvature-EFT
scale separation. At T3 the normal tidal ratios are approximately `4.05e6`,
`4.69e6` and `8.70e6` for uniform, pattern and coherent, respectively. Their
lowest positive proper-shell ranges reach hundreds to thousands of m2. Tiny
R4 and order-one-or-smaller Weyl ratios do not contradict those measurements.

Balanced-v2 changes only the uniform preparation's declared initial auxiliary
sector and independently constrained radius. Its balanced root satisfies the
conditional initial chi-double-dot criterion; higher adiabatic jets are not
matched. The balanced branch has physical component ratios `0.00248`, `0.00255`,
`0.00847`, then `3.31e4` at T0/T0.3/T1/T3. Its upper perturbation reaches
`7.20e4` at T3. Thus the initial criterion delays growth in this finite family
but does not keep the measured T3 state in a separated local-curvature regime.
The lower half-step branch stops at T about `1.482`; its ordinal-three checkpoint
is that last positive-chart state, not a T3 observation. The report retains the
actual time, chart-exit status and event record without extrapolation.

These are saved nf128 half-step diagnostics, not spatially certified continuum
statements or a global theory verdict. Loss of separation for a local EFT
assumption does not declare the full nonlocal theory invalid.

## Creation and checking

Default is pure read-only postprocessing and prints the full report. Both
datasets must be present; the consumer does not silently replace a missing
family or unfinished station. Write and check modes are exclusive. Writing
requires a full immutable producing commit containing the new consumer,
mathematical note, tests and its transitive local Python source closure.
Historical saved-producer closures are authenticated separately. Creation
refuses an existing exact output and makes the new JSON read-only. Its canonical
record digest binds the tables and metadata, and input hashes bind every consumed
manifest, checkpoint, payload and initial preparation. Check authenticates those
closures, requires matching current measurement sources, and remeasures the
saved table with an explicit numerical replay tolerance. It writes zero bytes.

After root freezes and commits the producer:

```sh
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_eft_domain.py
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_eft_domain.py \
  --write results/development/nsc-discovery-eft-domain-v1.json \
  --producer-commit FULL_COMMIT_SHA
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_eft_domain.py \
  --check results/development/nsc-discovery-eft-domain-v1.json
```

Focused tests include the exact constant-radius control (zero scalar and Weyl
but nonzero physical tides), uniform proper Dirac dispersion, the saved lower
branch's chart stop, no-preparation/no-evolution dispatch, creation-only output,
and integrity/closure/replay dispatch. Mocked provenance in temporary test records
is labelled as a dispatch test, not scientific producer evidence. Root creates
the actual committed diagnostic record after the consumer freeze.
