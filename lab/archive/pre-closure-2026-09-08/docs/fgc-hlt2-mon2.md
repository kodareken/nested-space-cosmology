# FGC-1-HLT2-MON2: branch-owned and convergent PROTO4 admission

## Decision

`FGC-1-HLT2-MON2` implements the premise-only successor that
[`FGC-1-CAL0-PREF2`](fgc-cal0-pref2.md) required and
[`FGC-1-PRO4-FRZ1`](fgc-pro4-frz1.md) froze. Its positive machine decision is

```text
PROTO4_successor_admission_monitor_implemented = true
```

This means that branch applicability, weighted spectral admission, and
common-event constraint convergence now have executable definitions. It does not authorize a fresh GR-0 calibration run, an SGB-L comparison, or the
FGC-QR holdout. No trajectory or collapse outcome is read.

The distinction matters. HLT1 remains the correct implementation of the old
PROTO3 contract. CAL0 remains the correct diagnosis that two of that
contract's uses were scientifically invalid: a largest-occupied-bin Boolean
treated an arbitrarily small Fourier coefficient like unresolved power, and
FGC-QR-owned limits were applied wholesale to GR-0. HLT2 does not edit either
result. It provides the prospectively frozen PROTO4 interpretation.

## Exact branch ownership

PROTO4 partitions all 26 historical HLT1 stop identifiers exactly once:

```text
9 universal runtime stops
11 candidate-branch-only stops
6 replaced numerical-admission decisions
```

The universal stops apply to GR-0, SGB-L, and FGC-QR. They protect positive
metric factors, the auxiliary cone, solver residual/iteration behavior, the
kinetic solve, and the remaining causal boundary budget.

GR-0 cannot be rejected by any of the eleven candidate-action stops. Its
calibration role has no nonlinear candidate root, Horndeski deformation,
effective-Planck deformation, candidate characteristic frame, regulator
amplitude, or candidate curvature proxy to which those limits could honestly
refer.

FGC-QR activates the eleven stops only under the exact owner bundle

```text
FGC-1-SRC1-NL1 + FGC-1-DOM4-RUN1 + FGC-1-HYP2-MD1.
```

Those artifacts share the unredefined ACT1/VAR1 equations and an explicit
FGC-QR scope. HLT2 refuses an unnamed or different owner.

SGB-L remains fail-closed. The present HYP2 result explicitly targets FGC-QR,
and SRC1's branch-preserving nonlinear source likewise certifies FGC-QR only.
Therefore HLT2 does not copy their numerical health values onto SGB-L. The
matched SGB-L control required by the frozen analysis protocol cannot run
until its own source/branch and multidirectional-health definition exists.
This is a newly exposed missing premise, not a failure of the SGB-L dynamics.

The six replaced HLT1 IDs are retained as historical controls but are never
active PROTO4 decisions. Their new owners are the spectral and constraint
admission rules below.

## Proper, dimensionless spectral fields

The spatial estimator acts on the declared six Minkowski-reference
deviations

```text
alpha - 1
shift/r
lambda - 1
R/r - 1
phi/Lambda
chi/Lambda.
```

At the exact regular centre, division by a numerical epsilon is forbidden.
HLT2 uses the analytic limits

```text
(shift/r)|r=0 = partial_r shift|r=0,
(R/r)|r=0     = partial_r R|r=0.
```

The physical radial coordinate is constructed from

```text
s(r) = integral_0^r lambda(r') dr'.
```

The implementation uses a deterministic trapezoidal integral and resamples
the fields linearly onto a uniform proper-distance grid. It serializes the
round-trip interpolation mismatch by field. That mismatch is a diagnostic
and a future error input; it is not represented as a rigorous continuum
interpolation bound.

Spatial data use a two-edge compact `C-infinity` window inside the retained
measurement interval. Temporal data require at least 64 proper-time samples
and use only the causal past history available at the event being evaluated.
Both are finite-window estimates, not exact continuum band-limit theorems.

## Weighted tail and RMS decision

For real samples `f_j`, HLT2 computes a one-sided real FFT with the correct
Parseval multiplicities. If `P_n` is the field power in angular-frequency or
angular-wavenumber bin `n`, it evaluates

```text
field tail       = sum_top(P_n) / sum_all(P_n),
derivative tail  = sum_top(k_n^2 P_n) / sum_all(k_n^2 P_n),
k_RMS            = sqrt(sum_all(k_n^2 P_n) / sum_all(P_n)).
```

The top band is the highest one eighth of resolved bins. A nonzero signal is
admitted only strictly below

```text
field tail                 < 2^-20,
derivative-weighted tail   < 2^-10,
k_RMS/Lambda or omega_RMS/Lambda < 1/4.
```

Every adjacent refinement must additionally reduce both tail fractions by a
factor strictly greater than four, and at least three nested resolutions are
required. A zero signal has zero tail and zero RMS scale. A tail that
reappears after an exact-zero denominator fails closed.

HLT1's last occupied bin (the largest occupied bin above the `2^-40` relative
floor) is still recorded
as a diagnostic. It never enters the decision. The canonical control applies
the new estimator to the declared amplitude-2 compact GR-0 initial profile at
`1025`, `2049`, and `4097` points. Several last-bin diagnostics remain true,
but every weighted field/derivative/RMS bound passes and every tail decreases
far faster than the required factor. This is the intended distinction:
resolved smooth compact structure may have nonzero coefficients arbitrarily
deep in a finite Fourier transform without carrying material unresolved
power.

The control is an initial-slice estimator check only. It is not one of the
fresh calibration trajectories and does not select a matter amplitude.

## Common-event constraint admission

This common-event constraint admission is deliberately stricter than a
single-grid small-residual test.

HLT2 orders the complete spherical constraints as

```text
Hamiltonian, radial momentum, gauge_t, gauge_r,
six kinematic reduction fields.
```

At every point and for every component, the normalization denominator is

```text
max(sum_j |unredefined term_j|, M_Pl^2/L0^2, 1e-30).
```

The absolute terms are summed before the maximum; cancellations in the field
equation cannot make the denominator artificially small. Every denominator
must be finite and positive.

Three or more nested resolutions must refer to one common event. GR-0 uses a
bitwise-equal binary64 coordinate time. DEF1 will additionally require the
same affine generator label and affine origin. The stored grid count must
equal the sample count from which each constraint norm was formed. HLT2 then
requires

```text
coarsest normalized infinity < 1/50,
finest normalized infinity   < 1/1000,
each adjacent observed order >= 3/2.
```

An exact-zero pair is recorded as exact zero instead of inventing an infinite
order. Reappearance of a nonzero residual after an exact-zero coarse value
fails. A single small-resolution value can never establish admission, and
numerical smallness is not promoted into a continuum constraint-propagation
theorem.

Manufactured second-order residuals pass the gate. Constant residuals and a
deliberately misaligned common event fail. Exact-zero controls pass under the
explicit zero convention.

## Richardson and observable-error boundary

For arrays in the same observable units, HLT2 supplies the usual asymptotic
estimate only when it also receives a passed common-event admission record:

```text
epsilon_R = max |O_fine - O_coarse->fine| / (rho^p - 1).
```

Every named error contribution must be nonnegative, a constraint contribution
is mandatory, and the final operation is addition only. No API permits a
constraint or discretization error to be subtracted from a claimed physical
margin.

One prerequisite intentionally remains false: a normalized constraint
residual has no automatic conversion into units of the Raychaudhuri
observable. Such a conversion needs a stability/sensitivity map for the
complete evolved system. HLT2 therefore records

```text
constraint_to_DEF1_observable_stability_map_supplied = false.
```

This does not block the next static ledger or GR-0 calibration work. It does
block any later attempt to claim a positive DEF1 margin by simply adding a
dimensionless constraint norm to a dimensionful observable.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt2_mon2.py \
  --output results/fgc-1-hlt2-mon2.json
```

The canonical result binds PROTO4/PRO4-FRZ1, historical HLT1/CAL0, and the
three FGC-QR owner artifacts. It writes no calibration or holdout run path.

## Nonclaims

HLT2 does not complete a fresh GR-0 calibration, construct PRO4-HLD1, evolve
SGB-L or FGC-QR, form a trapped interval, demonstrate regulator activation,
measure a positive complete Raychaudhuri margin, reject any action branch, or
say anything about gradients in general. It supplies neither retained-EFT
validity nor a physical transition, singularity resolution, child domain,
dark-sector mechanism, or varying locally measured speed of light.
