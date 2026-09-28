# FGC-1-HLT1-MON1: transactional classical health and typed stops

## Decision

`FGC-1-HLT1-MON1` turns PROTO3's health, constraint, scale, alias, and boundary
rules into executable pre-acceptance transactions. It passes the monitor gate,
not a trajectory gate: synthetic controls exercise every stop, while no
FGC-QR evolution or holdout outcome is opened.

The machine decision is

```text
classical_health_and_typed_stop_monitoring_verified = true
```

within that synthetic, pre-trajectory scope.

For each proposed Runge--Kutta stage, the numerical engine must construct one
complete `StageHealthSnapshot`. HLT1 evaluates every premise before committing
the stage. If one or more fail, the fixed scientific priority selects a typed
first failure, the last accepted state is unchanged, and later attempts return
the same first failure. There is no warning-only mode.

`stage_index` is a strictly increasing global monitor-transaction serial, not
a method-local Butcher-stage label. The physical `time` stored in a snapshot
is its actual stage abscissa and need not increase between internal stages:
the SSPRK3 comparator, for example, evaluates states at `0`, `1`, and `1/2`
of a step. Ordering by the serial prevents stale or duplicate transactions
without falsely rejecting that valid non-monotone stage sequence.

## Canonical dimensionless norms

The six fields are ordered

```text
(alpha, shift, lambda, R, phi, chi)
```

with natural-unit mass dimensions `(0,0,0,-1,1,1)`. For field dimension
`d_i`, HLT1 scales

```text
u_i     -> L0^d_i u_i,
p_i,q_i -> L0^(d_i+1) p_i,q_i.
```

Matrix monitors use the induced infinity norm after explicit row and column
scaling. Branch displacement is measured from the continuation root belonging
to the preceding accepted stage, never an arbitrary root at the same trial
state. Constraint components use

```text
max_i |C_i| / max(|physical source scale_i|, 1e-30).
```

The exact-zero convention is retained: zero divided by a genuinely zero
declared denominator maps to zero; any nonzero numerator with zero denominator
fails closed.

## Proper spectral and curvature proxies

The temporal estimator accepts only a causal past-time series sampled in
proper time; the radial estimator accepts a field sampled in proper radial
distance. Nonuniform data must be explicitly resampled before the estimator.
Spatial data use a compact C-infinity window inside the exact vacuum buffer.
The largest real-FFT bin above both an absolute `1e-30` floor and a relative
`2^-40` amplitude floor defines support. Occupation of the top one eighth of
the resolved spectrum is a distinct alias stop.

These choices are frozen before the holdout because the word “support” is not
operational without a numerical amplitude floor. They do not make a finite
Fourier window an exact continuum band-limit theorem.

The later [`FGC-1-CAL0-PREF2`](fgc-cal0-pref2.md) composition test shows that
this distinction is decisive: the largest-bin Boolean rejects PROTO3's own
smooth compact input on every frozen grid even while the weighted top-band
power converges rapidly toward zero. HLT1 remains a correct certificate of the
frozen estimator and transaction; CAL0 closes that estimator as a physical
admission rule for PROTO3. A successor must use a convergence- and error-
budgeted tail criterion rather than silently weakening this historical record.

PROTO3 defines the linear curvature scale

```text
kappa = max(|Kretschmann|^(1/4), sqrt(|Ricci eigenvalue|)).
```

Since `kappa` and the cutoff `Lambda` both have mass dimension one, HLT1 gives
the protocol key `curvature_scale_over_Lambda_squared` its dimensionally
consistent executable meaning

```text
(kappa/Lambda)^2.
```

The frozen maximum is `1/4`. This resolves an operational ambiguity; it does
not supply the missing omitted-operator or remainder estimate.

## Strict stage contract

Continuous upper limits are accepted only below their threshold, and
continuous lower limits only above it. Equality therefore stops before a
crossing can be accepted. The discrete Newton iteration count is accepted at
`16` and rejected above it. Required booleans must be true, and either temporal
or radial alias occupation stops the run.

The monitor composes:

- positive lapse, radial metric, noncentral areal radius, and hat-cone margin;
- SRC1 residual, iteration, monotonicity, condition, branch, and acceleration
  limits;
- HYP2 coefficient, kinetic, Planck, characteristic, eigenframe, and
  symmetrizer margins;
- CON4's absolute normalized constraint monitor and the cross-resolution
  convergence floor;
- field-amplitude, proper-frequency, proper-wavenumber, alias, and curvature
  scale proxies; and
- BND2's remaining all-cone causal margin.

All 26 stop reasons are independently injected. A simultaneous-failure test
also verifies the frozen priority and immutable first-failure receipt.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt1_mon1.py \
  --output results/fgc-1-hlt1-mon1.json
```

The canonical result is hash-bound to SRC1, CON4, DOM4, HYP2, BND2, ACT1,
VAR1, and the shared PROTO3 semantic run envelope.

## Nonclaims

HLT1 proves that the monitor formulas and stop transaction behave as frozen on
direct and injected controls. It does not prove that a future solver computes
the underlying fields correctly; NUM1 owns that integration and validation.
It does not prove an evolution exists through the desired time, that any stage
will be accepted, that collapse or defocusing occurs, or that the declared
scale proxies establish retained-EFT validity. A scale stop is a classical
reason to cease interpretation, not evidence of a UV completion.
