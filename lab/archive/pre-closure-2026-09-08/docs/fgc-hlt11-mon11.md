# FGC-1-HLT11-MON11 — PROTO13 restart-runtime authorization

`FGC-1-HLT11-MON11` is the pre-trajectory authorization for the PROTO13 GR-0
calibration restart. It is a numerical-premise artifact. It is not a collapse,
trapped-region, candidate-action, defocusing, retained-EFT, or physical result.

## Question owned by this gate

HLT11 asks only whether the six hash-bound states frozen by
`FGC-1-PRO13-FRZ1` can enter one fail-closed continuation engine without
changing the equations or pretending that unequal finest grids are the same
grid. The exact restart time is `t=23/16`; the first new event would be
`t=3/2`.

The frozen numerical runtime is CPython `3.14.3`, NumPy `2.5.1`, and
Accelerate on Darwin `arm64`. This is not selected for convenience: under that
backend the immutable RSP2 raw bundle and PREF14 endpoint reproduce exactly,
whereas the separately used OpenBLAS analytic runtime does not reproduce the
last bits of that evolved state. HLT11 fails before source evaluation or state
advance if the inherited RSP2 arithmetic contract differs. The eventual
manifest and terminal result must also serialize the complete observed Python,
NumPy, BLAS, and platform metadata so a resumed calculation cannot silently
change arithmetic backend.

The method-owned ladders are:

```text
RK4:    2049 -> 4097 -> 8193
SSPRK3: 4097 -> 8193 -> 16385
```

Every complete restart payload includes `u`, `p`, `q`, tracer positions,
tracer proper times, all 24 existing event samples, the runtime monitor, the
causal debit, accepted-step and transaction counters, source-retry count, and
CFL-retry count. No state may be initialized again, reprojected, interpolated,
fitted, or reset.

## Unequal-grid observable contract

The two methods independently evaluate the complete reference-balanced
constraint and PROTO12 pairwise spatial-spectral gates on their own ladders.
For the trapped-sign observable, all six states are restricted to the exact 2,049-node physical mesh
common to both ladders. Each method then produces its
own finest-pair Richardson interval using the unchanged minimum order `3/2`.

A trapped-sign event can pass only if:

1. both complete method-owned common-event admissions pass;
2. the two Richardson intervals overlap;
3. the smaller finest-grid sign margin exceeds four times the conservative
   combined error ledger; and
4. the later temporal-spectrum gate passes after at least 64 samples.

The restart event itself cannot count toward the eight new consecutive
qualified events. Constraint and spectral failures remain hard vetoes; they
are not converted into arbitrary expansion-error units.

## Authorization and temporal provenance

HLT11 verifies the immutable PROTO13 freeze and both raw checkpoint hashes,
recomputes all six SRC4 source prechecks without state mutation, recomputes
both complete method-owned common events, and proves the PROTO13 calibration
and holdout namespaces absent. It advances no trajectory.

The canonical result may authorize exactly one fresh GR-0 continuation in
`runs/fgc-2-sf1/proto13/calibration`. The authorization is committed before
the runner may create that directory. After launch, successor verification
uses the immutable authorization commit and raw hashes rather than replaying
the historical namespace-absence observation against a changed filesystem.

## Failure semantics

- Failure to restore any byte-bound payload is a provenance failure.
- Failure of an SRC4 restart precheck is a numerical/source premise failure.
- Failure of either complete restart common event is an admission failure.
- A runtime crash is invalid and cannot become a scientific negative.
- A later typed health, constraint, source, or boundary stop belongs to the
  calibration result, not to HLT11.

## Explicit nonclaims

HLT11 does not establish a fresh or eligible GR-0 calibration, a trapped
interval, SGB-L health, FGC-QR activation, metric-null defocusing, singularity
resolution, a child domain, retained-EFT validity, dark-sector emergence,
varying locally measured light speed, or a rejection of the gradient route.

Reproduce the authorization before launch with:

```bash
python3 scripts/reproduce_fgc_hlt11_mon11.py
make verify-fgc-pro13-prelaunch
```

The long trajectory is a separate, recoverable command and is not part of the
ordinary repository verifier.
