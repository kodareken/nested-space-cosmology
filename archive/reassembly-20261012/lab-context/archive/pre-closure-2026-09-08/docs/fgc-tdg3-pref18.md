# FGC-1-TDG3-PREF18: native-grid history diagnosis

`FGC-1-TDG3-PREF18` binds the outcome-neutral calculation frozen at immutable
commit `5c8d3cc…`. It restores the exact CAL10 terminal checkpoint read-only,
applies both TDG3 estimators to all six method-owned resolution histories,
advances no state, resumes no trajectory, and reads no SGB-L or FGC-QR data.

## Result in one sentence

Removing common-grid interpolation resolves a substantial subset of the old
temporal failures, but it does not create a valid replacement gate: the
native-grid outcome remains mixed and more than half of the formerly
enclosure-dominated finest-failure classifications remain unresolved or
estimator-dependent.

## Complete native-grid count

TDG3 evaluates `576` method/tracer/field ladders. Each ladder has two power
measures and two independent native estimators, giving `1,152` combined and
`2,304` estimator-specific classifications.

| Method and measure | Zero-contracting | Nonzero | Arithmetic floor | Unresolved | Estimators disagree |
|---|---:|---:|---:|---:|---:|
| RK4 phase field | 59 | 34 | 40 | 147 | 8 |
| RK4 phase derivative | 59 | 35 | 40 | 145 | 9 |
| SSPRK3 phase field | 19 | 47 | 81 | 132 | 9 |
| SSPRK3 phase derivative | 19 | 51 | 81 | 130 | 7 |

The two estimators agree on `1,119/1,152` classifications. That high agreement
is not equivalent to resolution: `554` of those agreements are jointly
unresolved, while `242` are jointly below the declared arithmetic floor.

Every per-estimator nonzero result is a convergent-nonzero native-surrogate
classification. None is produced merely by overlapping resolution-independent
arithmetic intervals.

## What happened to the historical failures

TDG2 had `152` finest-resolution individual failures. Counting its field and
derivative diagnoses separately gives three source groups:

| TDG2 source group | Count | TDG3 zero | TDG3 nonzero | TDG3 floor | TDG3 unresolved | TDG3 disagreement |
|---|---:|---:|---:|---:|---:|---:|
| Zero-contracting | 92 | 92 | 0 | 0 | 0 | 0 |
| Tentatively nonzero | 9 | 6 | 0 | 0 | 3 | 0 |
| Enclosure dominated | 203 | 56 | 15 | 28 | 101 | 3 |

The first row is the strongest surviving signal: all `92/92` TDG2
zero-contracting finest-failure classifications remain zero-contracting under
both no-resampling estimators.

The second row removes an apparent signal rather than preserving it. None of
the nine TDG2 classifications provisionally labeled nonzero under the
conditional interpolation enclosure remains nonzero in TDG3; six become
zero-contracting and three remain unresolved.

The third row isolates the remaining wall. Of the `203` TDG2
enclosure-dominated field/derivative classifications, TDG3 resolves `99` into
zero, nonzero, or floor categories. The remaining `104` are unresolved or
estimator-dependent. Common-grid interpolation therefore obscured real
classification structure, but it was not the sole owner of the historical
failure.

Within the narrower `169` cases whose TDG2 debit was explicitly owned by
common-grid interpolation, TDG3 gives:

```text
50  zero-contracting
15  nonzero native-surrogate power
101 unresolved
 3  estimator disagreement
```

No arithmetic-floor result occurs in this subset because those cases were
selected specifically by an interpolation-owned rather than DFT-floor-owned
TDG2 enclosure.

## Method agreement

RK4 and SSPRK3 return the same combined class for `170/288` field signals and
`167/288` derivative signals. Their shared classes include `19` zero,
`40` floor, `8` field (`10` derivative) nonzero, and `102` field (`98`
derivative) unresolved signals. One field signal is estimator-disagreement in
both methods.

The cross-method mismatch is material. Both methods consume the same 64-event
history shape but have different spatial and temporal discretizations. A
replacement admission cannot simply vote between them or count shared
unresolved labels as evidence of continuum control.

## The deeper limitation exposed by TDG3

TDG3 deliberately measures bins `28..32` from only `64` samples. Those bins
sit at the edge of the sampled phase bandwidth. The exact piecewise-linear
integral and native quadrature are independent calculations of two declared
sampled surrogates, but neither constrains all smooth functions that pass
through the same 64 values.

This is now the smallest unresolved mathematical question:

> Can any finite set of 64 history samples bound the continuum top-band power
> without an explicit regularity, bandlimit, or evolution-residual theorem?

The likely answer is no. Smooth between-sample perturbations can vanish at
every sampled point while carrying independently adjustable high-frequency
power. TDG4 is therefore authorized only as a prospective sampling-
identifiability theorem. It must prove the exact non-identifiability statement,
state the minimal additional assumptions needed for a bound, and attack those
assumptions. It is not another estimator and it does not define a new gate.

If the theorem closes as expected, the present spectral-tail rule should be
retired as a continuum admission. A successor calibration would then need a
directly convergence-controlled proper-time observable or denser stage-level
history with an a priori evolution-error bound.

## Immutable and clean-clone boundary

PREF18 verifies the TDG3 freeze commit and every frozen tracked blob, hashes
the four-file CAL10 raw bundle, restores the terminal histories through the
independent CAL10 loader, recomputes all `576` ladders with the frozen TDG3
module, and cross-tabulates them against the exact PREF17 signal records.

When the ignored raw bundle is present, verification reproduces the complete
diagnosis. A clean clone may verify the compact canonical record and all
tracked hashes without the raw bundle. A partially present raw bundle fails
closed.

## Claim boundary

```text
TDG3_actual_histories_diagnosed = true
TDG3_native_grid_outcome_is_mixed = true
common_grid_interpolation_is_the_sole_temporal_failure_owner = false
piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure = false
replacement_temporal_admission_defined = false
TDG4_sampling_identifiability_theorem_design_may_begin = true
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

Reproduce or verify the result with:

```bash
make fgc-tdg3-pref18
```
