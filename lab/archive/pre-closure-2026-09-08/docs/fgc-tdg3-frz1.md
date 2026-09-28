# FGC-1-TDG3-FRZ1: prospective native-grid discriminator

`FGC-1-TDG3-FRZ1` freezes the next numerical question after the mixed TDG2
result at immutable commit `8585dbe…`. It verifies the compact PREF17 evidence
and may hash the ignored terminal checkpoint, but it does not load a history
array, resume a trajectory, inspect SGB-L or FGC-QR, or define a replacement
temporal admission.

## The question TDG3 isolates

TDG2 independently cleared binary64 DFT arithmetic as the owner of every
finest historical failure. Its remaining uncertainty is narrower: most
unresolved cases are dominated by the conditional debit assigned when three
slightly different proper-time grids are interpolated onto one common uniform
grid.

TDG3 asks:

> If each 64-sample history is measured directly on its own native proper-time
> grid, without interpolating a value onto a common grid, do two independently
> defined native estimators give the same zero-contracting, nonzero, floor, or
> unresolved classification?

This is an instrument-ownership question. A favorable answer would not make
the historical PROTO13 gate valid, and it would not turn a sampled surrogate
into the unknown continuum history.

## Native phase and retained band

Each history maps its own increasing proper times to

```text
x_j = (tau_j - tau_0) / (tau_63 - tau_0),    0 <= x_j <= 1.
```

The exact native spacings and proper-time span remain public. No common overlap
or common-grid value resampling is performed. The same `1/8` compact window is
evaluated at the native `x_j`. TDG3 measures bins `28..32` in normalized phase,
with separate field-shape and phase-derivative power:

```text
P_tail = 2 sum_k |C_k|^2,
D_tail = 2 sum_k (2 pi k)^2 |C_k|^2.
```

These are dimensionless shape diagnostics. They are deliberately not called a
physical proper-frequency admission. A future temporal gate would still need
an operational proper-time observable and its own error theorem.

## Two no-resampling estimators

The primary estimator defines `y_PL(x)` as the piecewise-linear interpolant of
the native compact-windowed samples and evaluates every segment integral in
closed form:

```text
C_k^PL = integral_0^1 y_PL(x) exp(-2 pi i k x) dx.
```

The comparator uses positive native Voronoi/trapezoid weights `q_j`, whose sum
is one:

```text
C_k^Q = sum_j q_j y_j exp(-2 pi i k x_j).
```

Neither route constructs a common grid. The first chooses a piecewise-linear
surrogate between samples; the second chooses a native quadrature rule. Their
spread is recorded but is not mislabeled as an error enclosure. Agreement is
required before TDG3 emits a combined zero, nonzero, or floor label.

## Arithmetic contract

Both coefficient sets are recomputed from exact Decimal conversions of the
binary64 phase and weighted samples at `72` and `96` decimal digits. Separate
binary64 implementations use closed-form segment sums and native weighted
sums. Low/high precision disagreement, binary64/reference disagreement, and
outward rounding define the only interval debit.

This arithmetic interval encloses calculation of the declared sampled
surrogate. It does not enclose the unknown continuum signal or all possible
between-sample reconstructions.

## Outcome-neutral classes

Each estimator first uses five power classes:

```text
contracts_to_zero_at_required_power_order
nonzero_resolution_independent_within_arithmetic_enclosure
converges_to_nonzero_native_power
arithmetic_floor_dominated
unresolved_native_power_behavior
```

The two routes then return exactly one combined class:

```text
native_estimators_agree_zero_contracting
native_estimators_agree_nonzero
native_estimators_agree_floor_dominated
native_estimators_agree_unresolved
native_estimators_disagree
```

The two nonzero estimator classes form one combined nonzero family. Any other
cross-estimator mismatch remains disagreement; it is never rounded, voted, or
tuned away. No combined class is a replacement admission.

## Synthetic preflight

Seven no-history controls must pass for both power measures:

```text
exact zero                         -> both estimators floor dominated
nonzero below the normalized floor -> both estimators floor dominated
amplitude contracting by eight     -> both estimators zero contracting
identical nonzero native history   -> both estimators nonzero
nonzero limit plus shrinking error -> both estimators nonzero
nonmonotone amplitude              -> both estimators unresolved
jittered native grids              -> TDG2 interpolation owner, TDG3 nonzero
```

The last fixture is a discriminator control, not campaign evidence: it proves
that a known continuous signal can be rendered enclosure dominated by TDG2's
common-grid round trip while both no-resampling TDG3 estimators retain its
nonzero classification.

Mutations of a signal bit, estimator, band, minimum orders, Decimal precision,
source hash, no-resampling/noncontinuum statements, or claim booleans must fail
closed.

## What the freeze establishes

The freeze establishes only that two independently defined native-grid
estimators can distinguish the declared synthetic alternatives without
common-grid resampling. It does not establish which class occurs in any real
PROTO13 history.

After this freeze is committed, a separate TDG3 binder may load the immutable
six-history ladder read-only and cross-tabulate TDG2's enclosure-dominated
cases against both native estimators. It may not resume the checkpoint or
advance a state.

## Claim boundary

```text
TDG2_actual_histories_diagnosed = true
TDG3_interpolation_ownership_discriminator_frozen = true
TDG3_actual_histories_diagnosed = false
common_grid_interpolation_removed_from_TDG3_measurement = true
piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure = false
individual_temporal_budget_failure_cause_fully_derived = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

Reproduce or verify the prospective freeze with:

```bash
make fgc-tdg3-frz1
```
