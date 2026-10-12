# FGC-1-TDG2-FRZ1: prospective absolute-tail discriminator

`FGC-1-TDG2-FRZ1` freezes the next read-only numerical question before the
immutable PROTO13 terminal histories can influence its design. It consumes the
compact `FGC-1-TDG1-PREF16` record at commit `873e13b…`, verifies the ignored
checkpoint's hash when the file is present, but does not load a checkpoint
array, resume a trajectory, or inspect SGB-L or FGC-QR.

## Why another discriminator is necessary

TDG1 established two different facts that must not be merged:

1. Dividing each resolution's temporal tail by its own total power cancels
   amplitude convergence. The inherited normalized cross-resolution ratio is
   therefore not a general spatial-convergence observable for a nonzero
   continuum history.
2. The historical individual temporal-budget failures remain real outputs of
   PROTO13. Their counts decrease with refinement, but TDG1 did not establish
   whether their **absolute** high-frequency power vanishes, approaches a
   nonzero resolved value, or is smaller than the numerical instrument can
   distinguish.

TDG2 addresses only the second question. It neither rewrites PROTO13 nor
defines a successor admission rule.

## Frozen measurement

For every method, tracer, field, and resolution, TDG2 will align the three
`64`-sample histories on their common proper-time intersection and apply the
same compact `1/8` edge taper. It measures the unnormalized one-sided power in
DFT bins `28` through `32`:

```text
P_tail = sum_{k=28}^{32} m_k |F_k|^2 / N^2,

D_tail = sum_{k=28}^{32} omega_k^2 m_k |F_k|^2 / N^2,

N = 64,
m_32 = 1,
m_k = 2 otherwise.
```

No division by total field or derivative power is allowed. Field-tail and
derivative-tail power are classified separately.

The five complex coefficients are independently recomputed with deterministic
Decimal direct DFTs at `72` and `96` decimal digits. The existing NumPy
binary64 rFFT remains an audit implementation. Their disagreement, the
low/high-precision disagreement, and outward binary64 rounding remain public
arithmetic debits; none may be silently discarded.

## Conditional interpolation enclosure

The three histories generally occupy slightly different proper-time grids.
After piecewise-linear alignment, TDG2 reverses that map to the source samples
and records the maximum round-trip discrepancy `delta`.

The analysis may use `delta` as a declared pointwise perturbation radius. For a
window `w`, each top-band coefficient then has the conditional bound

```text
E_k <= ||w||_1 delta.
```

Consequently the public power debits are

```text
Delta P_tail <= sum m_k (2 |F_k| E_k + E_k^2) / N^2,

Delta D_tail <= sum omega_k^2 m_k (2 |F_k| E_k + E_k^2) / N^2.
```

This implication is exact **conditional on** `delta` being used as the
pointwise perturbation radius. The round-trip discrepancy is not claimed to be
a theorem bounding the unknown continuum interpolation error. That limitation
is part of the frozen claim boundary.

## Outcome-neutral classes

Each absolute-power ladder receives exactly one of five descriptive classes:

| Class | Frozen meaning |
|---|---|
| `contracts_to_zero_at_required_power_order` | Both conservative power ratios decrease with order at least `3`, twice the minimum `3/2` history-amplitude order. |
| `nonzero_resolution_independent_within_enclosure` | All three power intervals have a strictly positive common intersection. |
| `converges_to_nonzero_resolved_power` | Same-direction power differences contract with conservative order at least `3/2`, and the declared geometric remainder leaves a strictly positive continuum lower bound. |
| `instrument_floor_or_interpolation_enclosure_dominated` | All coefficient intervals lie below the declared `1e-30` DFT-amplitude floor, the finest power interval reaches zero, or the conditional enclosures prevent a directional distinction. |
| `unresolved_absolute_tail_behavior` | The power is nonmonotone, below the required order, or otherwise not separated by the preceding classes. |

These are diagnosis classes, not pass/fail decisions. In particular, nonzero
resolved temporal content may be legitimate continuum dynamics, a window
effect, or evidence that a future admission rule needs a different physical
observable. A zero-contracting result may identify truncation contamination.
Neither interpretation is selected at the freeze.

## Synthetic preflight

Seven no-history controls must reproduce before the actual diagnosis can run:

```text
exact zero                         -> floor/enclosure dominated
nonzero below 1e-30                -> floor/enclosure dominated
amplitude contracting by 8         -> power contracts toward zero
identical nonzero top-band signal  -> nonzero resolution-independent
nonzero limit plus shrinking error -> converges to nonzero power
nonmonotone amplitude              -> unresolved
jittered-grid interpolation        -> interpolation-enclosure dominated
```

Both field and derivative power must return the expected class. A one-bit
signal mutation must change the absolute evidence. Mutations of the method
ladder, minimum orders, Decimal precision, source hashes, interpolation
nonclaim, or promotion booleans fail closed.

## What the freeze establishes

The preflight establishes only that the absolute discriminator preserves the
amplitude information destroyed by TDG1's normalized ratio and separates the
five declared synthetic cases. It does not establish which case occurs in the
PROTO13 histories.

After this freeze is committed, a separate `FGC-1-TDG2-PREF17` binder may load
the six immutable history ladders read-only, classify every
method/tracer/field/measure record, and preserve the exact evidence. It may not
resume the checkpoint or advance a state.

## Claim boundary

```text
TDG1_terminal_histories_diagnosed = true
TDG2_absolute_tail_discriminator_frozen = true
TDG2_actual_histories_diagnosed = false
source_checkpoint_history_arrays_loaded = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

Reproduce or verify the prospective freeze with:

```bash
make fgc-tdg2-frz1
```
