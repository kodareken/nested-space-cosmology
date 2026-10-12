# FGC-1-TDG1-PREF16: terminal-history diagnosis

`FGC-1-TDG1-PREF16` binds the read-only diagnosis authorized by
FGC-1-TDG1-FRZ1 at immutable commit `e8310bd…`. It consumes the exact CAL10
terminal checkpoint, advances no state, resumes no trajectory, and reads no
SGB-L or FGC-QR data.

## What the actual histories establish

The synthetic result survives contact with the campaign data: the normalized
raw cross-resolution tail ratio is structurally invalid as a general spatial-
convergence observable. That conclusion is mathematical and does not depend
on selecting a favorable subset of the terminal histories.

The remaining individual temporal-budget layer cannot yet be dismissed. Of
the `1,728` signal-resolution records (`2` methods times `3` resolutions times
`48` tracers times `6` fields), `544` fail the inherited individual budget:

```text
149  have FFT peaks at or below the declared absolute 1e-30 floor
395  remain above that floor
```

The failure count contracts with refinement on both method-owned ladders:

```text
RK4:    110 -> 106 -> 93
SSPRK3:  96 ->  80 -> 59
```

That is evidence against interpreting every failure as a fixed continuum
high-frequency obstruction, but it is not enough to pass the layer. Some
signals change classification after mean removal or affine detrending; others
fail every frozen view. Downsampling changes only a small subset. These views
separate possible window, trend, sampling, arithmetic-floor, and resolved-
dynamics effects; none is a replacement admission.

## Raw history-difference evidence

On the common proper-time intersection, the unnormalized history differences
give the following `288` classifications per method:

| Method | Exact identical | Fine pair exact | Order at least `3/2` | Below `3/2` | Nonmonotone | Interpolation-limited |
|---|---:|---:|---:|---:|---:|---:|
| RK4 | 15 | 8 | 237 | 28 | 0 | 83 |
| SSPRK3 | 33 | 1 | 233 | 21 | 0 | 51 |

Most histories therefore contract spatially, but complete admission does not
follow: `49` signals remain below the diagnostic order and `134` comparisons
are no larger than their public interpolation debit. The history-difference
view was prospectively declared diagnostic only, so even a perfect count could
not retroactively pass PROTO13.

## Scientific interpretation

PREF16 separates a proved instrument defect from a still-open numerical
question:

1. The normalized nested-ratio layer cannot be retained as a general spatial-
   convergence test because its normalization cancels amplitude convergence.
2. The actual individual-budget failures are real outputs of the historical
   rule. Their count decreases with resolution, but `395` above-floor failures
   and `49` below-order history differences prevent an outcome-neutral pass.
3. The missing discriminator is absolute, unnormalized tail convergence with
   explicit arithmetic and interpolation enclosure. It must determine whether
   the offending tail power contracts toward zero, converges to resolved
   physical power, or remains below the instrument floor.

This authorizes design of **TDG2**, a no-trajectory absolute-tail discriminator.
It does not yet authorize a replacement temporal gate, PROTO14, GR-0
eligibility, SGB-L, FGC-QR, COL1, DEF1, retained-EFT evolution, or a physical
claim.

## Immutable and clean-clone boundary

The exact raw bundle is verified under CPython `3.14.3`, NumPy `2.5.1`, and
Accelerate on Darwin `arm64`. A clean clone may verify the compact tracked
PREF16 record without the ignored raw bundle. A partially present raw bundle
or any locally present file with the wrong hash fails closed.

The terminal checkpoint may not resume. A later successor may consume it only
under a new prospective contract, and no TDG2 calculation may be read before
that contract is committed.

## Claim boundary

```text
TDG1_terminal_histories_diagnosed = true
normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test = false
individual_temporal_budget_failure_cause_fully_derived = false
TDG2_absolute_tail_discriminator_design_may_begin = true
TDG2_frozen = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
GR0_case_eligible = false
FGCQR_holdout_execution_authorized = false
```

Reproduce or verify the result with:

```bash
make fgc-tdg1-pref16
```
