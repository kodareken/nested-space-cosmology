# FGC-1-TDG9-UR1-FRZ1: retry-3 `u:R` envelope-owner diagnostic

## Result

**FGC-1-TDG9-UR1-FRZ1 is prospective and unexecuted.** It authorizes only a
future direct-successor commit to run one read-only retry-3 SSPRK3-on-inherited-
SBP4 diagnostic for the single `u:R` channel. The prelaunch certificate binds
sealed HEAD `22f7dc5e181af45eb43f6f19dee23fbc29f89c0c`, the completed AC1-PREF1
compact result, both canonical AC1 raw leaves, the exact replay receipt, the
115-leaf campaign store, and the existing TI2 retry-3 `u:R` row identity.

Prelaunch restores and fingerprints retry 3 without preparing a compositor or
constructing a shadow proposal. It requires the future UR1 raw namespace to be
absent and leaves both that namespace and the campaign store untouched.

## Frozen question

After the exact authority image is committed as the direct successor,
the separate RUN1 may build the same seven shadow paths and 28 SSPRK3 stage and
endpoint records used by AC1. It publishes only the `u:R` diagnostic:

- the exact 2,044-row binary64 Hermite stream is hashed in domain
  `TDG9-AR1-BINARY64-HERMITE-ROWS-v1` and compared with
  `b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969`;
- D01 must replay as
  `[0x0.0p+0, 0x1.2d198e246e459p-38]`, D12 as
  `[0x0.0p+0, 0x1.2d1c31bb91376p-38]`, and their design-only classifier as
  `order_inconclusive`;
- every public `Binary64CubicEnvelope` field is serialized for D01 and D12;
- the lower bound is independently reconstructed as
  `max(0, Fraction(raw) - Fraction(coefficient debit) - Fraction(arithmetic debit))`
  followed by downward binary64 rounding;
- Bernstein certification slack is recorded but is not included in that lower
  clip; and
- each zero lower is assigned exactly one closed owner:
  `raw_candidate_maximum_is_zero`,
  `coefficient_plus_arithmetic_debit_clips_positive_raw_maximum`, or
  `downward_binary64_rounding_of_positive_exact_lower`.

Only after the row hash matches may the raw terminal record the related TI2
radius-free complete-C class `sufficient_contraction_pass`. That historical
class is explicitly not a replacement production admission and does not alter
AC1's `order_inconclusive` production-classifier result.

## Bounded terminals and raw contract

RUN1 may publish only canonical `manifest.json` and `terminal.json` under
`runs/fgc-2-sf1/tdg9-ur1/retry3-ur-envelope-owner`, both with schema
`UR1-raw-v1`. Its closed terminal classes are:

- `completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity`;
- `completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ`;
- `row_stream_identity_mismatch`;
- `replayed_production_intervals_differ_from_sealed_AC1`;
- `shadow_proposal_premise_stop`; or
- `invalid_provenance_or_implementation`.

`invalid_provenance_or_implementation` is fail-closed console output and is not
publishable raw evidence. Status authenticates an already published canonical
terminal; a run requires the namespace to be absent.

No temporal-admission call, fine-path commit, PDE-state write, campaign-store
write, continuation, calibration, candidate branch, mechanism result, EFT
retention, transition, or physical claim is authorized. A matching completed
owner result licenses only design of a later, separately prospectively frozen
envelope or admission change. It does not itself authorize such a change.

## Verification

Compact verification is store-, raw-, and shadow-blind:

```bash
python scripts/reproduce_fgc_tdg9_ur1_frz1.py --verify-compact
```

The bounded live prelaunch writes only the tracked compact certificate while
restoring no proposal:

```bash
python scripts/reproduce_fgc_tdg9_ur1_frz1.py --write
```

Neither path invokes RUN1. A later committed authority requires an explicit
`--authority-commit`; the runner remains deliberately separate from prelaunch.
