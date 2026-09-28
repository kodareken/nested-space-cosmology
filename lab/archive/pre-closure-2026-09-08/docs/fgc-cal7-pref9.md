# FGC-1-CAL7-PREF9: PROTO10 centre-roundoff diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration numerical-source diagnosis;
no mechanism result

```text
PROTO10_campaign_terminated_normally = true
PROTO10_GR0_case_eligible = false
PROTO10_center_adjacent_binary64_source_floor_localized = true
PROTO11_well_balanced_reference_map_required = true
FGCQR_holdout_execution_authorized = false
```

## What the fresh campaign established

HLT8 authorized exactly one fresh GR-0 calibration under PROTO10. The campaign
ran from immutable commit `9c8aac20954c47b1c844c847fda91e2c423d6c85`,
terminated normally after about 3,158 seconds, and wrote 264 canonical events
under campaign identity
`3601547df421a6920d4b9b95d48a9d8ce8f785465f30e9fba9291cc93a510e69`.
The event log contains two guarded `t=0` common events and 262 wholly
unaccepted source-only proposal records.

Both amplitudes pass the complete guarded PROTO10 source, constraint,
absolute-spectrum, direct-coarse, conditioning, and whole-profile gate set at
`t=0`. Both raw PROTO9 spectral failures remain public. Neither amplitude is
initially trapped.

Both amplitudes then stop in the same member before the first evolved common
event:

```text
member                 RK4-8193
classification         scientific_source_retry_exhausted
reason                 newton_residual_limit
exhaustion             minimum_step_size
last common event      t = 0
target common event    t = 1/16
accepted step index    15
last rejected step     about 1.5522e-9
next forbidden step    about 7.7610e-10
frozen minimum step    about 9.3132e-10
```

The `5/2` and `3` terminal accepted-member times are respectively
`0.004389551467835599` and `0.004389551467678107`, an absolute difference of
about `1.57e-13`. Their last candidate source residuals are respectively
`1.3592093373108497e-12` and `1.3592093372686369e-12`, against the unchanged
raw ceiling `1e-12`. All rejected proposals preserve fields, time, accepted
stage count, causal debit, and external transaction state. No non-source stop
is relabelled as retryable, and the final checkpoint verifies cross-member
rollback to the last complete common event.

The campaign result is therefore exactly:

```text
classification = calibration_failed_no_eligible_GR0_case
selected_amplitude = null
holdout_execution_authorized = false
SGBL_outcome_read = false
FGCQR_outcome_read = false
mechanism_question_answered = false
```

## Where the source floor lives

CAL7/PREF9 reconstructs all twelve frozen HLT8 initial states and exactly
recovers every accepted-state source residual. For the primary fourth-order
RK4 representation, the maximum is always the `alpha` equation at the first
positive grid point. Both amplitudes give the same values:

```text
points   radius       raw residual
2049     1/16         4.9145872556740196e-14
4097     1/32         1.9658349022696078e-13
8193     1/64         7.8633396090784959e-13
```

The residual multiplies by four when the spacing halves. More specifically,

```text
abs(residual) * radius^2 / binary64_epsilon = 83/96
```

to the serialized arithmetic accuracy on all three grids. The maximum is in
the exact centre vacuum buffer, not in the matter pulse, and it is identical
for the two matter amplitudes. The scaling is therefore the signature of a
centre-adjacent binary64 cancellation amplified by the spherical `1/r^2`
terms. It is not evidence that the continuum GR-0 acceleration root, the
FGC-QR restoring branch, or a gradient mechanism fails.

The numerical origin is explicit. PROTO10 projects the auxiliary first-order
field with the frozen SBP derivative of the full ADM state. The physical
reference state is

```text
u_ref = (1, 0, 1, r, 0, 0)
q_ref = (0, 0, 0, 1, 0, 0).
```

Although these fields are analytically constant or linear, floating SBP
boundary-coefficient cancellation leaves tiny nonzero first and second
derivatives in the first centre-adjacent rows. The annular spherical source
then amplifies that noise as the grid is refined. At `8193` points the initial
raw residual is already about `78.6%` of the absolute source ceiling, leaving
too little representational headroom for a finite evolution step.

## Outcome-neutral well-balanced control

The certificate recomputes the same initial data after differentiating only
the departure from the exact reference:

```text
q   = D_h (u - u_ref) + q_ref
q_r = D_h (q - q_ref).
```

No equation, physical input, solver, residual, or threshold is changed in this
control. On the primary RK4 grids the source residuals become approximately

```text
2049     2.9532e-14
4097     3.2752e-14
8193     7.5791e-14
```

and every method, amplitude, and resolution remains below the unchanged
`1e-12` gate. The fine primary residual falls by a factor of about `10.375`.
This is a `t=0` diagnostic control only. Its projected state hashes differ
from PROTO10, so it cannot be retroactively substituted into the completed
campaign or described as a successful trajectory.

## Required successor boundary

The smallest evidence-justified successor is a separately frozen
**FGC-2-SF1-PROTO11** numerical-map repair. It may change only the
semidiscrete reference-state differentiation:

- differentiate `u-u_ref` and restore `q_ref` analytically;
- differentiate `q-q_ref` when constructing `q_r`;
- rebuild every projected input hash under the new map;
- prove exact Minkowski preservation on every frozen grid;
- retain nontrivial generic REF1/source controls so the reference subtraction
  cannot hide a physical perturbation;
- retain the equations, source linear solve, raw `1e-12` residual gate,
  physical data, methods, CFL rule, half-step retry rule, 32-retry budget,
  minimum step, constraint/spectrum/health/boundary gates, and observables;
- use a fresh namespace and a fresh unread campaign.

PROTO11 is required here but is not yet frozen. CAL7/PREF9 does not authorize
another run.

## Scientific boundary

This is a stronger result than an unexplained failed run: it identifies a
specific numerical premise, its exact grid location, its `epsilon/r^2`
scaling law, its amplitude independence, and an outcome-neutral control that
removes the scaling without relaxing the raw threshold. It is still not a
physical or mechanism result.

CAL7/PREF9 does not select an amplitude, complete GR-0 calibration, observe a
dynamically trapped sphere, read SGB-L or FGC-QR, test regulator activation,
derive a Raychaudhuri margin, reject FGC-QR or gradients, validate a retained
EFT, resolve a singularity, derive a child domain or dark sector, or vary a
locally measured speed of light. A new campaign must earn every one of those
later steps.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal7_pref9.py --check
```
