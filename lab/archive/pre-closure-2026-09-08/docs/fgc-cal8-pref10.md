# FGC-1-CAL8-PREF10: PROTO11 evolved-run diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration numerical diagnosis; no
mechanism result

```text
PROTO11_campaign_terminated_normally = true
PROTO11_both_amplitudes_reached_evolved_common_event = true
PROTO11_GR0_case_eligible = false
PROTO11_evolved_affine_source_arithmetic_obstruction_localized = true
PROTO11_direct_coarse_phi_derivative_veto_localized = true
FGCQR_holdout_execution_authorized = false
```

## What the fresh campaign established

HLT9 authorized exactly one fresh GR-0 calibration under PROTO11. The campaign
ran from immutable commit
`199c86e55a7504e780ccce5a840648528482f819`, terminated normally after about
5,217 seconds, and wrote 168 canonical records under campaign identity
`63ea00c3a781d3a024c69d9841cadc74bf994cd7ed3750bc12cd29bbf24c6548`.
The ledger contains four synchronized common events and 164 wholly unaccepted
source-only proposals. Both amplitudes reached the first genuinely evolved
event at `t=1/16`; neither produced a qualified trapped event.

The terminal campaign classification is exactly:

```text
classification = calibration_failed_no_eligible_GR0_case
selected_amplitude = null
holdout_execution_authorized = false
SGBL_outcome_read = false
FGCQR_outcome_read = false
mechanism_question_answered = false
```

This is not a failed physical prediction. The candidate branch was never opened.
It is a failed GR-0 calibration contract with two separately localized causes.

## Amplitude 5/2: evolved affine source-arithmetic wall

At `t=1/16`, all six amplitude-`5/2` members passed the complete constraint and
spatial admissions. The finest trapped scores from RK4 and SSPRK3 were about
`-0.0771108564` and `-0.0771108579`, so the event was not trapped. The temporal
spectrum was honestly unavailable because the frozen 64-sample minimum had not
yet been reached.

While advancing toward `t=1/8`, only `RK4-8193` encountered source failures.
CAL8 proves directly from the immutable event ledger that:

- all 164 rejected proposals contain only `newton_residual_limit`;
- fields, time, accepted stages, causal debit, and external transaction state
  are preserved on every rejection;
- no non-source veto is relabelled as retryable;
- every retry group is an exact binary half-step sequence;
- between groups, the inferred accepted increment is the next half-step within
  two floating-time ULPs; and
- cross-member rollback returns to the last complete event.

The rejections occur in 16 accepted-time groups with retry counts
`1,2,3,4,5,6,7,8,12,13,14,15,17,18,19,20`. The last accepted state is at
coordinate time `0.0706352245695038`, step 97, transaction 485. Its last legal
failed step is about `1.55220413173635e-9`; the next required half-step is
about `7.76102065868175e-10`, below the unchanged minimum
`2^-30`. The last complete raw source residual is about
`1.0659145473163184e-12` against the unchanged `1e-12` ceiling.

This localizes an evolved affine source-evaluation arithmetic problem. It does
not demonstrate a continuum source inconsistency, branch loss, or failure of
FGC-QR. The present campaign does not contain the complete affine matrix and
constant term at the wall, so it cannot yet distinguish conditioning from the
specific binary64 solve/residual evaluation path.

## Amplitude 3: independent direct spectral veto

Amplitude `3` reached `t=1/16` with no source retry. Both methods passed their
constraint gates, finest-pair individual spectrum budgets, and whole-profile
contraction requirements. All direct coarse-to-medium field tails and every
derivative tail except `phi/Lambda` also passed.

The unchanged direct ceiling is `<1/4`. The serialized derivative-tail ratios
are:

```text
method    coarse-to-medium     medium-to-fine
RK4       0.326717281113685    0.22812049822366257
SSPRK3    0.26675704013240803  0.2142551901915683
```

Thus both finest pairs pass while both direct coarse pairs fail. PROTO11 makes
the coarse pair a public veto, so neither apparent convergence nor a diagnostic
map round trip can saturate or erase it. This could be a pre-asymptotic coarse
grid or a persistent resolution problem; one event cannot decide which.

## Required successor boundary

The smallest justified next gate is **FGC-1-SRC2-PREF11**, not PROTO12 and not
a threshold change. SRC2 must replay and capture the evolved amplitude-`5/2`
source wall, binding the complete affine matrix, constant term, state, raw
residual authority, and branch identity. It may compare stable arithmetic
evaluations only while preserving:

- the unredefined continuum equations and PROTO11 reference map;
- the complete raw residual and its `1e-12` ceiling;
- amplitudes, methods, grids, CFL rule, retry factor and budget, and minimum
  step;
- all constraint, health, boundary, spectrum, and observable rules; and
- the direct coarse-to-medium `<1/4` veto.

Only a nontrivial arithmetic control that passes those unchanged conditions
could justify freezing a fresh PROTO12 namespace. The amplitude-`3` spectral
failure remains binding and may not be renamed or tuned away.

## Scientific boundary

CAL8/PREF10 converts a long run into useful, falsifiable information: PROTO11
can evolve both controls beyond initial data, but no GR-0 amplitude satisfies
the complete calibration contract. It identifies one reproducible source
arithmetic wall and one independent resolution-spectrum veto. It changes no
threshold and authorizes no trajectory.

It does not select a GR-0 amplitude, open SGB-L or FGC-QR, test regulator
activation, derive a complete Raychaudhuri margin, reject FGC-QR or the general
gradient route, validate a retained EFT, resolve a singularity, derive a child
domain or dark sector, or vary a locally measured speed of light.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal8_pref10.py --check
```
