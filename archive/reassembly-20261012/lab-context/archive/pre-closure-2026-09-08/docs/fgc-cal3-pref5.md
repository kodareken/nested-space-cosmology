# FGC-1-CAL3-PREF5: PROTO6 terminal-proposal diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration transaction diagnosis; no
mechanism result

```text
PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed = true
PROTO7_general_transaction_revision_required = true
fresh_GR0_dynamic_calibration_completed = false
FGCQR_holdout_execution_authorized = false
```

## What the frozen campaign established

HLT4 authorized exactly one fresh GR-0 calibration under PROTO6. The run was
launched from immutable commit
`f8412dd202886dc60650acb887721be3fc62a737` and terminated normally after
about 581 seconds. Both amplitudes passed the complete primary/comparator
common event at `t=0`; neither completed the first nonzero common event at
`t=1/16`. The runner returned its frozen classification,
`calibration_failed_no_eligible_GR0_case`, and selected no amplitude.

Both amplitudes produced one public, exactly rolled-back source retry in the
finest RK4 member near `t=0.0594075`. The retry record names two failed
internal stages, `rk4_k2` and `rk4_k3`, with residuals close to
`1.364e-12` against the unchanged raw `1e-12` limit. A later
`newton_residual_limit` remained terminal, so all six members were restored to
the last complete cross-member event at `t=0`.

The manifest, event log, terminal checkpoint, and result are bound here by
SHA-256 before interpretation. No PROTO6 output is relabelled or overwritten.

## Exact terminal replay

The terminal checkpoint intentionally contains the synchronized rollback,
not a mixed-time state near the stop. CAL3 therefore reconstructs the frozen
amplitude-`3`, `RK4-4097` member from the exact HLT4 `t=0` state and replays it
deterministically to the first stop. The public retry record is reproduced
field for field.

At the last accepted state,

```text
t                         = 0.058593708693194996
raw source residual       = 5.370340537498194e-13
raw limit                 = 1.000000000000000e-12
accepted-state gate       = pass
```

The next source retry halves the proposal. In that smaller—but still wholly
unaccepted—proposal, the source gate passes at the accepted-state precheck,
`rk4_k1`, `rk4_k2`, and `rk4_k3`. It misses only at:

```text
rk4_k4                    = 1.3639071551646946e-12
candidate_endpoint        = 1.3639179325880144e-12
```

No non-source premise fails. The lapse remains positive, the kinetic
condition number remains far below `1e10`, the unresolved relative
acceleration correction is far below `1e-12`, and the maximum residual is
again localized at the first annular point, `r=1/32`.

PROTO6 asks whether *all* failures in a proposal are source-only failures at
internal RK stages. Because the candidate endpoint also misses, that answer
is false. The unchanged PROTO5 transaction is then invoked and latches the
earlier `rk4_k4` source miss as terminal. The candidate endpoint has not been
accepted, but its stage label prevents a second reduction.

CAL3 restores the bitwise last accepted state and evaluates one additional
factor-`1/2` proposal. Every source call and the complete transaction pass the
same unchanged gates. The experiment therefore distinguishes a rejectable
unaccepted proposal from an accepted-state continuum failure.

## Rethinking decision

PROTO5 treated every source miss as terminal. PROTO6 added a list of internal
stage names that may retry but excluded the candidate endpoint. Both rules
classify numerical validity by integrator labels. The stable distinction is
instead transactional:

```text
accepted state
    source miss -> terminal premise failure

unaccepted proposal (all RK stages plus candidate endpoint)
    source-only miss -> rollback, reduce, retry
    any non-source miss -> terminal under the original priority
```

A candidate endpoint is not physical or accepted merely because it is the
last evaluation in a proposal. It becomes accepted only after the complete
transaction commits. This rule removes the whole stage-specific exception
class rather than adding another endpoint patch.

The alternatives were rejected for explicit reasons:

- raising the raw tolerance or replacing it with the already-small relative
  correction after reading the outcome would change the numerical premise;
- changing precision, solver, resolution, or physical inputs would confound
  ownership with a new experiment;
- declaring GR-0 or FGC-QR rejected would assign physical meaning to states
  the integrator never accepted;
- retrying only candidate endpoints would preserve the same brittle
  stage-label design.

The raw residual and every physical, constraint, spectral, trapped-sign,
boundary, and retry-budget rule therefore remain unchanged. Relative
correction and conditioning are serialized as diagnostics, not substituted as
admission gates.

## Required successor boundary

CAL3 permits only a separately frozen PROTO7 design:

- evaluate and retain the accepted-state source precheck as terminal;
- preview the complete proposal before classification;
- allow any source-only miss anywhere in an unaccepted proposal to restore the
  exact accepted state, halve, and retry;
- forbid retry when any non-source failure occurs;
- retain the existing factor `1/2`, 32-retry limit, minimum step `2^-30`, raw
  `1e-12` gate, equations, data, methods, grids, and observables;
- serialize terminal and retry-exhaustion evidence with member, stage,
  residual, accepted-state hash, and rollback facts;
- use a new protocol version and fresh namespace.

PROTO7 must be frozen and independently implemented before another campaign.
Nothing in CAL3 authorizes SGB-L, FGC-QR, DEF1, retained-EFT promotion, or a
physical-transition claim.

## Scientific boundary

This is a stronger test-contract result than CAL2 because it identifies the
general ownership invariant and closes the stage-label formulation, but it is
still not the tangible physical result sought by FGC-2-SF1. It establishes no
collapse, trapped sphere, regulator activation, defocusing interval,
singularity resolution, child domain, dark-sector mechanism, or varying local
speed of light. The mechanism question remains unanswered.

## Reproduction

Fast validation of the immutable campaign and stored result:

```bash
python3 scripts/reproduce_fgc_cal3_pref5.py --check
```

Full reconstruction and deterministic finest-member replay:

```bash
python3 scripts/reproduce_fgc_cal3_pref5.py --check --replay
```

Generating the tracked certificate without an existing result also performs
the full replay:

```bash
python3 scripts/reproduce_fgc_cal3_pref5.py \
  --config configs/fgc/fgc-1-cal3-pref5.toml \
  --output results/fgc-1-cal3-pref5.json
```

Raw checkpoints remain outside Git under the repository's declared large-run
data boundary; their exact hashes are part of the tracked certificate.
