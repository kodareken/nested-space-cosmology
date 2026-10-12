# FGC-1-HLT5-MON5: PROTO7 runtime and one-run authorization

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-HLT5-MON5` implements the transaction contract frozen by
`FGC-1-PRO7-FRZ1` and authorizes exactly one fresh GR-0 calibration campaign
under `FGC-1-CAL4-RUN1-PLAN`.

```text
PROTO7_successor_runtime_compositor_implemented = true
PROTO7_fresh_GR0_dynamic_calibration_authorized = true
PROTO7_resolved_holdout_manifest_authorized = false
classical_spherical_diagnostic_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

The authorization is tied to immutable PROTO7 checkpoint
`ec8ac628528b822c600de35eb2ca6b2b49506051`. It reads the protocol, freeze,
HLT4 input certificate, and ID2 static ledger from that commit and verifies
their exact hashes. It reads no PROTO7 trajectory and creates no run
directory.

## Implemented transaction

Each proposed step has one accepted boundary and one unaccepted transaction.
The raw source gate is first evaluated on the bitwise last accepted state. A
failure there remains terminal and is latched as `newton_residual_limit`.

If that check passes, all Runge--Kutta evaluations and the independently
evaluated candidate endpoint are constructed as a proposal. The complete
nine-stop transaction is previewed before classification:

- if every failure is `newton_residual_limit`, the proposal is unaccepted,
  durably serialized, discarded, and retried at half the step size;
- if any non-source premise fails, retry is forbidden and the unchanged
  chronological and within-stage priority selects the terminal reason;
- if no premise fails, external tracer state is previewed and the proposal is
committed once through the inherited transaction guard.

Any non-source failure therefore remains a hard retry veto.

The retry factor `1/2`, maximum of 32 source retries per accepted step,
minimum step `1/1073741824`, raw residual limit `1e-12`, affine source solver,
kinetic condition gate, CFL rule, and all physical/numerical stops are
unchanged from PROTO6.

## Complete failure observability

Every rejected proposal binds:

- amplitude, method, grid size, member, target common event, initial step,
  attempted step, local retry count, and cumulative member retry count;
- accepted-state field hash, accepted time, accepted step index, and accepted
  transaction serial;
- every failed evaluation's stage name, time, transaction serial, raw source
  residual, raw threshold, and complete failure tuple;
- the union of all failures, retry classification, selected terminal reason,
  and whether a non-source failure vetoed retry;
- exact preservation of fields, time, accepted-stage count, transaction
  serial, causal debit, tracer state, and preaccept state.

The append-and-fsync event sink must return before a retry is proposed.
Retry exhaustion retains the last complete rejected-proposal record and its
frozen exhaustion condition. A terminal stop restores all six grid members to
the last fully completed common event, verifies that restoration, and embeds
the complete stop evidence in the campaign result.

## Adversarial controls

The machine certificate injects a source-only failure independently at RK4
`k1`, `k2`, `k3`, `k4`, and `candidate_endpoint`; all five must be retryable
without changing the transaction or causal state. It then injects a
simultaneous lapse failure at every one of those positions; all five must veto
retry and preserve `nonpositive_lapse` as the inherited priority.

A mixed control places a source miss earlier than a later non-source miss. The
later failure vetoes retry, while the earlier source miss remains the selected
terminal reason under chronological ownership. A healthy control confirms
that external preview occurs before exactly one commit. Runner-level controls
exercise durable endpoint retry, typed retry exhaustion, and cross-member
snapshot restoration.

## Frozen inputs and fresh namespace

The twelve GR-0 inputs retain amplitudes `[5/2, 3]`, methods `[RK4, SSPRK3]`,
resolutions `[1025, 2049, 4097]`, and the exact projected state hashes from
immutable HLT4. All twelve accepted-state source prechecks and all four
method/amplitude `t=0` common-event admissions must pass before authorization.
Only the PROTO7 ownership label, complete-observability bit, run-plan hash,
protocol hash, and fresh expanded configuration hash change.

The calibration and holdout roots are respectively
`runs/fgc-2-sf1/proto7/calibration` and
`runs/fgc-2-sf1/proto7/holdout`. Both must be absent or empty when HLT5 is
created. HLT5 creates neither. The authorized runner refuses a pre-existing
calibration root, tracked implementation drift, a dirty tracked worktree, or
a manifest/checkpoint identity mismatch. PROTO6 output is immutable diagnostic
history and cannot be relabelled as PROTO7 evidence.

## Epistemic boundary

This certificate proves that the repaired numerical transaction is
implemented and that the unchanged GR-0 inputs are eligible to run once. It
is not the run. It establishes no dynamically formed trapped sphere, no
SGB-L comparison, no FGC-QR health or evolution, no positive Raychaudhuri
margin, and no physical transition. It does not derive singularity
resolution, a child domain, a dark-sector mechanism, or a varying locally
measured speed of light, and it rejects neither FGC-QR nor the general
gradient route.

## Reproduction and authorized launch

Before launch:

```bash
python3 scripts/reproduce_fgc_hlt5_mon5.py \
  --output results/fgc-1-hlt5-mon5.json
python3 -m unittest tests.test_fgc_proto7_runtime \
  tests.test_fgc_gr0_campaign_runner_v7 \
  tests.test_fgc_hlt5_mon5_reproduction -v
python3 scripts/run_fgc_gr0_calibration_v7.py \
  --check-authorization-only
```

After the complete HLT5 unit is committed, the one authorized campaign is:

```bash
python3 scripts/run_fgc_gr0_calibration_v7.py
```
