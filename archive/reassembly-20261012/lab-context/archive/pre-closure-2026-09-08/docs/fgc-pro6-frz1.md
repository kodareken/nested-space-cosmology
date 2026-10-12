# FGC-1-PRO6-FRZ1: outcome-neutral PROTO6 freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-PRO6-FRZ1` freezes **FGC-2-SF1-PROTO6** as the smallest successor to
PROTO5's obstructed source-stop ownership. It is tied to immutable checkpoint
`47f4c83bae2df26c51c3544eed8c06507a0cd9eb`, which contains both the exact
PROTO5 protocol and the CAL2-PREF4 diagnosis.

```text
PROTO6_outcome_neutral_protocol_frozen = true
PROTO6_immutable_lineage_verified = true
PROTO6_accepted_state_source_terminal_contract_frozen = true
PROTO6_internal_trial_source_retry_contract_frozen = true
PROTO6_successor_runtime_compositor_implemented = false
PROTO6_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

This is a numerical premise repair, not a calibration, collapse, or mechanism
result.

## The one changed ownership rule

Before every proposal, the source gate must be evaluated on the last accepted
state. If that accepted state misses the unchanged raw `1e-12` residual gate,
the failure remains terminal under the original `newton_residual_limit` label.

If the accepted state passes but an internal stage of the still-unaccepted
proposal misses that same gate, the proposal may be rejected and retried from
the bitwise last accepted state. The retry uses only PROTO5's already frozen
factor `1/2`, maximum of 32 retries, and minimum step size `2^-30`. The failed
trial advances no field, time, constraint ledger, causal debit, or scientific
classification. Every rejected step size, stage name, residual, and retry
count remains public evidence.

No other failure becomes retryable. In particular, a singular kinetic block,
condition-limit failure, nonpositive metric field, lost characteristic cone,
constraint or spectral failure, trapped-sign failure, boundary contamination,
candidate-action health stop, or scale-control stop keeps its prior ownership
and fail-closed behavior.

## What remains identical

PROTO6 inherits PROTO5 by exact SHA-256 except for that source-stop ownership
and the required output names. It does not change:

- the action, equations, couplings, initial fields, pulse profiles, or
  amplitude order `[5/2, 3]`;
- the direct affine source solver, refinement rule, raw residual threshold, or
  kinetic condition threshold;
- native `q=D_hu` initialization, method-owned constraint guards, roundoff
  classification, common-event alignment, weighted spectra, or convergence
  rules;
- physical-null orientations, trapped-sphere classification, Raychaudhuri
  observable, outcome labels, error accounting, or robustness burden;
- any SGB-L or FGC-QR health requirement.

The new output roots are `runs/fgc-2-sf1/proto6/calibration` and
`runs/fgc-2-sf1/proto6/holdout`; the future resolved manifest is
`configs/fgc/fgc-1-pro6-hld1.toml`. PROTO5 output is immutable diagnostic
history and cannot be relabeled as PROTO6 evidence.

## Why execution is still closed

CAL2 proves the proposed distinction on the saved failure, but CAL2's replay
is design evidence rather than a completed calibration. `FGC-1-HLT4-MON4`
must implement the accepted-state precheck, trial rejection, rollback, retry
ledger, namespace precondition, and injected-failure controls before it may
authorize one fresh PROTO6 GR-0 campaign.

Consequently PRO6-FRZ1 does not answer the mechanism question. It authorizes
no trajectory and establishes no trapped sphere, regulator activation,
Raychaudhuri margin, transition, singularity resolution, child domain, dark
sector, varying locally measured light speed, or rejection of FGC-QR or the
general gradient route.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro6_frz1.py \
  --output results/fgc-1-pro6-frz1.json
python3 -m unittest tests.test_fgc_protocol_v6 -v
```
