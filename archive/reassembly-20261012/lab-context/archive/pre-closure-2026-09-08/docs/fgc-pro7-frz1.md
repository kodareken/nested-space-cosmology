# FGC-1-PRO7-FRZ1: outcome-neutral PROTO7 freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-PRO7-FRZ1` freezes **FGC-2-SF1-PROTO7** as the smallest general
successor to PROTO6's stage-labelled source-retry contract. It is tied to
immutable checkpoint `27c63e7aca1439e75e18b576c6da414e458280ff`, containing
the exact PROTO6 protocol and CAL3-PREF5 diagnosis.

```text
PROTO7_outcome_neutral_protocol_frozen = true
PROTO7_immutable_lineage_verified = true
PROTO7_accepted_state_source_terminal_contract_frozen = true
PROTO7_unaccepted_proposal_source_retry_contract_frozen = true
PROTO7_non_source_retry_veto_contract_frozen = true
PROTO7_complete_failure_observability_contract_frozen = true
PROTO7_successor_runtime_compositor_implemented = false
PROTO7_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

This is a numerical-premise revision, not a calibration, collapse, or
mechanism result.

## The transaction boundary

PROTO7 uses two ownership states only.

Before each proposal, the source gate is evaluated on the bitwise last
accepted state. A miss there remains the terminal scientific premise failure
`newton_residual_limit`.

The proposal then consists of every internal integrator evaluation plus the
independently evaluated candidate endpoint. None of those states is accepted
until the complete transaction commits. If every failure in that unaccepted
proposal is source-only, the entire proposal may be discarded and retried
from the exact accepted state with the inherited factor `1/2`. If any
non-source failure appears, retry is forbidden and the original failure
priority applies.

This closes the stage-label exception class. `rk4_k2`, `rk4_k4`,
`ssprk3_k3`, and `candidate_endpoint` do not have different physical status:
before commit they are all trial evaluations. The rule does not make an
accepted invalid state retryable.

## Complete failure evidence

Every retry, terminal stop, and retry-exhaustion record must retain enough
information to reconstruct the decision:

- member, method, grid size, initial step, attempted step, and retry count;
- accepted-state content hash, time, and transaction serial;
- every failed source call's stage name, stage time, transaction serial, raw
  residual, and raw threshold;
- the complete failure set, selected terminal reason, and whether a
  non-source failure vetoed retry;
- proof that fields, time, accepted-stage count, causal debit, and external
  transaction state were restored before retry or terminal return.

The runner may not collapse a failed transaction to only its first stage
label. This observability rule changes no admissibility threshold; it makes
the existing decision auditable.

## What remains identical

PROTO7 inherits PROTO6 by exact SHA-256 except for transaction ownership,
failure-record completeness, and fresh names. It does not change:

- the action, equations, couplings, initial profiles, amplitudes, or their
  frozen order `[5/2, 3]`;
- the affine source solver, refinement rule, raw `1e-12` residual threshold,
  kinetic condition threshold, retry factor, retry maximum, or minimum step;
- grids, methods, CFL rule, constraints, spectra, convergence, trapped-sign,
  boundary, health, scale, affine-null, Raychaudhuri, outcome, or robustness
  rules;
- any SGB-L or FGC-QR applicability or promotion requirement.

The new output roots are `runs/fgc-2-sf1/proto7/calibration` and
`runs/fgc-2-sf1/proto7/holdout`. PROTO6 output remains immutable diagnostic
history and cannot be relabelled as PROTO7 evidence.

## Why execution remains closed

CAL3 supplies design evidence for this transaction rule, not a PROTO7
calibration. `FGC-1-HLT5-MON5` must independently implement and attack the
new compositor, prove exact rollback and non-source vetoes at every proposal
position, freeze unchanged inputs, verify fresh namespaces, and authorize one
new GR-0 campaign.

Consequently PRO7-FRZ1 does not answer the mechanism question. It authorizes
no trajectory and establishes no trapped sphere, regulator activation,
Raychaudhuri margin, transition, singularity resolution, child domain, dark
sector, varying locally measured speed of light, or rejection of FGC-QR or
the general gradient route.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro7_frz1.py \
  --output results/fgc-1-pro7-frz1.json
python3 -m unittest tests.test_fgc_protocol_v7 -v
```
