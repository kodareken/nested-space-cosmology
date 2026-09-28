# FGC-1-PRO9-FRZ1: outcome-neutral PROTO9 freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced protocol and immutable-lineage freeze; no
runtime or mechanism result

```text
PROTO9_outcome_neutral_protocol_frozen = true
PROTO9_resolution_ladder_contract_frozen = true
PROTO9_successor_runtime_implemented = false
FGCQR_holdout_execution_authorized = false
```

## Frozen change

PROTO9 consumes FGC-2-SF1-PROTO8 and FGC-1-CAL5-PREF7 from immutable commit
`a2c929e33ec8f2ee5464b086a57db8ed7f9277f0`. It changes exactly one
numerical axis:

```text
PROTO8: 1025 -> 2049 -> 4097
PROTO9: 2049 -> 4097 -> 8193
```

The entire nested ladder moves one level. The coarsest grid remains a public
convergence witness, the medium and finest grids retain every absolute
spectral budget, both adjacent field and derivative tails must remain below
the unchanged `1/4` ceiling, and the same coarse/fine constraint magnitude
guards plus `3/2` minimum component order remain in force.

This is not a threshold repair. CAL5 showed that the existing `2049` states
already pass the unchanged prospective coarse constraint guards and the
existing `4097` states already pass all absolute spectral budgets. It also
showed that amplitude `3` still fails the `2049 -> 4097` `phi/Lambda`
derivative-tail contraction and that amplitude `5/2` SSPRK3 still misses its
order gate at the terminal PROTO8 event. PROTO9 preserves those observations
and requires genuinely new `8193` data. Conditional Richardson projections
remain non-admissive planning evidence.

## Exact inheritance

The following remain byte- or contract-inherited from PROTO8 and its
predecessors:

- the ACT1 action, branches, couplings, profiles, and physical initial data;
- amplitudes `[5/2, 3]` in the same order;
- RK4 and SSPRK3, their spatial orders, CFL maximum, dissipation, retry factor,
  retry budgets, minimum step, event interval, and final time;
- the raw `1e-12` source gate, solver, accepted/unaccepted proposal ownership,
  non-source veto, and complete rollback evidence;
- evolution-owned constraint rows, public full residuals, accumulated-roundoff
  order-only classification, magnitude guards, and order threshold;
- spatial and temporal spectrum definitions, every absolute budget, every
  nested-tail ceiling, and the trapped-sign/error test;
- causal boundary, health, scale, affine-observable, outcome, and robustness
  rules; and
- every candidate, retained-EFT, physical, global, and cosmological nonclaim.

The only other syntactic change is provenance: future calibration and holdout
outputs must live under fresh `proto9` namespaces, and no PROTO8 state may be
relabelled as a PROTO9 outcome.

## Immutable evidence and fail-closed boundary

PRO9-FRZ1 verifies the predecessor protocol and CAL5 config/result directly
from the checkpoint commit, validates their declared SHA-256 hashes, and
checks that the checkpoint is an ancestor of the current tree. Both new output
roots are absent when frozen. Mutation controls reject a different resolution,
a relaxed gate, old-output reuse, or any promoted claim.

The freeze performs no evolution and does not implement HLT7. It does not
complete GR-0 calibration, select an amplitude, resolve a holdout manifest,
read SGB-L or FGC-QR, classify a trapped sphere, test regulator activation,
answer the mechanism question, validate the retained EFT, reject gradients,
resolve a singularity, derive a child domain or dark sector, or vary locally
measured `c`.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro9_frz1.py
```
