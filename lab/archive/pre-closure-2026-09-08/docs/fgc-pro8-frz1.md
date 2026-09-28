# FGC-1-PRO8-FRZ1: outcome-neutral PROTO8 freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-PRO8-FRZ1` freezes **FGC-2-SF1-PROTO8** as the smallest successor to
PROTO7's common-event diagnostic ownership contract. It is tied to immutable
checkpoint `219ba9493899029bc51bd313e3c5f689f3809c2c`, containing the exact
PROTO7 protocol and CAL4-PREF6 diagnosis.

```text
PROTO8_outcome_neutral_protocol_frozen = true
PROTO8_immutable_lineage_verified = true
PROTO8_evolution_owned_constraint_contract_frozen = true
PROTO8_accumulated_roundoff_order_contract_frozen = true
PROTO8_finest_pair_nested_spectral_contract_frozen = true
PROTO8_successor_runtime_compositor_implemented = false
PROTO8_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

This is a premise-only diagnostic revision. It is not a calibration, collapse,
or mechanism result.

## Constraint ownership at a common event

The outer projector intentionally fixes four rows. A derivative stencil that
touches those projector-fixed rows is also outside the reduction diagnostic's
evolution-owned domain. PROTO8 therefore excludes exactly:

```text
RK4 / fourth-order SBP:       4 fixed + 3 stencil-reach rows = 7
SSPRK3 / second-order SBP:    4 fixed + 1 stencil-reach row  = 5
```

This exclusion applies only to the common-event convergence diagnostic. The
full-domain and owned-domain raw residuals, their maxima, and their locations
must both remain public. The owned raw norm decides the unchanged coarsest and
finest magnitude guards. An interior reduction defect cannot hide behind the
ownership boundary.

The former `4096*epsilon_64` enclosure described one diagnostic evaluation,
while `q` and `D_h u` accumulate independently through accepted stages.
PROTO8 freezes the deterministic order-classification envelope

```text
4096 * binary64_epsilon * (1 + accepted_stage_count).
```

It may classify an order as indistinguishable from accumulated roundoff. It
may not change a raw magnitude guard, subtract error from a physical or
Raychaudhuri margin, or excuse a nonconvergent interior component. The
unchanged minimum finite finest-pair order remains `3/2`.

## Spectral roles across the nested grids

PROTO8 gives the three grids nonconflicting jobs:

- the coarsest grid remains a public convergence witness;
- the medium and finest grids must each pass every unchanged absolute field,
  derivative, and RMS-scale budget; and
- every coarse-to-medium and medium-to-fine field and derivative tail ratio
  must remain below the unchanged `1/4` ceiling.

A coarse-only absolute-budget miss therefore cannot pass by itself. It passes
only when the full nested sequence contracts and the finest pair is already
resolved. Any noncontracting tail remains terminal. All raw per-grid budgets
and ratios remain serialized.

## What remains identical

PROTO8 inherits PROTO7 by exact SHA-256 except for these two diagnostic
ownership adapters and fresh names. It does not change:

- the action, equations, couplings, initial profiles, amplitudes, or frozen
  order `[5/2, 3]`;
- the affine source solver, raw `1e-12` residual gate, kinetic condition,
  transaction boundary, retry factor, 32-retry limit, or minimum step;
- methods, grids, CFL rule, constraint magnitude/order thresholds, spectral
  thresholds, trapped-sign test, boundary, health, scale, affine-null,
  Raychaudhuri, outcome, or robustness rules; or
- any SGB-L or FGC-QR applicability or promotion requirement.

The fresh output roots are `runs/fgc-2-sf1/proto8/calibration` and
`runs/fgc-2-sf1/proto8/holdout`. PROTO7 output remains immutable diagnostic
history and cannot be relabelled as PROTO8 evidence.

## Why execution remains closed

CAL4 supplies design evidence for these ownership rules, not a PROTO8
calibration. `FGC-1-HLT6-MON6` must independently implement the two adapters,
attack their exact boundaries, prove the unchanged transaction/source path,
reconstruct the twelve frozen inputs, verify fresh namespaces, and authorize
one new GR-0 campaign.

Consequently PRO8-FRZ1 does not answer the mechanism question. It authorizes no
trajectory and establishes no selected amplitude, trapped sphere, regulator
activation, Raychaudhuri margin, transition, singularity resolution, child
domain, dark sector, varying locally measured speed of light, or rejection of
FGC-QR or the general gradient route.

## Reproduction

```bash
python3 scripts/reproduce_fgc_pro8_frz1.py \
  --output results/fgc-1-pro8-frz1.json
python3 -m unittest tests.test_fgc_protocol_v8 -v
```
