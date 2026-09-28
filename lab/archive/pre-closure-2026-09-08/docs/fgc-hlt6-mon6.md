# FGC-1-HLT6-MON6: PROTO8 runtime and one-run authorization

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-HLT6-MON6` implements the two common-event diagnostic ownership
adapters frozen by `FGC-2-SF1-PROTO8` and authorizes exactly one fresh GR-0
calibration in the new `proto8` namespace.

```text
PROTO8_successor_runtime_compositor_implemented = true
PROTO8_fresh_GR0_dynamic_calibration_authorized = true
PROTO8_resolved_holdout_manifest_authorized = false
classical_spherical_diagnostic_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

This is a pre-trajectory implementation and authorization certificate. It is
not a calibration, collapse, trapped-sphere, regulator, defocusing, mechanism,
or physical-model result.

## Immutable boundary

The authorization is tied to immutable checkpoint
`b4c3cdcd216e305ef55e2f2700fbee4111afc95c`, which contains the exact PROTO8
protocol and PRO8-FRZ1 freeze. It also binds the prior HLT5 runtime, ID2 static
inputs, and CAL4 plan by their exact hashes. Every implementation file used by
the inherited PROTO7 source and transaction is checked against HLT5 before the
new authorization can reproduce.

The disclosed PROTO7 campaign supplied the diagnosis that motivated PROTO8.
HLT6 reads no PROTO8 trajectory and creates no run directory.

## Runtime delta

The evolution path is still the exact PROTO7 path:

- the accepted-state source gate remains terminal;
- every source-only failure in a wholly unaccepted proposal may retry;
- any non-source failure vetoes retry;
- the raw `1e-12` source gate, affine root solver, monotonicity and kinetic
  guards, retry factor, retry budget, minimum step, grids, methods, CFL rule,
  dissipation, boundary ledger, and checkpoint transaction are unchanged; and
- cross-member failure still restores the last atomic common event.

HLT6 adds no stage or trajectory logic. Only after all six members have
committed to one bitwise-identical common time does it replace the shared
diagnostic compositor.

## Evolution-owned constraint diagnostic

Both full-domain and evolution-owned raw residuals are serialized. The owned
domain excludes the four projector-fixed rows and every row whose derivative
stencil touches them:

```text
RK4 / fourth-order SBP:       4 + 3 = 7 excluded outer rows
SSPRK3 / second-order SBP:    4 + 1 = 5 excluded outer rows
```

Only the owned raw norm decides the unchanged coarsest and finest magnitude
guards. The deterministic envelope

```text
4096 * binary64_epsilon * (1 + accepted_stage_count)
```

may classify constraint convergence order as indistinguishable from
accumulated roundoff. It never changes a raw magnitude guard and cannot be
subtracted from a trapped-sign or Raychaudhuri margin.

Injected controls prove that a defect confined to the projector/stencil
interface remains public but is not attributed to evolution, while the same
defect moved into the owned interior remains terminal.

## Nested spatial-spectrum diagnostic

The coarsest grid is a public convergence witness. It may miss an absolute
budget only if both the medium and finest grids pass every unchanged field,
derivative, and RMS-scale budget and every coarse-to-medium and
medium-to-fine tail ratio remains below the unchanged `1/4` ceiling.

Injected controls establish all three directions:

- a coarse-only miss with a resolved finest pair and contracting tails passes;
- a medium or finest absolute-budget miss fails; and
- a noncontracting adjacent tail fails even when the finest pair passes its
  absolute budgets.

The causal-past temporal spectrum is inherited unchanged. The trapped-sign
observable and its additive two-method error ledger are also unchanged.

## Frozen inputs and namespace

The twelve inputs are the Cartesian product

```text
amplitudes:  [5/2, 3]
methods:     [RK4, SSPRK3]
resolutions: [1025, 2049, 4097]
```

Every projected `u,p,q` state is reconstructed and must match HLT5 bitwise.
All twelve accepted-state source prechecks pass the unchanged raw gate, and
all four amplitude/method groups pass the PROTO8 common-event compositor at
`t=0`. The new expanded-input and campaign hashes contain the PROTO8 plan and
diagnostic ownership but no trajectory result.

The fresh output roots are:

```text
runs/fgc-2-sf1/proto8/calibration
runs/fgc-2-sf1/proto8/holdout
```

Both are absent or empty at authorization, and HLT6 creates neither. The
runner refuses a dirty tracked worktree, implementation drift, authorization
drift, pre-existing output, overwrite, or a mismatched resume manifest.
PROTO7 output remains immutable diagnostic history and cannot be relabelled as
PROTO8 evidence.

## What execution may establish

The authorization permits one fresh outcome-neutral GR-0 amplitude
calibration. It does not authorize SGB-L, FGC-QR, a holdout manifest, retained
EFT evolution, or a physical transition claim. A later calibration result
must be reduced and independently classified before any successor gate can
open.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt6_mon6.py \
  --output results/fgc-1-hlt6-mon6.json
python3 -m unittest \
  tests.test_fgc_proto8_runtime \
  tests.test_fgc_gr0_campaign_runner_v8 \
  tests.test_fgc_hlt6_mon6_reproduction -v
```

After the complete HLT6 unit is committed, authorization can be checked
without creating output:

```bash
python3 scripts/run_fgc_gr0_calibration_v8.py \
  --check-authorization-only
```

The fresh campaign is a separate, post-commit action:

```bash
python3 scripts/run_fgc_gr0_calibration_v8.py
```
