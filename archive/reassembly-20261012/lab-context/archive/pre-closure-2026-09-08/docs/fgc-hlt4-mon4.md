# FGC-1-HLT4-MON4: PROTO6 GR-0 runtime authorization

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

## Decision

`FGC-1-HLT4-MON4` implements the one numerical ownership change frozen by
`FGC-2-SF1-PROTO6` and authorizes one fresh GR-0 recalibration in the `proto6`
namespace. It does not authorize SGB-L, FGC-QR, DEF1, or a retained-EFT or
physical-transition claim.

```text
PROTO6_successor_runtime_compositor_implemented = true
PROTO6_fresh_GR0_dynamic_calibration_authorized = true
PROTO6_resolved_holdout_manifest_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

The certificate reads PROTO6, PRO6-FRZ1, HLT3, and ID2 from immutable commit
`f3a32ff9faddfa4dba0e0f5e192d58ef94d5e094`. It reconstructs all twelve
amplitude/method/resolution inputs, requires their projected state hashes to
match HLT3 exactly, and evaluates the unchanged raw source gate on every
accepted `t=0` state. The four amplitude/method common events must also retain
the inherited constraint and proper-spectrum admission.

## Exact transaction boundary

Before every proposal, the same affine GR-0 source is evaluated on the
bitwise last accepted state. A raw residual at or above `1e-12` is immediately
latched as the terminal `newton_residual_limit`; no step-size change can rescue
an invalid accepted state.

For a proposed step, a shadow copy of the complete PROTO5 transaction checks
the source, metric, cone, kinetic, causal-boundary, and CFL premises without
mutating the real transaction. A retry is permitted only when the complete
failure tuple contains `newton_residual_limit` alone and the failing record is
an internal RK4 or SSPRK3 stage. The explicit candidate endpoint remains
terminal because PROTO6 froze only internal-stage ownership. A simultaneous
metric, kinetic, boundary, or other non-source failure also remains terminal.

A rejected source trial:

- preserves the accepted arrays bitwise;
- advances no time, step index, transaction serial, causal debit, or tracer;
- is fsynced to the append-only event ledger with its member, stage, stage
  time, proposed step, observed residual, raw threshold, and retry count;
- halves the step using the frozen factor `1/2`;
- stops on the frozen count of 32 retries or minimum step `2^-30`; and
- restarts only from the last accepted state.

The runner additionally snapshots every resolution/method member at each
common-event boundary. If any member stops before the shared event is
complete, all members are restored to that boundary before the terminal
checkpoint is written. The attempted retry records remain evidence, while a
partially advanced cross-resolution state cannot masquerade as a common
event.

## What was not changed

The accepted-state RHS produced by the PROTO6 adapter is checked bitwise
against the original GR-0 operator on both numerical methods. The action,
initial data, amplitudes, grids, affine root and refinement, raw and kinetic
thresholds, CFL rule, constraint guards, spectra, trapped-sign observable,
causal boundary, outcome classification, and robustness burden are inherited
unchanged. The new run plan is mechanically normalized back to CAL2 and must
then match HLT3's exact validated plan; only the declared retry ownership and
fresh namespace may differ.

## Epistemic boundary

HLT4 uses synthetic failures to test transaction ownership and reconstructs
physical inputs only at `t=0`. It reads no PROTO6 trajectory and classifies no
collapse or trapped-sphere outcome. The prior PROTO5 campaign remains
diagnostic history; it cannot be relabeled as PROTO6 evidence.

The next allowed operation is therefore exactly one fresh campaign:

```bash
python3 scripts/run_fgc_gr0_calibration_v6.py
```

The campaign result must be reduced independently before any holdout manifest
or candidate evolution can be considered.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt4_mon4.py \
  --output results/fgc-1-hlt4-mon4.json
python3 -m unittest tests.test_fgc_proto6_runtime -v
python3 -m unittest tests.test_fgc_hlt4_mon4_reproduction -v
```
