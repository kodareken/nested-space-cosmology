# FGC-1-PRO5-FRZ1: outcome-neutral PROTO5 freeze

## Decision

`FGC-1-PRO5-FRZ1` freezes **FGC-2-SF1-PROTO5** as the premise-only successor
to the semidiscretely obstructed PROTO4 contract. It is tied to immutable
checkpoint `0781d7d218e9c129858b134c3181fad07ac692b0`, which contains the exact
PROTO4 protocol and CAL1-PREF3 diagnosis named by the amendment.

```text
PROTO5_outcome_neutral_protocol_frozen = true
PROTO5_immutable_lineage_verified = true
PROTO5_native_semidiscrete_initialization_contract_frozen = true
PROTO5_method_owned_constraint_admission_frozen = true
PROTO5_roundoff_zero_classification_contract_frozen = true
PROTO5_successor_runtime_compositor_implemented = false
PROTO5_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

This is a repaired input-to-discretization contract, not a trajectory or a mechanism result.

## What changed

PROTO5 changes only the three premises diagnosed before trajectory access:

1. ID2 continues to own the physical initial fields `u,p`; both must remain
   bitwise unchanged. Before the first stage, only the kinematic auxiliary
   variable is initialized as `q=D_hu` with the native SBP operator of each
   method. This is discrete reduction-constraint initialization, not a change
   to the physical data or field equations.
2. Raw normalized residuals still decide every magnitude guard. The primary
   D4-2/RK4 path retains PROTO4's `1/50` coarsest and `1/1000` finest limits.
   The independently lower-order D2-1/SSPRK3 comparator uses the prospectively
   frozen `1/10` and `1/200` limits. Both remain subject to common-event
   alignment, three-grid monotone refinement, and finest-pair order at least
   `3/2`.
3. A public `4096 epsilon_64` enclosure may classify an analytically zero
   component as roundoff-enclosed for convergence order only. It never changes
   a raw residual, never relaxes a magnitude guard, and can enter later error
   accounting only as a nonnegative contribution—never as a subtraction.

The initial-event values used to diagnose PROTO4 were seen before this freeze
and are disclosed in CAL1. No calibration trajectory was seen. The two
eligible amplitudes remain `[5/2, 3]` in that order.

## What did not change

PROTO5 inherits the complete PROTO4 overlay by exact SHA-256. The action,
couplings, profiles, held-out cases, physical equations, null observable,
outcome labels, robustness burden, branch-owned stop partition, weighted
spectral admission, and trapped-sign classification are unchanged. PROTO1
through PROTO4 remain immutable history.

The output roots are newly versioned under `runs/fgc-2-sf1/proto5/`, and the
future manifest is `configs/fgc/fgc-1-pro5-hld1.toml`. Earlier output, if any,
is development-only and cannot be promoted into PROTO5 evidence.

## Why calibration is still closed

CAL1 demonstrated the repair on the initial event, but a demonstration is not
the runtime owner. **FGC-1-HLT3-MON3** must bind the exact PROTO5 rules into
the stage transaction, re-run the initial composition, verify the untouched
namespace and immutable candidate inputs, and issue a separate calibration
authorization. PROTO5 itself therefore keeps fresh calibration, the resolved
manifest, SGB-L, and every FGC-QR execution gate false.

## Reproduction and nonclaims

```bash
python3 scripts/reproduce_fgc_pro5_frz1.py \
  --output results/fgc-1-pro5-frz1.json
python3 -m unittest tests.test_fgc_protocol_v5 -v
```

The certificate reads neither output namespace. It derives no collapse,
trapped interval, regulator activation, defocusing, transition, singularity
resolution, child domain, dark sector, varying locally measured light speed,
or rejection of FGC-QR or the general gradient route.
