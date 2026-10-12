# FGC-1-DEF1-STAB1-PREF1: candidate-blind pre-holdout error-map readiness

**FGC-1-DEF1-STAB1-PREF1** independently binds the candidate-blind
pre-holdout DEF1 conversion/error-map as **map readiness only**. It is the
direct successor of sealed FRZ1 commit
`32438805126266999a089396f2025e94f8e6bad4` (parent
`c4feb941e408a7c72913b49071801b7531183c1a`).

This is not a trajectory result and not Wave-5 `FGC-1-DEF1-PREF1`. It does
not evaluate the nine saved DEF1 booleans, consume measured `Q` or actual
error radii, authorize holdout or ROB1, or claim a mechanism or physical
result. `DEF1_error_map_passed=true` here means `map_readiness_only=true`.

The partial instrument owner remains [FGC-1-DEF1-STAB1](fgc-def1-stab1.md).
The freeze owner remains [FGC-1-DEF1-STAB1-FRZ1](fgc-def1-stab1-frz1.md).
Current campaign routing stays in [PLAN.md](../PLAN.md). PRO20-PREF1 remains
the critical frontier; RSRC1 remains the critical-path next design. This
binder is a parallel nonexecuting gate artifact.

## Independent reconstruction

PREF1 does not import or call `def1_stab1_frz1_certificate.py` or
`def1_stab1_freeze_contract.py`. Explicit live construction authenticates
the exact FRZ1 commit, parent, and delta, reads the tracked freeze config
SHA-256 `c88676a7…` and result SHA-256 `c2e5c8c2…`, and treats those bytes
as evidence to challenge. It then independently reconstructs the
conversion/error-map contract from the low-level DEF1 owners:

1. Exact RED1/MHG2 six-row identity and FO1/ADM 26-input slot orders.
2. All 21 provider routes: each positive control and each injected-failure
   refusal.
3. Componentwise `Q`-error assembly independent of measured signal size; a
   `measured_q` argument is refused.
4. Misner–Sharp flux/conservation enclosure, including a typed
   `mass_flux_inconsistency` veto.
5. Activation/control comparison under canonical GR-0 `phi = 0`.
6. Trappedness and complete-`Q` margin helpers, including strict-boundary
   failures.
7. Exact BASE-to-ADM two-jet roundtrip packing all 26 geometry inputs.
8. Conditional IMP1 18-channel-to-`Q` debit that cannot itself set a global
   PDE theorem.
9. Live owner SHA-256 values of `def1_stab1.py`, `def1_geometry_error.py`,
   `def1_stab1_providers.py`, and `def1_stab1_qualification.py`.
10. The complete remaining-work taxonomy from the provider owner.

The remaining-work items are reclassified, not discarded. Conditional
Richardson/IMP1, caller-supplied C-to-jet or complete residuals,
inverse-J/Gronwall/derivative premises, and typed refusal on absent inputs
are the operational premises of the map. COL1 trajectory binding and the
nine-boolean DEF1 classification remain later application. A universal
nonlinear PDE theorem is not a prerequisite and is not claimed.

If and only if every independent predicate closes, PREF1 sets

```text
DEF1_error_map_passed = true
map_readiness_only = true
```

and keeps all nine trajectory DEF1 booleans null.

## Ordinary verification

```bash
python3.14 -I -B scripts/reproduce_fgc_def1_stab1_pref1.py --check
make fgc-def1-stab1-pref1
```

`--check` is Git/raw/`runs/`/source/trajectory/reconstruction blind. It never
executes a runner, never reconstructs from live owner modules, and fails
closed on altered canonical config or result bytes. Initial construction may
use `--emit-config`, `--emit-result`, or `--write-result`; those paths inspect
Git and low-level DEF1 owners and must not read `runs/` or candidate/control
trajectories. The binder module itself never writes state.

## Claims

True:

- `DEF1_error_map_passed` with `map_readiness_only`
- independent reconstruction of conversion identities, all 21 provider
  routes, `Q`-error assembly, Misner–Sharp enclosure, activation,
  trappedness and complete-`Q` margins, BASE-to-ADM roundtrip, conditional
  IMP1 debit, owner hashes, and remaining-work taxonomy
- FRZ1 commit/delta and tracked freeze hashes authenticated and challenged
- operational premises and typed refusal make the map ready before holdout

False:

- `def1_booleans_evaluated` (all nine names remain unevaluated)
- candidate/control trajectory read
- actual error bounds evaluated
- global PDE error, measured-`Q` dependence, fitted threshold
- ROB1, holdout, mechanism, physics
- `GR0_case_eligible` and `SGBL_branch_owned_and_healthy` (separate gates)

Richardson and IMP1 remain conditional premises.
