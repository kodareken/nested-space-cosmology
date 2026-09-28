# FGC-1-DEF1-STAB1-FRZ1: candidate-blind conversion/error-map instrument freeze

**FGC-1-DEF1-STAB1-FRZ1** prospectively freezes the candidate-blind
pre-holdout DEF1 stability/error-map **instrument contract**. It binds the
current DEF1 owner bytes on implementation image
`c4feb941e408a7c72913b49071801b7531183c1a` and directly exercises the
already-owned conversion, coverage, algebra, and margin controls.

This is not [FGC-1-DEF1-STAB1-PREF1](fgc-def1-stab1.md). PREF1 is not
implemented here. FRZ1 licenses only a separately implemented independent
PREF1 binder. It does not set `DEF1_error_map_passed`, evaluate the nine
trajectory DEF1 booleans, read a candidate or control trajectory, open
holdout or ROB1, or claim a mechanism or physical result.

The partial instrument owner remains [FGC-1-DEF1-STAB1](fgc-def1-stab1.md).
Current campaign routing stays in [PLAN.md](../PLAN.md). PRO20-PREF1 remains
the critical frontier; RSRC1 remains the critical-path next design. This
freeze is a parallel nonexecuting gate artifact.

## What is frozen

The compact config and result prove that the conversion/error-map instrument
is complete as a candidate-blind contract:

1. Exact RED1/MHG2 six-row identity and FO1/ADM 26-input slot orders.
2. All fourteen-component provider routes: each positive control and each
   injected-failure refusal.
3. Componentwise `Q`-error assembly from sensitivity times input-error
   enclosures, independent of measured signal size; a `measured_q` argument
   is refused.
4. Misner–Sharp flux/conservation enclosure, including a typed
   `mass_flux_inconsistency` veto.
5. Activation/control comparison under canonical GR-0 `phi = 0`.
6. Trappedness and complete-`Q` margin helpers, including strict-boundary
   failures.
7. Exact BASE-to-ADM two-jet roundtrip packing all 26 geometry inputs.
8. Conditional IMP1 18-channel-to-`Q` debit that cannot set
   `DEF1_error_map_passed` or a global PDE theorem.

Owner SHA-256 values of `def1_stab1.py`, `def1_geometry_error.py`,
`def1_stab1_providers.py`, and `def1_stab1_qualification.py` are bound
together with the non-promoting freeze-contract payload.

## Ordinary verification

```bash
python3.14 -I -B scripts/reproduce_fgc_def1_stab1_frz1.py --check
make fgc-def1-stab1-frz1
```

`--check` is Git/raw/`runs/`/source/trajectory blind. It never executes a
runner, never reconstructs from live owner modules, and fails closed on
altered canonical config or result bytes. Initial construction may use
`--emit-config`, `--emit-result`, or `--write-result`; those paths reconstruct
from live owners and must not read `runs/` or candidate/control trajectories.

## Claims

True:

- `conversion_instrument_contract_complete`
- direct exercise of conversion identities, provider-route coverage,
  `Q`-error assembly, Misner–Sharp enclosure, activation, trappedness and
  complete-`Q` margins, BASE-to-ADM roundtrip, and IMP1-to-`Q` debit
- `licenses_separate_independent_pref1_binder`

False:

- `DEF1_error_map_passed`
- `def1_booleans_evaluated` (all nine names remain unevaluated)
- candidate/control trajectory read
- global PDE error, measured-`Q` dependence, fitted threshold
- ROB1, holdout, mechanism, physics
- PREF1 implemented

Richardson and IMP1 remain conditional premises.
