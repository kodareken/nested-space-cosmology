# FGC-1-TDG5-IMP1: atomic temporal-refinement runtime

`FGC-1-TDG5-IMP1` implements and synthetically validates the bounded successor
authorized by `FGC-1-TDG5-PREF20` at immutable commit `76da687…`. It reads no
PROTO13 history or checkpoint, resumes no campaign, and evaluates no SGB-L or
FGC-QR state.

## Result

The frozen TDG5 same-grid runtime pair now exists for both inherited methods.
From one bitwise-identical accepted state it constructs

```text
coarse shadow: one full step of width Delta t
fine shadow:   two consecutive steps of width Delta t / 2.
```

Each shadow uses the unchanged numerical engine and PROTO7 source, health,
CFL, causal-boundary, and transaction machinery. The coarse and fine paths
have independent monitor and causal ledgers. Preparation cannot mutate the
real accepted state or its real ledgers. A commit token rechecks the state
hash, coordinate time, step index, member transaction serial, monitor state,
and causal state before adopting anything. It then adopts only the fine
two-half-step endpoint and the fine shadow ledgers.

The synthetic controls show the resulting accounting explicitly:

```text
RK4:    two fine proposals, 10 accepted stage records, member serial 10
SSPRK3: two fine proposals,  8 accepted stage records, member serial  8
```

The corresponding coarse proposal is retained as evidence and is never
committed. A deliberately injected non-source failure at `t=3/32`, which is
visited only by the second RK4 half-step in the control, leaves the real state,
monitor, and causal ledger exactly unchanged. A post-prepare ledger mutation
also makes commit fail closed.

This is a two-phase core transaction, not yet the production campaign
compositor. Normal-flow tracer preview/commit, durable rejection records,
retry policy, and checkpoint integration remain obligations of the later
prospectively frozen runtime gate.

## Continuous complete-state comparison

For every dimensionless component of `(u,p,q)` on noncentre rows excluding the
four projector-owned outer rows, IMP1 constructs a cubic Hermite record from

```text
left endpoint value,
Delta t times the fresh left-endpoint RHS,
right endpoint value,
Delta t times the fresh right-endpoint RHS.
```

The coarse cubic is restricted to each normalized half interval and subtracted
from the corresponding fine cubic. For every resulting cubic, the runtime
evaluates both endpoints and every numerically real interior root returned by
the scaled quadratic derivative solve. Near-degenerate discriminants are
counted explicitly rather than silently classified.

The implementation keeps four quantities distinct:

1. the largest value found at the theorem's endpoint/stationary-point
   candidates;
2. coefficient-construction roundoff;
3. binary64 root and value arithmetic debit; and
4. any additional gap required by an independently computed cubic Bernstein
   convex-hull envelope.

The final continuous upper bound dominates both the stationary-point route and
the inflated Bernstein route. This makes the runtime fail closed even when a
derivative discriminant is numerically delicate. It is intentionally more
conservative than pretending a rounded root is exact. The Bernstein gap is
reported as certification slack, not mislabeled as arithmetic roundoff.

The controls include the endpoint-only bump

```text
4 theta (1-theta),
```

whose raw candidate maximum is exactly one at its single interior stationary
point, a genuine cubic with two distinct interior stationary points, an exact
zero control that remains zero without an invented numerical floor, four
dense independent audits with `131073` points each, and a one-bit coefficient
mutation.

## Error semantics

The continuous upper bound is divided by the method-owned conditional
Richardson denominator:

```text
RK4:    15
SSPRK3:  7.
```

That quantity is a prospective fine-path numerical debit. It is not an
acceptance threshold, a proof that a production trajectory is asymptotic, or a
rigorous global PDE error estimate. Spatial ladders, cross-method agreement,
constraints, source/health/boundary controls, and the later complete DEF1
error ledger remain independent requirements.

## What the implementation earns

The runtime pair and its synthetic validation are complete. The next bounded
task is prospective threshold and admission design using independent controls.
No threshold may be selected from a trapped-sphere or Raychaudhuri outcome,
and no old CAL10 event is reclassified.

```text
TDG5_stage_complete_refinement_theorem_completed = true
runtime_refinement_pair_implementation_authorized = true
runtime_refinement_pair_implemented = true
runtime_refinement_pair_synthetic_validation_passed = true
runtime_threshold_design_authorized = true
runtime_thresholds_frozen = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

The historical PROTO13/CAL10 stop remains valid and unreclassified. IMP1 is an
implemented numerical instrument, not a collapse, defocusing, transition, or
physical result.
