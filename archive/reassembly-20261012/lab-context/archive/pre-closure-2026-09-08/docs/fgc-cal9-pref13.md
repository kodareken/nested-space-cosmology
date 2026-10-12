# FGC-1-CAL9-PREF13: PROTO12 late-event constraint-order diagnosis

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced post-calibration numerical diagnosis; no
mechanism result

```text
PROTO12_campaign_terminated_normally = true
PROTO12_both_amplitudes_reached_late_evolved_common_events = true
PROTO12_GR0_case_eligible = false
SRC4_arithmetic_wall_cleared_for_both_amplitudes = true
PROTO12_dual_amplitude_SSPRK3_radial_momentum_order_veto_localized = true
RSP2_high_ladder_constraint_preflight_required = true
FGCQR_holdout_execution_authorized = false
```

## What the fresh campaign established

HLT10 authorized exactly one fresh GR-0 calibration under PROTO12. The
campaign ran from immutable commit
`2f9bddbd4967e4a8998845513f345df3475a4bea`, terminated normally after about
16,924 seconds, and wrote 49 canonical common-event records under campaign
identity
`32ac798e8bd7ccb984dd66262439d55ec01ad7c5cabfea785fe30d1fc60eff36`.
There are no other event types in the ledger.

The terminal campaign classification is exactly:

```text
classification = calibration_failed_no_eligible_GR0_case
selected_amplitude = null
holdout_execution_authorized = false
SGBL_outcome_read = false
FGCQR_outcome_read = false
mechanism_question_answered = false
retained_EFT_evolution_authorized = false
```

This is neither a candidate trajectory nor a failed physical prediction. The
SGB-L and FGC-QR branches were never opened. It is a completed GR-0 numerical
calibration whose two declared amplitudes failed the same frozen comparator
admission.

## What changed relative to PROTO11

PROTO11 stopped near the beginning of the amplitude-`5/2` evolution because
the prior source-evaluation graph exhausted its binary-halving source-retry
contract. SRC2, SRC3, and SRC4 isolated and replaced only that numerical
evaluation graph while preserving the unredefined GR-0 equations and strict
raw `1e-12` residual gate. PROTO12 then evolved amplitude `5/2` through 25
common events to `t=3/2` and amplitude `3` through 24 common events to
`t=23/16`.

Every serialized member source-retry count is zero. Thus the earlier affine
source-arithmetic wall is cleared for both amplitudes under SRC4. This does not
prove the continuum equations healthy in every state; it proves only that the
specific PROTO11 arithmetic obstruction did not recur over these two frozen
trajectories.

## Source retries are not CFL step reductions

The terminal member diagnostics also contain nonzero cumulative
`CFL_retry_count` values:

| amplitude | RK4-2049 | RK4-4097 | RK4-8193 | SSPRK3-2049 | SSPRK3-4097 | SSPRK3-8193 |
|---|---:|---:|---:|---:|---:|---:|
| `5/2` | 227 | 455 | 149 | 189 | 377 | 147 |
| `3` | 217 | 443 | 161 | 185 | 358 | 187 |

These counters record permitted adaptive reductions of proposed timesteps to
satisfy the frozen CFL rule. They are not source-solver retries, and neither
amplitude terminated by exhausting the CFL-retry budget. The tracked result
therefore does not say “zero rejected numerical proposals.” It says, more
precisely:

```text
serialized_source_retry_count = 0
adaptive_CFL_step_reductions_occurred = true
campaign_stop_was_not_CFL_retry_exhaustion = true
```

## The exact terminal veto

At every preterminal event, RK4 and SSPRK3 both passed the complete PROTO12
common-event admission. At each terminal event:

- RK4 passed its constraint and pairwise spatial-spectrum gates;
- SSPRK3 passed alignment, coarse and fine absolute constraint guards, and
  monotone three-grid refinement;
- SSPRK3 passed the complete PROTO12 pairwise spatial-spectrum admission;
- only SSPRK3's minimum finest-pair constraint order failed; and
- `radial_momentum` was the only finite-order constraint component below the
  unchanged minimum `3/2`.

The endpoint values are:

| amplitude | last passing radial order | terminal radial order | shortfall below `3/2` | terminal RK4 minimum order |
|---|---:|---:|---:|---:|
| `5/2` | 1.5024552112430836 | 1.494957481780717 | 0.005042518219283032 | 3.2606752232454888 |
| `3` | 1.522307748775793 | 1.4990947384363291 | 0.0009052615636708783 | 3.2663994781424783 |

The final amplitude-`5/2` SSPRK3 component orders were Hamiltonian
`1.732974569918137`, radial momentum `1.494957481780717`, gauge-t
`2.0168292274237456`, and gauge-r `1.9676995306631015`. The corresponding
amplitude-`3` values were `1.690211681563169`, `1.4990947384363291`,
`1.9919927231725072`, and `1.9644243027261274`. Every reduction-constraint
pair was enclosed by the accumulated binary64 roundoff allowance.

The radial-momentum orders over the last six events decrease as follows:

```text
amplitude 5/2:
1.6043050187351453
1.5756441352509614
1.5489849605973036
1.5241421047640726
1.5024552112430836
1.494957481780717

amplitude 3:
1.6276453690940582
1.6025218682802533
1.5738425370783709
1.547166329439042
1.522307748775793
1.4990947384363291
```

This is why neither shortfall may be rounded into a pass. The gate is
`p >= 3/2`, not “approximately second order” or “close enough.”

## Independent order reconstruction

The frozen admission does not compute order from the raw medium/fine norms
alone. For each grid and constraint component it first forms

```text
effective norm = max(raw owned-domain norm - accumulated roundoff enclosure, 0)
```

and then computes the finest-pair order from the effective norms. PREF13
independently repeats that operation from the serialized raw norms and
enclosures and requires binary64 equality with every stored component order.

For example, the terminal amplitude-`3` SSPRK3 radial-momentum raw norms are

```text
0.01782165048810245
0.0060603980047813765
0.0021440219481765766
```

and the respective enclosures are

```text
1.179614628199488e-9
2.296474121976644e-9
3.5915945773012936e-9.
```

After those subtractions, the reconstructed order is exactly
`1.4990947384363291`. A naïve raw-norm ratio is not the declared evaluator.

## What two matching vetoes do and do not imply

The same configured gate failed in the same method and component for both
amplitudes. That makes the next question narrower and better instrumented. It
does **not** establish that the amplitudes share one physical cause, one
continuum defect, or even one numerical cause. The present three-resolution
records cannot decide whether the late-time SSPRK3 behavior is:

- a pre-asymptotic effect that clears on a genuinely finer nested pair;
- a persistent order reduction in this constraint/method combination;
- an interaction between the owned-domain boundary and the late pulse;
- or another reproducible numerical mechanism.

Those alternatives require new data fixed in advance. Inferring the answer
from the near-threshold value would be post-outcome fitting.

The finest trapped scores remain negative at both endpoints:

| amplitude | RK4 | SSPRK3 |
|---|---:|---:|
| `5/2` | -0.0769817625296822 | -0.07698176427149771 |
| `3` | -0.07414441277408346 | -0.07414441636828581 |

No qualified trapped event was observed. These values are GR-0 calibration
observables, not FGC-QR regulator or defocusing outcomes.

## Required successor boundary

The smallest justified next gate is **FGC-1-RSP2-FRZ1**, a prospective
higher-ladder constraint study, not PROTO13 and not SGB-L/FGC-QR execution.
RSP2 should freeze amplitude `3` as the nearest-threshold, lower-cost witness
and test the late SSPRK3 radial-momentum constraint on at least one finer
nested resolution. Before reading that outcome it must fix:

- the unchanged `3/2` minimum order and its pass/fail rule;
- the exact late common-event interval and endpoint;
- the unchanged equations, physical input, reference map, SRC4 source backend,
  CFL rule, and owned-domain constraint definition;
- the new resolution ladder and memory/runtime budget;
- the classification of a cleared, persistent, or invalid/nonconverged result;
  and
- the rule that no threshold reduction, fitted asymptote, or post-outcome
  parameter change is allowed.

Only that outcome-neutral study can determine whether a PROTO13 calibration
design is justified. Until then, no GR-0 case is eligible.

## Scientific boundary

CAL9/PREF13 is tangible numerical progress: it demonstrates that SRC4 removed
the prior source-arithmetic wall over both long calibration trajectories, then
localizes the new obstruction to one late SSPRK3 radial-momentum convergence
gate. It preserves the raw evidence, threshold, and failure.

It does not select a calibration amplitude, open SGB-L or FGC-QR, test
regulator activation, derive affine metric-null defocusing, reject FGC-QR or
the general gradient programme, validate a retained EFT, resolve a singularity,
derive a child domain or dark sector, or vary a locally measured speed of
light.

## Reproduction

```bash
python3 scripts/reproduce_fgc_cal9_pref13.py --check
```
