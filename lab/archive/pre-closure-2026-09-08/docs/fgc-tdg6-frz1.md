# FGC-1-TDG6-FRZ1: temporal-admission threshold and production-contract freeze

`FGC-1-TDG6-FRZ1` spends only the threshold-design permission earned by
`FGC-1-TDG5-IMP1` at immutable commit `986ceb1…`. It reads the compact tracked
TDG5 result and the inherited protocol/runtime sources. It does not load a
PROTO13 history or checkpoint, resume CAL10, advance a production state, or
read SGB-L or FGC-QR data.

## Decision

The replacement temporal admission will not use a fitted absolute state
tolerance. It will ask two separate questions:

1. Does another same-grid refinement reduce the complete continuous numerical
   path difference with at least the already-declared order `3/2`?
2. How large is the remaining finest-pair uncertainty, channel by channel,
   when it is carried forward without cancellation?

The first is a numerical admission gate. The second is an error debit. A large
debit is never hidden by declaring a sufficiently generous state tolerance; it
must later be paid by the four-times-combined-error trapped or Raychaudhuri
sign margin. Conversely, an arbitrarily tiny difference does not pass if it
fails the convergence gate.

This removes the last outcome-dependent choice from the threshold design. No
number is selected from CAL10, a trapped sphere, or a desired defocusing
margin.

## Three-level same-grid discriminator

For one macro interval of width `Delta t`, all levels start from one bitwise-
identical accepted state and use the same method, source, projector, spatial
operator, and physical equations:

```text
level 0:  one step of width Delta t
level 1:  two steps of width Delta t/2
level 2:  four steps of width Delta t/4.
```

Level 0 and level 1 are evidence only. Only the complete four-quarter-step
level 2 path may commit, and all four fine proposals commit atomically or none
does. The proposal counts are therefore `1 + 2 + 4 = 7`. With the inherited
fresh candidate-endpoint source record, the complete shadow and committed
fine accounting is:

| method | all shadow stage records | committed fine stage records |
|---|---:|---:|
| RK4 | 35 | 20 |
| SSPRK3 | 28 | 16 |

Every path retains the unchanged source, health, CFL, boundary, scale, and
transaction gates. The future production compositor must also shadow the
normal-flow tracers. Preparation may mutate none of the accepted member,
monitor, causal, tracer, retry, or temporal-debit ledgers. Commit revalidates
all of them and adopts only the level-2 state.

The macro step is still proposed from the inherited CFL and target-event
distance. A failed temporal-order or order-inconclusive assessment may halve
that macro step under a separate counter, the inherited factor `1/2`, the
prospectively fixed maximum `32`, and minimum step `2^-30`. This is a numerical
recovery path, never a physical classification. Source-only and non-source
failures keep their existing ownership and priority.

## Continuous channel intervals

TDG5 already compares a coarse Hermite cubic with two fine-half cubics on the
complete owned `(u,p,q)` state. TDG6 adds one more level. It compares:

```text
D01: level 0 against level 1,
D12: level 1 against level 2.
```

The comparison remains continuous inside the declared piecewise-cubic
numerical class. The future binder must derive the quarter-interval restriction
maps rather than infer maxima from more samples. The production implementation
must retain a certified nonnegative interval for every one of the 18 channels

```text
(u,p,q) x (alpha,v,lambda,R,phi,chi),
```

after maximizing over every owned row and subinterval. Write those intervals
as

```text
D01 in [L01,U01],
D12 in [L12,U12].
```

The regular centre and four projector-owned outer rows remain under their
existing exact centre/projector, constraint, and causal-boundary contracts.

## Exact `3/2` threshold

For resolvable nonzero `D01`, an observed order of at least `3/2` requires

```text
D12 <= D01 / 2^(3/2).
```

The outward interval test uses the fine upper and coarse lower bounds. Squaring
the nonnegative quantities removes an irrational floating threshold:

```text
8 U12^2 <= L01^2.
```

Equality passes. The exact rational control at the boundary passes, while a
rational perturbation of `2^-52` below it fails. The gate therefore cannot be
changed by decimal rendering or by rounding a near miss.

Each channel is classified as exactly one of:

- `exact_zero`: both certified upper bounds are zero; no infinite order is
  invented and the debit is zero;
- `enclosure_dominated_debit_only`: both intervals contain zero and the fine
  upper bound does not exceed the coarse upper bound; no order is claimed and
  the complete fine upper bound remains a debit;
- `resolved_order_pass`: the exact squared inequality passes;
- `resolved_order_failure`: the inequality fails and the macro step may be
  retried at half width; or
- `order_inconclusive`: the outward intervals do not support either the
  enclosure-dominated class or a resolved lower-bound comparison, so the macro
  step may be retried and cannot commit.

This gate is deliberately sensitive to convergence rather than magnitude. An
exact control with a tiny but nonconvergent difference fails. A large but
convergent control passes the order gate while retaining the large debit that
will prevent an unsupported sign claim.

## Why there is no absolute temporal pass number

A fixed `atol + rtol*scale` controller would require choosing the scale and
tolerance before the constraint-to-Raychaudhuri stability map exists. Choosing
it from the observed collapse or sign margin would be outcome-driven. Choosing
it from the old CAL10 values would merely refit the failed instrument.

TDG6 therefore defines the per-channel accepted-step debit as

```text
E_temporal,channel = U12.
```

Accepted macro-step debits accumulate by nonnegative addition without
cancellation. The formal RK4 and SSPRK3 Richardson divisions by `15` and `7`
remain public conditional diagnostics, but they are not used to shrink the
admission debit. At a common event, both methods must still pass their own
spatial and constraint ladders and agree through their separately enclosed
observables. The future trapped and DEF1 sign ledgers must include the temporal
vector through the independently derived stability map.

This is conservative, but it is not advertised as a rigorous global PDE error
bound. Propagation can amplify local perturbations, so the later stability map
remains mandatory. A large accumulated debit makes a result inconclusive; it
does not become evidence for or against FGC-QR.

## Alternatives considered

| route | decision | reason |
|---|---|---|
| Fixed absolute/relative state tolerance | Rejected | It introduces an arbitrary scale before the observable stability map and can be tuned to the desired outcome. |
| One full step versus two half steps only | Rejected as the final gate | It supplies a difference but cannot test the asymptotic assumption used by its Richardson denominator. |
| More temporal-history samples | Rejected | TDG4 proves that finite point samples retain an unrestricted smooth sampling kernel. |
| Full nonlinear a posteriori PDE energy estimate | Deferred | It would be stronger, but no useful nonexploding semidiscrete stability estimate presently exists for this system. |
| Three same-grid levels plus an additive debit | Selected | It tests contraction prospectively, keeps magnitude public, changes no physical equation, and connects directly to the existing four-times-error decision rule. |

The selected route is computationally expensive: seven shadow proposals are
evaluated per accepted macro interval and the four-quarter path is the only
one retained. That cost is explicit rather than concealed. It buys a
trajectory-local convergence discriminator instead of assuming that a method's
formal order automatically applies to the evolving solution.

## Freeze boundary

The threshold and replacement-admission semantics are now prospectively
defined. They have not yet been independently rebound or implemented in the
production compositor.

```text
TDG6_threshold_and_admission_design_frozen = true
runtime_thresholds_frozen = true
replacement_temporal_admission_defined = true
TDG6_independent_binder_completed = false
production_compositor_implementation_authorized = false
production_compositor_implemented = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

After this freeze is committed, only a separately hash-bound TDG6 binder may
independently derive the quarter restriction, interval-order, channel/debit,
transaction, and control results and decide whether production implementation
may begin. CAL10 remains historical and unreclassified. No collapse,
defocusing, transition, retained-EFT, or physical claim follows.
