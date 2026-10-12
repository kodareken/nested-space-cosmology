# FGC-1-TDG6-PREF21: independent temporal-admission binder

`FGC-1-TDG6-PREF21` independently checks the temporal-admission design frozen
by `FGC-1-TDG6-FRZ1` at immutable commit `dff13cf…`. It consumes the compact
freeze certificate and immutable runtime source only. It does not import the
freeze's design module as a theorem, load a PROTO13 history or checkpoint,
resume CAL10, advance a state, or read SGB-L or FGC-QR data.

## Result

All independent controls pass. The frozen three-level rule is algebraically
well defined inside its declared piecewise-cubic numerical class, its
conservative interval test really implies at least three-halves contraction,
its 18-channel debit is an all-of noncancelling ledger, and the immutable
runtime exposes the extension points needed to implement the compositor.

This authorizes one bounded successor:

> Implement and synthetically qualify the frozen TDG6 three-level GR-0
> production compositor, including independent monitor, causal, and tracer
> shadows; complete continuous `D01` and `D12` intervals; fine-only atomic
> commit; typed temporal retry; and checkpoint round trips of every new ledger.

It does not authorize a production trajectory or PROTO14.

## Direct quarter restriction

Let

```text
P(theta) = sum_(n=0)^3 a_n theta^n.
```

On quarter `j` with local coordinate `s`, substitute

```text
theta = (j+s)/4,  j in {0,1,2,3}.
```

Expanding the binomial independently gives the coefficient map

```text
b_k = sum_(n>=k) a_n C(n,k) j^(n-k) / 4^n.
```

Every resulting matrix is upper triangular with diagonal

```text
(1, 1/4, 1/16, 1/64),
```

so each determinant is `1/4096`. Fraction-preserving elimination independently
returns the same determinant. All four maps also agree exactly with composing
the corresponding two half-interval maps, and all 16 cubic basis/quarter
controls agree at five rational probe points.

The determinant is nonzero, so restriction loses no coefficient information
inside the declared cubic class. The continuous comparisons therefore have a
finite exact organization:

```text
D01: outer cubic against two medium half cubics;
D12: two medium half cubics against four fine quarter cubics.
```

The existing outward continuous cubic-envelope evaluator can be applied on
each subinterval. No extra time samples are being treated as a continuum
proof.

## Conservative interval theorem

For nonnegative true differences satisfying

```text
d01 in [L01,U01],
d12 in [L12,U12],
L01 > 0,
```

the frozen test

```text
8 U12^2 <= L01^2
```

implies

```text
8 d12^2 <= 8 U12^2 <= L01^2 <= d01^2.
```

Thus `d12/d01 <= 2^(-3/2)` for every value admitted by the outward intervals.
The factor eight is independently recovered from
`2^(2p)` at `p=3/2`. Equality passes exactly; reducing the right-hand squared
boundary by the rational amount `2^-52` fails.

The special classifications also retain their frozen meanings:

- exact zero passes with zero debit and no invented infinite order;
- an enclosure-dominated pair passes only when both intervals include zero
  and the fine upper bound does not grow;
- a resolved order failure or inconclusive interval requests numerical retry;
- no absolute magnitude or physical signal enters the classification.

## Complete channel and debit ownership

The binder reconstructs the exact ordered 18-channel state

```text
(u,p,q) x (alpha,v,lambda,R,phi,chi).
```

Admission is all-of: mutating one otherwise passing channel to a resolved
failure vetoes the complete step. The retained debit is the full fine-pair
upper bound in every channel. Accepted-step vectors accumulate by componentwise
nonnegative addition, so a later positive and negative perturbation cannot
cancel numerically and conceal uncertainty.

The debit is not divided by `15` or `7` for admission. It is also not promoted
to a global PDE error bound. A future independently derived stability map must
propagate it into the trapped-sphere and complete Raychaudhuri error ledgers.

## Immutable runtime feasibility

The binder parses the immutable source syntax and verifies:

- TDG5's prepared pair retains the initial state hash, time, step and serial,
  original monitor and causal state, and fine monitor and causal state;
- TDG5 preparation creates two transaction clones and three shadow attempts,
  then proves that the accepted state and real transaction did not move;
- TDG5 commit revalidates the state/time/step/serial/monitor/causal boundary
  and adopts only the fine monitor and causal state;
- PROTO7 exposes tracer preview, commit, and rollback inputs; and
- PROTO13's atomic checkpoint already owns the field, monitor, causal, tracer,
  event, and state-hash surfaces that the successor must extend.

One missing piece is deliberately visible: TDG5's present shadow attempt uses
`preaccept=None`, while tracer previews live in the runner. PROTO13 also has no
TDG6 temporal-debit vector or temporal-retry ledger. The production successor
must integrate those states; the binder does not pretend they already exist.

Source-shape feasibility means that the extension can be made without changing
the physical equations or discarding the existing transaction model. It is
not a runtime proof and not an implementation.

## Mutation controls

The canonical certificate must fail closed under mutations of:

- the immutable freeze commit or any bound hash;
- the quarter determinant or direct/composed-map requirement;
- the squared multiplier, equality rule, or physical-signal nonnormalization;
- the complete channel order or all-of semantics;
- fine-only commit ownership;
- tracer or checkpoint extension obligations;
- production-trajectory authorization;
- historical CAL10 preservation; or
- any candidate, retained-EFT, or physical promotion.

Runtime-source attacks that remove tracer preview, redirect the fine monitor
commit, or remove the checkpoint writer also make the independent feasibility
decision false.

## Scientific boundary

The theorem applies to the declared numerical extension, not the exact PDE
history. A three-level contraction result does not by itself prove a complete
trajectory is asymptotic. A local accumulated debit is not a rigorous global
error estimate. Spatial convergence, independent-method agreement, constraints,
causal isolation, physical health, and the later DEF1 stability ledger remain
mandatory.

The historical PROTO13/CAL10 stop remains unchanged.

```text
TDG6_threshold_and_admission_design_frozen = true
TDG6_independent_binder_completed = true
production_compositor_implementation_authorized = true
production_compositor_implemented = false
production_trajectory_authorized = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

The next artifact must implement and synthetically qualify the compositor. It
may not mutate the threshold, inspect the old terminal trajectory as a new
result, freeze PROTO14, or answer the collapse/defocusing question.
