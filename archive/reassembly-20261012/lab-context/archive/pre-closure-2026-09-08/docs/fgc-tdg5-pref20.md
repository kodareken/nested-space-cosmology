# FGC-1-TDG5-PREF20: stage-complete temporal-refinement theorem

`FGC-1-TDG5-PREF20` independently executes the theorem frozen by
`FGC-1-TDG5-FRZ1` at immutable commit `1716cc4…`. It consumes the compact
freeze record and tracked numerical-engine source only. It does not load a
PROTO13 history or checkpoint, resume the terminal campaign, advance a state,
or read an SGB-L or FGC-QR trajectory.

## Result

The frozen finite-dimensional construction is algebraically sufficient to
remove the sampling-kernel ambiguity identified by TDG4 **inside its declared
numerical class**, and the existing step proposal exposes the endpoint data
needed to implement it.

For a normalized accepted interval, write

```text
H(theta) = a0 + a1 theta + a2 theta^2 + a3 theta^3,
theta in [0,1].
```

The independent binder expands the four Hermite basis functions rather than
reusing the freeze's matrix inverse. Fraction-preserving Gaussian elimination
then gives

```text
det(S) = 1,
det(S^-1) = 1,
||S||_infinity = 6,
||S^-1||_infinity = 9,
kappa_infinity(S) <= 54.
```

All four rational unit-data records round-trip exactly. Thus the endpoint
values and endpoint scaled slopes determine one and only one member of the
declared cubic class.

## Exact restriction and continuous maximization

If `P(theta)=a+b theta+c theta^2+d theta^3`, restriction to the left and
right normalized half intervals is linear in `(a,b,c,d)`. The independently
derived restriction matrices both have determinant `1/64`. Exact basis and
composition controls verify

```text
P_left(s)  = P(s/2),
P_right(s) = P((1+s)/2)
```

for `s in [0,1]`. The difference between a coarse cubic and either fine-half
cubic is therefore itself cubic.

On a compact half interval, the absolute value of a cubic attains a maximum.
Unless the polynomial is identically zero, every nonzero interior maximizer
is a root of its derivative. Because that derivative has degree at most two,
the complete candidate set contains at most four locations:

```text
the two endpoints + every real interior derivative root.
```

The exact endpoint-only attack

```text
4 theta (1-theta)
```

has equal zero endpoint values and an interior absolute maximum of one at
`theta=1/2`; its endpoint slopes expose it. A genuine cubic control has two
distinct interior stationary points, demonstrating that both roots must be
retained.

This proves the finite candidate theorem. It does **not** implement the future
binary64 root solver or its outward roundoff enclosure. That implementation is
an explicit successor obligation.

## Method-owned refinement debit

Suppose a method of formal order `p` is in an asymptotic local-error regime
with one shared leading coefficient. One step of width `Delta t` carries
leading coefficient `1`, while two half steps carry `1/2^p`. Their difference
therefore carries `(2^p-1)/2^p`, and the leading fine-path error is the
difference divided by `2^p-1`.

The method-owned factors are

```text
RK4:    p=4, denominator 15
SSPRK3: p=3, denominator 7.
```

This is a numerical refinement debit under stated assumptions. PREF20 does
not prove that a future trajectory is already asymptotic, and it does not
turn the factor into a rigorous global PDE error bound.

## Existing engine compatibility

The binder parses the immutable numerical engine's syntax tree and verifies:

- the complete `StageRecord(stage_name,time,state,rhs)` shape;
- the complete `StepProposal` shape;
- all RK4 and SSPRK3 internal-stage evaluations;
- an independent `candidate_endpoint` evaluation at `final_time` on the
  candidate state;
- a fresh right-hand-side call inside that evaluation; and
- retention of the resulting `StageRecord`.

The frozen runtime pair can therefore be built around the existing proposal
shape. This source-level compatibility result is not the implementation of an
atomic coarse/fine transaction.

## What is now authorized

PREF20 authorizes one bounded successor:

> Implement and synthetically validate the frozen same-grid one-full-step and
> two-half-step refinement pair, including the complete continuous cubic
> difference and an outward binary64 arithmetic debit.

That implementation must preserve the bitwise initial state, source,
projector, spatial operator, and every inherited stage-level source, health,
CFL, boundary, scale, and transaction guard. The coarse path remains evidence
only. The two fine half steps must eventually commit atomically or neither may
commit.

Threshold selection remains a separate prospective freeze. No production
trajectory is authorized by this binder.

## Scientific boundary

The declared cubic is a numerical continuous extension, not a claim that the
exact PDE history is cubic. Finite-dimensional injectivity prevents hidden
freedom inside that class; it does not establish PDE accuracy. Spatial
convergence, cross-method agreement, source and constraint control, and the
complete later DEF1 error ledger remain mandatory.

The historical PROTO13/CAL10 stop is preserved and unreclassified.

```text
TDG5_stage_complete_refinement_theorem_completed = true
runtime_refinement_pair_implementation_authorized = true
runtime_refinement_pair_implemented = false
runtime_thresholds_frozen = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

The next result must be an implementation and synthetic-validation artifact.
It may not select a threshold, reopen the old terminal checkpoint, define
PROTO14, or answer the candidate mechanism question.
