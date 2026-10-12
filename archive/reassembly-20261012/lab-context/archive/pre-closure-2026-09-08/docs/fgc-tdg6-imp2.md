# FGC-1-TDG6-IMP2: production-shaped temporal compositor

`FGC-1-TDG6-IMP2` implements and synthetically qualifies the bounded runtime
object authorized by `FGC-1-TDG6-PREF21` at immutable commit `b42f954…`. It
reads no PROTO13 trajectory array or campaign checkpoint, resumes no campaign,
and evaluates no SGB-L or FGC-QR state.

## Result

The frozen three-level temporal admission now exists as executable runtime code
for both inherited methods. From one bitwise-identical accepted state, one
macro interval of width `Delta t` is represented by seven independent shadow
proposals:

```text
outer:   1 x Delta t
medium:  2 x Delta t / 2
fine:    4 x Delta t / 4.
```

Every level has an independent PROTO7 monitor ledger, causal ledger, and
normal-flow tracer clone. Tracer preview is part of every shadow proposal's
preaccept boundary. Shadow acceptance mutates only the corresponding clone.
Preparation therefore leaves the real accepted state, monitor, causal,
tracer, and temporal-debit ledgers unchanged.

The synthetic controls execute all seven proposals with both numerical
methods:

```text
RK4:    35 total shadow stage records; 20 fine records committable
SSPRK3: 28 total shadow stage records; 16 fine records committable.
```

The outer and medium paths remain evidence only. After revalidating coordinate
time, state hash, step index, transaction serial, monitor state, causal state,
complete tracer snapshot, and temporal ledger, the commit boundary adopts only
the four-quarter-step fine endpoint and its ledgers. A late exception during
tracer commit restores the prior monitor, causal, and complete tracer state.

## Continuous admission

For each of the 18 ordered dimensionless channels

```text
(u,p,q) x (alpha,v,lambda,R,phi,chi),
```

the runtime constructs continuous cubic-Hermite differences on the owned
radial rows. `D01` compares the outer history against the two medium halves;
`D12` compares the medium history against the four fine quarters. The existing
TDG5 envelope evaluates endpoints and every numerically real stationary-point
candidate and retains its independent Bernstein convex-hull certification.

The lower bound subtracts the complete coefficient-construction and evaluation
debits with exact rational arithmetic. The upper bound is the full certified
continuous envelope. The frozen classifier consumes the exact rational images
of those binary64 bounds and applies

```text
8 U12^2 <= L01^2.
```

Admission is all-of: one failed channel vetoes the macro step. An admitted
channel contributes its complete `D12` upper bound to a componentwise,
nonnegative, outward-accumulated temporal debit. There is no Richardson
division, cancellation, absolute state tolerance, or normalization by an
observed physical signal.

The synthetic controls distinguish magnitude from convergence. An exact-zero
channel passes without an order claim; a deliberately tiny but resolved
nonconvergent channel vetoes the whole step; and a large convergent difference
passes while retaining a correspondingly large debit. A dense polynomial
audit independently lies between the runtime lower and upper bounds.

## Retry and checkpoint semantics

A failed temporal admission requests a prospective half-step retry. This
counter is separate from PROTO7 source refinement and CFL ownership. Before a
retry ledger can exist, a complete rejection record containing all 18 channel
decisions must be written through the caller's durable sink. A failed sink
therefore leaves the immutable input ledger untouched.

The numerical recovery is bounded by 32 retries per macro step and a minimum
step of `2^-30`. Exhaustion is a typed numerical stop, never a physical
classification. A later accepted macro step resets only the current-step
counter; cumulative counts and the complete canonical rejection history remain.

The checkpoint extension round-trips:

- the ordered 18-channel accumulated-debit vector and its content hash;
- accepted-step and current/cumulative retry counts; and
- every canonical serialized temporal rejection.

Hash, channel-order, schema, missing-field, and claim-promotion mutations fail
closed.

## Failure ownership

A source-only failure on any shadow remains owned by the unchanged PROTO7
source retry. An inherited non-source health, CFL, boundary, scale, transaction,
or tracer failure remains terminal evidence at its named shadow path. Neither
class is relabeled as TDG6 order failure. A temporal retry is created only after
all seven paths succeed and the complete continuous admission itself fails.

This separation matters: the compositor classifies its own numerical question
without taking ownership of failures belonging to another instrument layer.

## What the implementation earns

The production-shaped compositor and its synthetic qualification are complete.
The next bounded task is a prospective `PROTO14` freeze that integrates this
runtime into the GR-0 calibration protocol, defines raw/checkpoint namespaces,
and binds restart and terminal classification before any trajectory is read.

```text
TDG6_independent_binder_completed = true
production_compositor_implementation_authorized = true
production_compositor_implemented = true
production_compositor_synthetic_qualification_passed = true
PROTO14_freeze_design_authorized = true
production_trajectory_authorized = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
DEF1_execution_authorized = false
```

The historical PROTO13/CAL10 result remains valid and unreclassified. The
piecewise cubic is a numerical extension rather than the exact PDE history; a
three-level pass does not prove trajectory-wide asymptotics; and the local
debit is not a global PDE or Raychaudhuri error theorem. Spatial ladders,
cross-method agreement, the constraint-to-Raychaudhuri stability map, candidate
execution, DEF1, transition, and every cosmological interpretation remain
separate future work.
