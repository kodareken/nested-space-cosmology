# FGC-1-PRO14-FRZ1: replacement-temporal calibration freeze

**Status:** prospective, machine-reproduced protocol freeze. It advances no
trajectory and authorizes no calibration or candidate execution.

```text
PROTO14_frozen = true
PROTO14_replacement_temporal_admission_frozen = true
PROTO14_successor_runtime_implemented = false
PROTO14_fresh_GR0_dynamic_calibration_authorized = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

## Question owned by this artifact

The historical PROTO13 run reached `t=63/16` with both method-owned spatial
and constraint admissions passing, then stopped at the first 64-sample
temporal-spectrum gate. CAL10/PREF15 bound that result without identifying its
cause. TDG4 later proved that finite samples cannot bound the declared
continuum top-band functional on the unrestricted smooth class. TDG5 and TDG6
then defined, independently bound, implemented, and synthetically qualified a
finite-dimensional replacement.

PROTO14 answers the next prospective design question:

> What exact GR-0 calibration may apply the TDG6 one/two/four same-grid
> admission to every accepted post-restart macro step without reclassifying
> PROTO13, weakening any spatial or physical premise, or opening a candidate
> branch?

The answer is frozen before the successor runtime and before either new output
namespace exists.

## What changes—and what does not

Only the temporal admission and output namespace change. Every unlisted
PROTO13 premise is inherited by the exact predecessor hash.

The old 64-sample normalized temporal spectrum remains serialized as a public
diagnostic. It can neither admit nor reject PROTO14. The replacement admission
uses the implemented `FGC-1-TDG6-IMP2` compositor on both RK4 and SSPRK3:

```text
outer:   1 x Delta t
medium:  2 x Delta t / 2
fine:    4 x Delta t / 4
```

All seven paths start from one bitwise-identical accepted state. Every one of
the 18 ordered `(u,p,q) x (alpha,v,lambda,R,phi,chi)` continuous cubic-Hermite
channels must pass the exact outward rule

```text
8 U12^2 <= L01^2,
```

or belong to the already frozen exact-zero or enclosure-dominated debit-only
class. Admission is all-of. Only the four-quarter fine path can commit.
Every finest-pair upper bound is accumulated componentwise and outward with no
cancellation, absolute state tolerance, or normalization by a physical signal.

A failed TDG6 admission is durably serialized before a prospective half-step
retry. The retry counter is independent of PROTO7 source refinement and the
ordinary CFL counter. It is bounded by 32 retries and the unchanged minimum
macro step `2^-30`. Exhaustion is an invalid or nonconverged numerical run,
never a physical classification.

## Restart boundary and its deliberate limitation

The same six method-owned members are restored at `t=23/16`:

```text
RK4/D4-2:      2049 -> 4097 -> 8193
SSPRK3/D2-1:   4097 -> 8193 -> 16385.
```

Their exact state, tracer, proper-time history, monitor, causal, step, stage,
transaction, source-retry, and CFL-retry records remain bound to the immutable
PROTO12 and RSP2 checkpoints. The terminal PROTO13 checkpoint at `t=63/16`
cannot resume.

PROTO14 initializes a new TDG6 ledger at zero at `t=23/16`. That ledger bounds
only the admitted post-restart numerical composition. This artifact does not
prove a temporal error bound from `t=0` to the restart. The restart arrays are
therefore fixed **discrete calibration inputs**, not a physical trajectory
certificate, and they may not seed an FGC-QR or DEF1 claim. If amplitude `3`
becomes eligible, it is eligible only for the later prerequisite work; the
eventual candidate experiment must close its own complete trajectory and error
contract.

## Common-event and terminal contract

The common-event schedule remains `1/16`, the final coordinate time remains
`32`, and recoverable checkpoints remain spaced by `1/4`. A selected GR-0 case
requires all of the following:

- every prior accepted macro step passed TDG6;
- both method-owned spatial and constraint admissions pass independently;
- each method supplies its own finest-pair Richardson interval;
- the two trapped-observable intervals overlap on the exact common physical
  nodes with the fixed future-null orientation;
- the trapped sign exceeds four times the combined spatial error; and
- eight consecutive common events meet the complete condition.

The accumulated temporal-debit vector is serialized at every common event,
but it is not yet a Raychaudhuri or global PDE error map. A selected result is
therefore a calibration input to later premise gates, not a physical result.

The only valid campaign classifications are:

```text
GR0_calibration_completed_selected_amplitude
calibration_failed_no_eligible_GR0_case
invalid_implementation_or_nonconverged_run
```

No trapped interval by `t=32` yields no eligible GR-0 case. A constraint,
spatial, health, boundary, scale, or typed numerical stop remains owned by its
premise and cannot be relabeled as physics. An unexpected error or exhausted
TDG6 retry yields an invalid/nonconverged run.

## Immutable lineage and temporal provenance

The freeze consumes immutable commit
`27badf9ad3c0bc3aa1c8a4153a42a26bdaf4c00f`. It verifies the exact PROTO13
protocol and freeze, CAL10/PREF15 result, TDG6-IMP2 configuration/result/runtime,
and both ignored restart checkpoints. It restores all six restart payloads
without interpolation or refitting.

Before generating the canonical result, it proves both new roots absent:

```text
runs/fgc-2-sf1/proto14/calibration
runs/fgc-2-sf1/proto14/holdout
```

The freeze creates neither path. That absence is temporal prelaunch evidence.
After an authorized successor creates a root, a later binder must verify this
freeze from its immutable commit rather than rerun the historical absence
observation.

## Claim boundary and successor

PROTO14 advances no production state, selects no GR-0 case, and answers no
mechanism question. It does not authorize SGB-L, FGC-QR, COL1, DEF1,
retained-EFT evolution, singularity resolution, a child domain, a dark-sector
mechanism, or varying locally measured light speed.

The successor runtime owner is **FGC-1-HLT12-MON12**. It must integrate the
immutable compositor into the six-member campaign, extend atomic checkpoints,
prove recovery and every typed failure route, validate the fresh namespace,
and separately authorize exactly one GR-0 calibration. Only a later
outcome-neutral post-run binder may classify that trajectory.

Before the successor namespace exists, reproduce the freeze with:

```bash
make fgc-pro14-frz1
```

The canonical result is `results/fgc-1-pro14-frz1.json`. The command advances
no state and creates no run directory.
