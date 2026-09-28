# FGC-1-TDG5-FRZ1: stage-complete temporal refinement design freeze

`FGC-1-TDG5-FRZ1` freezes the first prospective replacement architecture
authorized by TDG4-PREF19. It consumes only tracked evidence at immutable
commit `b866d24…`. It does not load the PROTO13 checkpoint or temporal arrays,
resume the terminal campaign, advance a state, or read SGB-L or FGC-QR data.

## Why this is a different measurement

TDG4 proved that values at any finite collection of times leave an
infinite-dimensional sampling kernel on an unrestricted smooth function
space. Another 64, 128, or 16,384 point-value spectrum would therefore repeat
the same logical mistake if it were promoted into a continuum enclosure.

TDG5 does not attempt to reconstruct an arbitrary smooth history. It declares
the continuous numerical object carried by one accepted time step. On the
normalized interval `theta in [0,1]`, every monitored scalar component is a
cubic Hermite polynomial

```text
H(theta) = a0 + a1 theta + a2 theta^2 + a3 theta^3,
```

determined by

```text
(H(0), H'(0), H(1), H'(1)).
```

Here “stage-complete” does not mean that the interpolant treats the internal
Runge--Kutta stage states as independent interpolation nodes. It means that
every inherited internal-stage and candidate-endpoint source/health premise
still has to pass, while the continuous complete-state record itself is fixed
by accepted endpoint values and freshly evaluated endpoint right-hand sides.
This distinction is part of the frozen runtime contract.

The physical endpoint derivatives are multiplied by the step width before
entering this normalized data vector. For the monomial coefficient vector, the
sampling matrix and its inverse are

```text
S = [[1, 0, 0, 0],       S^-1 = [[ 1,  0,  0,  0],
     [0, 1, 0, 0],               [ 0,  1,  0,  0],
     [1, 1, 1, 1],               [-3, -2,  3, -1],
     [0, 1, 2, 3]]               [ 2,  1, -2,  1]].
```

The exact controls give

```text
det(S) = 1,
||S||_infinity = 6,
||S^-1||_infinity = 9,
kappa_infinity(S) <= 54.
```

There is therefore no invisible between-record direction inside the declared
piecewise-cubic class. This is the finite-dimensional sufficient route named
by TDG4. It is not a statement that the exact PDE history is cubic.

## Why endpoint values are still insufficient

The exact adversarial control compares zero with

```text
H(theta) = 4 theta (1-theta).
```

Both have value zero at `theta=0` and `theta=1`, while the second reaches one
at the midpoint. Its normalized endpoint slopes are `+4` and `-4`, so the
stage-complete Hermite record distinguishes it exactly. A replacement that
stored only accepted endpoint values would retain the defect TDG4 exposed.

## The same-grid temporal discriminator

For every proposed production interval of width `Delta t`, the successor
runtime must construct two paths from one bitwise-identical accepted state and
one identical spatial operator:

```text
coarse:  one step of width Delta t
fine:    two steps of width Delta t/2.
```

Both paths must independently pass the unchanged source, health, CFL,
boundary, scale, and transaction premises. The coarse path is evidence only.
Only the two-half-step fine path may commit, and its two transactions commit
atomically or neither does.

The coarse Hermite polynomial is restricted to each half interval and
subtracted from the corresponding fine polynomial. Each difference is cubic,
so its continuous maximum occurs at an endpoint or at a real root of its
quadratic derivative. The runtime must evaluate all such locations and add an
outward binary64 arithmetic debit. It may not infer the maximum from another
finite sample grid.

The monitored vector is the complete dimensionless `(u,p,q)` state on the
owned interior rows. The regular centre and four projector-owned outer rows
remain covered by their existing exact projector, regular-centre, constraint,
and causal-boundary gates. The old normal-flow tracer spectrum remains public
as a diagnostic but is no longer an admission veto.

## Method order and error semantics

The no-history preflight checks the exact Butcher order conditions for the two
existing methods:

```text
RK4:    formal order 4, step-doubling denominator 2^4-1 = 15
SSPRK3: formal order 3, step-doubling denominator 2^3-1 = 7.
```

Exact `y'=y` controls at rational step sizes confirm that both one-step versus
two-half-step differences contract under refinement. These facts justify the
method-owned Richardson debit used by a future runtime design. They do not
turn that debit into a rigorous global PDE error bound or prove that an
unknown physical trajectory is already in its asymptotic regime.

The future admission must therefore retain the independent protections that
already exist:

- each method's spatial ladder and constraint convergence;
- RK4 versus SSPRK3 agreement;
- every stage-level source, branch, health, CFL, and causal stop;
- a threshold frozen from independent controls rather than from the observed
  trapped or Raychaudhuri margin;
- and the complete later DEF1 error ledger.

## Alternatives considered

A full a posteriori Gronwall estimate for the semidiscrete system would require
a useful nonlinear energy bound for a state with roughly hundreds of thousands
of degrees of freedom. A naive infinity-norm Lipschitz constant grows with
inverse grid spacing and is expected to exponentiate into a nondiscriminating
bound. That route is deferred until a sharper energy estimate exists.

An independently bounded continuum derivative would also satisfy TDG4, but
the present GR-0 health contract does not bound every required derivative
between steps. Additional finite sampling is forbidden by TDG4, and
endpoint-only interpolation fails the exact bump control. The selected design
reuses the existing RK stages and transaction machinery, introduces no runtime
dependency, and isolates time error without changing the spatial operator.

## Freeze boundary

At this checkpoint the replacement architecture and its exact no-history
controls are frozen, but its independent theorem/result binder has not yet
executed. No runtime pair or threshold exists.

```text
TDG5_replacement_design_frozen = true
TDG5_exact_no_history_controls_pass = true
TDG5_stage_complete_refinement_theorem_completed = false
runtime_refinement_pair_implementation_authorized = false
runtime_refinement_pair_implemented = false
runtime_thresholds_frozen = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

After this freeze is committed, only a separately hash-bound TDG5 result may
independently rederive the algebra and decide whether implementation of the
paired runtime may begin. It may not reclassify CAL10, resume the old terminal
campaign, define PROTO14, or answer a candidate or physical question.
