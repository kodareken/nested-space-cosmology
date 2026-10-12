# FGC-1-BND2-CP1: causal isolation of the spherical measurement region

## Decision

`FGC-1-BND2-CP1` passes the boundary-or-domain-of-dependence prerequisite for
the scoped classical diagnostic. It does so by the second route in that name:
the measured region is kept outside the outer boundary's continuum domain of
dependence, with a declared finite-difference stencil pad. It is explicitly
not a nonlinear constraint-preserving boundary theorem.

The machine decision is

```text
spherical_boundary_or_domain_of_dependence_control_passed = true
```

within that exact causal-isolation scope.

The certificate consumes BND1's exact frozen incoming-mode count, CON4's
physical/gauge/reduction subsidiary classification, ID1's exact outer vacuum
buffer, DOM4's run geometry, and HYP2's all-covector characteristic clusters.
No collapse or holdout trajectory is evaluated.

## From local-frame roots to coordinate speeds

HYP2 works in a physical orthonormal frame. Its separated real root clusters
imply the conservative local bounds

```text
physical: |z| <= 1 + 1/5 = 6/5
tilde:    |z| <= 1/2 + 1/20 = 11/20
hat:      |z| <= 1/3 + 1/20 = 23/60.
```

For the spherical ADM line element

```text
ds^2 = -N^2 dt^2 + Lambda^2 (dr + shift dt)^2 + R^2 dOmega^2,
```

the coordinate speed is `dr/dt=-shift+(N/Lambda)z`. The runtime all-cone
envelope is therefore

```text
v_coord,max = |shift| + (N/Lambda) (6/5).
```

This conversion matters. HYP2 alone does not bound coordinate speeds because
the lapse, shift, and radial metric can evolve even while the local principal
system remains healthy.

## Causal ledger and fail-closed step rule

At each proposed Runge--Kutta step, before accepting its final state, the solver must supply the
maximum coordinate-speed envelope over the complete radial slice at the
previous endpoint, every internal stage used by that step, and the proposed
final Runge--Kutta combination. The final candidate endpoint is evaluated
explicitly because the last internal RK stage need not equal the accepted
state. The ledger debits

```text
Delta d = Delta t * max(v_previous, v_stage_1, ..., v_stage_s, v_endpoint).
```

The proposed state is admissible only when

```text
r_outer - r_measure - d_accumulated - stencil_reach > 16.
```

The inequality is strict. A failed proposal raises the typed stop
`boundary_causal_buffer`; the immutable last accepted ledger is not advanced.
HLT1 owns the first-failed-premise transaction that composes this check with
the other health and scale stops.

The stencil term is the largest spatial offset used by the declared SBP
operator: three grid intervals for the primary fourth-order/interior,
second-order/boundary operator and one for the comparator. This is a
conservative extraction pad around a continuum causal calculation. It is not
a claim that a method-of-lines matrix exponential has an exactly compact
semidiscrete graph domain of dependence.

## Outer-boundary move control

The robustness radii are `96`, `128`, and `160`. Rather than holding the point
count fixed and silently changing the measurement resolution, BND2 preserves
each nominal interior spacing. At the coarsest nominal spacing `1/8`, for
example, the three grids contain `769`, `1025`, and `1281` points. The
measurement radius `24` is the same grid node in every variant.

On the flat reference envelope, the all-cone speed is exactly `6/5`. Even the
smallest boundary and widest primary stencil retain

```text
remaining buffer at t=32 = 1329/40,
margin above the required 16 = 689/40.
```

Every grid/boundary combination passes this exact reference calculation. This
is the pre-evolution outer-boundary move check required here: the variants are
grid-aligned and independently capable of protecting the frozen measurement
window. Agreement of actual evolved observables after moving the boundary is
still a NUM1/ROB1 requirement.

## Constraint scope

BND1 has six incoming main-system modes at each annular end. CON4 finds two
incoming gauge-subsidiary components at each end and zero-speed kinematic
reduction constraints, while explicitly leaving a nonlinear
constraint-preserving main-variable boundary map open.

BND2 does not erase that open problem. Instead, the causal ledger guarantees
that any incoming physical, gauge, reduction, reflection, or truncation error
originating at the outer boundary cannot enter a retained measurement
interval. Thus incoming constraint characteristics are controlled **inside
the certified measurement domain**, not on the entire numerical domain.

## Reproduction

```bash
python3 scripts/reproduce_fgc_bnd2_cp1.py \
  --output results/fgc-1-bnd2-cp1.json
```

The result is canonical unique-key JSON and binds the same PROTO3 semantic
run-envelope hash as RUN1.

## Nonclaims

BND2 is not a quasilinear existence theorem, nonlinear boundary-stability
proof, asymptotic boundary treatment, collapse evolution, affine-null result,
retained-EFT authorization, or physical-transition claim. In particular, a
run that exhausts the causal budget stops; it is not repaired by weakening the
buffer or ignoring a fast auxiliary cone.
