# FGC-1-SGB1-CTL1 constraint continuity and conditional propagation

This is a branch-owned algebraic continuity map plus one *conditional*
boundary-free constraint-propagation theorem. It is not a CON4/HYP1
certificate, not an existence or IBVP result, and not
`SGBL_branch_owned_and_healthy`. Final FRZ1/PREF1 is a later owner.

## Claim

On one exact SGB-L source point with `F = Mpl^2 > 0` and the locked hat
factor `q = 9`, the following are re-evaluated independently:

1. The MHG extension does not modify the scalar equations `E_phi` and `E_chi`.
2. The complete six-row residual remains affine in the ADM accelerations,
   `R(a) = R0 + J a`.
3. The unredefined physical projections `H` and `M` are independent of those
   accelerations.
4. At `F = Mpl^2`, `F' = 0`, the hat-projector contraction on the scalar
   shell is the hat-wave operator `(F/2) hatq`.
5. On the MHG metric shell, the map from unredefined `(H, M)` to
   `(nabla_t C^t, nabla_t C^r)` is the inverse of the exact `2 x 2` ADM
   matrix obtained by projecting `E = -X` on two unit normal-derivative
   columns.
6. The locked-hat radial null polynomial has two real roots of opposite
   sign, so the homogeneous gauge subsidiary system has both radial
   hat-cone branches.
7. The FO1 kinematic identity `partial_t C + partial_r D - K = 0` holds
   field-by-field, including on an off-shell probe where `D`, `C`, and `K`
   themselves need not vanish.

If those identities hold, the immutable
`SGBLConstraintPropagationContract` may set

```text
conditional_boundary_free_constraint_propagation_proved = true
```

That boolean is the implication below. It is not existence of a solution
and it does not open aggregate health.

## Conditional statement

If a sufficiently smooth compatible SGB-L solution exists, `C = H = M =`
radial-reduction constraints `= 0` on a spacelike initial slice, and the
complete MHG metric equations and both unmodified scalar equations hold,
then the exact invertible normal map, the homogeneous hat-wave subsidiary
system, linear normally-hyperbolic uniqueness, and the FO1 kinematic
identity preserve those constraints throughout the boundary-free domain of
dependence of the slice.

### Logical chain

```text
C = 0 on the initial slice
        -> spatial nabla C = 0
H = M = 0 and the MHG metric equations
        -> normal nabla C = 0 by the invertible 2x2 map
unmodified scalars plus MHG
        -> homogeneous hat-wave subsidiary system
named linear normally-hyperbolic uniqueness
        -> C = 0 in the boundary-free domain of dependence
C = 0
        -> X = 0
        -> MHG recovers the unredefined metric equation
        -> H = M = 0 remains true
FO1 identity partial_t C + partial_r D - K = 0
        -> the six radial reduction constraints remain zero
```

Uniqueness is a *named premise*. This module does not import a CON3
certificate, does not derive an energy estimate, and does not construct
the solution to which the uniqueness statement would apply.

## Exact physical-to-normal map

Write the spherical base as

```text
(h_tt, h_tr; h_tr, h_rr),   h_rr > 0,   ell^2 = -h_tt + h_tr^2 / h_rr > 0.
```

On the MHG shell `E = -X` with `C = 0` as a function on the slice (so
spatial covariant derivatives of `C` vanish), the unredefined projections
satisfy

```text
(H, M)^T = (F q / 2) [[ell^2, 0], [-h_tr, -h_rr]] (nabla_t C^t, nabla_t C^r)^T
```

with `F = Mpl^2` and locked `q = 9`. The determinant is

```text
det = - F^2 q^2 ell^2 h_rr / 4 < 0.
```

The inverse is therefore the map from unredefined `(H, M)` to the normal
time derivatives of `(C^t, C^r)`. On fixture A this evaluates to

```text
forward matrix = [[72, 0], [-27/2, -81/2]]
det            = -2916
inverse        = [[1/72, 0], [-1/216, -2/81]]
```

The two columns of the forward matrix are obtained by inserting exact
unit `nabla_t C` data into the MHG extension residual and projecting
`E = -X`. The analytic ADM matrix is constructed independently. Entrywise
equality is required. No CON4 aggregate certificate is imported.

## Hat-cone characteristics

The homogeneous subsidiary principal part is `(F/2) hat_g^{ab} nabla_a nabla_b C^mu`.
For the radial covector `xi_a = (-c, 1)` the locked-hat null polynomial on
fixture A is

```text
(7/36) - (3/2) c - (9/4) c^2 = 0
```

with exact roots `-7/9` and `1/9`. The signs are `(-1, +1)`: one incoming
and one outgoing radial branch. The product of the roots is negative
because the constant/leading coefficient ratio is negative.

## FO1 kinematic reduction

For each of the six FO1 fields,

```text
D = partial_t u - p,    C = q - partial_r u,    K = partial_t q - partial_r p.
```

Commutation of exact `t, r` partials gives `partial_t C + partial_r D - K = 0`
identically, including off the kinematic shell. The subsidiary principal
speed is zero, so this constraint supplies no incoming boundary field.
The source-point two-jet used here is the commuting embedding of the
SGB-L spherical jet; it is not a claim that the source residual already
solves the reduction constraints.

## Explicit assumptions

`True` in the stored assumption map means the implication *uses* that
premise. It does not mean the premise was proved in this module.

| Assumption | Role |
|---|---|
| sufficiently smooth compatible solution assumed | existence is not supplied |
| `C = H = M =` radial-reduction `= 0` on the initial slice assumed | no compatible hypersurface is constructed |
| complete MHG metric equations assumed | the solution, if it exists, must satisfy the extended metric rows |
| both unmodified scalar equations hold along the solution assumed | needed to homogenize the hat-wave operator |
| linear normally-hyperbolic uniqueness for the hat wave assumed | named PDE premise, not an imported CON3 pass bit |
| FO1 kinematic equations `D` and `K` hold as differentiable identities assumed | reduction constraints then remain zero |
| boundary-free domain of dependence | the theorem stops where a boundary enters |

## Typed stops and adversarial controls

| Outcome | Meaning |
|---|---|
| `constraint_continuity_compatible` | Pointwise identities (1)--(4) hold. |
| `conditional_boundary_free_constraint_propagation_proved` | Identities (1)--(7) were re-evaluated and the implication is recorded. |
| `broken_map` | Direct MHG columns and the locked-hat analytic matrix disagree, or the determinant is not `-F^2 q^2 ell^2 h_rr / 4`. |
| `broken_sign` | Analytic sign mutation, or a nonnegative map determinant. |
| `broken_hat` | Hat factor is not the locked `q = 9`, or `q <= 1`. |
| `broken_root` | Hat-cone roots are not one incoming and one outgoing radial branch. |
| `broken_gauge` | Scalar rows were modified, or the hat-wave contraction failed. |
| `broken_assumption` | FO1 identity failed field-by-field. |
| `action_identity_error` | Not the linear SGB-L action with `F > 0`. |
| `source_identity_failed` | Affine reconstruction missed the residual, or `H, M` saw accelerations. |

Mutation of derived health, existence, IBVP, FRZ1/PREF1/holdout, the frozen
conditional statement, the named uniqueness premise, or the logical chain
is refused. Replacing
`conditional_boundary_free_constraint_propagation_proved` independently of
the re-evaluated identities is refused.

## What remains false

```text
smooth_solution_exists = false
compatible_initial_data_supplied = false
regular_center_supplied = false
boundary_condition_supplied = false
IBVP_proved = false
unconditional_domain_propagation = false
SGBL_branch_owned_and_healthy = false
FRZ1 = false
PREF1 = false
holdout = false
constraint_propagation_on_a_domain_proven = false
constraint_preserving_boundary_map_supplied = false
```

The remaining gaps are therefore a regular centre, a compatible
nonzero-width initial slice, an incoming constraint-preserving boundary
map, a quasilinear existence theorem, and any unconditional (as opposed
to conditional and boundary-free) domain statement.

## What this slice is not

It is not a neighboring HYP1/CON1/CON2/CON3/CON4 result, not an interval
cone, not a source runtime, and not permission to open holdout or
physical execution. Branch-neutral tensor and map primitives may be
called, but every returned identity is re-checked field-by-field on SGB-L
data.
