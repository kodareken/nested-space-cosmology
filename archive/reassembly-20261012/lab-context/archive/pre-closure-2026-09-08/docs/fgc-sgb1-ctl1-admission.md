# FGC-1-SGB1-CTL1 — parametric source admission (prospective)

This owner is an implementation map for a declared-box uniqueness
certificate of the affine SGB-L source. It is not a frozen production
admission, not `FRZ1`/`PREF1`, and not `SGBL_branch_owned_and_healthy`.

## Claim

At fixed lower jets `z` the complete six-row MHG residual of the linear
branch is the affine map already owned by `sgb1_ctl1_source.py`:

```text
R(a; z) = R0(z) + J(z) a,
a = (alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt).
```

On a *caller-declared* lower-jet box `Z` and a *caller-declared*
acceleration box `A` whose interior contains the exact algebraic root at
the box centre, a parametric interval-Newton/Krawczyk certificate proves
a unique acceleration root for every `z` in `Z` when it passes.

## Certificate conditions (not fitted thresholds)

The declared boxes are inputs. A failed inclusion does not enlarge,
shrink, or otherwise fit them after the outcome.

The certificate passes only if all of the following hold:

1. The exact centre Jacobian is nonsingular. An exact `det J(z0)=0` remains
   the existing `singular_source_jacobian` stop.
2. A Neumann inverse of the interval Jacobian `J(Z)` proves `rho_∞<1`.
   Wide wrapping with `rho>=1` is typed `interval_inconclusive`. That
   inequality is the Banach contraction hypothesis, not a numerical floor.
3. Every displacement interval of `A` about the declared acceleration
   centre contains zero strictly in its interior.
4. The parametric Krawczyk image
   `-C F(a0,Z) + (I-C J(Z)) (A-a0)`
   lies *strictly inside* the displacement box, with a strictly positive
   componentwise inclusion margin.
5. The affine residual enclosure of that image, `F(a0,Z)+J(Z)(image-a0)`,
   contains the origin in every row.

Condition 5 is a validated residual enclosure. It is not PROTO4
`newton_residual_limit`, not `1e-12`, and not “the residual is finite”.
A merely finite residual with failed contraction or failed inclusion is
inconclusive.

An interval-Newton image `-J(Z)^{-1} R0(Z)` is reported when a Neumann
inverse exists. Strict Newton inclusion is additional evidence. The
uniqueness bit requires the Krawczyk contraction/inclusion pair and the
residual enclosure, because those are the hypotheses of the existing
exact parametric Krawczyk theorem.

## Local versus aggregate flags

`unique_acceleration_root_for_every_declared_parameter_point` may be
true as a local theorem on the declared box.

The following remain false:

- `source_solve_admission_qualified`
- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- execution and holdout authorization

Final production conjunction is a later owner.

## Typed stops

| Reason | Meaning |
|---|---|
| `singular_source_jacobian` | Exact centre `det J=0`. Not wrapping. |
| `acceleration_box_not_interior` | Declared acceleration half-width is zero. |
| `neumann_rho_not_below_one` | Wrapping; contraction not proved. |
| `krawczyk_contraction_not_below_one` | Same contraction failure from the Krawczyk operator. |
| `krawczyk_image_not_strictly_inside` | Contraction may hold; inclusion failed. The box is not refit. |
| `residual_enclosure_misses_origin` | Image residual box does not contain 0. |
| `resource_limit` | Residual-evaluation or rational-bit cap. |
| `interval_chart_domain` / `zero_in_interval_reciprocal` | Declared box left the ADM chart. |

## Non-claims

This owner does not enclose the cone, qualify the evolving centre, run a
trajectory, copy FGC-QR/GR-0 health evidence, or freeze a physical
production width.
