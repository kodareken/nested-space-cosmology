# FGC-1-SGB1-CTL1 — whole-cell source-acceleration admission (prospective)

This owner is an implementation map for a *nonuniform family product-box*
uniqueness certificate of the affine SGB-L source. It is not a frozen
production admission, not `FRZ1`/`PREF1`, and not
`SGBL_branch_owned_and_healthy`. It does not edit the existing uniform
parametric admission or the matched-family principal feeder.

## Claim

At fixed lower jets `z` the complete six-row MHG residual of the linear
branch is the affine map already owned by `sgb1_ctl1_source.py`:

```text
R(a; z) = R0(z) + J(z) a,
a = (alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt).
```

On an *explicit nonuniform* lower-jet product box `Z` whose thirty ADM
slots each carry an exact interval and an owner/formula, together with a
*predeclared* acceleration box `A` whose interior contains the exact
algebraic root at the declared centre, a parametric interval-Newton/Krawczyk
certificate proves a unique acceleration root for every `z` in `Z` **only
when that inclusion closes**.

The existing uniform owner applies one `parameter_half_width` to every
non-acceleration slot. The family feeder therefore certifies only a
declared witness and names

```text
parametric_Krawczyk_uniqueness_on_the_non_uniform_family_cell_product_box
```

as a missing owner. This module is that owner, implemented as a new
fail-closed instrument. It never replaces the reconstructed family radii
by one uniform half-width and never promotes a singleton witness to a
covering theorem.

## Slot provenance

Every lower-jet slot is an input with an owner. There is no
`parameter_half_width` field.

| Slot | Family owner | Formula |
|---|---|---|
| `alpha.*` | unit-lapse polar-areal gauge | `alpha=1`, other components `0` |
| `shift.value, shift.dr, shift.drr` | zero-shift plus `C^t=C^r=0` | `shift=0` |
| `shift.dt, shift.dtr` | zero-shift plus `C^t=C^r=0` | `shift_t=(2 L^3-2 L+L_r r)/(4 L^3 r)` |
| `lambda.value` | Picard `lambda_box` | authenticated cell interval |
| `lambda.dt` | polar-areal `K^r_r=-2k` | `2 L k` |
| `lambda.dr` | affine Hamiltonian | `-H0 / H_L` |
| `lambda.dtr` | polar-areal | `2(L_r k + L k_r)` |
| `lambda.drr` | interval differentiation of that RHS | needs `phi_rrr` |
| `areal_radius.value` | authenticated cell radius | polar-areal `R=r` |
| `areal_radius.dt, dr, dtr, drr` | polar-areal identity | `R_t=-r k`, `R_r=1`, `R_tr=-k-r k_r`, `R_rr=0` |
| `phi.value, dr, drr` | `sgbl_compact_bump_enclosure` | `A_phi (B, B_r, B_rr)` |
| `phi.dt, dtr` | declared compact family | `0` |
| `chi.*` | complete interval family fields | `A_chi` bump / `r` two-jet, including `chi_tr` |
| `coordinate_radius` | cell radius or exact source-point radius | MHG reference connection at the declared exact centre |
| ADM `dtt` | predeclared acceleration box | Krawczyk image of `R(a;z)=R0+J a`; never a hidden zero |

The MHG reference connection remains the frozen annulus chart at the
declared exact radius centre. Polar-areal `R` is the independent
`areal_radius.value` product slot. Whole-cell uniqueness is uniqueness
on that declared ADM jet product box.

A family cell box is reconstructed from
`sgbl_family_adm_lower_jets`. Tampering those intervals, replacing them
by one shared half-width, or collapsing them to the witness singleton is
a typed refusal. Omitted slots are refused.

## Certificate conditions (not fitted thresholds)

The declared product box and acceleration box are inputs. A failed
inclusion does not enlarge, shrink, uniformize, or otherwise fit them
after the outcome.

The certificate passes only if all of the following hold:

1. The exact centre Jacobian is nonsingular. An exact `det J(z0)=0`
   remains the existing `singular_source_jacobian` stop.
2. A Neumann inverse of the interval Jacobian `J(Z)` proves `rho_∞<1`.
   Wide wrapping with `rho>=1` is typed `neumann_rho_not_below_one`.
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
residual enclosure.

Non-singleton residual and Jacobian endpoints may be outwardly rounded
onto the declared `2^{-64}` dyadic grid. That is a valid (wider)
enclosure. Singleton exact boxes are not rounded, so fixture A
reproduces the frozen source Jacobian.

## Local versus aggregate flags

`unique_acceleration_root_for_every_declared_parameter_point` may be
true as a local theorem on the declared product box.

The following remain false:

- `source_solve_admission_qualified`
- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- `promoted_singleton_witness`
- `uniformized_nonuniform_radii`
- execution and holdout authorization

The existing feeder is untouched: it may still pass a witness uniqueness
certificate and retain the named missing covering owner. Witness
uniqueness is not copied into this owner.

## Controls

- **Flat nonuniform box.** Minkowski vacuum with distinct declared
  `lambda.value` and `areal_radius.value` radii, other slots exact.
  Inclusion is required to close.
- **Exact fixture A singleton.** All thirty slots are exact rationals.
  The centre determinant is `6519803496633/6553600000`. Inclusion
  closes, and the Jacobian reduces to the frozen source matrix.
- **`A_chi=1/8` completed family cells.** Reconstruct the nonuniform
  cell product box. Either whole-cell inclusion closes, or the exact
  wrapping/resource/chart obstruction is preserved. The box is not
  refit and the feeder witness is not promoted.
- **Nominal incomplete graph.** `A_chi=3` remains
  `picard_strict_self_map_failed`. The graph entry returns one typed
  `family_geometry_incomplete` record: no cell certificates, no health
  pass, missing theorem
  `complete_validated_constraint_ODE_graph_covering_the_compact_support`.
- **Attacks.** Omitted slots, uniform-half-width aliases, correlated
  replacement of distinct family radii by one shared width, singleton
  witness promotion, zero-containing denominators, residual-evaluation
  resource caps, and forged uniqueness without a Krawczyk payload are
  typed refusals.

## Typed stops

| Reason | Meaning |
|---|---|
| `singular_source_jacobian` | Exact centre `det J=0`. Not wrapping. |
| `acceleration_box_not_interior` | Declared acceleration box has empty interior. |
| `neumann_rho_not_below_one` | Wrapping; contraction not proved. |
| `krawczyk_contraction_not_below_one` | Same contraction failure from the Krawczyk operator. |
| `krawczyk_image_not_strictly_inside` | Contraction may hold; inclusion failed. The box is not refit. |
| `residual_enclosure_misses_origin` | Image residual box does not contain 0. |
| `resource_limit` | Residual-evaluation or rational-bit cap. |
| `interval_chart_domain` / `zero_in_interval_reciprocal` | Declared box left the ADM chart. |
| `omitted_slot` | A required lower-jet slot is missing. |
| `uniformization_alias` | A single uniform half-width was offered as the product box. |
| `correlated_uniformization` | Distinct family radii were replaced by one shared width. |
| `promoted_singleton_witness` | A non-singleton family cell was collapsed to its witness. |
| `tampered_slot_map` | Family slots no longer match the reconstructed owners. |
| `family_geometry_incomplete` | Compactness graph does not cover the support. |
| `unauthenticated_cell` | Picard self-map / residual / compactness invariant failed. |
| `float_enclosure_refused` | Binary64 was offered as an enclosure. |

## Non-claims

This owner does not enclose the cone, qualify the evolving centre, run a
trajectory, copy FGC-QR/GR-0 health evidence, freeze a physical
production width, or rewrite the existing feeder/admission modules.
Aggregate health stays false even after a local whole-cell uniqueness
proof.
