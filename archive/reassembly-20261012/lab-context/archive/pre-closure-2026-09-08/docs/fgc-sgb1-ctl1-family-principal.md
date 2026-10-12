# FGC-1-SGB1-CTL1 matched-family principal-box feeder (prospective)

This owner is an implementation map for feeding authenticated matched-family
geometry into the existing cone instrument. It is not
`SGBL_branch_owned_and_healthy`, not FRZ1/PREF1, and not a holdout pass.

## Claim

Given an authenticated completed initial-health ODE cell or compactness
record together with the exact family specification, enclose complete
interval lower/full two-jets and emit an orthonormal
`SGBLPrincipalBackgroundBox` for `sgbl_enclose_principal_cone`.

Binary64 values are not treated as enclosures. A required derivative or
source-box theorem that cannot be implemented rigorously is named in an
immutable missing-owner record. It is not filled by zero.

## Slot ownership

| Slot | Owner | Formula |
|---|---|---|
| `r` | authenticated ODE cell or buffer/exterior radius | interval or singleton |
| `lambda`, `k` | Picard graph `lambda_box`, `k_box` | Minkowski buffer `(1,0)` |
| `phi, phi_r, phi_rr` | `sgbl_compact_bump_enclosure` | `A_phi (B, B_r, B_rr)` |
| `phi_rrr` | `sgbl_compact_bump_third_derivative_enclosure` | `A_phi B_rrr` |
| `phi_pi, phi_pi_r` | declared compact family | `0` |
| `chi` | complete interval family fields | `A_chi B / r` |
| `chi_r` | complete interval family fields | `A_chi (B_r/r - B/r^2)` |
| `chi_rr` | complete interval family fields | `A_chi (B_rr/r - 2 B_r/r^2 + 2 B/r^3)` |
| `chi_pi` | complete interval family fields | `A_chi B_r / r = chi_t` |
| `chi_pi_r = chi_tr` | complete interval family fields | `A_chi (B_rr/r - B_r/r^2)` |
| `lambda_r` | affine Hamiltonian | `-H0 / H_L` |
| `k_r` | affine momentum | `-M0 / M_k` |
| `lambda_rr, k_rr` | interval differentiation of that RHS | `IntervalFirstTangent` along `(r, lambda(r), k(r), fields(r))`; needs `phi_rrr` |
| `alpha` | unit lapse | `1` |
| `shift, shift_t` | zero shift plus `C^t=C^r=0` | `shift_t=(2 L^3-2 L+L_r r)/(4 L^3 r)` |
| `lambda_t` | polar-areal `K^r_r=-2k` | `2 L k` |
| `R=r, R_t=-r k, R_r=1, R_tr=-k-r k_r, R_rr=0` | polar-areal identity | |
| ADM `dtt` | `sgbl_parametric_source_admission` | Krawczyk image of `R(a;z)=R0+J a` |
| `Riemann_coord` | `direct_4d_curvature` | existing RED1 owner |
| `Hess(phi)_coord` | spherical covariant Hessian | `partial_a partial_b phi - Gamma^c_ab partial_c phi` |
| orthonormal frame | checked positive frame | `n=-lapse g^{-1} dt`, `e_r=(0,1/sqrt(h_rr))`, `e_theta=e_phi=1/R` with outward rational square roots |
| `Riemann_orth`, `Hess(phi)_orth` | frame pullback | `R_abcd=e^mu_a e^nu_b e^rho_c e^sigma_d R_munurhosigma` |
| box | `SGBLPrincipalBackgroundBox` | `F=Mpl^2`, `F'=0`, `f'=alpha_gb`, `Hess(f)=alpha_gb Hess(phi)` |

The third-profile bound `|B'''(x)|<39312` follows from `e<3` and
`e^{-u} u^n <= n!`. It is not a sampled float maximum.

## Coverage

- **Exact Minkowski buffer.** `lambda=1`, `k=0`, scalars vanish, `dtt=0`.
  The orthonormal box is curvature-free and the existing cone instrument
  proves the continuum cone. Aggregate health stays false.
- **Completed support cells.** Lower jets are the cell product box. The
  principal box is built at a declared exact witness inside that product
  box, with `dtt` from a parametric uniqueness certificate at
  `parameter_half_width=0`. Uniform covering uniqueness on the whole
  cell product box is a named missing owner, not a hidden zero.
- **Analytic exterior.** Only when the compactness graph covers the
  compact support. The support-boundary state is the authenticated input.
  An incomplete graph does not fabricate `C=2M/r`.

## Typed incomplete

Nominal `A_chi=3` remains `picard_strict_self_map_failed` on the
declared tree. The feeder returns one
`family_geometry_incomplete` record: no exterior box, no health pass,
and the missing theorem
`complete_validated_constraint_ODE_graph_covering_the_compact_support`.

## Nonflat fixture A

The locked nonzero-shift two-scalar source fixture is converted by exact
interval curvature and a checked rational-square-root frame, not by
`float.as_integer_ratio` as the definition of the tensor. The box
contains the independent NumPy adapter image. Its cone remains
`interval_inconclusive` / `resolvent_margin_not_strictly_positive`.

## What this slice is not

It is not a trajectory, holdout, production width, FGC-QR/HYP2 pass,
or `SGBL_branch_owned_and_healthy`. FRZ1, PREF1, execution and holdout
stay false. Existing cone and nonflat outcomes are preserved because
this module does not edit those owners.
