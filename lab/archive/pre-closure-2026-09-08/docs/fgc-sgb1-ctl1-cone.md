# FGC-1-SGB1-CTL1 cone / symmetrizer slice (prospective)

This is a branch-owned implementation map, not a frozen certificate and not
`SGBL_branch_owned_and_healthy`. Final FRZ1/PREF1 is a later owner.

## Claim

On a caller-declared orthonormal principal-background box for the linear
branch `f(phi)=alpha_gb phi`, `beta=eta=0`, enclose the complete ACT1
coefficient tensors `(A, B_i, C_ij)` and decide a continuum all-unit-spatial-
covector kinetic / cone / symmetrizer statement from n-independent operator
norms plus an **exact unnormalized radial Einstein--two-scalar eigenframe**.

A real complete 24-characteristic basis and a positive Kovacs--Reall
symmetrizer are reported only when every strict margin is positive. Finite
angular samples and binary64 SVD/eig are regressions. They are never the
proof.

## SGB-L coefficient box

The box is declared in the physical orthonormal Minkowski frame. The linear
branch is part of the contract, not an FGC-QR or HYP2 pass bit:

```text
F = Mpl^2          (interval, strictly positive)
F' = 0             (identically)
f' = alpha_gb      (interval; zero only as a named control)
Hess(f) = alpha_gb * Hess(phi)
```

Public enclosures of `A`, `B_i`, and `C_ij` remain the Frobenius-normalized
quadratic coefficients of

```text
P(xi) = A xi_0^2 + B_i xi_0 n_i + C_ij n_i n_j
```

so they can be compared with the branch-owned NumPy principal adapter. The
**proof** uses the unnormalized symmetric-tensor basis, where every ESF
entry is an exact rational singleton.

The adapter `sgbl_principal_box_from_state` and the locked nonzero-shift
fixture `sgbl_nonflat_source_principal_box` ground the same derivative
contract on the branch-owned principal image. Float64 orthonormal curvature
is stored as exact dyadics. That is not a matched-family curvature box.

## Exact ESF reference

In the unnormalized basis, extract singleton rational radial blocks `A`,
`B_r`, `C_rr` at `n=(1,0,0)` and form the exact 24-by-24 companion

```text
M = [[ 0, I ], [ -A^{-1} C_rr, -A^{-1} B_r ]].
```

For each known real speed

```text
{ -1, -1/2, -1/3, 1/3, 1/2, 1 }
```

compute an exact rational kernel of `M - λ I`. Each kernel has dimension
four. Concatenate the six kernels to a 24-by-24 frame `V`. Exact rank is
24 and `V` is inverted over `Q`. Each column satisfies `M v = λ v`
exactly.

A rational 2-norm upper bound is `||X||_2 ≤ sqrt(||X||_1 ||X||_∞)`. Then

```text
κ_2(V) ≤ ||V||_2^upper ||V^{-1}||_2^upper.
```

Declared contour radii (`1/8` physical, `1/16` auxiliary) are required to
be strictly less than half the nearest spectral gap. On each contour the
distance to the spectrum is the radius, so Bauer--Fike gives

```text
σ_min(zI - M) ≥ radius / κ_2(V).
```

No floating SVD and no ulp allowance enter this lower bound. Arc samples
may exist as regressions; they do not certify.

## Physical energy

On the unnormalized action blocks `(A_act, B_act)` the first-order form is

```text
H_s = sign(s) [[B_act, A_act], [A_act, 0]]
```

with `sign(+1) = -1` and `sign(-1) = +1`. Restricting to the four exact
physical eigenvectors at each of `±1` yields a 4-by-4 rational matrix.
Sylvester's criterion (strictly positive leading minors `4, 16, 32, 64`)
and a positive Gershgorin lower bound (`2`) prove coercivity. The operator
2-norm of `H_s` is bounded by `4`.

## Continuum all-direction argument

The companion `M(n)` is quadratic in the unit spatial covector `n`. For
`|n|=1`,

```text
||M(n)-M_ESF(n)||_2
  ≤ sum_i ||Δ(A^{-1} B_i)||_2 + sum_ij ||Δ(A^{-1} C_ij)||_2
```

in the **unnormalized** basis. The right-hand side does not depend on `n`.
Homogeneity extends the bound to every nonzero spatial covector. Pass when
this deformation is strictly below every certified resolvent lower bound,
the kinetic Neumann `σ_min(A)` is strictly positive, and the transported
physical energy remains strictly positive.

## Structural identities versus wrapping

Ungauged Hessian symmetry and the homogeneous cubic pure-gauge identity
are algebraic statements of the ACT1 formula. Linear-GB principal
dependence is linear in the independent background `Hess(phi)` slots and
the algebraic Riemann slots, so the identities are executed on a
**complete exact basis**, not on a symmetry representative:

- all 10 symmetric `4 x 4` Hessian generators;
- 20 linearly independent algebraic-curvature generators, constructed as
  Kulkarni--Nomizu products of those symmetric tensors with exact
  rank-20 selection.

Each curvature generator is independently required to obey pair
antisymmetry, pair exchange, and the first Bianchi identity. The 55 KN
candidates span a rank-20 space; a selected 20-generator subset has exact
rank 20. `representative_generator_sufficient` is therefore false: Hess(00)
and R(0101) lie in the span but do not replace the basis. Mutation or
omission of a generator fails closed as `broken_gauge`.

A valid nonflat box is **not** rejected because independently intervalized
identities wrapped. `broken_gauge` is reserved for an actual algebraic
failure, an incomplete basis, or an injected mutation.

## Typed stops versus wrapping

| Outcome | Meaning |
|---|---|
| `continuum_cone_proved` | All strict exact margins positive and structural identities hold. |
| `interval_inconclusive` | Wrapping, `ρ_∞ ≥ 1`, or a nonpositive deformation/resolvent/energy margin. Not a kinetic failure. |
| `lost_kinetic` | Center `A` is exactly singular, or the ESF eigenframe/energy form fails exactly. |
| `broken_gauge` | Algebraic Hessian or pure-gauge identity failed on the complete 10+20 generator basis, a rank/Bianchi failure, or an injected mutation/omission. |
| `cone_chart_error` / `sign_chart_error` | Wrong frame, `F ≤ 0`, or unlocked auxiliary cones. |
| `resource_limit` | Declared bit, evaluation, width, contour, or copied-`512` cap. |

The locked nonflat nonzero-shift source fixture contains the independent
NumPy principal adapter. Its cone is `interval_inconclusive` with exact
obstruction `resolvent_margin_not_strictly_positive` (`δ ≈ 31.8` against
an ESF resolvent floor `≈ 9.9e-5`). Structural identities still hold.
That is not a matched-family box and does not open branch health.

No post-outcome margin fitting. Resource policies are prospective.

## What this slice is not

It is not a trajectory, holdout, production width, retained-EFT statement,
matched-family curvature/Hessian box, or branch-health aggregate. It does
not freeze FRZ1/PREF1. A later owner must feed authenticated family boxes
before `SGBL_branch_owned_and_healthy` can be considered.
