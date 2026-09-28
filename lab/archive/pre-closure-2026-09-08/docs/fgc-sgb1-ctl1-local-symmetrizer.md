# FGC-1-SGB1-CTL1 point-local / nonflat cone certificate (prospective)

This is a branch-owned implementation map, not a frozen certificate and not
`SGBL_branch_owned_and_healthy`. It does not change the ESF-perturbation
continuum owner in `sgb1_ctl1_cone`. Final FRZ1/PREF1 is a later owner.

## Claim

On a caller-declared orthonormal `SGBLPrincipalBackgroundBox` for the linear
branch `f(phi)=alpha_gb phi`, `beta=eta=0`, that is an exact singleton or a
declared small interval, enclose the complete ACT1 coefficient tensors and
decide a *point-local* kinetic / cone / symmetrizer statement.

Floating binary64 spectra seed eigenvalue candidates only. Real clusters,
eigenvectors, and Riesz projectors are certified by exact rational kernels
together with `V Λ V^{-1} = M`, or by rational interval Krawczyk/Neumann on
a reduced pencil. A positive symmetrizer is `H = (V^{-1})^T V^{-1}` with
Sylvester minors, or a coefficientwise construction. Complete 24-mode
multiplicity, basis, and coercivity are retained only when every strict
margin is positive.

All-direction is claimed only from a coefficientwise verified construction.
Otherwise that slot is typed incomplete. Finite angular samples, binary64
eig/SVD, fitted gaps, HYP2 cluster certificates, and FGC-QR health bits are
never the proof.

### Local Interpretation

From inside one orthonormal frame, at one declared background box, the
radial companion is a 24-by-24 matrix. When it is an exact rational
singleton, every candidate speed that survives an exact kernel check is a
true eigenvalue. Six four-dimensional kernels at

```text
{ -1, -1/2, -1/3, 1/3, 1/2, 1 }
```

recover the Einstein--two-scalar reference. Concatenating them yields a
rank-24 frame `V` inverted over `Q`. Exact Riesz projectors are
`P_J = V I_J V^{-1}`: idempotent, trace four, commuting with `M`, and pure
eigen-projectors. The symmetrizer `H = (V^{-1})^T V^{-1}` is positive by
Sylvester, and `H M` is exactly symmetric. Physical action-energy
restrictions at `±1` keep the leading minors `4, 16, 32, 64` and a
Gershgorin lower bound `2`.

### Non-Local Interpretation

The same collapse/expansion pair is not re-proved by transporting an ESF
resolvent across a deformation `||M(n)-M_ESF(n)||`. That argument remains
the existing cone owner. This owner either certifies the local companion
directly or, for all spatial covectors, exhibits a coefficientwise
polynomial identity or symmetrizer. A parent/child change of chart is not
inferred from a local pass.

### Simplest Relational Analogue

A calculator that returns an endless `9` or `0` has hit a chart wall. Here
a floating eigenvalue is only a seed. The certified object is an exact
kernel, an exact projector, or a typed obstruction that the local chart
could not close.

### Candidate Mathematics

Exact rational Gaussian kernels; Neumann inverse of the 12-by-12 kinetic
block `A`; parametric Krawczyk on a simple eigenpair of dimension at most
16; exact Riesz projectors `V I_J V^{-1}`; Sylvester minors; quadratic
polynomial identities on the unit sphere for the two physical scalar modes.

### Inherited Invariant

The linear-branch derivative contract `F=Mpl^2`, `F'=0`, `f'=alpha_gb`,
`Hess(f)=alpha_gb Hess(phi)`, the complete 10+20 generator identities, and
the locked auxiliary cones `tilde=4`, `hat=9`.

### Observable or Logical Consequence

The exact ESF/Minkowski reference proves the point-local 24-mode cone and
keeps aggregate health false. The locked nonflat source fixture A remains
one typed obstruction: auxiliary kernels at `±1/2` and `±1/3` are exact of
dimension four, but the split physical clusters are not a complete
rational 24-mode basis. A later family-box entrypoint accepts a
`SGBLPrincipalBackgroundBox` and runs the same local certificate; uniform
covering of a family cell product box is a named later owner.

### Open Variable

A coefficientwise 24-mode all-direction symmetrizer `H(n)` for non-isotropic
backgrounds, and a Krawczyk/Riesz close of the six nonrational nonflat
modes, remain open. They are not filled by floating SVD or sampled
directions.

## SGB-L coefficient box

The box is the existing orthonormal principal-background product box. This
module does not edit that type. Widths are declared by the caller. A
non-singleton companion is `interval_inconclusive` /
`interval_companion_not_singleton`, not a kinetic failure and not a fitted
floor.

## Exact local reference

Unnormalized radial blocks `A`, `B_r`, `C_rr` at `n=(1,0,0)` give the
companion

```text
M = [[ 0, I ], [ -A^{-1} C_rr, -A^{-1} B_r ]].
```

Floating `eig(M)` may suggest the six ESF speeds. The proof is the exact
kernel of `M - λ I` at each surviving rational candidate, the identity
`M v = λ v`, rank 24, `V Λ V^{-1} = M`, the six exact Riesz projectors,
and the Sylvester symmetrizer. No ulp allowance and no contour sample
count enter those identities.

## Continuum all-direction

All-direction is not inherited from the ESF-perturbation deformation bound.
It is claimed only when a coefficientwise construction is verified. The
physical scalar modes satisfy the quadratic identities

```text
(M(n) - s I) v_{phi,chi}(s) = (|n|^2 - 1) w
```

for `s = ±1` on a determining quadratic grid. Those identities do not
supply the remaining twenty-two modes. The all-direction slot therefore
stays `typed_incomplete` in this slice.

## Family-box API

```text
entrypoint     = sgbl_local_symmetrizer_from_family_box
accepted slot  = SGBLPrincipalBackgroundBox
later owner    = uniform_family_cell_product_box_local_symmetrizer
covering       = false
health         = false
```

A family principal box may be fed here as a declared local box. Uniform
covering uniqueness on a cell product box is named, not replaced by a
hidden zero.

## Typed stops versus wrapping

| Outcome | Meaning |
|---|---|
| `point_local_cone_proved` | Exact 24-mode real basis, positive symmetrizer, positive kinetic and energy margins, algebraic identities. |
| `typed_incomplete` | Exact kernels did not complete 24 modes. One obstruction is retained. Not a kinetic failure. |
| `interval_inconclusive` | Non-singleton companion or `ρ_∞ ≥ 1`. Not a kinetic failure. |
| `lost_kinetic` | Center `A` exactly singular, or `V Λ V^{-1} = M` failed. |
| `broken_gauge` | Algebraic Hessian or pure-gauge identity failed on the complete 10+20 basis. |
| `defective_cluster` | Geometric multiplicity strictly below algebraic multiplicity. |
| `complex_spectrum` | Isolated real 2-block with negative discriminant. |
| `cluster_coalescence` | Declared centres are not strictly separated, or ESF six-speed multiplicities merged. |
| `eigenframe_incomplete` | Rank below 24 after an injected omission or a singular frame. |
| `coercivity_failed` | Sylvester or Gershgorin energy/symmetrizer test failed. |
| `resource_limit` | Declared bit, evaluation, candidate, or copied-`512` cap. |
| `cone_chart_error` / `sign_chart_error` | Wrong frame, `F ≤ 0`, or unlocked auxiliary cones. |

The locked nonflat nonzero-shift source fixture remains
`typed_incomplete` with obstruction `complete_24_mode_basis_not_certified`.
Structural identities still hold. That is not a matched-family box and
does not open branch health.

No post-outcome margin fitting. Resource policies are prospective.

## Attacks

The test surface injects defective Jordan blocks, complex 2-blocks,
coalesced centres, eigenframe omissions, Hessian/gauge mutations, energy
sign flips, resource caps, and digest tampering. Each fails closed with
the matching typed reason. None of those attacks may set
`SGBL_branch_owned_and_healthy`.

## What this slice is not

It is not the ESF-perturbation continuum cone owner, not a trajectory,
holdout, production width, retained-EFT statement, matched-family
curvature box, or branch-health aggregate. It does not freeze FRZ1/PREF1.
It does not import HYP2 pass bits or an FGC-QR builder. A later owner
must feed authenticated family boxes and close all-direction before
`SGBL_branch_owned_and_healthy` can be considered.
