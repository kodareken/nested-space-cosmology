# FGC-1-HYP2-MD1: all-covector classical early-kill gate

## Decision

`FGC-1-HYP2-MD1` supplies a quantitative sufficient weak-coupling health ball
for the complete four-dimensional principal system of the frozen ACT1
FGC-QR branch in modified harmonic gauge. It retains all ten metric
components, the regulator `phi`, and the independent collapsing scalar `chi`.
The corresponding first-order characteristic matrix has dimension 24.

The gate passes:

```text
quantitative_all_covector_weak_coupling_health_envelope_passed = true
```

This is an all-spatial-direction **classical principal-health** statement. It
is not retained-EFT authorization and it does not say that a future nonlinear
trajectory stays inside the health ball. HLT1 must evaluate these same
inequalities at every accepted stage and stop before the first exit.

## The theorem being specialized

Kovács and Reall establish strong hyperbolicity for weakly coupled Horndeski
theories in a modified harmonic formulation. Their result is a continuity
theorem: sufficiently small deformations of the Einstein--scalar reference
retain the required spectral and physical-subspace structure. It does not
publish one universal numerical epsilon for every action and normalization.
See [Well-posed formulation of scalar-tensor effective field theory](https://arxiv.org/abs/2003.08398).

HYP2 therefore does not quote an epsilon. It specializes the continuity
argument to the exact ACT1/VAR1 row normalization, the `1<tilde<hat` auxiliary
cones `4` and `9`, and a Frobenius-normalized symmetric-tensor basis. It then
imposes deliberately conservative executable inequalities.

## Covariant principal system

At one regular background, the fields are

```text
(g_00,g_01,g_02,g_03,g_11,g_12,g_13,g_22,g_23,g_33,phi,chi).
```

For a covector `xi`, the unredefined metric and scalar equations define the
quadratic symbol

```text
P(xi)=A xi_0^2 + B_i xi_0 n_i + C_ij n_i n_j.
```

The metric variables are covariant components `g_ab`. The action-conjugate
metric row is consequently `-E^ab/2`, not `E_ab/2`. Raising both equation
indices and retaining this sign makes the ungauged metric--`phi`--`chi`
matrix the symmetric Hessian of the frozen action. The modified-harmonic block
is added in the same row normalization and changes only the ten-metric block.

The independent `chi` equation remains its canonical scalar principal block.
It therefore contributes one physical characteristic per sign as a direct
sum; it is not silently absorbed into the regulator sector.

## Every spatial covector, not an angular scan

An off-diagonal symmetric metric component is multiplied by `sqrt(2)`. In
this basis its Euclidean norm is the tensor Frobenius norm, and spatial
rotations act orthogonally. HYP2 extracts every coefficient matrix `A`, `B_i`,
and `C_ij` from the complete symbol.

For `|n|=1`, triangle and operator-norm inequalities then bound the entire
companion deformation from flat Einstein plus two canonical scalars:

```text
||M(n)-M_ESF(n)||_2
 <= sum_i ||Delta M_i||_2 + sum_ij ||Delta M_ij||_2.
```

The right-hand side is independent of `n`. Homogeneity extends the result from
unit directions to every nonzero spatial covector. Directional samples are
used only in unit tests as weaker regressions; no finite angular scan supplies
the certificate claim.

The diffeomorphism pure-gauge identity is checked in the same way. Since
`P(xi)` is quadratic and

```text
h_ab = xi_a X_b + xi_b X_a
```

is linear, `P(xi)h(xi,X)` is a homogeneous cubic. The implementation extracts
and symmetrizes every cubic coefficient and bounds their complete operator-norm
sum. Three or any other finite number of chosen directions cannot substitute
for that coefficient test.

## Spectral separation

The flat reference has six clusters,

```text
physical:  -1, +1
tilde:     -1/2, +1/2
hat:       -1/3, +1/3,
```

each with multiplicity four. Physical contours have radius `1/5`; auxiliary
contours have radius `1/20`. The locked lower singular-value bounds for the
reference resolvent are `1/32` and `1/1024`, respectively.

The reference contour calculation uses 512 binary64 SVD samples. The maximum
unsampled arc is covered by the exact chord distance and the one-Lipschitz
property of singular values; an explicit enlarged roundoff allowance is then
subtracted. This is a conservative machine certificate, not a rational
interval or formal proof of the SVD implementation.

If `delta=||M-M_ESF||_2` is below a reference resolvent lower bound, the
Neumann estimate keeps the contour in the resolvent set. The Riesz projectors
therefore retain their dimensions and the six contours continue to partition
all 24 characteristics.

## Physical modes and real diagonalizability

Spectral separation alone would allow complex eigenvalues within a contour.
The physical groups require the second part of the Kovács--Reall argument: the
symmetric ungauged action block defines a positive form on each reference
physical subspace.

HYP2 bounds the candidate/reference Riesz-projector displacement. For a unit
candidate vector `v`, with reference projection `w=P0 v`,

```text
||v-w|| <= ||P-P0||,
||w|| >= 1-||P-P0||.
```

These inequalities are applied directly to the reference form, including its
cross terms; no orthogonality of the generally non-normal Riesz projectors is
assumed. The independently bounded deformation of the action energy form is
then subtracted. A strictly positive remainder supplies the physical
coercivity margin that forces real, diagonalizable physical characteristics
inside the declared ball.

The frozen sufficient limits are

```text
companion deformation              <= 1/4096
action-energy deformation          <= 1/8192
coordinate-time kinetic deformation<= 1/64
effective Planck coefficient       >= 2
action-Hessian symmetry defect     <= 1e-10
pure-gauge cubic coefficient bound <= 1e-10.
```

They are study bounds, not inferred physical constants and not a universal
Horndeski threshold.

## Independent checks and witnesses

The covariant 12-field symbol is restricted to the six spherical
perturbations at the exact activated COMP1 fixture. All twelve resulting roots
agree with the independently implemented exact spherical modified-harmonic
symbol to the frozen binary64 tolerance. This checks field ordering, equation
row normalization, sign, and auxiliary-cone convention.

Four complete ID1-to-REF1 witnesses are then evaluated:

- central amplitude/width and the base regulator seed;
- the same matter data with the regulator seed doubled;
- the high-amplitude, widest-pulse stress slice at its profile centre; and
- that stress slice at its peak-compactness radius.

Each witness first solves the physical initial constraints, reconstructs the
complete spatial second jet, and solves the six complete REF1 accelerations.
Only then is it mapped to a physical orthonormal covariant background. Every
strict inequality passes. By continuity, each witness lies in a nonzero open
local health neighborhood.

The witnesses prove the ball is nonempty and intersects relevant initial
data. They do not enclose every point of a complete hypersurface or any future
trajectory. This distinction is why HLT1 remains a separate required RUN1
artifact.

## What changed relative to DOM3

DOM3 proved a compact nonflat **radial** eigenframe and symmetrizer on its
implicit branch graph. HYP2 addresses the missing angular/covector question
with the full local four-dimensional principal symbol. It does not convert a
spherical evolution into a nonspherical stability result: backgrounds and
later dynamics remain symmetry reduced even though the local principal
perturbation test covers every covector direction.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hyp2_md1.py \
  --output results/fgc-1-hyp2-md1.json
```

The result is canonical unique-key JSON. It consumes DOM4's exact secondary
domain-definition hash and retains PROTO3's semantic hash as the common RUN1
run-envelope binding.

## Nonclaims

HYP2 supplies no Wilsonian operator-completeness or remainder estimate, no
frequency/cutoff authorization, no nonlinear existence time, no boundary
theorem, no trajectory, no collapse or affine-null result, and no
nonspherical robustness result. It does not derive a physical transition,
singularity resolution, child topology, dark-sector mechanism, or variation
of locally measured `c`.
