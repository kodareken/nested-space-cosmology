# FGC-1-HYP1-CON1-COMP1: activated local compatible-constraint witness

**FGC-1-HYP1-CON1-COMP1** composes an exact local compatible datum with the
existing QIFT1 implicit-branch theorem.  At the fixed flat-reference coordinate
radius \(r=4\), it constructs one activated local two-jet with

\[
\phi=\frac1{131072},
\qquad
\partial_r\phi=\frac1{131072},
\]

and proves that the unique QIFT1 coordinate-time acceleration root at that
datum also solves the *unredefined* ACT1 metric equation and both unmodified
scalar equations at that point.

This is one local, off-hypersurface two-jet witness.  It is neither an explicit
list of the root accelerations nor a constraint-compatible initial
hypersurface, much less a constraint-propagation or evolution theorem.

## 1. Fixed local variables and the activated datum

As in [FO1-RC1](fgc-hyp1-fo1-rc1.md) and
[QIFT1](fgc-hyp1-dom1-qift1.md), write the complete REF1 relation as

\[
\mathcal R(a;z)=0,
\qquad
a=p_t,
\qquad
z=(u,p,q,p_r,q_r),
\]

in the frozen base-field order

\[
u=(h_{tt},h_{tr},h_{rr},R,\phi,\chi).
\]

The coordinate radius remains fixed at \(r=4\); it is not one of the thirty
QIFT1 parameter axes.  Start with the flat FGC-QR reference datum and activate
both the regulator value and its radial first jet,
\(u.\phi=q.\phi=1/131072\), retaining the declared nonzero action couplings.
Since

\[
\left|\frac1{131072}\right|<\frac1{65536},
\]

both activated coordinates are strictly inside QIFT1's uniform parameter
half-width before the two spatial second-jet entries below are selected.  The
exact certificate checks that the resulting full \(z_\star\), including those
two entries, lies strictly inside the same QIFT1 component box.  Thus it is a
local activated datum near the flat reference root, not the flat inactive root
itself.  The cone control pulls this complete compatible state, including both
solved spatial second jets, back to the generalized ADM variables and rejects
the calculation unless the round trip reproduces every `SphericalState` jet
exactly.  At that same point, its nonzero regulator first jet makes the
regulator characteristic quadratic distinct from the physical metric
quadratic.  This is a pointwise property only; CON1-COMP1 does not claim a
nonflat uniform-hyperbolicity box or a full-system symmetrizer.

For a rational, normalization-free physical normal direction in the two
dimensional \((t,r)\) base, use

\[
w^a=\left(1,-\frac{h_{tr}}{h_{rr}}\right).
\]

It obeys \(h_{ra}w^a=0\).  No square root, sign choice, or unit normalization
is introduced: \(w\) is deliberately unnormalized and rational on the
declared chart where \(h_{rr}\ne0\).  Define the corresponding physical
projection by

\[
\mathcal H_{ww}:=E_{ab}w^aw^b,
\]

and denote the radial physical momentum projection by \(\mathcal M\).  Here
\(E_{ab}\) is the unredefined ACT1 metric residual, with the project’s fixed
metric-variation convention.

## 2. Exact compatible local construction

All local data other than two selected spatial second jets are frozen by the
CON1-COMP1 configuration.  Those two rational spatial second jets are solved
simultaneously from

\[
\mathcal H_{ww}=0,
\qquad
\nabla_rC^r=0.
\tag{COMP-1}
\]

The resulting exact values are

\[
\partial_r^2h_{tt}
=-\frac{584115552257}{18889465931341141901312},
\qquad
\partial_r^2R
=-\frac{584115552257}{4722366482835285475328}.
\tag{COMP-1a}
\]

The same exact evaluation verifies, at the constructed datum,

\[
C^\mu=0,
\qquad
\nabla_r C^t=0,
\qquad
\mathcal M=0.
\tag{COMP-2}
\]

The gauge quantity is the metric-defined reference-gauge quantity

\[
C^\mu=-\widetilde g^{\rho\sigma}
\left(\Gamma^\mu{}_{\rho\sigma}-\bar\Gamma^\mu{}_{\rho\sigma}\right),
\]

not an independently prescribed \(C^\mu\) probe.  Consequently (COMP-1) and
(COMP-2) are local compatibility conditions on the same metric/scalar two-jet
that enters the complete REF1 evaluator.

The certificate uses exact rational arithmetic.  It checks every denominator,
the nonsingularity of the two-equation solve, the exact residual equalities,
and strict QIFT1-box membership.  Failure of any check rejects the witness
rather than replacing it by a nearby floating-point fit.

## 3. Why the physical projections vanish for every acceleration

Hold \(z=z_\star\) fixed and leave the six coordinate-time accelerations
\(a^I\) formal.  The artifact extracts the exact multivariate coefficients of
\(\mathcal H_{ww}(a;z_\star)\) and \(\mathcal M(a;z_\star)\) through total
degree two.  Every such coefficient is zero.  Therefore

\[
\boxed{
\mathcal H_{ww}(a;z_\star)=0,
\quad
\mathcal M(a;z_\star)=0
\quad\text{for every }a\in\mathbb R^6.
}
\tag{COMP-3}
\]

This is a polynomial identity at one fixed local datum, not a sampled
statement about a finite list of accelerations.  The degree bound has a direct
term-by-term origin:

- the connection and curvature are affine in the selected second jets, hence
  affine in \(a\) at fixed \(z_\star\);
- the Einstein term and the \((g_{ab}\Box-\nabla_a\nabla_b)F\) term are at
  most affine in \(a\), because \(F= M_{\rm Pl}^2+\beta\phi^2\) and its
  Hessian is affine at a fixed first jet;
- the Gauss--Bonnet metric contribution is a contraction of a curvature-linear
  double-dual tensor with a Hessian of \(f(\phi)\), so it is at most the
  product of two affine factors and therefore has total degree at most two;
- the scalar Gauss--Bonnet source likewise contains the curvature-squared
  invariant with fixed \(f'(\phi)\), hence has degree at most two; and
- potentials, first-derivative scalar stresses, the independent \(\chi\)
  stress, and fixed-chart inverse-metric coefficients contribute degree zero
  or one, never a higher acceleration degree.

There are no acceleration-dependent inverse-metric denominators: inverse
metrics and the unnormalized \(w\) depend only on the fixed component data.
Thus exact total-degree-\(\le2\) extraction is exhaustive for these two
physical projections.  A nonzero coefficient of any degree \(0\), \(1\), or
\(2\), or the discovery of a degree above two, falsifies (COMP-3).

## 4. Composition with QIFT1

QIFT1 proves, for every \(z\) in its declared continuous box, the existence
and uniqueness of an acceleration vector in its declared acceleration box:

\[
\forall z\in Z,\quad \exists!\,a(z)\in A
\quad\mathcal R(a(z);z)=0.
\]

Because \(z_\star\in\operatorname{int}Z\), there is therefore one and only
one acceleration vector

\[
a_\star:=a(z_\star)\in A.
\tag{COMP-4}
\]

The construction **defines** \(a_\star\) by this exact existence-and-uniqueness
statement.  It does not print or claim closed-form rational coordinates for
\(a_\star\).  The QIFT1 fixed-preconditioner iteration is an enclosure and
existence tool; it is not being relabelled as an explicit solution formula.

Applying (COMP-3) at \(a=a_\star\) gives the physical Hamiltonian and radial
momentum projections at the unique root without evaluating any numerical root
coordinates.  The datum is exactly shift-free,
\(h_{tr}/h_{rr}=0\), so at this point the two projection identities are also
\(\mathcal H_{ww}=E_{tt}=0\) and \(\mathcal M=E_{tr}=0\) for every formal
acceleration.  This zero-shift fact is an explicit serialized premise of the
composition below, not an implicit specialization of the displayed general
projection formula.

## 5. From the modified root to the unredefined local equations

The complete REF1 relation contains the two modified-harmonic metric root
equations in the \(tt\) and \(tr\) slots, together with the two unmodified
scalar equations.  At \(a_\star\), those root equations hold by (COMP-4),
while (COMP-3) supplies the zero physical projections.

At the constructed datum, \(C^\mu=0\).  Direct tensor evaluation verifies that
every non-normal component of \(\nabla_aC^\mu\) vanishes, including the radial
pair in (COMP-2); the only components not fixed by the parameter datum are its
normal derivatives at the QIFT1 root.  This stronger all-non-normal-row fact is
an explicit serialized premise.  The artifact forms the exact linear map from
the two remaining normal \(\nabla C\) components to the \(tt,tr\) components of
the full gauge extension

\[
X_{ab}=F\,\widehat P_{\alpha}{}^{\beta}{}_{ab}\nabla_\beta C^\alpha.
\]

It verifies that this map is invertible at \(z_\star\).  Because zero shift and
(COMP-3) give \(E_{tt}=E_{tr}=0\) for every acceleration, the modified
\(tt,tr\) root equations force the corresponding extension components to
vanish; invertibility gives

\[
\nabla_{\perp}C^\mu=0.
\tag{COMP-5}
\]

Together with (COMP-2), this proves the complete first derivative vanishes,

\[
\nabla_aC^\mu=0,
\qquad X_{ab}=0.
\tag{COMP-6}
\]

Hence the unique modified-harmonic root has zero extension and satisfies the
unredefined ACT1 metric equation locally.  The complete residual evaluator
independently checks and serializes that the two REF1 scalar rows equal the two
unmodified ACT1 scalar residuals exactly, so the same root satisfies both
scalar equations locally as well.

This is a local algebraic composition theorem only.  It does not assert that
\(C^\mu\) remains zero at another point, at another time, or under a numerical
evolution.

## 6. Relation to the primary literature

The structural precedent is Ková​cs and Reall,
[*Well-posed formulation of Lovelock and Horndeski theories*](https://arxiv.org/abs/2003.08398),
especially Sec. 2.2 Eq. (12) and Sec. 4 Eqs. (86)--(88): a
modified-harmonic metric equation, the unmodified scalar equation, and the
diffeomorphism identity can compose into a gauge-constraint argument when
compatible physical data are available.  Their result is a one-scalar
weak-coupling formulation theorem; it does not establish this dynamic-
\(F(\phi)\), Gauss--Bonnet, independent-\(\chi\) witness or its propagation.

Gundlach and Martín-García,
[*Hyperbolicity of second-order in space systems of evolution equations*](https://arxiv.org/abs/gr-qc/0506037),
Sec. III Eqs. (15)--(21), provide the standard distinction between defining
first-order reduction constraints and deriving a closed homogeneous
reduction-constraint system.  CON1-COMP1 does neither the latter derivation
nor a hyperbolicity proof; its use of local jets must not be described as
reduction-constraint propagation.

East and Ripley,
[*Evolution of Einstein-scalar-Gauss-Bonnet gravity using a modified harmonic formulation*](https://arxiv.org/abs/2011.03547),
Eqs. (1)--(4) and (13)--(15), are the closest compatible-data comparator.
They explicitly use a special zero-scalar initial-data branch rather than a
general scalar--Gauss--Bonnet constraint solve.  Their result is therefore a
caution against importing a general activated ACT1 initial-data family from a
GR-like special branch.

## 7. What this does not establish

CON1-COMP1 does **not** establish any of the following:

- a full initial hypersurface or a nonzero-width family of compatible data;
- regular-centre or outer-boundary conditions, an IBVP, or a finite-mass
  collapse datum;
- metric-defined gauge, reduction, Hamiltonian, momentum, or other physical
  constraint propagation;
- a complete nonlinear first-order solved evolution system;
- a nonflat uniform strong-hyperbolicity box, complete real characteristic
  structure, or full-system symmetrizer;
- an EFT cutoff, retained-operator ledger, or omitted-operator envelope;
- evolution, independent-\(\chi\) collapse, a restoring interval, metric-null
  affine defocusing, singularity resolution, or a child spacetime.

The witness is falsified if the exact two-jet solve is singular or leaves the
QIFT1 box; if any declared physical/gauge residual is nonzero; if the
degree-\(\le2\) coefficient ledger contains a nonzero entry; if QIFT1 fails
to provide its unique root; if the normal-gauge-extension map is singular; or
if the scalar rows are found not to be the unmodified ACT1 scalar residuals.

## 8. Reproduction

Run

```bash
python3 scripts/reproduce_fgc_hyp1_con1_comp1.py \\
  --config configs/fgc/fgc-1-hyp1-con1-comp1.toml \\
  --output results/fgc-1-hyp1-con1-comp1.json
```

The generated record is the authoritative machine-readable certificate.  It
binds the frozen configuration, the exact compatible-constraint construction,
the QIFT1 predecessor and box-membership proof, the coefficient extraction,
the normal-extension invertibility check, the source implementation, this
owner document, and the fail-closed nonclaims.
