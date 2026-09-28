# FGC-1-HYP1-MHG3-IMP1: local implicit acceleration branch

**FGC-1-HYP1-MHG3-IMP1** closes one narrow algebraic prerequisite left open
by REF1. At one exact flat FGC-QR vacuum/reference-gauge two-jet, it
differentiates the complete six-equation REF1 residual with respect to the six
base-field coordinate-time accelerations, obtains an exact nonsingular
Jacobian, and applies the finite-dimensional real implicit function theorem.

The result is one unique smooth local acceleration branch through that root.
Its neighborhood exists but is not quantified. This is not yet a first-order
PDE formulation, a propagation theorem, an open hyperbolicity or retained-EFT
domain, constraint-compatible initial data, or permission to evolve.

## Residual map and frozen variable order

The complete REF1 equation order is

\[
\mathcal R_A=
\left(
E^{\rm MHG}_{tt},
E^{\rm MHG}_{tr},
E^{\rm MHG}_{rr},
E^{\rm MHG}_{\theta\theta},
E_\phi,
E_\chi
\right).
\]

IMP1 uses the same base-field order as RED1, MHG1, and REF1,

\[
u_B=(h_{tt},h_{tr},h_{rr},R,\phi,\chi),
\]

and freezes the acceleration vector

\[
\boxed{
a_B=
\left(
\partial_t^2h_{tt},
\partial_t^2h_{tr},
\partial_t^2h_{rr},
\partial_t^2R,
\partial_t^2\phi,
\partial_t^2\chi
\right).
}
\]

These are base-metric coordinate-time accelerations. They are not silently
relabeled generalized-ADM accelerations. Away from the flat root, a change to
`(alpha,shift,lambda,R,phi,chi)` requires its own nonlinear chain rule.

Let `z` denote the remaining local data: the six field values, their first
time and radial derivatives, their mixed `tr` and spatial `rr` second
derivatives, together with the fixed action parameters, reference connection,
and reference coordinate radius. The finite-dimensional residual map is

\[
\mathfrak F(a;z)=\mathcal R(u,\partial u,a,\partial_t\partial_r u,
\partial_r^2u).
\]

## Exact forward tangent, not finite differencing

IMP1 adds a dependency-free exact scalar

\[
x+\varepsilon\dot x,
\qquad \varepsilon^2=0.
\]

For example,

\[
(x+\varepsilon\dot x)(y+\varepsilon\dot y)
=xy+\varepsilon(\dot x y+x\dot y),
\]

and

\[
\frac{1}{x+\varepsilon\dot x}
=\frac1x-\varepsilon\frac{\dot x}{x^2}.
\]

One acceleration slot is seeded with unit tangent and the complete REF1
residual evaluator is run unchanged. The tangent parts of all six output rows
are therefore the exact column

\[
K_A{}^B=
\frac{\partial\mathfrak F_A}{\partial a_B}.
\]

The primal part of every seeded run must exactly reproduce the unseeded root
residual. This guards against replacing the nonlinear residual evaluation by a
precomputed principal matrix. No floating-point perturbation, finite step,
tolerance, or computer-algebra dependency enters the derivative.

## Declared FGC-QR root

At reference radius

\[
r_0=4>r_{\min}=\frac12,
\]

the root has

\[
h_{tt}=-1,\qquad h_{tr}=0,\qquad h_{rr}=1,\qquad
R=r_0,qquad \partial_rR=1,
\]

with every other metric and areal-radius first or second jet zero, and

\[
\phi=\chi=0,qquad
\partial\phi=\partial\chi=0,qquad
\partial^2\phi=\partial^2\chi=0.
\]

The action is not switched to GR-0. IMP1 retains the exact FGC-QR parameters

\[
M_{\rm Pl}=2,qquad
\beta=-\frac14,qquad
\mu=3,qquad
g_4=\frac12,qquad
\alpha=0,qquad
\eta=\frac12.
\]

Thus

\[
F(0)=M_{\rm Pl}^2=4,qquad
F'(0)=0,qquad f'(0)=0.
\]

The nonzero FGC-QR couplings are present but inactive at this vacuum root.
Flat spherical physical and reference connections agree, so

\[
C^\mu=0,qquad \nabla_\nu C^\mu=0,qquad
F\widehat P\nabla C=0.
\]

The unredefined metric and both scalar equations also vanish. Therefore the
complete gauge-fixed residual obeys

\[
\boxed{\mathfrak F(0;z_0)=0.}
\]

This is an exact local equation root. It is not a constructed spacetime or a
nontrivial regulator configuration.

## Exact acceleration Jacobian

The exact forward-tangent calculation returns

\[
K=
\begin{pmatrix}
36&0&9&9&0&0\\
0&72&0&0&0&0\\
4&0&1&-1&0&0\\
64&0&-16&0&0&0\\
0&0&0&0&-1&0\\
0&0&0&0&0&-1
\end{pmatrix}.
\]

It verifies exact rank six, two-sided inverse identities over the rationals,
and

\[
\boxed{\det K=-165888\ne0.}
\]

The calculation also independently constructs the complete MHG1 radial
symbol for the same FGC-QR state. Entry by entry, `K` equals the coefficient
of `c^2` in

\[
\mathcal P_{\rm MHG}(-c,1),
\]

namely MHG1's complete coordinate-time kinetic block. This comparison is
performed after differentiating the complete REF1 residual; it is a
regression, not the derivative's source.

## Smooth finite-dimensional domain at the root

The direct ACT1/REF1 evaluator is rational in its finite local jet coordinates
on the regular chart, with polynomial curvature and Hessian contractions in
the numerator. IMP1 checks at the root that

\[
\det h_{AB}=-1<0,qquad
g^{00}=-1<0,qquad
R=4>0,qquad
F=4>0,
\]

that the reference radius is strictly inside its annulus, and that both
auxiliary inverse metrics remain Lorentzian. All denominators used by the
physical inverse metric, the normal outer product, the auxiliary metrics, and
the reference connection are therefore nonzero there. By continuity they
remain nonzero on some unspecified neighborhood, and the residual map is
smooth on that neighborhood.

## Implicit-function conclusion

At `(a_0,z_0)=(0,z_0)`, IMP1 has established

\[
\mathfrak F(a_0;z_0)=0,
\qquad
\det\left(
\frac{\partial\mathfrak F}{\partial a}
\right)_{(a_0,z_0)}\ne0.
\]

The ordinary finite-dimensional real implicit function theorem therefore
gives neighborhoods `U_a` of `a_0` and `U_z` of `z_0`, and one unique smooth
map

\[
\boxed{a=\mathcal A(z),\qquad \mathcal A(z_0)=0,}
\]

such that

\[
\mathfrak F(\mathcal A(z);z)=0
\]

for `z` in `U_z` and `\mathcal A(z)` in `U_a`.

The neighborhoods are existential. IMP1 publishes no radius, operator-norm
bound, condition-number estimate, retained-EFT margin, or uniform
hyperbolicity margin. A later domain gate must derive those quantities rather
than treating continuity as a numerical certificate.

Kovacs and Reall likewise separate invertibility of a coordinate-time
second-derivative block from their subsequent weak-coupling
strong-hyperbolicity and symmetrizer analysis. Their one-scalar Horndeski
result motivates the architecture but does not prove this ACT1 calculation;
IMP1 derives the complete ACT1/REF1 Jacobian independently
([arXiv:2003.08398](https://arxiv.org/abs/2003.08398)).

## What this changes

REF1 had an exact complete residual but no proof that its six equations could
be solved locally for its six coordinate-time accelerations. IMP1 closes that
one-point algebraic obstruction for the FGC-QR branch. The project no longer
needs to assume a regular local acceleration branch at the declared flat
vacuum root.

The downstream [FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate now
derives and cross-checks the conditional lower-order gauge operator on an
independent spherical `C^mu` two-jet. The project still does not have a
complete nonlinear first-order PDE. The downstream
[FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) gate now introduces the exact
eighteen-variable implicit lift, derives the complete flat-root branch
derivative and linearized matrices, and closes the chosen kinematic radial
reduction identity. Constructing a usable system still requires a controlled
nonlinear source map, metric-defined homogeneous gauge and physical-constraint
propagation, and proof that the resulting
characteristic system is uniformly controlled on a quantified retained-EFT
domain. The downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md)
certificate quantifies the local REF1 branch on a full-dimensional rational
box. The later [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result selects
one activated interior point and proves that its unique root is locally
metric-defined compatible with the unredefined equations. Neither downstream
result closes propagation, an initial hypersurface, or the retained-domain
requirements.

## Fail-closed nonclaims

IMP1 does not establish:

- an IMP1-native quantified implicit-branch neighborhood (the downstream
  QIFT1 record supplies one for the frozen REF1 formulation);
- a nonzero activated solution, constraint-compatible initial data, or a
  nontrivial spacetime;
- complete lower-order first-order sources;
- gauge-constraint, reduction-constraint, or physical-initial-constraint
  propagation;
- a uniform eigenbasis, symmetrizer, open hyperbolicity domain, or retained
  EFT box;
- a regular centre or boundary system;
- any evolution, collapse, metric-null affine defocusing, singularity
  resolution, child spacetime, dark sector, varying local light speed,
  spectrum, or theory of everything.

The exact result is narrower and useful: the complete REF1 residual has one
regular local FGC-QR coordinate-time acceleration branch through its declared
flat vacuum/reference-gauge root.
