# FGC-1-HYP1-MHG1: modified-harmonic principal formulation gate

## Status and result

**HYP1-MHG1** adds a complete gauge-fixed spherical principal system to the
unredefined ACT1/VAR1 equations.  It is the next exact step after RED1's raw
second-jet matrix and SYM1's physical gauge quotient.  The selected
formulation is the modified-harmonic construction of Kovacs and Reall, not
ordinary generalized harmonic gauge.

On each frozen exact fixture, MHG1:

1. constructs two Lorentzian auxiliary inverse metrics whose unphysical
   gauge and gauge-constraint cones are distinct from each other and from all
   physical scalar characteristics;
2. adds the exact spherical restriction of the modified-harmonic
   gauge-fixing principal term to the unredefined VAR1 metric equations;
3. derives the complete six-by-six second-order radial symbol and the standard
   twelve-variable first-order principal reduction;
4. factorizes the complete determinant exactly as

   \[
   \det\mathcal P_{\rm MHG}(c)
   =K\,\widetilde N(c)^2\widehat N(c)^2N_g(c)Z_\phi(c),
   \qquad K\ne0;
   \]

5. verifies the required two-dimensional kernels at every repeated rational
   gauge root and, where the physical factors coincide in a control, at every
   repeated physical root; and
6. proves from the exact multiplicities, positive scalar discriminants,
   nonzero cross-cone resultants, and kernel dimensions that the frozen radial
   first-order matrices have complete real pointwise eigenbases.

The activated FGC-QR fixture retains SYM1's distinct regulator quadratic.  Its
two roots do not intersect either auxiliary cone or the physical metric cone.
The three `X=0` controls remain regular because the implementation works
directly with the smooth VAR1 equations and never divides by
`X=-nabla(phi)^2/2` or introduces `log X`.

This is **not** the completed HYP1 theorem; it is not yet the open-domain theorem.
MHG1 has not serialized the
complete lower-order gauge-fixed source, the lower-order gauge-constraint
wave coefficients, or the propagation of all first-order reduction
constraints.  It has not produced a uniform symmetrizer or an explicit
retained EFT box.  No evolution, collapse, or defocusing claim is authorized.

## Why modified harmonic rather than conventional generalized harmonic

The conventional generalized-harmonic construction leaves physical,
pure-gauge, and gauge-condition-violating polarizations on coincident cones.
That degeneracy is harmless in GR but is generically unstable under Horndeski
principal corrections. Papallo's analysis shows that
Einstein--dilaton--Gauss--Bonnet gravity is not strongly hyperbolic in any
generalized-harmonic gauge of the class analyzed there on generic weak-field
backgrounds
([arXiv:1710.10155](https://arxiv.org/abs/1710.10155)).

Kovacs and Reall instead define a distinct modified-harmonic formulation and
assign the pure-gauge and gauge-constraint sectors two separate auxiliary
cones. They prove strong hyperbolicity for their weakly coupled one-scalar
Horndeski system in that formulation, subject to their auxiliary-cone
hypotheses
([arXiv:2003.08398](https://arxiv.org/abs/2003.08398)).  East and Ripley give a
full nonlinear ESGB implementation of the same architecture
([arXiv:2011.03547](https://arxiv.org/abs/2011.03547)).  Those sources motivate
and constrain MHG1; they are not substituted for an ACT1 calculation.

ACT1 differs in two material ways.  Its Einstein coefficient

\[
F(\phi)=M_{\rm Pl}^2+\beta\phi^2
\]

is dynamic, and it has an independent canonical matter field `chi`.  MHG1
therefore adds the gauge term to ACT1's own VAR1 tensor and retains `chi` in
the full principal system.  Positivity of `F` is a hard local condition.

## Physical equations and gauge-surface equivalence

Write the unredefined equations as

\[
E_{\mu\nu}=0,
\qquad E_\phi=0,
\qquad E_\chi=0,
\]

where `E_mn` is exactly the VAR1 bulk equation,

\[
E_{\mu\nu}
=F G_{\mu\nu}
 +(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)F
 -T^{(\phi)}_{\mu\nu}-T^{(\chi)}_{\mu\nu}
 -T^{({\rm GB})}_{\mu\nu}.
\]

Introduce an auxiliary inverse metric `tilde g` and the gauge constraint

\[
H^\mu=-\widetilde g^{\rho\sigma}
\Gamma^\mu{}_{\rho\sigma}.
\]

For a production spherical-coordinate implementation, `Gamma` is replaced by
the difference from a declared smooth reference connection or an equivalent
gauge source.  That replacement changes lower-order terms, not the radial
principal symbol certified here.

Define

\[
\widehat P_\alpha{}^{\beta\mu\nu}
=\delta_\alpha^{(\mu}\widehat g^{\nu)\beta}
 -\frac12\delta_\alpha^\beta\widehat g^{\mu\nu}.
\]

The selected ACT1 gauge-fixed metric equation is

\[
\boxed{
E_{\rm MHG}^{\mu\nu}
=E^{\mu\nu}
 +F\,\widehat P_\alpha{}^{\beta\mu\nu}
   \partial_\beta H^\alpha=0.
}
\]

The positive factor `F` matches the local Einstein principal normalization;
its derivatives contribute only lower-order terms to gauge-constraint
propagation.  The scalar equations are not gauge modified.

If `H` and its first derivative vanish, the added term vanishes exactly and

\[
E_{\rm MHG}^{\mu\nu}=E^{\mu\nu}.
\]

Thus MHG is an off-gauge-surface extension, not a redefinition of an ACT1
solution.  Physical equivalence still requires compatible initial ACT1
constraints and homogeneous propagation of `H`; the principal part of that
propagation is derived below, while its complete lower-order serialization is
the next gate.

## Auxiliary cone convention

Let `n^mu` be the physical unit normal to the coordinate-time slices.  MHG1
freezes

\[
\widetilde g^{\mu\nu}
=g^{\mu\nu}-(q_{\sim}-1)n^\mu n^\nu,
\qquad q_{\sim}=4,
\]

\[
\widehat g^{\mu\nu}
=g^{\mu\nu}-(q_{\wedge}-1)n^\mu n^\nu,
\qquad q_{\wedge}=9.
\]

No irrational normalization enters the executable.  It uses the identity

\[
n^\mu n^\nu
=-\frac{g^{\mu0}g^{\nu0}}{g^{00}}.
\]

For the generalized radial ADM chart

\[
ds^2=-\alpha^2dt^2+\Lambda^2(dr+vdt)^2+R^2d\Omega^2,
\]

the physical radial coordinate speeds are

\[
c_{g,\pm}=-v\pm\frac{\alpha}{\Lambda},
\]

while the unphysical auxiliary roots are

\[
c_{\sim,\pm}=-v\pm\frac{\alpha}{2\Lambda},
\qquad
c_{\wedge,\pm}=-v\pm\frac{\alpha}{3\Lambda}.
\]

Their tangent-space cones lie strictly inside the physical tangent cone; the
corresponding cotangent cones lie outside it, in the convention used by
Kovacs--Reall.  These speeds govern coordinates and gauge-constraint errors,
not matter, photons, or information.  They imply no variable local speed of
light and no superluminal signal.

## Exact gauge-fixed radial symbol

Use the RED1/SYM1 perturbation order

\[
t=(\delta h_{tt},\delta h_{tr},\delta h_{rr},
   \delta R,\delta\phi,\delta\chi)
\]

and radial covector `xi_A=(-c,1)`.  At principal order the linearized
gauge-fixing response is evaluated directly from Kovacs--Reall's projector
formula,

\[
\mathcal P_{\rm GF}(\xi)^{\mu\nu\rho\sigma}t_{\rho\sigma}
=-F\,\widehat P_\alpha{}^{\gamma\mu\nu}\xi_\gamma
g^{\alpha\beta}
\widetilde P_\beta{}^{\delta\rho\sigma}\xi_\delta
t_{\rho\sigma}.
\]

The executable lowers the response with the physical metric and maps
`delta R` to

\[
t_{\theta\theta}=2R\,\delta R,
\qquad
t_{\varphi\varphi}=2R\,\delta R
\]

at the regular equatorial evaluation point.  It then adds that four-row block
to RED1's unredefined six-by-six polynomial symbol.  Neither scalar row is
modified.

The complete determinant obeys, coefficient by coefficient,

\[
\boxed{
\det\mathcal P_{\rm MHG}(c)
=K\,\widetilde N(c)^2\widehat N(c)^2
 N_g(c)Z_\phi(c),\qquad K\ne0.
}
\]

`N_g Z_phi` is not recomputed by assumption: it is the exact SYM1 scalar Schur
determinant, independently obtained after quotienting the covariant
diffeomorphism/Bianchi degeneracy.  Agreement therefore checks that gauge
fixing lifts only the unphysical degeneracy and leaves the physical factors
unchanged.

## Standard first-order principal reduction

Write the complete second-order spherical system as

\[
A\,u_{tt}+B\,u_{tr}+C\,u_{rr}=\text{lower order},
\qquad
u=(h_{tt},h_{tr},h_{rr},R,\phi,\chi)^T.
\]

With `p=u_t` and `q=u_r`, MHG1 constructs the twelve-variable principal
system

\[
\begin{pmatrix}A&0\\0&I\end{pmatrix}
\partial_t\begin{pmatrix}p\\q\end{pmatrix}
+
\begin{pmatrix}B&C\\-I&0\end{pmatrix}
\partial_r\begin{pmatrix}p\\q\end{pmatrix}
=0.
\]

Every frozen coordinate-time kinetic matrix `A` is exactly invertible. A
block determinant identity predicts

\[
\det[-cA^t+A^r]=\det\mathcal P_{\rm MHG}(c),
\]

The executable constructs the full twelve-by-twelve polynomial pencil
`-c A^t+A^r`, evaluates its determinant with fraction-free exact elimination,
and verifies coefficient by coefficient that it equals the independently
constructed six-by-six second-order determinant. It separately verifies that
`det A` is the leading coefficient of that complete polynomial.

The two roots of each auxiliary quadratic have algebraic multiplicity two.
Exact rational rank calculations on both the six-by-six second-order symbol
and the twelve-by-twelve first-order pencil give kernel dimension two at each
of those roots; the two dimensions agree through the standard `p=-c q`
eigenvector map. In the three controls, `Z_phi` is proportional to `N_g`;
the same two-rank calculation gives the required physical kernel dimension
two. In the activated fixture, `N_g` and `Z_phi` have nonzero resultant, and
each physical or regulator root is simple. All four quadratics have positive
exact discriminants, and every auxiliary-to-physical resultant is nonzero.
A simple determinant root has a one-dimensional eigenspace, so no algebraic
extension rank calculation is needed at the irrational regulator roots.

Consequently, the frozen twelve-by-twelve radial first-order matrices have
complete real eigenbases.  This is a **pointwise** statement.  It does not yet
supply a smooth uniformly bounded eigenbasis or positive symmetrizer over a
declared open state domain.

## Activated exact speed controls

For `FGCQR_activated_generic`, the exact roots are

\[
c_{\sim}=\left\{-\frac{17}{84},\frac{73}{84}\right\},
\qquad
c_{\wedge}=\left\{-\frac1{42},\frac{29}{42}\right\},
\]

\[
c_g=\left\{-\frac{31}{42},\frac{59}{42}\right\}.
\]

The regulator roots retain SYM1's exact rational isolating intervals rather
than being rounded into a proof.  Their resultant with every auxiliary and
physical quadratic is nonzero.

These fixture values are algebra controls, not physical parameter estimates
or solution data.

## Gauge-constraint propagation boundary

The covariant ACT1 Noether identity has the schematic form

\[
\nabla_\mu E^{\mu\nu}
=E_\phi\nabla^\nu\phi+E_\chi\nabla^\nu\chi
\]

up to the fixed equation-normalization convention.  On the two scalar
equations, taking the divergence of the MHG metric equation yields a
homogeneous equation for `H`.  Its principal part is

\[
\boxed{
\frac{F}{2}\widehat g^{\rho\sigma}
\partial_\rho\partial_\sigma H^\mu.
}
\]

MHG1 does not insert this block as an expected answer. The executable performs
the full four-index projector contraction

\[
Q^\nu{}_\alpha(\xi)
=F\,\xi_\mu\widehat P_\alpha{}^{\beta\mu\nu}\xi_\beta
\]

and verifies every entry of the four-by-four identity

\[
Q^\nu{}_\alpha(\xi)
=\frac{F}{2}\widehat g^{\rho\sigma}\xi_\rho\xi_\sigma
\delta^\nu_\alpha.
\]

Its spherical `t,r` restriction is therefore two copies of `hat N(c)` after
monic normalization. Because `F>0` and `hat g` is Lorentzian, this is a
hyperbolic principal propagation channel for gauge-constraint error.

The omitted terms are linear and homogeneous in `H` and its first derivatives
for the smooth auxiliary-metric construction, as in the source theorem.  The
project nevertheless keeps the lower-order propagation gate false until
those ACT1-specific coefficients, the reference-connection choice, boundary
conditions, and the first-order reduction constraints are all written and
checked in one system.

## Direct regularity at `X=0`

At `X=0`, a change to a standard Horndeski `G_i(phi,X)` representation would
require a separate check that its individual coefficient functions and their
cancellations remain regular. MHG1 does not invoke such a rewriting or any
`log X` coefficient.

MHG1 instead uses VAR1's smooth tensor

\[
T^{({\rm GB})}_{\mu\nu}
=-8P_{\mu\rho\nu\sigma}\nabla^\rho\nabla^\sigma f
\]

and exact rational jets.  No coefficient divides by `X`, `phi`, or a scalar
gradient.  GR-0, the constant-field SGB-L control, and the zero-regulator
Schwarzschild control all have zero local scalar gradient and pass the same
complete factorization and kinetic checks.  This proves regularity of the
implemented principal coefficients at those `X=0` points; it is not a global
theorem about every lower-order chart or boundary.

## Weak coupling and the remaining HYP1 burden

For their one-scalar Horndeski system, Kovacs--Reall prove an open
strong-hyperbolicity neighborhood when the Horndeski correction is weak
relative to the Einstein--canonical-scalar principal system and the auxiliary
cones remain separated. MHG1's exact
fixtures demonstrate the needed pointwise structure, including one activated
point, but finite fixtures do not define a retained domain.

The downstream
[FGC-1-HYP1-MHG2-REF1](fgc-hyp1-mhg-reference.md) gate now freezes the flat
spherical reference connection on a positive-radius annulus, evaluates the
complete `F hat_P nabla C` gauge extension, verifies exact
gauge-surface/off-gauge controls, and recovers this MHG1 block entry by entry.
It also exposes that ACT1 is generally nonlinear in second derivatives, so an
implicit acceleration-branch gate must precede the domain theorem. The
downstream [FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate now proves one
such branch at an exact flat FGC-QR vacuum/reference-gauge root, but its
neighborhood is existential and unquantified. The downstream
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate now derives the
conditional lower-order reference-gauge operator on independent spherical
`C^mu`, cross-checks two exact assemblies, and recovers this MHG1 principal
block. It does not prove that the metric-defined gauge and reduction
constraints propagate.

QIFT1 now supplies the compact rational box and strict regular-chart margins
needed for the local REF1 implicit branch, including a uniform acceleration
inverse bound. It does not supply the remaining characteristic, constraint, or
validity quantities below. The eventual retained-domain health gate must
publish strict lower bounds for at least:

- `F` and the coordinate-time kinetic determinant;
- the auxiliary-to-physical cone resultants;
- the physical and regulator discriminants;
- every denominator introduced by the first-order solve;
- the smallest eigenvalue of a full-system symmetrizer, or an equivalent
  uniform eigenbasis condition;
- the EFT/weak-coupling operator margin relative to the declared GR control.

The later [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result selects one
activated point strictly inside QIFT1, verifies a distinct physical/regulator
cone there, and composes one metric-defined compatible local root. It is a
point control, not the uniform eigenbasis/symmetrizer required above.

Before that domain gate, the project must derive metric-defined homogeneous gauge and
reduction-constraint propagation and initial/boundary compatibility. If no
useful open box survives, or if the
hoped-for restoring regime lies beyond it, HYP1 fails there even though MHG1's
pointwise fixtures remain algebraically correct.

## Fail-closed nonclaims

MHG1 does not establish:

- complete lower-order first-order evolution equations;
- serialized lower-order gauge-constraint coefficients;
- propagation of all first-order reduction constraints;
- a uniform open-domain symmetrizer or retained EFT region;
- a regular centre, horizon boundary system, or constrained initial data;
- a nonlinear evolution or collapse solution;
- an affine Raychaudhuri sign or magnitude;
- singularity resolution, a child spacetime, a dark sector, varying local
  light speed, a particle spectrum, or a theory of everything.

The exact result is narrower and useful: the selected modified-harmonic lift
removes the spherical gauge degeneracy without altering the independently
derived physical characteristic factors, and the entire frozen radial
first-order principal matrix is real and diagonalizable at every declared
fixture.

The downstream [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) artifact now embeds
that principal block in an exact local eighteen-variable implicit relation and
publishes the complete flat-root linearized matrices. Its downstream
[FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) certificate now provides a
quantified full-dimensional local nonlinear branch box, but not complete
constraint propagation, uniform nonflat hyperbolicity, or retained EFT
validity required to promote MHG1 into an evolution theorem. COMP1's one
activated compatible local root does not change that open-domain boundary.

## Primary sources

- Kovacs and Reall,
  [*Well-posed formulation of Lovelock and Horndeski theories*](https://arxiv.org/abs/2003.08398):
  modified-harmonic gauge, auxiliary-cone separation, gauge propagation, and
  weak-coupling strong-hyperbolicity theorem.
- East and Ripley,
  [*Evolution of Einstein-scalar-Gauss-Bonnet gravity using a modified harmonic formulation*](https://arxiv.org/abs/2011.03547):
  full nonlinear ESGB implementation and constraint-damping comparator.
- Papallo,
  [*On the hyperbolicity of the most general Horndeski theory*](https://arxiv.org/abs/1710.10155):
  obstruction for EdGB in the generalized-harmonic class analyzed there.
- Thaalba et al.,
  [*Hyperbolicity in scalar-Gauss-Bonnet gravity: a gauge invariant study for spherical evolution*](https://arxiv.org/abs/2410.16264):
  independent physical-symbol and gauge-quotient comparator.
