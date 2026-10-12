# FGC-1-HYP1-MHG2-REF1: full reference-gauge residual

**FGC-1-HYP1-MHG2-REF1** completes one equation-level prerequisite left open
by MHG1.  It freezes a genuine reference connection on a positive-radius
spherical annulus, evaluates the full contravariant gauge constraint and its
covariant derivative on exact rational two-jets, adds the resulting full
modified-harmonic gauge extension to the unredefined ACT1/VAR1 metric equation,
and differentiates that complete extension back to MHG1's principal gauge
block.
The extension itself is a second-order term.  Its reference-, connection-, and
source-dependent remainder supplies the lower-differential-order information
that was absent from MHG1.

This is not yet a first-order evolution system. ACT1's Gauss--Bonnet source
is generally nonlinear in second derivatives, so REF1 itself does not identify
or certify a regular local implicit acceleration branch. The downstream
[FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate now closes that algebraic
obstruction at one exact flat FGC-QR vacuum/reference-gauge root only. REF1
is now followed by [FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md), which
derives the conditional lower-order operator for independent spherical
`C^mu`. REF1 and PROP1 still do not prove the metric-derived gauge-propagation equation, reduction-constraint
propagation, an open retained-EFT domain, a regular centre, collapse, or
defocusing.

## Why MHG1 needed an equation-level completion

MHG1 used the principal gauge condition

\[
H^\mu=-\widetilde g^{\rho\sigma}
\Gamma^\mu{}_{\rho\sigma}
\]

and correctly noted that a production spherical implementation must replace
the bare Christoffel symbols by a connection difference or an equivalent
declared gauge source.  Bare `Gamma` is not a tensor, and its spherical
coordinate terms are nonzero even in flat spacetime.  That omission does not
alter MHG1's highest-derivative matrix, but it prevents MHG1 alone from being a
complete gauge-fixed equation.

Kovacs and Reall define the modified-harmonic principal architecture and also
allow prescribed coordinate gauge sources.  Their construction supplies the
formulation pattern, not an ACT1-specific gauge-fixed equation
([arXiv:2003.08398](https://arxiv.org/abs/2003.08398)).  East and Ripley show
that the same architecture can be implemented in a neighboring one-scalar
Einstein--scalar--Gauss--Bonnet system; their source choices and numerical
health region are not transferred to ACT1
([arXiv:2011.03547](https://arxiv.org/abs/2011.03547)).

## Frozen reference connection and domain

REF1 selects the Levi-Civita connection of the fixed reference metric

\[
d\bar s^2=-dt^2+dr^2+r^2d\Omega^2
\]

in the coordinate order `(t,r,theta,phi)`.  It is evaluated at
`theta=pi/2` on the declared annulus

\[
r>r_{\min}=\frac12.
\]

The centre is not included.  That boundary is explicit: no spherical
coordinate cancellation is promoted into a regular-centre theorem.

At the equator, the nonzero reference components include

\[
\bar\Gamma^r{}_{\theta\theta}
=\bar\Gamma^r{}_{\varphi\varphi}=-r,
\]

\[
\bar\Gamma^\theta{}_{r\theta}
=\bar\Gamma^\theta{}_{\theta r}
=\bar\Gamma^\varphi{}_{r\varphi}
=\bar\Gamma^\varphi{}_{\varphi r}=\frac1r.
\]

The intrinsic two-sphere Christoffels that are zero at the equator still have
nonzero angular derivatives.  REF1 retains, in particular,

\[
\partial_\theta\bar\Gamma^\theta{}_{\varphi\varphi}=1,
\qquad
\partial_\theta\bar\Gamma^\varphi{}_{\theta\varphi}
=\partial_\theta\bar\Gamma^\varphi{}_{\varphi\theta}=-1.
\]

Those terms prevent a radial-only surrogate from passing as the full
four-dimensional spherical connection.

The reference coordinate radius is fixed input data.  It is not the evolving
areal-radius field `R`.  The four algebra fixtures evaluate the reference at
`r=4,5,8,9/2`, respectively.

## Tensorial gauge constraint

Let

\[
\Delta\Gamma^\mu{}_{\rho\sigma}
=\Gamma^\mu{}_{\rho\sigma}
-\bar\Gamma^\mu{}_{\rho\sigma}.
\]

The difference of two connections is a tensor.  REF1 therefore freezes the
contravariant gauge constraint

\[
\boxed{
C^\mu=-\widetilde g^{\rho\sigma}
\Delta\Gamma^\mu{}_{\rho\sigma}.
}
\]

It uses MHG1's same exact auxiliary inverse metrics,

\[
\widetilde g^{\mu\nu}
=g^{\mu\nu}-3n^\mu n^\nu,
\qquad
\widehat g^{\mu\nu}
=g^{\mu\nu}-8n^\mu n^\nu,
\]

corresponding to `q_tilde=4` and `q_hat=9`.  The normal outer product is
evaluated without a square root,

\[
n^\mu n^\nu
=-\frac{g^{\mu0}g^{\nu0}}{g^{00}}.
\]

The executable differentiates this expression exactly.  With `A` a
coordinate index,

\[
\partial_A C^\mu
=-
(\partial_A\widetilde g^{\rho\sigma})
\Delta\Gamma^\mu{}_{\rho\sigma}
-\widetilde g^{\rho\sigma}
\partial_A\Delta\Gamma^\mu{}_{\rho\sigma},
\]

and

\[
\boxed{
\nabla_A C^\mu
=\partial_A C^\mu
+\Gamma^\mu{}_{A\nu}C^\nu.
}
\]

Thus the complete gauge data depend on the metric through its second jet, but
not on any undeclared third derivative.

## Full reference-gauge modified-harmonic equation

The unredefined ACT1 equations remain

\[
E_{\mu\nu}=0,
\qquad E_\phi=0,
\qquad E_\chi=0.
\]

With

\[
F(\phi)=M_{\rm Pl}^2+\beta\phi^2>0
\]

and

\[
\widehat P_\alpha{}^{\beta\mu\nu}
=\delta_\alpha^{(\mu}\widehat g^{\nu)\beta}
-\frac12\delta_\alpha^\beta\widehat g^{\mu\nu},
\]

REF1 defines

\[
\boxed{
E_{\rm MHG}^{\mu\nu}
=E^{\mu\nu}
+F\widehat P_\alpha{}^{\beta\mu\nu}
\nabla_\beta C^\alpha=0.
}
\]

The extension is first assembled with upper metric-equation indices and then
lowered using the physical metric before it is added to VAR1's lower-index
residual.  The scalar equations are copied without modification.  No damping
term is selected in REF1.

The calculation retains the full four-dimensional connection and projector
at the regular equatorial point.  It verifies symmetry of both index forms and
equality of the `theta-theta` and `phi-phi` lower components there.

## Exact gauge-surface controls

The GR-0 and constant-coupling SGB-L fixtures use flat spherical geometry with
`R=r`.  Their physical and reference connections and first derivatives agree
exactly.  Consequently,

\[
C^\mu=0,
\qquad
\nabla_\beta C^\mu=0,
\qquad
E_{\rm MHG}^{\mu\nu}=E^{\mu\nu}
\]

at the complete declared two-jets.  The SGB-L fixture is an equation-form
control, not a claim that its nonzero constant scalar and potential solve all
unredefined equations.

The Schwarzschild Painleve--Gullstrand and activated FGC-QR fixtures are
deliberately off this reference gauge.  Their `C^t,C^r` and gauge extensions
are exact, finite, and nonzero, while `C^theta=C^phi=0` by spherical symmetry.
The Schwarzschild fixture separately retains its exact zero unredefined VAR1
residual.  It is not called a solution of the gauge-fixed equation because it
does not satisfy the selected gauge at that jet.

This separation matters: gauge fixing changes the off-gauge extension of the
PDE, not the physical solution set obtained after compatible gauge data and
homogeneous propagation are established.

## Independent principal regression

For every field

\[
u=(h_{tt},h_{tr},h_{rr},R,\phi,\chi)
\]

and every second-jet slot `(dtt,dtr,drr)`, REF1 perturbs the complete gauge
extension by exact rational `+1` and `-1` directions.  The central
difference is exact because the gauge extension is affine in those second
jets.  The resulting radial polynomial is assembled as

\[
P(-c,1)=A_{rr}-cA_{tr}+c^2A_{tt}.
\]

Coefficient by coefficient, all four fixture matrices equal the independently
implemented MHG1 principal projector block.  The reference connection, covariant
`Gamma*C` term, and derivatives of the auxiliary inverse metric contribute
only lower-order terms, exactly as required.

This regression is stronger than observing the same determinant: it compares
every entry of the complete six-by-six polynomial gauge block.  Adding that
block to RED1's unredefined principal matrix is the MHG1 full metric-symbol
construction; REF1 does not relabel the gauge block as the entire symbol.

## Why this is not yet a first-order system

The ACT1 metric source contains

\[
P_{\mu\rho\nu\sigma}
\nabla^\rho\nabla^\sigma f(\phi),
\]

while the regulator equation contains `f'(phi) G`.  Curvature and the scalar
Hessian both contain second derivatives.  Their products make the full
equations generally nonlinear in the second-jet variables even though their
linearization has the well-defined MHG1 principal matrix.

The next gate must therefore treat

\[
\mathcal R_A
(u,\partial u,\partial_t^2u,
 \partial_t\partial_ru,\partial_r^2u)=0
\]

as an implicit algebraic system for the six normal accelerations, identify a
declared exact solution jet, and prove that its acceleration Jacobian is
nonsingular.  Only a local implicit-function branch may then be called the
first-order evolution representation.  REF1 does not replace that proof by
multiplying with a frozen kinetic matrix.

## Remaining propagation and domain burden

On the scalar equations, diffeomorphism invariance should turn the divergence
of the gauge-fixed metric equation into a homogeneous linear equation for
`C`.  Its principal term is the MHG1 hat-metric wave operator.  The complete
operator also contains derivatives of `F`, the hat projector and auxiliary
metric, connection-curvature commutator terms, and the fixed reference data.
Those coefficients are now evaluated by PROP1 for independent spherical
`C^mu`; differentiating the metric-defined gauge vector and proving the full
constraint theorem remain separate gates.

After IMP1's one-point implicit branch and PROP1's conditional operator, HYP1
now also has [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md)'s exact local eighteen-variable implicit lift, complete
flat-root linearized branch map, and chosen kinematic radial reduction
identity. Its downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md)
gate certifies the nonlinear branch on a nonzero-width thirty-dimensional
state/jet box and its regular-chart margins. HYP1 still needs direct
metric-derived propagation equations plus all characteristic discriminants and
resultants, a declared EFT correction measure, and a uniformly bounded full
eigenbasis or symmetrizer. Finite fixtures and continuity language do not
supply that theorem.

The later [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result uses this
metric-defined `C` and full extension—not PROP1's independent probe—to close
one activated compatible local QIFT1 root. It does not differentiate the
metric-defined gauge constraint through a third jet or prove propagation.

## Fail-closed nonclaims

REF1 by itself does not establish:

- a regular implicit second-time-derivative branch;
- a complete lower-order first-order evolution system;
- complete lower-order gauge-constraint propagation;
- reduction-constraint or physical-initial-constraint propagation;
- a uniform open-domain symmetrizer or retained EFT region;
- a regular centre or boundary system;
- constraint-compatible initial data or any evolution;
- collapse, metric-null affine defocusing, singularity resolution, a child
  spacetime, a dark sector, varying local light speed, a particle spectrum, or
  a theory of everything.

The result is narrower: the project now has one explicit, tensorial, full
reference-gauge spherical residual whose gauge extension has the
already-certified MHG1 principal gauge block as its exact directional
derivative.
