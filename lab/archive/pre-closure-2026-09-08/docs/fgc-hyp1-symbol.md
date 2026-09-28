# FGC-1-HYP1-SYM1: physical radial symbol and formulation gate

## Status and narrow result

**HYP1-SYM1** is the exact covariant-symbol step between the unredefined
spherical equations in HYP1-RED1 and a future gauge/constraint-complete
evolution system.  It does four things that RED1 deliberately did not:

1. verifies the two radial diffeomorphism right-null directions and the two
   principal Bianchi left-null identities as polynomial identities;
2. constructs a declared gauge-quotient chart and performs an exact metric
   Schur elimination at every candidate physical characteristic;
3. factorizes the remaining scalar characteristic determinant into the
   physical-metric null factor for ordinary collapsing matter `chi` and a
   separately derived quadratic factor for the regulator `phi`; and
4. translates the general ADM second jets into normal/radial first-order
   variables, proving the candidate Hamiltonian and momentum projections are
   free of normal accelerations and that the pointwise metric--`phi` kinetic
   block is nonsingular on the frozen fixtures.

For the activated FGC-QR fixture the `phi` quadratic is not proportional to
the matter metric-null quadratic.  Its exact discriminant is positive, and an
exact positive companion-block symmetrizer exists at that point.  This is the
first project artifact to derive a distinct pointwise real regulator
characteristic quadratic from the frozen action rather than assume one.

That pointwise result is **not** full strong hyperbolicity.  SYM1 does not yet
provide a single nonlinear first-order evolution system, a gauge driver, the
lower-order radial constraint solve, homogeneous constraint propagation, a
uniform open state domain, boundary/centre conditions, or an evolution.  It
therefore authorizes neither collapse nor affine-defocusing claims.

## Frozen equations and radial symbol

The action, physical metric, and equations remain exactly those of ACT1,
VAR1, and RED1.  No equation is fixed, excised, disformally redefined, or used
to hide a second derivative.  On the two-dimensional orbit space, write the
base-field perturbation order as

\[
 (\delta h_{tt},\delta h_{tr},\delta h_{rr},
   \delta R,\delta\phi,\delta\chi),
\]

and set the nonzero radial covector to

\[
 \xi_A=(-c,1).
\]

If `M[q.dtt]`, `M[q.dtr]`, and `M[q.drr]` are the exact RED1
second-jet coefficients, the polynomial symbol column for `q` is

\[
 \mathcal P_q(c)
 =M[q.dtt]c^2-M[q.dtr]c+M[q.drr].
\]

The sign in the mixed term follows directly from
`partial_t -> -c` and `partial_r -> 1`; no factor of two is inserted because
RED1 stores the coefficient of the full mixed derivative itself.

## Exact diffeomorphism quotient

At principal order a radial gauge vector has the two metric polarizations

\[
 g_{(t)}=(-2c,1,0,0,0,0)^T,
 \qquad
 g_{(r)}=(0,-c,2,0,0,0)^T.
\]

SYM1 requires the unredefined six-by-six polynomial symbol to satisfy

\[
 \mathcal P(c)g_{(t)}(c)=0,
 \qquad
 \mathcal P(c)g_{(r)}(c)=0
\]

coefficient by coefficient on every exact fixture.  It independently raises
the covector with the RED1 base inverse metric and verifies the principal
Bianchi identities

\[
 \xi^tE_{tt}+\xi^rE_{tr}=0,
 \qquad
 \xi^tE_{tr}+\xi^rE_{rr}=0
\]

as left-null polynomial identities.  The scalar rows are not discarded, and
the lower-order Noether identity involving background scalar gradients is not
misrepresented as a completed constraint-propagation proof.

For `xi_r=1`, `delta h_tt` is a nonsingular representative of the metric
right-gauge quotient: the removed-column gauge-generator minor has constant
determinant two.  SYM1 uses the quotient variable order

\[
 q_Q=(\delta h_{tt},\delta R,\delta\phi,\delta\chi).
\]

The primary Bianchi-row patch of this quotient atlas uses
`(E_rr,E_theta_theta,E_phi,E_chi)`.  This is a declared coordinate patch on
the quotient, not a claim that `h_tt` is a physical graviton.  Partition the
resulting four-by-four symbol
as

\[
 \mathcal P_Q=
 \begin{pmatrix}
   \mathcal P_{gg}&\mathcal P_{gs}\\
   \mathcal P_{sg}&\mathcal P_{ss}
 \end{pmatrix},
 \qquad g=(\delta h_{tt},\delta R),\quad s=(\delta\phi,\delta\chi).
\]

The exact quotient metric-block determinant has the form

\[
 \det\mathcal P_{gg}(c)=K_g(c+v)^2,
 \qquad K_g\ne0,
\]

where `v=h_tr/h_rr` is the ADM shift.  The factor is the Bianchi-row atlas
Jacobian: the retained `E_rr` row patch is singular at the normal-advection
speed `c=-v`.  It is not a universal kinetic determinant or a physical
characteristic.  SYM1 therefore constructs a second patch using
`(E_tt,E_theta_theta,E_phi,E_chi)`, verifies its metric block is nonzero at
`c=-v`, and requires the two patches to yield exactly the same scalar
determinant coefficient by coefficient.  On either valid patch it forms the
covector-level Schur
complement

\[
 \mathcal P_{\rm scalar}
 =\mathcal P_{ss}-\mathcal P_{sg}
  \mathcal P_{gg}^{-1}\mathcal P_{gs}
\]

and requires that at least one atlas block be invertible at every root retained
from the scalar determinant.  The executable record checks exactly that
neither the physical-metric factor nor the derived `phi` factor vanishes at
the primary patch boundary and that the alternate block covers that boundary.

## Matter cone, regulator cone, and exact pointwise symmetrizer

In the generalized radial ADM chart

\[
 ds^2=-\alpha^2dt^2+\Lambda^2(dr+vdt)^2+R^2d\Omega^2,
\]

the normalized physical-metric radial null polynomial is

\[
 N_g(c)=(c+v)^2-\frac{\alpha^2}{\Lambda^2}.
\]

Its roots are the physical coordinate speeds

\[
 c_{g,\pm}=-v\pm\frac{\alpha}{\Lambda}.
\]

Exact Schur factorization gives, up to a nonzero rational normalization,

\[
 \det\mathcal P_{\rm scalar}(c)=N_g(c)Z_\phi(c).
\]

The first factor belongs to the independent canonical matter field `chi`.
The second is derived from the coupled metric--`phi` block; it is not imported
from constant-Einstein-coefficient scalar--Gauss--Bonnet theory.  GR-0, the
constant-`f` SGB-L control, and the zero-regulator Schwarzschild control give
`Z_phi proportional to N_g`.  The activated FGC-QR fixture gives a distinct
quadratic

\[
 Z_\phi(c)=a c^2+b c+d,
 \qquad a\ne0,
 \qquad \Delta_\phi=b^2-4ad>0.
\]

SYM1 publishes the exact coefficients, discriminant, rational isolating
intervals for both roots, and their separation from `c=-v`.  No decimal root
is used as proof.

For either quadratic `p(c)=a c^2+b c+d` with positive discriminant, define

\[
 C_p=
 \begin{pmatrix}
 0&-d/a\\1&-b/a
 \end{pmatrix},
 \qquad
 H_p=
 \begin{pmatrix}
 1&-b/(2a)\\
 -b/(2a)&-d/a+b^2/(2a^2)
 \end{pmatrix}.
\]

Then `H_p C_p` is exactly symmetric and

\[
 H_{p,11}=1>0,
 \qquad
 \det H_p=\frac{b^2-4ad}{4a^2}>0.
\]

The direct sum of the `N_g` and `Z_phi` companion blocks therefore has an
exact positive pointwise scalar-quotient symmetrizer.  This is stronger than
printing real roots of a determinant.  It remains narrower than a symmetrizer
for one complete nonlinear gauge/constraint evolution system.

## ADM normal-variable preflight

SYM1 also translates RED1's coordinate second jets into the action-independent
normal variables

\[
 \partial_\perp=\partial_t-v\partial_r,
 \qquad
 A=K^r{}_r=-\frac1\alpha
 \left(\frac{\partial_\perp\Lambda}{\Lambda}-\partial_rv\right),
 \qquad
 B=K^\theta{}_\theta=-\frac{\partial_\perp R}{\alpha R},
\]

\[
 P_I=\alpha^{-1}\partial_\perp I,
 \qquad
 Q_I=\Lambda^{-1}\partial_r I,
 \qquad I\in\{\phi,\chi\},
\]

with `L=partial_r ln Lambda` and
`S=(partial_r R)/(Lambda R)`.  At principal order,

\[
 \begin{aligned}
 \Lambda_{rr}&=\Lambda L_r,&
 \Lambda_{tr}&=v\Lambda L_r-\alpha\Lambda A_r+\Lambda v_{rr},&
 \Lambda_{tt}&=v^2\Lambda L_r-2\alpha v\Lambda A_r
                 -\alpha\Lambda A_\perp
                 +\Lambda(v_{tr}+v v_{rr}),\\
 R_{rr}&=\Lambda R S_r,&
 R_{tr}&=v\Lambda R S_r-\alpha R B_r,&
 R_{tt}&=v^2\Lambda R S_r-2\alpha vR B_r-\alpha R B_\perp,\\
 I_{rr}&=\Lambda Q_{I,r},&
 I_{tr}&=v\Lambda Q_{I,r}+\alpha P_{I,r},&
 I_{tt}&=v^2\Lambda Q_{I,r}+2\alpha vP_{I,r}
            +\alpha P_{I,\perp}.
 \end{aligned}
\]

Lapse and shift second derivatives are declared gauge-source terms in this
preflight; they are not silently relabeled as physical variables. The
executable Jacobian retains all six `alpha`/`v` second-jet source columns and
the displayed correlated `Lambda` terms required to hold `A` fixed. With
`n^a=alpha^-1(partial_t-v partial_r)^a` and
`s^a=Lambda^-1(partial_r)^a`, the four metric projections are

\[
 E_{nn}=\frac{E_{tt}-2vE_{tr}+v^2E_{rr}}{\alpha^2},
 \quad
 E_{ns}=\frac{E_{tr}-vE_{rr}}{\alpha\Lambda},
 \quad
 E_{ss}=\frac{E_{rr}}{\Lambda^2},
 \quad
 E_{\hat\theta\hat\theta}=\frac{E_{\theta\theta}}{R^2}.
\]

The executable matrix substitution proves that `E_nn` and `E_ns` contain no
`A_perp`, `B_perp`, `P_phi_perp`, or `P_chi_perp` coefficient, while the
same rows also contain no prescribed-gauge second-jet coefficient. The
four-by-four block of `(E_ss,E_theta_theta,E_phi,E_chi)` against those normal
accelerations has nonzero determinant on every frozen fixture.  This
establishes the **principal classification** of the candidate constraints and
the pointwise solvability of the evolution accelerations.  The complete
lower-order constraints and their propagation remain to be derived.

## A rejected shortcut

The same exact data expose a formulation failure that must not be hidden.  If
one merely prescribes lapse and shift, retains dynamic `R`, appends the four
kinematic compatibility equations, and treats the four projected equations
as an unconstrained eight-variable first-order system, the GR-0 control has a
zero-speed generalized eigenspace of dimension four but an ordinary
zero-speed eigenspace of dimension one.  That comparator is defective and
therefore cannot serve as the HYP1
strongly hyperbolic evolution system.

This rejection does not diagnose the physical scalar quotient as elliptic.
It says that gauge and radial constraints must be incorporated correctly.
The downstream [FGC-1-HYP1-MHG1](fgc-hyp1-modified-harmonic.md) gate selects
the modified-harmonic route and constructs the complete radial principal
system. Flat Painleve--Gullstrand remains a conditional comparator because
`R=r` and `Lambda=1` must both be preserved by any constrained alternative.

## Relation to primary literature

The covariant metric/scalar principal-block method and the need to quotient
gauge polarizations follow the gauge-invariant spherical analysis of Thaalba
et al., [arXiv:2410.16264](https://arxiv.org/abs/2410.16264), especially its
principal blocks and spherical physical-mode reduction.  That work also
includes an `h(phi) R` Ricci coupling, but its published one-scalar,
no-independent-`chi` reduction is not substituted for ACT1: the present
quotient is rederived in ACT1's own normalization and field content.

The generalized first-order block and Schur-complement logic is compared with
the flat-PG formulation and Appendix-A hyperbolicity definitions in Thaalba
et al., [arXiv:2409.11398](https://arxiv.org/abs/2409.11398).  That source uses
a specialized flat slicing and different action normalization, so it is a
method and regression source rather than an ACT1 equation source.

Hegade, Ripley, and Yunes,
[arXiv:2211.08477](https://arxiv.org/abs/2211.08477), provide an independent
covariant 2+2 comparator and the warning that an effective scalar metric can
change character in strong fields.  Their constant Einstein coefficient means
their final effective tensor cannot be imported unchanged.  The polar-areal
Ricci-coupled comparison in
[arXiv:2306.01695](https://arxiv.org/abs/2306.01695) likewise motivates the
coupling choice but does not prove the present formulation.

## Fail-closed boundary and next proof

SYM1 establishes exact local polynomial identities and pointwise fixture
properties only.  The fixtures need not satisfy the nonlinear constraints or
field equations.  It does not establish:

- a gauge-preserved nonlinear first-order evolution system;
- the complete lower-order Hamiltonian and momentum constraints;
- a homogeneous constraint-propagation system from the Noether identity;
- regular-centre or trapped-region boundary conditions;
- a uniformly positive symmetrizer or bounded characteristic basis on an open
  retained EFT domain;
- nonlinear predictivity through a restoring/defocusing interval;
- a collapse solution, finite affine interval, or metric-null Raychaudhuri
  sign change.

MHG1 promotes this quotient into one explicit gauge-fixed first-order
principal system and closes pointwise diagonalizability on the declared
fixtures. The downstream [FGC-1-HYP1-MHG2-REF1](fgc-hyp1-mhg-reference.md)
gate now serializes the complete reference-gauge residual and recovers MHG1's
principal gauge block by exact differentiation of the gauge extension. The
subsequent [FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate proves one
unquantified local coordinate-time acceleration branch at an exact flat
FGC-QR vacuum/reference-gauge root. The subsequent
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate derives the
conditional independent-`C^mu` lower-order gauge operator and matches its
principal columns to MHG1. The subsequent
[FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) gate supplies the exact local
first-order implicit lift, complete flat-root branch derivative and linearized
matrices, and the chosen kinematic radial reduction identity. The remaining
downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) gate controls the
nonlinear implicit branch on one nonzero-width thirty-dimensional rational box.
The later [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) gate composes one
activated point with a distinct physical/regulator cone and exact local
metric-defined compatibility, but does not provide a uniform eigenbasis or
constraint propagation. The remaining gates must certify complete metric-derived gauge/physical
constraint propagation, then prove the kinetic, constraint, hyperbolicity, and
EFT inequalities uniformly on an appropriate retained domain. If that
construction fails, HYP1 fails for that formulation;
pointwise Lorentzian quadratics and eigenbases are not permission to evolve
past the failure.
