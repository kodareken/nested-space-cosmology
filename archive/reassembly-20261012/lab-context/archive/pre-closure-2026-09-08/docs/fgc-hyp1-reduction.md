# FGC-1-HYP1-RED1: unredefined spherical reduction specification

## Status and owner boundary

**HYP1-RED1** is the source-bound, exact pointwise spherical-reduction
preflight for the FGC-1-ACT1/VAR1 action. It evaluates the unredefined
covariant equations in a general spherical chart and extracts their
uneliminated second-jet coefficient matrix. It is **not** a selected 1+1
evolution/constraint system, first-order principal symbol, hyperbolicity
certificate, solver, collapse result, or Raychaudhuri/defocusing result.

The retained action and physical metric are exactly those of
[FGC-1-ACT1](fgc-action-gate.md) and [FGC-1-VAR1](fgc-metric-variation.md):

\[
S=\int d^4x\sqrt{-g}\left[
 \frac12F(\phi)R-\frac12(\nabla\phi)^2-V(\phi)
 -\frac12(\nabla\chi)^2+f(\phi)\mathcal G\right],
\]
\[
F=M_{\rm Pl}^2+\beta\phi^2,\qquad f=\frac{\eta}{8}\phi^2.
\]

Matter, photons, and any later DEF1 null congruence couple to the same physical
metric \(g_{ab}\). An effective metric from a principal polynomial diagnoses a
field block; it is never silently substituted for \(g_{ab}\) or its null cone.

HYP1 must start from the unreduced VAR1 equations
\[
F G^{\rm Ein}_{ab}+(g_{ab}\Box-\nabla_a\nabla_b)F
=T^{(\phi)}_{ab}+T^{(\chi)}_{ab}-8P_{acbd}\nabla^c\nabla^d f,
\]
\[
\Box\phi-V'(\phi)+\frac12F'(\phi)R+f'(\phi)\mathcal G=0,\qquad\Box\chi=0.
\]
No equation of motion may hide metric second derivatives before construction of
the complete principal block.

## What the executable certificate establishes

The standard-library implementation represents every spherical field by the
exact local two-jet

\[
(q,\partial_tq,\partial_xq,
  \partial_t^2q,\partial_t\partial_xq,\partial_x^2q)
\]

with rational arithmetic. It constructs the curvature by two separate routes:

1. a direct four-dimensional Christoffel/Riemann calculation at the regular
   angular point \(\theta=\pi/2\), including the angular derivatives of
   \(\sin^2\theta\); and
2. a separately assembled warped-product metric inverse, base and mixed
   connection, scalar Hessians, and curvature for \(h_{AB}+R^2d\Omega^2\),
   including
   \(R_{AiBj}=-R\nabla_A\nabla_BR\,\gamma_{ij}\) and the spherical
   \(R_{ijkl}\) block.

The two geometry routes agree exactly for every connection component, scalar
Hessian, Riemann and Ricci component, \(R\), and \(\mathcal G\); when each is
fed into the same frozen VAR1 algebraic equation assembly, all six independent
metric/scalar projections also agree on every rational fixture. This is an
independent geometry/reduction check, not a second independent variation of the
action. Symmetric exact directional differentiation then records the full
six-by-eighteen coefficient matrix of the independent equations with respect to

\[
\partial_{tt},\partial_{tx},\partial_{xx}
\quad\text{of}\quad
(\alpha,v,\Lambda,R,\phi,\chi).
\]

This differentiation is exact rather than a finite-difference approximation:
the unredefined equations are at most quadratic in the second jets, so the
centered rational first tangent equals the analytic derivative. The record
contains the complete matrix and a digest for each fixture. It also verifies
that the independent \(\chi\) principal factor is precisely the physical
metric-null quadratic, that the activated FGC-QR fixture exercises principal
mixing in both metric-to-\(\phi\) directions, and that exact Schwarzschild PG
data with \(r_s=2\), \(R=8\) give

\[
R=0,\qquad \mathcal G=\frac{3}{16384},\qquad
E_{ab}=E_\phi=E_\chi=0.
\]
The record serializes the exact maximum absolute metric residual and both
absolute scalar residuals as zero rather than exposing only pass/fail labels.

An independent analytic spatially flat FLRW control fixes
\(a=2\), \(H=3/2\), and \(\dot H=-2/3\) at radial coordinate \(x=5\).
The spherical implementation returns exactly
\[
R=6(\dot H+2H^2)=23,\qquad
\mathcal G=24H^2(\dot H+H^2)=\frac{171}{2}.
\]

Those are reduction and coefficient-extraction results. The raw matrix is
diffeomorphism-degenerate and is deliberately not called
\(\mathsf A^t\), \(\mathsf A^x\), a kinetic matrix, or a characteristic
polynomial.

## Primary RED1 chart: generalized spherical ADM/PG

The primary chart is
\[
ds^2=-\alpha(t,x)^2dt^2+
\Lambda(t,x)^2[dx+v(t,x)dt]^2+R(t,x)^2d\Omega^2.
\]

Here \(x\) is a general radial coordinate, \(R>0\) is the dynamical areal
radius, \(\alpha>0\) is the lapse, \(v\) is the radial shift, and
\(\Lambda>0\) is the dynamical radial spatial-metric field. RED1 does not yet
fix the radial gauge. The areal-coordinate condition \(R=x\) and the
flat-spatial-slice condition \(\Lambda=1\), used together in some
scalar--Gauss--Bonnet evolutions, are comparator specializations rather than
assumptions of the reduction. HYP1 may impose either only after deriving the
constraints and showing the chosen condition can be maintained on the
retained branch.

The physical-metric radial null coordinate speeds are
\[
\left.\frac{dx}{dt}\right|_{g\text{-null}}=-v\pm\frac{\alpha}{\Lambda}.
\]
They are geometry/GR regressions, not modified scalar or tensor
characteristics. This chart is horizon-penetrating-compatible when
\(\alpha,\Lambda,v,R\) and the coordinate map stay finite and regular across
the marginal surface; the ansatz alone does not prove that a selected slicing
has this property. The flat-PG comparator has \(R=x\), \(\Lambda=1\), and
\(v=A\zeta\), giving \(A(\pm1-\zeta)\). Its \(\zeta=1\) horizon condition must
not be copied to the general chart; derive both physical null expansions from
\(g_{ab}\).

Retain independent first-order variables
\[
Q_\phi=\partial_x\phi,\quad
P_\phi=\alpha^{-1}(\partial_t\phi-v\partial_x\phi),\qquad
Q_\chi=\partial_x\chi,\quad
P_\chi=\alpha^{-1}(\partial_t\chi-v\partial_x\chi).
\]
Thus
\[
\partial_t\phi=\alpha P_\phi+vQ_\phi,\qquad
\partial_tQ_\phi=\partial_x(\alpha P_\phi+vQ_\phi),
\]
with the identical pair for \(\chi\). All remaining scalar, metric, and radial
constraint equations must be derived from the covariant equations above.

## Principal block HYP1 must compute

Let \(U\) contain these scalar variables and every metric/connection variable
of the first-order reduction. HYP1 must publish
\[
\mathsf A^t(U)\partial_tU+\mathsf A^x(U)\partial_xU=\mathsf S(U),
\]
with radial constraints separately stated. Before eliminating a constraint
block, show the relevant kinetic/time and constraint matrices are invertible.
For \(\xi_\mu=(\xi_t,\xi_x,0,0)\),
\[
\det[\mathsf A^t\xi_t+\mathsf A^x\xi_x]=0,\qquad c=-\xi_t/\xi_x.
\]

A retained point needs a nonsingular kinetic block, real characteristic roots,
and a complete eigenbasis or explicit symmetrizer. A positive discriminant of
one reduced two-by-two block is only a necessary diagnostic for that block, not
a strong-hyperbolicity certificate for the coupled metric--\(\phi\)--\(\chi\)
system.

The \(\chi\) sector is decoupled from curvature couplings in the action but not
from the full principal problem: its stress changes metric constraints and can
enter reduced \(\phi\)/metric characteristics after constraint elimination.
HYP1 must report both uneliminated and reduced blocks.

An independent covariant 2+2 route is mandatory:
\[
g_{ab}dx^adx^b=h_{AB}(x)dx^Adx^B+R(x)^2d\Omega^2,\qquad A,B\in\{0,1\}.
\]
Construct its orbit-space principal polynomial directly from ACT1. For a
single scalar block, a Lorentzian effective tensor is a necessary local
hyperbolicity diagnostic and a definite tensor signals ellipticity. In ACT1's
coupled system it remains only a block diagnostic. Agreement with the chart
polynomial is required, but neither replaces a full strong-hyperbolicity proof.

## Source normalization and formulation map

Thaalba et al. use
\[
S_T=\frac{1}{16\pi}\int\sqrt{-g}
\left[R+X-\left(\frac{\beta_T}{2}R-\alpha_T\mathcal G\right)
f_T(\phi)\right]d^4x.
\]
The raw coefficients of \(R\phi_T^2\) and \(\mathcal G\phi_T^2\) alone are
**not** an ACT1 parameter map because the source scalar has not yet been
canonically normalized relative to the ACT1 scalar. After an irrelevant
overall action rescaling and the canonical field map
\[
\phi_T=\frac{\sqrt{2}}{M_{\rm Pl}}\phi_{\rm FGC},
\]
matching the source action to
\[
\frac12F(\phi_{\rm FGC})R
-\frac12(\nabla\phi_{\rm FGC})^2
+\frac{\eta_{\rm FGC}}8\phi_{\rm FGC}^2\mathcal G
\]
gives
\[
\boxed{
\beta_{\rm FGC}=-\frac{\beta_T}{2},\qquad
\eta_{\rm FGC}=4\alpha_T.
}
\]
Equivalently, in source-regression units \(M_{\rm Pl}^2=2\), the scalar
variables coincide and the same map follows. Without this canonical
normalization, retain only the robust sign map: positive source \(\beta_T\)
corresponds to negative ACT1 \(\beta\), while positive source \(\alpha_T\)
corresponds to positive ACT1 \(\eta\). These maps are not parameter bounds;
ACT1 additionally has independent \(\chi\), a potential, and dynamic
\(\Lambda\).

The closest horizon-penetrating comparator is Thaalba et al.,
[*The dynamics of spherically symmetric black holes in scalar--Gauss--Bonnet
gravity with a Ricci coupling*](https://arxiv.org/abs/2409.11398), Sec. II
Eqs. (3) and (7)--(11), and Appendix A Eqs. (25)--(31). Its flat-PG line element
\[
ds^2=-A^2dt^2+[dr+A\zeta dt]^2+r^2d\Omega^2
\]
and variables \(Q=\partial_r\phi\), \(P=A^{-1}\partial_t\phi-\zeta Q\) are
comparators, not ACT1 equations.

Thaalba et al., [arXiv:2306.01695](https://arxiv.org/abs/2306.01695), Eq. (1)
and Eq. (7), provide the scalar-cloud/Ricci-coupling comparison: positive
published \(\beta_T\) mitigates naked elliptic regions for tested spherical
data. It uses polar-areal coordinates and is not the horizon-penetrating HYP1
chart authority.

Hegade, Ripley, and Yunes,
[*Where and why does Einstein--Scalar--Gauss--Bonnet theory break down?*](https://arxiv.org/abs/2211.08477),
Eq. (4), Eq. (11), and Appendix Eqs. (A14)--(A22), provide the independent
gauge-covariant spherical principal-symbol construction. Their Einstein
coefficient is constant and their matter is first-derivative minimally coupled.
Their final scalar effective principal tensor therefore **cannot be imported
unchanged**: \(F(\phi)R\) changes ACT1 metric--scalar principal mixing.

Thaalba et al.,
[*Hyperbolicity in scalar--Gauss--Bonnet gravity: a gauge invariant study for
spherical evolution*](https://arxiv.org/abs/2410.16264), supplies the essential
caveat: derivative-dependent field redefinitions, including the disformal
example studied there, can alter the initial-value formulation and apparent
equation character; this does not apply indiscriminately to every invertible
linear transformation. HYP1 must freeze unredefined ACT1 fields, equations,
constraint elimination, and physical metric before reporting a discriminant or
effective metric.

The downstream [FGC-1-HYP1-SYM1](fgc-hyp1-symbol.md) artifact now uses this
matrix without redefining the fields. It verifies the radial gauge and Bianchi
polynomial kernels, constructs two overlapping quotient-row patches, and
factorizes the pointwise scalar determinant into the metric-null `chi` factor
and a separately derived quadratic `phi` factor. That advances the physical
principal-symbol dependency but does not retroactively turn RED1 into a
first-order system or close SYM1's remaining constraint, open-domain, and
evolution nonclaims.

The next [FGC-1-HYP1-MHG1](fgc-hyp1-modified-harmonic.md) artifact now adds
the modified-harmonic gauge-fixing principal term directly to this same RED1
matrix. It preserves SYM1's physical factors, separates the two auxiliary
gauge sectors, and proves complete real pointwise eigenbases for the frozen
twelve-variable radial first-order matrices. That advances RED1 into a
pointwise principal formulation without changing the physical equations on
the gauge surface. It does not make RED1 a complete lower-order evolution
system or close the open-domain, retained-EFT, constraint-propagation, or
collapse gates.

The later REF1 and IMP1 artifacts supply the full reference-gauge residual and
one local implicit acceleration branch. The subsequent
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) artifact derives the
conditional lower-order gauge operator on independent spherical `C^mu` and
regresses its principal columns to MHG1. None of those later results turns
RED1's raw matrix into a complete metric-derived gauge/reduction-constraint
propagation system.

The subsequent [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) artifact does turn
the complete REF1 residual into an exact local eighteen-variable first-order
implicit relation and publishes its complete flat-root linearization. Its
radial constraint result is the kinematic identity for the declared coordinate
rows only. Its downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md)
certificate closes the local nonlinear implicit-branch box, but metric-gauge
and physical constraint propagation, nonflat hyperbolicity, retained validity,
and evolution remain open. The later
[FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result supplies one activated
metric-defined compatible local root only; it does not promote RED1 into an
initial-value or propagation system.

## Fail-closed RED1 boundary

RED1 authorizes neither numerical evolution nor use of excision,
fixing-the-equations drivers, artificial viscosity, or post-failure
continuation as evidence for the original FGC-1 action. It does not establish:

- explicit remaining evolution/constraint equations;
- invertibility of kinetic or constraint blocks;
- real complete characteristics or a symmetrizer;
- EFT validity, constraint propagation, regular-centre data, or convergence;
- a trapped region, metric-null defocusing, or affine Raychaudhuri margin;
- a nonsingular core, bounce, horizon transition, or child domain.

If the unredefined system loses hyperbolicity, a kinetic block becomes
singular, a constraint fails, \(F\) crosses its positive lower bound, or EFT
control fails before the proposed interval, that branch is rejected. Excision
may be an exterior-observable technique in a separately scoped study; it
cannot close HYP1, DEF1, or any interior FGC claim.

## Primary records

- [Hegade, Ripley & Yunes, arXiv:2211.08477](https://arxiv.org/abs/2211.08477)
  — covariant spherical principal-symbol route; [source](https://arxiv.org/e-print/2211.08477).
- [Thaalba et al., arXiv:2306.01695](https://arxiv.org/abs/2306.01695)
  — scalar-cloud Ricci-coupling comparison; [source](https://arxiv.org/e-print/2306.01695).
- [Thaalba et al., arXiv:2410.16264](https://arxiv.org/abs/2410.16264)
  — gauge-invariant/field-redefinition caveat; [source](https://arxiv.org/e-print/2410.16264).
- [Thaalba et al., arXiv:2409.11398](https://arxiv.org/abs/2409.11398)
  — horizon-penetrating PG comparator and characteristic diagnostic; [source](https://arxiv.org/e-print/2409.11398).
