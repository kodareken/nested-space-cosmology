# FGC-1-VAR1: covariant metric variation

## Status and exact boundary

**FGC-1-VAR1** derives the complete covariant bulk metric equation of the
action frozen by [FGC-1-ACT1](fgc-action-gate.md) and cross-checks the
Gauss--Bonnet sign and coefficient through four mutually constraining routes:

1. the compact divergence-free double-dual form;
2. a separately implemented seven-term curvature--Hessian expansion on exact
   generic algebraic-curvature fixtures;
3. the trace and algebraic curvature contraction required by the
   diffeomorphism/Noether balance; and
4. exact automatic differentiation with respect to the lapse in the spatially
   flat FLRW reduction.

After setting Thaalba et al.'s `f_T=alpha_T phi^2/2` and retaining their common
overall action normalization, the result maps exactly to their double-epsilon
metric source. This closes the **covariant bulk metric-variation gate** for the
declared conventions. It is not a tensor-CAS calculation or a second
independent functional variation, and it does not derive boundary or junction
terms.

Nothing in VAR1 supplies the spherical equations, constraints, principal
symbol, kinetic eigenvalues, nonlinear hyperbolicity, collapse, or affine
Raychaudhuri result. Those remain separate gates. In particular, “the field
equations contain no derivatives above second order” is not a health proof.

## Convention package

Use signature `(-+++)`, vary with respect to the inverse metric `g^{ab}`, and
define

\[
R^\rho{}_{\sigma\mu\nu}
=\partial_\mu\Gamma^\rho{}_{\nu\sigma}
-\partial_\nu\Gamma^\rho{}_{\mu\sigma}
+\Gamma^\rho{}_{\mu\lambda}\Gamma^\lambda{}_{\nu\sigma}
-\Gamma^\rho{}_{\nu\lambda}\Gamma^\lambda{}_{\mu\sigma},
\]

\[
R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu},
\qquad
\mathcal G
=R_{abcd}R^{abcd}-4R_{ab}R^{ab}+R^2.
\]

Antisymmetrization has weight one half. The action is

\[
S=\int d^4x\sqrt{-g}\left[
 \frac12F(\phi)R
 -\frac12(\nabla\phi)^2-V(\phi)
 -\frac12(\nabla\chi)^2
 +f(\phi)\mathcal G
\right],
\]

\[
F(\phi)=M_{\rm Pl}^2+\beta\phi^2,
\qquad
V(\phi)=\frac{\mu^2}{2}\phi^2+\frac{g_4}{4}\phi^4.
\]

The signs immediately before the four Lagrangian sectors are part of the
convention package. Importing a formula derived for `-f G/2`, varying
`g_ab` instead of `g^ab`, reversing the Riemann tensor, or silently changing
the order of the double-dual slots changes apparent signs and factors.

## Functional variation from the action

Define the standard four-dimensional double-dual/Lovelock tensor by its
explicit curvature expression:

\[
\boxed{
P_{abcd}
=R_{abcd}
-2g_{a[c}R_{d]b}
+2g_{b[c}R_{d]a}
+R g_{a[c}g_{d]b}.
}
\]

The expanded definition, rather than an epsilon formula, is the implementation
authority. It fixes orientation-independent component signs and gives

\[
\nabla^aP_{acbd}=0,
\qquad
g^{ab}P_{acbd}=-G^{\rm Ein}_{cd},
\]

where `G^Ein_ab` is the Einstein tensor, not the Gauss--Bonnet scalar.

Starting with

\[
\delta\sqrt{-g}=-\frac12\sqrt{-g}\,g_{ab}\delta g^{ab}
\]

and the Palatini curvature identity, varying the three curvature-squared
pieces of `f G` and integrating the derivatives of `delta g^{ab}` by parts
organizes the bulk coefficient as

\[
\delta\!\int d^4x\sqrt{-g}\,f\mathcal G
=\int d^4x\sqrt{-g}\left[
 f H^{(2)}_{ab}
 +4P_{acbd}\nabla^c\nabla^d f
\right]\delta g^{ab}
+\hbox{boundary},
\]

where

\[
H^{(2)}_{ab}
=2\left(
 RR_{ab}-2R_{ac}R^c{}_b-2R^{cd}R_{acbd}
 +R_a{}^{cde}R_{bcde}
\right)-\frac12g_{ab}\mathcal G.
\]

The four-dimensional Lanczos identity sets `H^(2)_ab=0` identically. With the
project's stress definition

\[
T^{(\mathrm{GB})}_{ab}
\equiv-\frac{2}{\sqrt{-g}}
\frac{\delta}{\delta g^{ab}}
\int d^4x\sqrt{-g}\,f(\phi)\mathcal G,
\]

the exact bulk contribution is therefore

\[
\boxed{
T^{(\mathrm{GB})}_{ab}
=-8P_{acbd}\nabla^c\nabla^d f.
}
\]

For constant `f`, its Hessian vanishes and so does the bulk source, as required
because the four-dimensional Gauss--Bonnet integral is topological up to
boundary data.

The nonminimal Ricci sector varies independently as

\[
\delta\!\int d^4x\sqrt{-g}\,\frac12FR
=\frac12\int d^4x\sqrt{-g}\left[
 FG_{ab}^{\rm Ein}
 +(g_{ab}\Box-\nabla_a\nabla_b)F
\right]\delta g^{ab}
+\hbox{boundary}.
\]

Consequently, with

\[
T^{(\phi)}_{ab}
=\nabla_a\phi\nabla_b\phi
-g_{ab}\left[\frac12(\nabla\phi)^2+V\right],
\]

\[
T^{(\chi)}_{ab}
=\nabla_a\chi\nabla_b\chi
-\frac12g_{ab}(\nabla\chi)^2,
\]

the complete covariant metric equation is

\[
\boxed{
F G^{\rm Ein}_{ab}
+(g_{ab}\Box-\nabla_a\nabla_b)F
=T^{(\phi)}_{ab}+T^{(\chi)}_{ab}
-8P_{acbd}\nabla^c\nabla^d f.
}
\]

The scalar equations remain those derived in ACT1:

\[
\Box\phi-V'(\phi)+\frac12F'(\phi)R+f'(\phi)\mathcal G=0,
\qquad
\Box\chi=0.
\]

## Seven-term expansion

Writing `H_ab=nabla_a nabla_b f` and `H=Box f`, the compact source expands to

\[
\boxed{
\begin{aligned}
T^{(\mathrm{GB})}_{ab}={}&
-8R_{acbd}H^{cd}
+8g_{ab}R_{cd}H^{cd}
-8R_{bc}H_a{}^c
-8R_{ac}H_b{}^c\\
&+8R_{ab}H
-4Rg_{ab}H
+4RH_{ab}.
\end{aligned}
}
\]

The executable expanded route assembles these seven basis tensors without
constructing `P`. The compact route constructs `P` and contracts it directly.
Their only shared curvature inputs are the metric and algebraic Riemann tensor;
neither implementation calls the other.

## Exact branch specializations

For a composite coupling,

\[
\nabla_a\nabla_b f(\phi)
=f'(\phi)\nabla_a\nabla_b\phi
+f''(\phi)\nabla_a\phi\nabla_b\phi.
\]

The frozen branches therefore give:

| Branch | Gauss--Bonnet metric source | Nonminimal derivative source |
|---|---|---|
| GR-0 | `0` | `0` |
| SGB-L | `-8 alpha P_acbd nabla^c nabla^d phi` | `0` |
| FGC-QR | `-eta P_acbd nabla^c nabla^d(phi^2)` | nonzero when the displayed `beta` source is nonzero |

For `F=M_Pl^2+beta phi^2`, the nonminimal derivative term has the separately
checked expansion

\[
\boxed{
\begin{aligned}
(g_{ab}\Box-\nabla_a\nabla_b)F
=2\beta\big[&g_{ab}\big(\phi\Box\phi+(\nabla\phi)^2\big)
-\phi\nabla_a\nabla_b\phi\\
&-\nabla_a\phi\nabla_b\phi\big].
\end{aligned}
}
\]

The factor-safe FGC-QR identity follows because
`f=eta phi^2/8`: multiplying its Hessian by `-8 P` leaves exactly
`-eta P nabla nabla(phi^2)`.

## Trace and Noether-contraction checks

Contracting the metric equation gives

\[
\boxed{
-FR+3\Box F
=-(\nabla\phi)^2-4V-(\nabla\chi)^2
+8G_{\rm Ein}^{cd}\nabla_c\nabla_d f.
}
\]

The exact fixtures verify both

\[
g^{ab}P_{acbd}=-G^{\rm Ein}_{cd}
\]

and the four-dimensional curvature identity

\[
P_{acbd}R^d{}_e{}^{ac}
=-\frac14\mathcal G\,g_{be}.
\]

Together with `nabla^a P_acbd=0`, the latter fixes the divergence:

\[
\boxed{
\nabla^aT^{(\mathrm{GB})}_{ab}
=\mathcal G\,\nabla_b f.
}
\]

Here `nabla^a P_acbd=0` follows analytically from the differential Bianchi
identity. The exact pointwise fixtures have curvature and scalar jets but no
derivative-of-curvature jet, so the executable check is the displayed
algebraic curvature contraction—not a numerical evaluation of `nabla P` or of
the full divergence.

The scalar stresses obey

\[
\nabla^aT^{(\phi)}_{ab}=(\Box\phi-V')\nabla_b\phi,
\qquad
\nabla^aT^{(\chi)}_{ab}=(\Box\chi)\nabla_b\chi.
\]

Taking the divergence of the metric equation and using the scalar equations
then gives `-R nabla_b F/2` on both sides. This simultaneously checks the
relative signs of `F R/2`, `f G`, and the scalar equation.

## Independent FLRW lapse check

For

\[
ds^2=-N(t)^2dt^2+a(t)^2d\mathbf x^2
\]

and homogeneous `f(t)`, the Gauss--Bonnet action per unit comoving volume can
be written

\[
N a^3 f\mathcal G
=\frac{d}{dt}\left(\frac{8f\dot a^3}{N^3}\right)
-\frac{8\dot f\dot a^3}{N^3}.
\]

The executable represents `N` as an exact value/tangent pair and differentiates
the same boundary-reduced expression by the product and integer-power rules.
It then applies
`T_00=-(N^2/a^3) partial L_reduced/partial N`. Setting `N=1` gives

\[
\boxed{T^{(\mathrm{GB})}_{00}=-24H^3\dot f.}
\]

The compact `P` contraction on a separately constructed FLRW Riemann/Hessian
fixture gives the same result, independent of `dot H` and `ddot f`. The full
homogeneous `00` equation is

\[
3FH^2+3H\dot F
=\frac12\dot\phi^2+V+\frac12\dot\chi^2-24H^3\dot f.
\]

This direct lapse route fixes the overall `-8` sign and coefficient without
using the compact-to-expanded identity.

## Executable certificate

The frozen [VAR1 configuration](../configs/fgc/fgc-1-metric-variation.toml)
uses a documented 64-bit linear congruential generator to build eight generic
algebraic Riemann tensors as sums of two Kulkarni--Nomizu products. Every
calculation uses `fractions.Fraction`; zero residuals are exact, not floating-
point tolerances.

The [implementation](../src/recursive_horizons/fgc/metric_variation.py) checks:

- both Riemann antisymmetries, pair exchange, and the first Bianchi identity;
- nonzero curvature, Gauss--Bonnet scalar, scalar jet, and source coverage;
- compact versus expanded Gauss--Bonnet sources component by component;
- the branch-specific composite Hessians for GR-0, SGB-L, and FGC-QR;
- compact versus expanded nonminimal `F` sources;
- stress symmetry, the `P` trace, full trace, and Noether curvature
  contraction identity;
- a constant-coupling topological zero; and
- the independent FLRW lapse sign and factor.

`python3 scripts/reproduce_fgc_metric_variation.py` writes the source-bound
[schema-1 result](../results/fgc-1-metric-variation.json). The record hashes
both FGC configurations, `action.py`, the VAR1 implementation, this derivation,
and the reproducer. A mutation test changes a single seven-term coefficient
and requires a nonzero exact residual, preventing an unexercised identity from
passing merely because both sources vanish.

VAR1 is classified as a covariant equation derivation and cross-check, not a
principal-symbol or collapse result. It explicitly records that no tensor CAS
or second independent functional variation was performed. The independent
checks are instead the primary-source epsilon mapping, exact double-dual
expansion, trace/Noether curvature-contraction identities, and automatic lapse
differentiation of the reduced FLRW action.

## Primary-source coefficient check

Equation (1) of Thaalba et al., [*Spherical collapse in scalar--Gauss--Bonnet
gravity: taming ill-posedness with a Ricci coupling*](https://arxiv.org/abs/2306.01695),
contains

\[
f_T(\phi)=\frac{\alpha_T}{2}\phi^2
\]

as the coefficient of `G`. Their metric Eq. (3) is expressed with two Levi--
Civita tensor densities. Expanding that double-epsilon term with its stated
symmetrization gives

\[
-4\alpha_T P_{acbd}\nabla^c\nabla^d(\phi^2)
=-8P_{acbd}\nabla^c\nabla^d f_T,
\]

which matches VAR1 exactly.

The printed compact normalization in Lara et al.,
[*Scalarization of isolated black holes in scalar Gauss--Bonnet theory in the
fixing-the-equations approach*](https://arxiv.org/abs/2403.08705), is not used
as a coefficient authority: its displayed action and its Eqs. (2), (6), and
(7) use a source normalization whose relation to the variational
`T^(GB)_ab` convention is not made explicit in those displayed equations. This
is recorded as a literature-convention warning rather than silently choosing
the formula that agrees with the project.

## HYP1 handoff

VAR1 exposes why HYP1 is unavoidable. `P_acbd nabla^c nabla^d f` mixes second
derivatives of the metric with second derivatives and gradient products of
`phi`; `(g_ab Box-nabla_a nabla_b)F` adds further metric--scalar principal
mixing. A spherical solver must therefore not be written by treating these
terms as prescribed lower-order sources.

[FGC-1-HYP1-RED1](fgc-hyp1-reduction.md) now completes the exact local
reduction preflight: it evaluates the unredefined equations in a general radial
ADM chart, cross-checks direct four-dimensional curvature against an
independent warped-product construction, and records the full uneliminated
second-jet coefficient matrix. That result makes the mixing executable but is
not the first-order principal symbol or a hyperbolicity certificate.

[FGC-1-HYP1-SYM1](fgc-hyp1-symbol.md) now performs the next local step: it
verifies the exact radial gauge/Bianchi quotient atlas, factorizes the scalar
Schur determinant into the canonical `chi` cone and a separately derived
quadratic `phi` factor, and constructs positive pointwise scalar companion
symmetrizers on the frozen fixtures. It also rejects one naive fixed-gauge
first-order shortcut in the GR control. That is a physical-symbol preflight,
not the complete nonlinear gauge/constraint system or an open-domain
strong-hyperbolicity result.

[FGC-1-HYP1-MHG1](fgc-hyp1-modified-harmonic.md) now performs the explicit
formulation step at principal order. It adds the modified-harmonic gauge term
directly to the unredefined VAR1 metric equations, retains both scalars,
factorizes the complete gauge-fixed determinant, and proves complete real
pointwise eigenbases for the standard radial first-order principal reduction
on the frozen fixtures. Its auxiliary cones govern gauge and
gauge-constraint error, not physical light. MHG1 is not yet a complete
lower-order evolution/constraint system or an open-domain strong-hyperbolicity
and retained-EFT theorem.

[FGC-1-HYP1-MHG2-REF1](fgc-hyp1-mhg-reference.md) now completes the narrower
full reference-gauge residual: it freezes a tensorial spherical
connection difference on a positive-radius annulus, evaluates
`F hat_P nabla C`, and recovers MHG1's principal gauge block entry by entry.
It is not yet a regular
implicit acceleration solve, propagation theorem, or retained-domain result.
The downstream [FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate now proves
one exact unquantified local coordinate-time acceleration branch at the flat
FGC-QR vacuum/reference-gauge root; it is not a complete first-order system or
an open-domain theorem. The subsequent
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate uses VAR1's exact
Noether sign to derive and cross-check the conditional lower-order gauge
operator for independent spherical `C^mu`; it does not directly differentiate
the complete metric-defined residual divergence or prove constraint
propagation.

The later [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) gate now performs that
exact local first-order recomposition, publishes the complete flat-root
linearized branch map, and derives the chosen kinematic radial reduction
identity. Its downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md)
certificate supplies a quantified full-dimensional local implicit branch box
with a unique enclosed acceleration root. The subsequent
[FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) certificate closes one
activated local metric-defined compatible root, but not a propagation system
or initial hypersurface. Direct metric-derived constraint propagation,
uniform nonflat hyperbolicity, retained EFT validity, and evolution remain
open.

The remaining admissible construction is the open-domain,
gauge/constraint-complete first-order unredefined spherical system:

1. retain REF1's frozen modified-harmonic reference-connection/gauge source;
2. use RED1's reduction of the complete covariant equations, not a
   gauge-fixed action;
3. use FO1-RC1's exact implicit lift, QIFT1's certified local nonlinear branch
   box, and COMP1's activated compatible local point as the base for a complete
   constraint-propagation and initial-data formulation;
4. construct the fully nonlinear linearized principal symbol before using equations of motion
   to hide metric second derivatives;
5. verify the complete metric-derived homogeneous gauge and reduction-constraint propagation identities;
6. require a real complete characteristic basis and nonsingular kinetic
   matrix throughout the intended domain; and
7. stop at loss of hyperbolicity or EFT control.

Until that sequence closes, `full_fgc1_health_gate_passed` and every collapse,
defocusing, singularity, and global-continuation claim remain false.
