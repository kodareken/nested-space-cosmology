# FGC-1: covariant action and health gate

## Status and scope

This document defines the first **candidate tournament** for Finite Gradient
Closure (FGC).  It is an action-level gate, not a collapse result.  The action
below has a curvature-sensitive regulator and a covariant low-gradient limit;
its covariant bulk metric variation is now closed separately by
[FGC-1-VAR1](fgc-metric-variation.md), but it has **not** yet been reduced to a
complete gauge/constraint-first-order spherical system or certified strongly
hyperbolic on a collapse background. The narrower
[FGC-1-HYP1-RED1](fgc-hyp1-reduction.md) preflight evaluates the exact
unredefined spherical equations and their uneliminated second-jet matrix; it
is followed by the pointwise physical-symbol
[FGC-1-HYP1-SYM1](fgc-hyp1-symbol.md) gate. SYM1 derives the local scalar
characteristic factors and pointwise symmetrizers. The downstream
[FGC-1-HYP1-MHG1](fgc-hyp1-modified-harmonic.md) gate supplies a direct
modified-harmonic radial principal system and pointwise complete real
eigenbases. [FGC-1-HYP1-MHG2-REF1](fgc-hyp1-mhg-reference.md) now supplies its
exact tensorial full reference-gauge residual. The downstream
[FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate closes one unquantified
local coordinate-time acceleration branch at an exact flat FGC-QR root. The
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate then derives the
conditional lower-order operator for independent spherical `C^mu`, but does
not close metric-derived constraint propagation or the open-domain HYP1
obligations. No `1+1` collapse
solver, bounce, singularity-avoidance, trapped-to-anti-trapped transition,
child spacetime, dark-sector, or observational claim is authorized from this
document alone.

The purpose of the tournament is narrower and falsifiable: determine whether
one explicitly stated, low-energy scalar--Gauss--Bonnet candidate can retain a
healthy predictive domain while providing an action-derived source that might
pass the local metric-null Raychaudhuri gate.  A failure rejects the candidate,
not Einstein gravity or the broader FGC question.

The required next gate is stated in [Finite Gradient Closure](finite-gradient-closure.md):
the complete source must overcome both expansion-squared and shear terms for a
finite affine interval, rather than merely produce a minimum of a coordinate
scale factor.

## Conventions and action

Use signature `(-+++)`, natural units `c=hbar=1`, and the curvature convention

\[
R^\rho{}_{\sigma\mu\nu}
=\partial_\mu\Gamma^\rho{}_{\nu\sigma}
-\partial_\nu\Gamma^\rho{}_{\mu\sigma}
+\Gamma^\rho{}_{\mu\lambda}\Gamma^\lambda{}_{\nu\sigma}
-\Gamma^\rho{}_{\nu\lambda}\Gamma^\lambda{}_{\mu\sigma}.
\]

Then

\[
\mathcal G
=R^2-4R_{\mu\nu}R^{\mu\nu}
 +R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}.
\]

The declared action is

\[
S=\int d^4x\,\sqrt{-g}\left[
 \frac{M_{\rm Pl}^2+\beta\phi^2}{2}R
 -\frac12\nabla_\mu\phi\nabla^\mu\phi-V(\phi)
 -\frac12\nabla_\mu\chi\nabla^\mu\chi
 +f(\phi)\mathcal G
\right],
\]

\[
V(\phi)=\frac{\mu^2}{2}\phi^2+\frac{g_4}{4}\phi^4,
\qquad
M_{\rm eff}^2(\phi)=M_{\rm Pl}^2+\beta\phi^2.
\]

`chi` is an ordinary collapsing matter control; `phi` is the proposed
regulator.  Keeping them distinct prevents a regulator seed from being
mistaken for ordinary matter collapse.  The three frozen branches are:

| Branch | Parameters | Purpose |
|---|---|---|
| GR-0 | `alpha=eta=beta=0`, with `phi=0` initial data | Einstein--`chi` collapse control.  The even potential preserves `phi=0`; a collapse implementation must reproduce the corresponding GR solution and its metric-null focusing. |
| SGB-L | `f(phi)=alpha phi`, `alpha!=0`, `beta=eta=0` | A linearly curvature-sourced regulator. |
| FGC-QR | `f(phi)=eta phi^2/8`, `eta>0`, signed `beta!=0`, `alpha=0` | The main quadratic/Ricci-coupled trigger candidate. |

The dimensions are `[phi]=[chi]=[mu]=mass`, `[g4]=[beta]=1`,
`[alpha]=mass^-1`, and `[eta]=mass^-2`.  A numerical scan should use a
declared pulse width `L0` and dimensionless coordinates
`alpha/L0`, `eta/L0^2`, `mu L0`, `g4`, and `beta`.

## Convention-locked field equations

Varying the stated action before imposing any symmetry gives

\[
\boxed{\quad
 \Box\phi-V'(\phi)+\beta\phi R+f'(\phi)\mathcal G=0,
 \qquad \Box\chi=0,
 \quad}
\]

where

\[
V'(\phi)=\mu^2\phi+g_4\phi^3,
\qquad
f'(\phi)=
\begin{cases}
\alpha,&\text{SGB-L},\\
\eta\phi/4,&\text{FGC-QR}.
\end{cases}
\]

Thus the FGC-QR linear scalar equation about `phi=0` is

\[
\Box\phi-m_{\rm lin}^2\phi+O(\phi^3)=0,
\qquad
m_{\rm lin}^2=\mu^2-\beta R-\frac{\eta}{4}\mathcal G.
\]

For the metric equation, define the ordinary scalar stresses

\[
T^{(\phi)}_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
-g_{\mu\nu}\left[\frac12(\nabla\phi)^2+V\right],
\]

\[
T^{(\chi)}_{\mu\nu}=\nabla_\mu\chi\nabla_\nu\chi
-\frac12g_{\mu\nu}(\nabla\chi)^2,
\]

and lock the Gauss--Bonnet stress definition to the action itself,

\[
T^{(\mathrm{GB})}_{\mu\nu}
\equiv-\frac{2}{\sqrt{-g}}
\frac{\delta}{\delta g^{\mu\nu}}
\int d^4x\sqrt{-g}\,f(\phi)\mathcal G.
\]

The exact tensor equation is then

\[
\boxed{\quad
M_{\rm eff}^2G_{\mu\nu}
+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)M_{\rm eff}^2
=T^{(\phi)}_{\mu\nu}+T^{(\chi)}_{\mu\nu}
 +T^{(\mathrm{GB})}_{\mu\nu}.
\quad}
\]

This variational definition is the sign authority for a future symbolic
reduction.  Equivalently, after fixing the double-dual curvature convention,
the Gauss--Bonnet contribution can be written as a tensor linear in
`P . nabla nabla f`, where `P` is the divergence-free double dual of the
Riemann tensor.  A code must generate that expression from the displayed
action and verify it against this functional definition; copying a component
formula with an unstated epsilon/sign convention is not permitted.

The equations contain at most second derivatives in this special scalar--
Gauss--Bonnet construction.  That fact is necessary, but it is **not** a proof
of a bounded Hamiltonian, strong hyperbolicity, causal characteristics, or
well-posed nonlinear evolution.

## Low-gradient GR limit

The exact GR-0 action is recovered at `alpha=eta=beta=0`.  For nonzero
couplings, the intended low-gradient regime requires both a weak effective
Planck-mass correction,

\[
\frac{|\beta|\phi^2}{M_{\rm Pl}^2}\ll1,
\]

and Gauss--Bonnet and nonminimal-curvature terms small compared with the
Einstein tensor on the solution under consideration.  A covariant operational
check is that the norm of

\[
(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)M_{\rm eff}^2
-T^{(\mathrm{GB})}_{\mu\nu}
\]

is subleading to `M_Pl^2 G_mn` away from components where that comparator
vanishes; the implementation must state its nonzero normalization for those
points.  It must also show that `phi`'s ordinary stress is negligible whenever
it claims recovery of a specified GR matter solution.  These are model and
solution-dependent inequalities, not a universal numerical threshold.

Every retained branch requires

\[
M_{\rm eff}^2\geq\epsilon_M M_{\rm Pl}^2>0
\]

throughout its stated domain.  Crossing this boundary is a hard failure, not
a phase transition.

## Invariant activation and Schwarzschild control

In a Ricci-flat Schwarzschild control with areal radius `r`,

\[
R=0,
\qquad
\mathcal G=K=\frac{12r_s^2}{r^6}.
\]

For positive `eta`, the FGC-QR linearized scalar becomes tachyonic when

\[
\frac{\eta}{4}\mathcal G>\mu^2,
\]

which gives the illustrative threshold

\[
r_*=
\left(\frac{3\eta r_s^2}{\mu^2}\right)^{1/6}
=3^{1/6}\left(r_s\sqrt{\eta}\,\mu^{-1}\right)^{1/3}.
\]

This is an invariant **linear-instability scale** for the stated control, not
a proved collapse transition, an EFT cutoff, or a universal inheritance law.
In nonvacuum collapse the `beta R` term, gradients, nonlinear potential, and
the metric equations all matter.  For SGB-L, \(\alpha\mathcal G\) is a source
rather than a tachyonic threshold, so it has no analogous sharp onset without
an additional, derived response criterion.

## Executable ACT1 certificate

**FGC-1-ACT1** makes only the preceding action/background algebra executable.
The frozen [TOML fixture](../configs/fgc/fgc-1-action-gate.toml) supplies all
three branches, the standard-library
[implementation](../src/recursive_horizons/fgc/action.py) evaluates the
coupling derivatives and activation control, and
`python3 scripts/reproduce_fgc_action.py` writes the versioned
[JSON certificate](../results/fgc-1-action-gate.json). The fixture values test
signs, dimensions, branches, and residuals; they are not a fit or a proposed
physical parameter point.

The certificate marks the declared action and scalar algebra complete while
retaining `full_fgc1_health_gate_passed=false` and
`publication_local_defocusing_gate_passed=false`. As a separately versioned
stage certificate, ACT1 also keeps its then-open expanded-metric-variation
field false.
[FGC-1-VAR1](fgc-metric-variation.md) is the later source-bound certificate that
closes that one dependency; it does not rewrite the historical ACT1 record.
Both records keep the spherical reduction, principal symbol, nonlinear ghost
freedom, constraints, collapse, singularity resolution, child spacetime,
dark-sector, varying-`c`, particle-spectrum, and theory-of-everything fields
false. Those flags are scientific boundaries, not placeholders to toggle.

## Degrees of freedom and health gates

On a branch with `M_eff^2>0`, the intended propagating content is the two
tensor polarizations plus the two scalar fields `phi` and `chi`.  The
Gauss--Bonnet density is topological only for constant `f`; the field-dependent
coupling changes the dynamics but does not by itself introduce the generic
massive spin-2 mode associated with arbitrary curvature-squared actions.

That degree-of-freedom count does not certify the candidate.  Before any
physical-collapse interpretation, each branch must pass all of the following:

1. **Covariant variation and reduction.** VAR1 derives and cross-checks the
   covariant bulk metric equation. RED1 evaluates that unredefined system in a
   declared general radial ADM spherical chart. SYM1 extracts its physical
   radial factors, and MHG1 supplies an explicit modified-harmonic principal
   first-order system. REF1 supplies the complete tensorial reference-gauge
   residual, IMP1 proves one local flat-root acceleration branch, and PROP1
   derives the conditional independent-`C^mu` lower-order gauge operator.
   [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) now additionally certifies
   a full-dimensional exact rational local REF1 implicit-acceleration box.
   HYP1 must still derive its reduced constraints and metric-derived
   Bianchi/gauge/reduction-constraint
   propagation identities and a nonflat hyperbolicity/validity domain. Gauge-fixing an action
   first is not a substitute
   for those missing equations.
2. **Weak-field kinetic gate.** Derive the quadratic action around the stated
   low-curvature branch and show positive physical kinetic eigenvalues.
3. **Nonlinear principal-symbol gate.** Choose a first-order reduction
   `A0(U) d_t U + Ar(U) d_r U = S(U)` and calculate
   `det[A0 xi_t + Ar xi_r]`.  Require invertible `A0`, real characteristics,
   and a complete eigenbasis or a symmetrizer on every retained background.
4. **Local source gate.** With affinely parametrized metric-null generators,
   calculate both sides of

   \[
   \frac{d\theta}{d\lambda}
   =-\frac12\theta^2-\sigma_{\mu\nu}\sigma^{\mu\nu}
    -R_{\mu\nu}k^\mu k^\nu.
   \]

   A candidate passes only if the complete, action-derived geometry has a
   resolved finite interval with the required sign **and magnitude**.  A
   coordinate minimum, a negative term in a rearranged effective stress, or a
   curvature threshold alone fails this gate.
5. **Finite/EFT gate.** Curvature, scalar gradients, stresses, and all active
   dimensionless expansion parameters must remain finite and inside an
   explicitly declared regime of validity.  Evolution must stop at loss of
   hyperbolicity or EFT control rather than be continued by numerical damping
   and relabeled predictive.
6. **Constraint and convergence gate.** A future spherical computation needs
   regular-centre conditions, convergent constraints, an affine-normalization
   check, finite Misner--Sharp mass/flux bookkeeping where asymptotic flatness
   is claimed, and resolution-independent sign margins.

For the radial null congruences used in a spherically symmetric calculation,
the null shear vanishes.  Such a calculation can test the zero-shear local
gate and radial stability only; it cannot establish the required robustness
to shear.  That requires at least a separately derived non-spherical
perturbation system and eventually a nonlinear non-spherical calculation.

## Explicit nonclaims and solver prohibition

ACT1 proves only that the displayed action is a defined candidate family and
that the scalar equations and invariant threshold control above follow from
it. VAR1 adds the covariant bulk metric equation and exact cross-checks. The
later [FGC-1-HYP1-FO1-RC1](fgc-hyp1-fo1-rc1.md) artifact now supplies an exact
local first-order implicit lift and complete flat-root linearization of the
REF1 residual. Its downstream
[FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) artifact certifies a unique
acceleration root throughout one declared thirty-dimensional local box, but
neither artifact supplies evolution. The downstream
[FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) artifact now composes one
activated interior parameter point with exact metric-defined gauge,
reduction, and physical normal-projection compatibility to obtain one local
unredefined-equation root. It is not a compatible initial hypersurface or a
propagation result. The combined records
still do **not** prove:

- strong hyperbolicity in the activated regime;
- absence of nonlinear ghosts or gradient instabilities;
- a finite-curvature collapse solution or a bounce;
- metric-null defocusing, horizon evolution, or a global continuation;
- an exterior matching, flux law, entropy law, or any cosmological result.

Accordingly, a collapse solver must not be represented as an FGC result until
the full nonlinear/spherical principal-symbol and constraint formulation has
been derived and independently checked.  If that derivation reveals an
elliptic region, singular kinetic matrix, failed constraint propagation, or
loss of EFT control before the local source gate, the affected branch is
rejected.

## Primary technical context

- A. Hegade K. R., J. L. Ripley, and N. Yunes, [*Where and why does
  Einstein--Scalar--Gauss--Bonnet theory break down?*](https://arxiv.org/abs/2211.08477)
  (first posted 2022; v3 2023): gauge-covariant hyperbolicity diagnostics and
  strong-field collapse
  failure regions in representative ESGB theories.
- F. Thaalba, M. Bezares, N. Franchini, and T. P. Sotiriou,
  [*Spherical collapse in scalar--Gauss--Bonnet gravity: taming ill-posedness
  with a Ricci coupling*](https://arxiv.org/abs/2306.01695) (2023): a concrete
  spherical-collapse study motivating the `beta phi^2 R` tournament branch;
  it is not evidence that this action passes the present gate.
- N. Franchini, M. Bezares, E. Barausse, and L. Lehner,
  [*Fixing the dynamical evolution in scalar--Gauss--Bonnet gravity*](https://arxiv.org/abs/2206.00014)
  (2022): an exploratory fixing procedure.  Such a modified evolution is not
  a certificate of well-posedness for the original action used here.
