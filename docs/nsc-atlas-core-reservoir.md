# A core supplied by its inherited parent

The preceding [core calculation](nsc-atlas-scalar-core.md) selected a
coherent representative and admitted it into a sourced initial atlas. That
construction added geometric mass to the parent. This connection instead
places the core, a depleted ambient annulus and an outward-current scalar
pulse in the same constrained geometry, with zero mass offset outside them.
Both original material walls and the surrounding parent data are retained.

The finite solution reduces the ambient scalar displacement by at most
18.5914%. Its positive core excess is 9.57140113e-5, the reservoir contributes
-9.64548506e-5 and the pulse contributes 7.40839324e-7 in geometric mass units.
The absolute canonical energy remains positive. A protected surrounding
annulus retains the reference homogeneous evolution and inherited vacuum;
it occupies at least 86.6472% of the actual initial parent volume.

This is an exact initial-constraint construction, evaluated numerically,
plus a conditional causal-domain theorem and exact mass-transport identities.
It does not assert that the isolated cooling run produced this initial state.

## Decision, attribution and finite scope

Compare a complete two-wall evolution, an averaged cosmological envelope,
and a nonlinear compensated initial state with covariant flux. The third
addresses the missing source of the core's mass without importing the
vacuum-limit cooling history into the excited parent. It introduces no
new coupling or field. Core mass/phase, shell supports and pulse preparation
are inputs; one depletion amplitude is determined by exterior matching.

The spherical mass law is borrowed from
[Hayward, Section II, equations 2.2--2.9](https://arxiv.org/html/gr-qc/9710089).
Einstein--scalar constraint data use
[Choquet-Bruhat, Isenberg and Pollack](https://arxiv.org/abs/gr-qc/0610045).
The separate evolution input is geometric uniqueness and finite propagation
for the hyperbolic Einstein--scalar Cauchy problem; the maximal-development
framework is due to
[Choquet-Bruhat and Geroch](https://doi.org/10.1007/BF01645389).
The additional work is their application to the same inherited scalar,
the retained core representative and the compact two-wall atlas. No general
mass law or generic compensation construction is claimed as new.

Four independent cases (pulse, tightened integration, no pulse, zero source)
run in four processes within 120 CPU seconds. No old P001 grids or raw field
histories are needed. A momentum-sign error would fail the independent
constraints; insufficient ambient excitation would prevent the compensation
root; a formation claim would require evolution beyond the initial theorem.

## The constraint which includes momentum

Use a local orthonormal future normal n and outward spatial normal s, with
v=n(R), beta=s(R)>0, Pi=n(phi), F=dphi/dR and s(phi)=beta F. The initial
radial chart covers the small ball on the increasing-radius side of the
parent pole. K=-Lie_n(h)/2 fixes the sign of extrinsic curvature. Set

\[
\rho=\frac{\Pi^2+\beta^2F^2}{2}+U,
\qquad p_r=\frac{\Pi^2+\beta^2F^2}{2}-U.
\]

The momentum constraint gives the two extrinsic eigenvalues

\[
K_t=-\frac{v}{R},\qquad
K_r=-v_R-\frac{\kappa R}{2}\Pi F.
\]

Write the geometric Misner--Sharp mass (G times physical mass, c=1) as

\[
m=\frac{R}{2}(1+v^2-\beta^2)
 =M+\frac{\kappa U_0 R^3}{6},\qquad U_0=\frac{21}{32}.
\]

Substituting these expressions into the Hamilton constraint yields

\[
M_R=\frac{\kappa R^2}{2}
       (\rho-U_0-v\Pi F).
\]

The momentum term cannot be dropped when the pulse carries current. The
construction fixes v=H0 R at the original parent expansion, rather than
solving a new wall speed. In the zero-current regions K_r=K_t=-H0; in the
pulse region the radial eigenvalue changes by the explicitly sourced term.

## One field: core, reservoir and radiation preparation

Retain omega=2000, the same inherited potential and fluxes, and the prior
mass-selected coherent amplitude zeta=.08266217132541569. The original
omega=20/P003 candidate remains unchanged. The scalar displacement above
the parent potential minimum is chi0=.00043301270189221924, with excitation
E0=omega^2 chi0^2/2=3/8 and vacuum density rho_v=9/32.

The previous tapered core occupies 0<=R<=Rc, Rc=.120974320413543. Define
B(x)=1024 x^5(1-x)^5 for 0<x<1, zero otherwise. It is C4 and has unit peak.
The annular field perturbation is -a chi0 B on [Rc,2Rc]. An additional
pulse on [3Rc,4Rc] has

\[
f_{\rm pulse}=.02\chi_0 B
  \sin\bigl(\omega(R-3R_c)\bigr),\qquad
\Pi=-\beta F.
\]

Here B uses the normalized coordinate of its own annulus. Pi=0 elsewhere.
One instantaneous outward characteristic is selected, giving nonnegative
outward current j=-Pi beta F=beta^2 F^2. The massive potential remains in
the equations: this is not an assertion that the pulse stays purely outgoing
at later times. Its shape, support and amplitude are supplied preparation.

Integrate M(0)=0 and determine a from M(4Rc)=0. For the finite benchmark,
M(4Rc)=9.64546314e-5 at a=0 and -2.57169724e-4 at a=1, so continuity gives
a compensating solution. Near zero core/pulse amplitude, the derivative
with respect to depletion is strictly negative:

\[
\left.\frac{\partial M(R_b)}{\partial a}\right|_0
 =-\frac{\kappa\omega^2\chi_0^2}{2}
   \int_{R_c}^{2R_c}R^2 B(R)\,dR<0.
\]

The implicit-function theorem therefore supplies a local compensation
family around the reference. This small-amplitude theorem does not alone
prove the displayed finite case; its bracket and numerical checks do.
No universal uniqueness or dynamically preferred reservoir shape is asserted.

| Quantity | Computed value, atlas units |
| --- | ---: |
| Solved depletion fraction a | .185913687401 |
| Core geometric excess | .0000957140112793 |
| Reservoir geometric excess | -.0000964548506037 |
| Outward-current pulse excess | .000000740839324389 |
| Total exterior mass offset | below 2e-11 checked tolerance |
| Maximum local compactness excess 2M/R | .00687441 |
| Minimum sampled absolute density | .52977708 |
| Exact canonical vacuum floor | .28125 |
| Outer-wall speed and acceleration changes | zero within 2e-12 |
| Parent proper length change | .000244008403 |

The finite radial sample minimum is a diagnostic; the analytic bound
rho>=rho_v follows from the canonical potential and nonnegative kinetic
and gradient terms. Negative *excess* mass in the reservoir is not negative
absolute energy. Its nonnegative energy floor also bounds capacity:

\[
M(R_2)-M(R_c)\geq
-\frac{\kappa E_0}{6}(R_2^3-R_c^3)
\quad (\Pi=0).
\]

Ignoring the small additional pulse already requires R2>=.14890848 to
supply this core. The chosen R2=.24194864 satisfies that necessary bound;
the actual constraint solve establishes sufficiency for its stated profile.

Proper energy excess, about 2.79814e-6, does not vanish when M(Rb)=0.
Geometric mass includes the slice geometry and momentum contribution; it
is not an unweighted integral of local energy. Independent integration
of R and beta in proper distance verifies the mass and radius, without
reusing the mass ODE. Symbolic Hamilton/momentum reductions, the no-pulse
control, zero-source recovery and tightened tolerances check the connection.

## The same evolving exterior, with a perturbed centre

At Rb=4Rc=.483897281654, f and its relevant derivatives vanish, Pi=0 and
M=0. Outside it the original scalar, spatial metric and extrinsic curvature
are recovered exactly. Continue through the original parent equator and
back to the outer material wall. The middle and child regions, both initial
Israel conditions and positive lapse collars are unchanged. The central
proper length and volume can differ; matching does not equate those integrals.

Let a0=6.62607707610 and use the reference parent angle theta. Its matching
sphere is theta_b=asin(Rb/a0)=.0730942965661. The original wall lies at
theta_w=2.85076435977. The same homogeneous reference solution has

\[
\eta_\infty\leq\eta_*=
\arcsin\left(\frac{1}{a_0H_\Lambda}\right)
=.515418396612.
\]

On the open domain theta_b+eta(t)<theta<theta_w-eta(t), smooth
Einstein--scalar domain-of-dependence uniqueness identifies any compatible
causal development with that reference. Consequently the fixed annulus

\[
.588512693178<\theta<2.335345963158
\]

remains inside the reference domain for all its future proper time. This
does **not** bound propagation inside the altered central core. A putative
singularity there cannot send a signal faster than the boundary of this
domain. Existence of a complete global wall/core continuation is not inferred.
The statement concerns the globally hyperbolic domain of dependence; it
does not constrain non-globally-hyperbolic extensions beyond that domain.

Integrating 4 pi a0^3 sin^2(theta) and using the **actual** modified initial
parent volume gives a protected fraction above .866472. The same-source
homogeneous equations therefore supply a genuine evolving surrounding
domain for the compensated preparation, even while its interacting centre
and wall future remain unresolved.

Take L_ref=a0, fixed by the retained initial geometry. Throughout this
reference domain the inherited vacuum part remains

\[
\Lambda_{\rm vac}L_{\rm ref}^2
=\kappa\rho_v a_0^2=12.3482523989285.
\]

The total instantaneous stress still includes the excited scalar. Local
rearrangement of that excitation neither changes the action's vacuum floor
nor creates an additional cosmological-constant term by mass bookkeeping.
This connects local structure and the same parent expansion; it does not
select the observed Lambda, a physical length, or the cosmic matter abundance.

## Flux and work in one history

Hayward's mass law in this frame is

\[
n(m)=\frac{\kappa R^2}{2}
       (\Pi\beta^2F-vp_r).
\]

An expanding normal sphere has pressure work as well as flux. For a fixed
areal-radius trajectory D=n-(v/beta)s, the outward pulse obeys

\[
D(m)=-\frac{\kappa R^2}{2}\beta F^2(\beta+v)^2\leq0.
\]

Where beta^2>v^2, division by sqrt(1-v^2/beta^2) gives the unit observer's
proper-time derivative. At the sampled peak current R=.422618949392,
j=.000298238586331 and that proper-time mass rate is -.0000381559218903.
The normal-sphere rate is positive there because expansion work dominates;
these different clocks/spheres must not be mistaken for contradictory fluxes.

This identity fixes what the next evolution must conserve. A cooling core,
its reservoir and escaping field can exchange mass while the enclosing
history keeps the appropriate flux and pressure work. Initial compensation
does not impose M=0 on arbitrary moving spheres for all later times.

## Reproduce and continue

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_core_reservoir.py --output results/summaries/nsc-atlas-core-reservoir-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_core_reservoir.py
```

The [compact record](../results/summaries/nsc-atlas-core-reservoir-v1.json)
contains the parameters, checks, controls and source hashes. Full radial
states are regenerated, not committed. The next calculation is coupled
evolution of this compensated preparation, with its ambient source and
outgoing field retained, followed by a local observer signal. It must not
reuse the isolated cooling endpoint as an already derived global history.

Two independent Grok reviews confirmed the signed mass/constraint laws
and the protected-annulus argument. The second checked the implementation
and actual-volume normalization against the compact record and found no
material defect. The constraint construction and the evolution-uniqueness
input are attributed separately above.

The [coupled evolution](nsc-atlas-core-evolution.md) now advances this source
and a same-mass diffuse preparation with their gravitational response and
received clock signals. Its finite concentration does not imply a unique
inherited cloud population or complete asymptotic formation.
