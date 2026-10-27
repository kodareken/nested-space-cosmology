# A sourced nested geometry and its gravitational boundary response

This is the next quantitative joint of the [flow dictionary](nsc-atlas-flow-dictionary.md)
and [atlas programme](nsc-cosmological-atlas.md). The question is whether a
compatible nested interior changes a gravitational signal while its exterior
mass, source and observer are fixed.

## Preparation, scope and decision

Completed inputs: Israel/Brown--York matching, the horizon-fluid shear limit,
the existing operator-composition method, and the previously completed scalar
formation/observer results in their own candidate. Those scalar histories are
not rescaled or identified with this static geometry.

Three routes were compared: a full rotating cosmological simulation, a Rindler
response with a freely supplied reflecting load, and a sourced four-dimensional
nested configuration with its actual gravitational perturbation equations.
The first requires a rotating state and formation law not yet selected; the
second leaves the crucial reflection coefficient as input. The third is chosen.
It connects compatible existing mathematics to a new finite response calculation.

Borrow [Visser--Wiltshire](https://arxiv.org/html/gr-qc/0310107v2) for material
shell dynamics, [Pani et al.](https://arxiv.org/html/0909.0287v2) for axial
gravitational perturbations and their general metric junctions, and
[Jampolski--Rezzolla](https://arxiv.org/abs/2310.13946) as the explicit
nested-gravastar antecedent. No novelty is claimed for those constructions.

Executed route:

1. Use a canonical massless scalar on each timelike shell, whose homogeneous
   stress has p=sigma, and a four-form for the regional vacuum stresses.
   Construct two shells with radii 2.5M and .25M and intermediate Kottler
   mass .1M, solving the curvatures and shared clocks from both junctions.
2. Verify signed matching, conservation and radial potential curvature.
   The initially considered R=2.2M one-shell member has V''<0 and is not
   used. The analytic R=2.5M one-shell member instead has
   M^2 V''=(4 sqrt(105)-28)/25>0. The two-shell member must pass its own check.
3. Derive the common-time axial master field and its shell transfer. Prove
   flux conservation, composability and the positive-energy stability domain.
4. At r_obs=10M calculate the canonical l=2 response for M omega=.4 and2.0.
   Compare a regular child with an ingoing black-hole termination while
   retaining exactly the same geometry and outer shell outside r=.25M.
   Fix these choices before evaluating either target response; this is a
   conditional synthetic calculation, not an astronomical forecast.
5. Check ODE tolerance/start-radius convergence, an independent analytic cap
   solution, horizon flux and common-clock scaling. Compute the earliest
   causal return interval for a disturbance from the retained observer.

New physical inputs are one exterior scale M, the supplied compactness/radius
and intermediate-mass ratios, the common shell/bulk actions and the no-incoming
or regular-center state. Solved curvatures, lapse factors, shell scalar amplitudes
and flux jumps are derived within that family. No free frequency-dependent
reflector, initial waveform fit or observational normalization is introduced.

The allowance is 180 CPU seconds and 500,000 RHS calls for the finite calculation,
including convergence controls and independent checks, with up to four useful
workers. Store scalars and limited plotting points; no field grid or raw history.
Stop on a failed junction, negative radial curvature, inconsistent axial matching
or noncontracting numerical indicators. Do not respond by extending a grid scan.

Completion means an actual compatible finite background, its retarded
gravitational response, an error indicator and a causal time/physical scale.
It does not mean a vortex cascade, cosmological formation or absolute Lambda
selection. The likely failure modes are an incorrect clock/normal in the wave
matching, an unstable material shell, or a nominal signal hidden below numerical
error. Each has a separate preflight check. A later direct time-domain evolution
that violates the predicted response or first-return bound would falsify this
finite calculation.

## Completed background and material source

The finite calculation below is complete. Write G=c=M=1. One shared time t
is normalized at infinity and each region has

$$
ds_i^2=-A_i^2 f_i dt^2+\frac{dr^2}{f_i}+r^2d\Omega^2,
\qquad f_i=1-\frac{2m_i}{r}-H_i^2r^2.
\tag{G1}
$$

The exterior has (m,H,A)=(1,0,1); the middle region has m=.1, and the
central de Sitter region has m=0. Walls lie at R_o=2.5 and R_i=.25.
Each is timelike with ordinary increasing-radius normals.

A compatible action is Einstein gravity plus a four-form with Lagrangian
-F_abcd F^abcd/(2*4!), and on each wall

$$
S_s=-\frac12\int_{\Sigma_s}\sqrt{-h}\,h^{ab}\partial_a\vartheta_s
\partial_b\vartheta_s\,d^3y+e_s\int_{\Sigma_s} A_3.
\tag{G2}
$$

The homogeneous scalar has sigma=p=(d theta/d tau)^2/2. The charge term
has no metric stress. Constant four-form states F4=q_i volume4 give
Lambda_i=4piG q_i^2=3H_i^2, with zero bare vacuum here. Shell charges
are their oriented flux jumps. Thus the vacuum difference has an actual
source law; it is not an unbalanced discontinuity in a cosmological constant.
The scalar amplitudes and charges are state data selected for this equilibrium
family. No quantization, primordial selection or universal charge ratio is assumed.
The gravitational and fixed-flux boundary terms follow the
[existing four-form action](nsc-atlas-flux-transition.md#the-action-determines-the-source-and-the-inheritance-law).
Vary that action before putting the four-form on shell; its on-shell bulk
Lagrangian alone would give the wrong vacuum sign.

The two static Israel equations are

$$
\sigma=\frac{\sqrt{f_-}-\sqrt{f_+}}{4\pi R},\qquad
p=\frac1{8\pi}\left[\frac{f'}{2\sqrt f}+\frac{\sqrt f}{R}\right].
\tag{G3}
$$

With b=sqrt(f_+(R)), a=sqrt(f_-(R)), p=sigma reduces to

$$
4a^2-Ca-\left(1-\frac{3m_-}{R}\right)=0,\qquad
C=\frac{1-m_+/R-2H_+^2R^2}{b}+2b,
\tag{G4}
$$

$$
a=\frac{C+\sqrt{C^2+16(1-3m_-/R)}}8,\qquad
H_-^2=\frac{1-2m_-/R-a^2}{R^2},\qquad A_-=A_+\frac ba.
\tag{G5}
$$

This is an actual finite composition rule: the next outer mass, vacuum and
clock supply the boundary data for the next inner interface. Accept a step
only if its lapses and density are positive and its stability requirement holds.
The local law is shared; radii, masses and four-form sectors need not be equal.
An unlimited sequence would additionally need control of accumulated curvature,
proper distances and its innermost state. No infinite spacetime is claimed here.

For the outer wall, H_parent^2=3(29-sqrt(265))/1000 exactly. The results are

| Quantity | Parent annulus / outer wall | Central region / inner wall |
| --- | --- | --- |
| H^2 M^2 | .0381635382117009 | 4.41709564870669 |
| A | .541737767915889 | .283041476224201 |
| sigma M=p M | .0120417630513390 | .129330203316059 |
| M^2 V'' | .442766777225571 | 51.2125679208436 |

The exterior lapse is positive for r>=2.5. The Kottler f is concave, so
positive values at .25 and2.5 guarantee positivity throughout that annulus.
The central f decreases from1 to .723931521956, remaining positive. There is
no central singularity or horizon in the retained nested body; the interfaces
carry the declared distributional matter.

There is no normal-tangential bulk flux. The surface scalar equation gives
sigma R^4=constant. Holding the regional masses/vacua fixed, put k=4pi sigma,
w=kR and

$$
\dot R^2+V(R)=0,\qquad
V=f_--\left(\frac{f_--f_++w^2}{2w}\right)^2.
\tag{G6}
$$

Both walls satisfy V=V'=0 and V''>0. Linearizing the acceleration gives
delta R''+(V''/2)delta R=0 in each wall's proper time. A finite radial
perturbation also needs constraint-compatible energy data at quadratic order.
The pressure Israel equation and the normal contraction
S^a_b mean(K^b_a)=[T_nn] agree independently. Polar and nonlinear stability
are not implied by this radial test.

The conditional vacuum products are

$$
\Lambda_{\rm parent}M^2=.114490614635103,\qquad
\Lambda_{\rm child}(.1M)^2=.132512869461201.
\tag{G7}
$$

They follow after the radius/mass ratios and material law are supplied. They
are not an absolute observed Lambda prediction or a vacuum generated solely
by changing resolution.

## A physical gravitational operator and its interface

Consider linear axial gravitational perturbations with ell>=2 and time factor
exp(-i omega t). These are coupled metric perturbations in vacuum/four-form
regions, not a test scalar or a phenomenological spin-two potential. In the
Regge--Wheeler gauge define the common-time master field and its traction by

$$
\psi=\frac{A f h_1}{r},\quad \frac{dx}{dr}=\frac1{Af},\quad
p_x=\frac{d\psi}{dx},\qquad
\psi_{xx}+\left[\omega^2-U(r)\right]\psi=0,
\tag{G8}
$$

$$
U=A^2f\left(\frac{\ell(\ell+1)}{r^2}-\frac{6m}{r^3}\right),\qquad
h_0=\frac{i}{\omega}(Af\psi+r p_x).
\tag{G9}
$$

The local master Q=f h1/r alone jumps when A changes. Pani et al.'s general
metric junctions, their equation43 and Appendix A, are [h0]=0 and
[sqrt(f) h1]=0. In the common-time variable these become

$$
[\psi]=0,\qquad [p_x]=q_s\psi,\qquad
q_s=-\frac{[Af]}R=4\pi\sigma\sqrt{F_s}>0,\quad F_s=A^2f|_{\Sigma_s}.
\tag{G10}
$$

Thus the gravitational wave feels a derived interface potential, not a fitted
reflectivity. For nonzero frequency, the surface scalar's axial covariant
velocity perturbation is zero. Its stress has delta S_tA=p delta h_tA.
The remaining temporal Israel equation is independently checked through
K_tA=sqrt(f)(h0'+i omega h1)/2; it reduces to
[sqrt(f) h0']=(2/R)[sqrt(f)]h0 once G10 is used. It holds by the bulk Einstein
wave equation and [f']-6[m]/R^2-2[f]/R=0. The four-form has no independent
propagating axial degree of freedom. These checks keep the material source
in the perturbation calculation.

### Composition and stability proof

Let T transport (psi,p_x) through one region. The first-order system has zero
trace, so det T=1; each shell has matrix [[1,0],[q_s,1]], also determinant1.
They preserve the Wronskian flux J=Im(conjugate(psi) p_x). Multiplication of
these physical matrices gives the full layered response. Its Robin form is

$$
Z_{\rm out}=\frac{T_{21}+T_{22}Z_{\rm in}}{T_{11}+T_{12}Z_{\rm in}}.
\tag{G11}
$$

Away from T12=0, the two-endpoint action has kernel

$$
K_T=\frac1{T_{12}}
\begin{pmatrix}T_{11}&-1\\-1&T_{22}\end{pmatrix}.
\tag{G12}
$$

Add adjacent kernels and the interface term q_s at the common endpoint, then
eliminate that endpoint by a Schur complement. The result is exactly the
kernel of the multiplied transfer matrices. This supplies a gravitational
realization of the composition method in this axial sector. No fermion kernel
or earlier abstract operator is silently substituted for the gravitational one.
At a Dirichlet pole retain the transfer representation; the physics is not
singular merely because that boundary parametrization fails.

For this candidate U>=0 in every retained region and q_s>0. Hence

$$
\mathcal E[\psi]=\int\left(|\psi_x|^2+U|\psi|^2\right)dx
+\sum_s q_s|\psi(x_s)|^2\geq0.
\tag{G13}
$$

With the regular-center condition and radiative exterior this nonnegative
operator has no finite-energy exponentially growing axial mode. For an assumed
growing outgoing solution, the endpoints decay and integration by parts would
give omega^2 norm(psi)^2=E>=0, contradicting Im omega>0. The horizon-terminated
control has the same positivity property with an ingoing horizon. This proof
covers the declared axial sector; it is not a proof of polar stability or a
nonlinear turbulence cascade.

## A resolved same-exterior response

Keep the outer wall, parent annulus and exterior exactly fixed. Compare the
regular central region plus inner wall against the same parent Kottler metric
extended to its black-hole horizon r_h=.200306715100933M. Impose an ingoing,
future-regular state there. The removed inner wall is not retained as an
additional reflector in that control.

Drive the finite worldtube r_o=10M in the canonical ell=2 channel. Its incoming
and returned characteristic amplitudes obey

$$
a_{\rm in}=\frac12\left(\psi-\frac{p_x}{i\omega}\right),\qquad
a_{\rm out}=\frac12\left(\psi+\frac{p_x}{i\omega}\right),\qquad
\mathcal R=\frac{a_{\rm out}}{a_{\rm in}}.
\tag{G14}
$$

These are the finite-boundary canonical energy channels, not plane-wave
amplitudes extrapolated from r_o to astronomical infinity. Their difference
of powers equals J/omega exactly, even where U is nonzero. The source protocol
and observing clock are the same for both geometries.

| M omega | abs(R_cap-R_horizon) | Horizon absorbed fraction | Indicator for absorbed fraction |
| --- | --- | --- | --- |
| .4 | .00000160391406 | .000000288534441 | 1.96e-10 |
| 2 | .258116878344 | .442777762923 | 3.69e-9 |

The regular body's returned fraction is1 in this lossless linear sector.
At M omega=2 the horizon control returns .557222237068. Thus the same exterior
static mass does not fix its gravitational response. At the lower frequency
the interior is strongly screened; its contrast is still resolved above the
numerical indicators. This does not distinguish nesting uniquely from all
other possible interiors.

The indicators use three tolerance/start-radius levels, another integrator,
an independent hypergeometric central solution and flux balance. They are
numerical comparison indicators, not rigorous enclosures or physical confidence
intervals. The medium-to-fine complex-response differences are 4.89e-11 at
the lower frequency and at most1.46e-9 at the higher one, down from2.79e-9
and5.47e-8 respectively. The largest analytic-cap comparison gap is1.92e-14.
The actual background parameters, matter law, loss mechanisms and exterior
measurement protocol remain physical assumptions.

## Causality, units and the next joint

For initially quiet perturbations, a new disturbance launched inward from the
worldtube cannot distinguish the two terminations before the round trip to
the first changed interface:

$$
\frac{\Delta\tau_{\rm first}}M
=2\sqrt{f_o}\left(\int_{R_o}^{r_o}\frac{dr}{f_{\rm ext}}
+\int_{R_i}^{R_o}\frac{dr}{A_{\rm parent}f_{\rm parent}}\right)
=34.3186793341.
\tag{G15}
$$

This is a causal lower bound, not a fitted echo peak. A null round trip to the
regular center is36.0742605915M. Independent null ODE and quadrature calculations
agree within2e-10M. At the supplied example mass10^6 nominal solar masses,
the first bound is169.036344s and the central null round trip177.683444s.
The proper frequencies of the two drives are .0144505908Hz and .0722529542Hz.
The conversion uses the [IAU nominal solar mass parameter](https://arxiv.org/abs/1510.07674),
GM_sun=1.3271244e20m^3/s^2 and c=299792458m/s. This is an independently
supplied ruler, not an inferred astronomical mass.

The [compact record](../results/summaries/nsc-atlas-gravitational-response-v1.json)
contains scalar results and23 plotting frequencies; it stores no wave history.
Four workers used145106 RHS calls and1.625 solver CPU seconds, inside the
180-second/500000-call allowance. Thirteen new controls plus31 retained
flow/parent controls pass. Grok independently checked background and axial
junction formulas; the calculations and proof remain independently inspectable.

Next: a direct retarded pulse to check G15 and the composed frequency response,
then a specified dissipative/rotating material state if pursuing actual vortex
transfer. Polar perturbations, formation, cosmological attachment and physical
selection of the input ratios remain separate connections. This completed
finite result is not reopened to refine an already resolved point.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_gravitational_response.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_gravitational_response.py
```
