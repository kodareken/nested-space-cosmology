# The unwrapped horizon domain and a Dirac source interface

The actual horizon benchmark has an unwrapped radial direction. Its
spacelike PG slice is R_rho times S2 times the declared compact interval.
The periodic axial cell used in v0.5.0 is a different global domain.
Consequently its flat axial Wilson-line minimum is not an automatic state
selection for this parent–child transition.

This calculation supplies a state-defined source in an admissible magnetic
sector of the existing U(1) candidate, and the equations that a corresponding
reduced metric solution would have to satisfy. The magnetic flux and the
common functional's coefficients remain unselected. No scalar stress field,
finite vacuum value or alternative metric is inserted to obtain a neck.

## Establish the domain before borrowing a state

For the specified R times S2 times I patch, the first homology is trivial.
There is no independent flat axial U(1) Wilson modulus. A locally finite
tree of simply connected rooms joined along simply connected sphere cuts
also creates no cycle unless additional identifications or return links are
specified. The recursive operator identity does not itself impose such an
identification. A path through one link and back through its adjoint has no
independent gauge-loop phase.

This statement does not remove the boundary self-energy or its amplitude.
It does not exclude curvature holonomy, non-flat gauge fields, fixed
endpoint phase references, other topology or nontrivial completion at
infinite depth. Nor does it erase the discrete compact-domain/parity choices.
The actual S2 factor has nontrivial second homology and can carry magnetic
flux even though its first homology is zero.

The existing causal Hamiltonian remains

\[
H_\kappa=-i(\sigma_2-\beta I)\partial_\rho
+\frac{i}{2}\beta'I+\frac{\kappa}{r}\sigma_1.
\]

Reuse its [tetrad and causal-domain result](nsc-lorentzian-transport.md).
Both characteristics inside the horizon travel toward decreasing rho.
The child truncation is outflow, not a reflecting endpoint or a periodic
return link. Norm-preserving coordinate-time evolution does not select a
regular horizon vacuum.

## A charged angular zero mode makes the transparent control physical

The existing unit-charge U(1) extension admits the standard monopole bundle

\[
F=\frac{q_{\rm mag}}2\sin\theta\,d\theta\wedge d\phi,\qquad
\frac1{2\pi}\int_{S^2}F=q_{\rm mag}\in\mathbb Z.
\]

This field is source-free on the stated spherical patch; its flux is
topological data, not a local inserted magnetic monopole. Its generation
or inheritance from the physical parent is not derived here.

Reuse [Maldacena–Milekhin–Popov, section2.3 and AppendixC](https://arxiv.org/html/1807.04726v3):
one four-dimensional massless Dirac field in nonzero flux gives |q_mag|
complex two-dimensional Dirac fields with zero angular eigenvalue.
The two compact bulk copies of the project supply one such 4D zero field;
they are not counted twice again. With no link mass and a fixed magnetic
background, these free radial channels have no angular barrier and unit
transmission. The previously formal kappa=0 transport control is therefore
physical in this magnetic sector. It is not a neutral sphere eigenmode.

Gauge dynamics, a link/fermion mass and higher modes can change the complete
response. In particular, angular index protection does not protect a
massless interacting 2D spectrum; the source discusses gauge-induced gaps.
No value of q_mag is chosen here, and no claim is made that the lowest
Landau level saturates the total radiation or vacuum stress.

## The state supplies flux without a return loop

Use the standard reduced conformal stress [Maldacena–Milekhin–Popov,
AppendixF](https://arxiv.org/html/1807.04726v3) and the fermionic
particle-plus-antiparticle normalization of
[Iso–Umetsu–Wilczek, equations25–28](https://arxiv.org/html/hep-th/0602146v2).
In signature +−, with ds2=A du dv and central charge c_eff=|q_mag|,

\[
T_{uu}=\frac{c_{\rm eff}}{192\pi}(2AA''-A'^2)+t_u,\qquad
T_{vv}=\frac{c_{\rm eff}}{192\pi}(2AA''-A'^2)+t_v,
\quad T_{uv}=\frac{c_{\rm eff}}{96\pi}AA''.
\]

Here u=t-r_star, v=t+r_star and dr_star/d rho=1/A. The full tensor is
conserved and has trace c_eff A''/(24 pi). Its spherical projection is
t_ab/(4 pi r²). This is the specified reduced-sector source, not the full
renormalized four-dimensional stress or a proof of finite-cutoff robustness.
The rest of the determinant, its matching terms and the complete state must
still be treated in the same functional.

The displayed tensor is the conformal/state part. An additional matched
local term proportional to g2 is not fixed by the anomaly; it is retained
in the source-potential interface below. Such a term has zero null
contraction and zero radial Killing flux. Higher-derivative terms are
outside this reduced approximation.

The projected free operator is normalized in the physical radial metric
g2. The inherited r-times-spinor identification removes the spherical
radial connection on the angular zero subspace. Its radius-independent
2D action has zero angular pressure in this restricted variation; that
statement is not extended to the complete four-dimensional determinant.
Other conformal-frame reductions require their corresponding cocycle.


Let P_H=c_eff kappa_h²/(48 pi). The standard stationary controls are:

| State | t_u | t_v | Interpretation |
|---|---|---|---|
| Boulware | 0 | 0 | Empty parent Killing modes; singular on the fixed horizon |
| Unruh | P_H | 0 | Regular outgoing horizon sector with no incoming parent particles |
| Hartle–Hawking | P_H | P_H | Equilibrium incoming bath as well as outgoing radiation |

These are reduced-state conditions. A Kerr Hadamard theorem is not a
theorem for this global black-universe continuation. The supplied
[horizon data](nsc-clock-horizon.md), a suitable extension and characteristic
state construction remain distinct requirements.

Transforming the tensor into the existing PG chart gives the exact current

\[
\boxed{4\pi r^2 T^\rho{}_\tau=t_u-t_v.}
\]

For the free massless magnetic channels in the Unruh branch, reuse the
recorded kappa_h=0.238325799634… rather than finding the horizon again.
It determines T_H=kappa_h/(2 pi) and

\[
\boxed{P_\infty=|q_{\rm mag}|\frac{\kappa_h^2}{48\pi}.}
\]

The standard source includes both particle and antiparticle occupation.
In ordinary units this is |q_mag| hbar c_light² kappa_h²/(48 pi), with
kappa_h in inverse length. It is parent Killing power. Killing time becomes
spacelike inside the horizon, so this is not automatically child proper
energy or the cosmological density rate Q. A static local counterterm has
no radial Killing-energy current; it cannot supply this state-dependent flux.

## Physical clock mapping, energy and radiation entropy

For a specified incoming conformal vacuum, let U(u) be the physical
ray-tracing relation to outgoing parent retarded time. Reuse
[Bianchi–Smerlak, equations6–9](https://arxiv.org/html/1404.0602v2):

\[
P(u)=-\frac{c_{\rm eff}}{24\pi}\{U,u\},\qquad
S_{\rm rad}(u)=-\frac{c_{\rm eff}}{12}\log U'(u).
\]

The entropy is a specified relative radiation entropy, not the entropy of
the whole parent or the necessary origin of existence. Its reference
normalization and the physical state matter. Constant peeling gives
P=c_eff kappa²/(48 pi) and S_rad'=c_eff kappa/12. These rates describe a
steady portion, not an eternal finite-energy evaporation or a complete
Page curve.

The Schwarzian composition law transports the state contribution under
successive actual ray maps. It is not a derived parent–child clock law.
Affine unit changes have zero Schwarzian. More generally, a passive
coordinate change does not create radiation: the metric polarization and
state terms must transform together. U(u) must come from physical
propagation and the stated vacuum, not be selected to fit an energy history.

## What the actual benchmark does with this source

The runner verifies conservation, the PG current and future-horizon
regularity. In the Unruh branch the horizon generator has a negative null
source, -P_H/(4 pi r_h²), while the outgoing part is regular.
The nonzero flux cannot source an exactly static metric unless another
flux cancels it or the geometry evolves.

At the original neck rho=0, however,

\[
A_0=1-\frac{3\pi}{2},\quad A'_0=6,\quad A''_0=-3\pi,\qquad
\frac{2A_0A''_0-A'^2_0}{192\pi}
=\frac{9\pi^2-6\pi-36}{192\pi}>0.
\]

The reduced Dirac source has positive contractions with both future radial
null directions there for all three listed controls. The imposed metric
requires negative contractions. Static magnetic and potential terms have
zero radial null contraction, so they cannot repair this sign mismatch.
The specified minimal sector therefore does not source this fixed neck.
This does not exclude higher-mode/nonlocal stress, the complete action,
another source-determined geometry or an evolving branch.

## The reduced equations define the next source-matching interface

In the static equilibrium branch set xi=G c_eff/(12 pi). Keep the source
potential (denoted U(r) in the implementation, distinct from the null ray
coordinate U(u)) equal to 8 pi G r² rho_pot(r) as an output to be obtained from the remaining
matched functional. For example, a constant-coupling potential has the
structure Lambda4 r²+rQ²/r²+u0. The finite projected vacuum term u0 is
not set to zero. U is not an inserted dark-fluid function and is not
reconstructed from a desired metric.

Reusing the spherical Einstein equations gives

\[
A''=\frac{-rU'-2A'rr'-\xi(A'^2-4\kappa^2)/(2A)}{r^2-\xi},
\qquad
r''=-\frac{\xi(2AA''-A'^2+4\kappa^2)}{4A^2r}.
\]
\[
{\cal C}=1-Ar'^2-A'rr'-U-\frac{\xi(A'^2-4\kappa^2)}{4A}=0.
\]

The runner checks that C'=0 under these equations and that the independent
angular and null-difference equations hold. Thus a satisfied constraint is
propagated; it is not replaced by an extremum along an imposed profile.
These are candidate two-derivative, spherical, free-LLL equations. Their
gravitational coefficients stand for a matched high-mode contribution,
with the explicit LLL counted once, not a second weighted gravitational
action. Higher derivative, interaction and domain contributions remain open.

At r'=0 define D=1-U0 and alpha=xi/r0². For nondegenerate data,

\[
A_0=\frac{\xi(A_1^2-4\kappa^2)}{4D},\qquad
A_2=\frac{-r_0U_1-2D}{r_0^2-\xi},\qquad
r_0r''_0=\frac{2D+\alpha r_0U_1}{2A_0(1-\alpha)}.
\]

For A0<0 and 0<alpha<1, a local expanding continuation requires
2D+alpha r0 U1<0. These are conditions on the independently computed
source and state. No values satisfying them have been chosen here.
The local ODE is regular away from A=0 and r²=xi; this does not prove
global horizon matching, stability, physical completeness or an
operational Planck resolution limit.

## Exact exterior reuse and the neutral branch

The parent exterior agrees with
[Heidari et al., equation6](https://arxiv.org/html/2608.23861v1) at G=M=a=1,
using atan(1/rho)=pi/2-atan(rho) for rho>0. The matter action is different,
and global continuation across rho=0 is not identified by that exterior
formula. Small-a/M expansions are not controlled at this parameter point.

The inspected direct trail includes Dirac QNMs and WKB transmission
estimates, but no exact neutral transmission/luminosity at this parameter
point was identified. The scalar partial waves are not Dirac results.
No priority claim follows from this bounded search. If the neutral branch
is pursued, the remaining calculation is its exact transmission factor,
not another Hawking or black-universe derivation.

For one complex 4D Dirac field the neutral flux sum has 8 kappa/(2 pi)
per positive kappa: the two angular signs, 2 kappa magnetic states and
particle plus antiparticle occupations. Angular signs are not
antiparticles. The magnetic |q_mag| zero-mode count is a different
reduction and must not multiply that neutral sum. The factor of two is checked against the two-dimensional particle-plus-antiparticle flux and the spin-averaged absorption convention; it is not supplied by the sign of kappa.

The [record](../results/development/horizon-source.json) authenticates the
prior geometry, charged-sector and source records. Reproduction checks
every field with 3e-13 absolute/relative float tolerance and exact
structure/source hashes:

    python scripts/check_nsc_horizon_source.py --check

The state tensor and reduced Einstein assembly were independently checked
through separate two-dimensional and four-dimensional symbolic calculations.
The checks include trace/conservation, PG current, the original neck
contractions, the classical RN limit and the source coefficient xi.

It reuses the recorded horizon and prior physics. No existing scientific
generator, collapse solution or greybody scan is rerun.
