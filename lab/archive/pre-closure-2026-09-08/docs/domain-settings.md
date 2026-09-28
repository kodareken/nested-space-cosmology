# Domain Settings and Scalar-Thermodynamic Synthesis

## Status

**DSF-1**, **DSF-2**, **DSF-3**, **DEB-1**, and **DST-1** are an explicit
research framework and controlled benchmark directions for Recursive Horizons
version 0.11.0. They are not a derivation of a black-hole-to-child transition, a detection of varying
constants, a dark-matter or dark-energy model, or evidence that our universe
is inside a black hole.

The goal is narrower and more useful: turn the language of “settings,”
gradients, resolution, and balancing forces into objects that can be defined
covariantly, compared operationally, and rejected if they do not reproduce
known physics. The active project boundaries remain the [claim ledger](claim-ledger.md),
[epistemic-status guide](epistemic-status.md), [controlled CCT-1 baseline](controlled-transition-model.md),
[information-map demonstrator](information-map.md), and
[perturbation/discriminator requirements](perturbations-and-discriminator.md). The
organizing law and its spectral specialization are
[Finite Gradient Closure](finite-gradient-closure.md) and the
[Nested Gradient Spectrum](nested-gradient-spectrum.md). NGS makes QFT modes,
EFT coarse-graining, bandwidth, thermodynamic occupation/dissipation, and
effective dark stress explicit H/Q calculations; it does not upgrade the
standard DSF/DST identities into an origin claim. The
two data-facing v0.6 checks are documented separately in the
[DESI DR2 BAO profile](desi-dr2-bao-likelihood.md) and
[GW170817 timing reconstruction](gw170817-relative-cone-data.md). The
[GMF-1B-PF1 preflight](global-matching-flux-closure.md#gmf-1b-pf1-two-preflight-exclusions-not-a-transition)
adds a separate matter-characteristic gate: an effective-fluid `dp/d rho`
failure is not evidence for varying local `c`, and a causal cone claim still
requires an action-derived principal symbol. Version 0.11's
[ECD source preflight](gmf-1b-ecd-symmetry.md) supplies one derivative-free
minimal-ECD contact source whose Dirac principal cone remains metric-null, and
the [interaction identity gate](gmf-1b-ecd-identity.md) translates that
interaction consistently between the canonical and numerical signature
conventions. Neither changes the principal cone or supplies a variable-`c`
result.

## DSF-1: domain settings are not one category of thing

### Operational thesis

> **DSF-1 — domain-settings framework.** A covariant state or order parameter
> `Phi` may select an effective causal domain through a derived metric or
> characteristic structure and a set of **dimensionless** effective settings.
> In the accessible domain, the construction must recover the common
> Lorentzian matter/photon cone and established effective physics.

This is a hypothesis and requirements statement, class **H/Q**. It becomes a
physical model only after it supplies an action, equations, a stable solution,
and a fixed observable. It is deliberately compatible with the project’s
existing RH-1 hypothesis that an order parameter may generate an effective
Lorentzian metric; DSF-1 specifies how such a claim must be articulated, not a
second mechanism.

A minimum mathematical target is

```text
S = integral d^4x sqrt(-g)
      [ M_Pl^2 R/2 + L(Phi, X) + L_matter(g, psi, ...) ],

X = -g^(mu nu) partial_mu Phi partial_nu Phi / 2,
G_a^(mu nu)(Phi, nabla Phi, state, ...) k_mu k_nu = 0.
```

The last equation is a characteristic equation for a mode `a`. It is not an
equation to insert by analogy: `G_a^(mu nu)` must be derived from the principal
part of the equations of motion. The domain’s settings must be specified as
operational, dimensionless quantities, for example

```text
s_i = { alpha_i, m_i/m_j, Lambda_i/Lambda_j,
        S_accessible/S_reference, transition ratios, ... }.
```

An alleged parent/child comparison additionally needs a physical transition
map and a common dimensionless convention. It cannot compare the bare
numerical value of a unit-dependent quantity across domains.

### A category table

| Object | What it is | What a DSF-1 claim may say | What it may not infer |
| --- | --- | --- | --- |
| `pi`, `e`, algebraic identities | Mathematical constants | Mathematics expresses reproducible relations compactly. | That a cosmological mechanism follows from mathematics being exact. |
| `c` in SI | A defined conversion factor; in Lorentzian physics, the local causal cone’s invariant speed | A pre-geometric phase may lack a Lorentzian cone, so local `c` can be **undefined** there. | That coordinate speed, a changed unit, or absence of photons makes locally measured `c` variable or infinite. |
| `l_P=sqrt(hbar G/c^3)` | A derived dimensional Planck combination | It marks the scale at which quantum-gravity effects are plausibly relevant. | That it is an experimentally established hard pixel or automatically creates a micro black hole. |
| `T=0` | A thermodynamic state/limit | It is a useful boundary for entropy and equilibrium questions. | That it means zero energy, no quantum fluctuations, or a proof that all cosmological dynamics must reverse. |
| `alpha`, mass ratios, entropy ratios | Dimensionless physical quantities | They can in principle be compared between regimes once measurement and dynamics are defined. | That they inherit, mutate, or vary without a distribution, action, or observable. |
| `alpha(mu)`, transport coefficients, sound speeds | Effective parameters/functions of state, scale, or medium | They may run or change across a derived phase. | That running with renormalization scale is automatically temporal variation of a fundamental law. |

This table preserves the distinction already required by
[epistemic status](epistemic-status.md): local invariants, coordinate
quantities, global causal statements, and cross-domain statements are not
interchangeable.

### Causal cones: when a speed claim is physical

In ordinary locally Lorentzian matter physics, a freely falling observer
measures nearby photon propagation on the local null cone. Coordinate speeds
can still vary. For a radial null curve in FLRW coordinates,

```text
d chi/dt = plus_or_minus c/a(t),
```

while the local causal speed remains `c`. Likewise, a Schwarzschild-coordinate
`dr/dt` is not a local measurement at a horizon. These boundaries are already
part of the repository’s [paper](../paper/recursive-horizons.md#3-causal-speed-photons-and-the-meaning-of-space).

The meaningful DSF-1 question is whether an action derives a *relative* cone.
If

```text
G_a^(mu nu) = Omega(Phi)^2 g^(mu nu),
```

then the null cone is conformally unchanged. A changed coefficient has not
created a physical variable speed of light. A genuine claim requires a
nonconformal characteristic structure and an observable such as

```text
delta_a(z, k, n) = c_a(z, k, n)/c_gamma(z, k, n) - 1,
```

with a stated reference sector, redshift/scale/direction dependence, stability
conditions, and a prediction fixed before the target data are examined.

For orientation, a noncanonical scalar can have a derived scalar sound speed,

```text
c_s^2 = P_X / (P_X + 2 X P_XX),
```

which illustrates the difference between a derived mode characteristic and a
verbal claim that `c` must vary. See Turner’s early analysis of coherent scalar
oscillations only for the oscillating-field result below, and Garriga and
Mukhanov’s [k-inflation perturbation analysis](https://doi.org/10.1016/S0370-2693(99)00602-4)
for this standard effective-sound-speed example. Any late-universe relative
gravitational/electromagnetic cone must also respect the multimessenger
GW170817 constraint, subject to its source-emission caveat:
[the collaboration’s primary report](https://arxiv.org/abs/1710.05834).

## DSF-2: a literal relative-cone gate

The first executable cone test is deliberately narrower than DSF-1. Define

```text
delta_T = c_T/c_gamma - 1.
```

It is a dimensionless comparison between tensor and photon characteristics,
not a statement that a locally measured photon speed or a dimensional `c` has
varied. A disformal effective metric can make distinct characteristics
mathematically possible, for example through a term proportional to
`nabla_mu Phi nabla_nu Phi`; a conformal rescaling cannot. The stationary
CCT-1 baseline has no such derived relative-cone effect and retains
`delta_T=0`. The covariant conformal-plus-gradient construction is described
in Bekenstein's primary [physical/gravitational geometry paper](https://doi.org/10.1103/PhysRevD.48.3641).

For the restricted subclass with constant `delta_T` over the accessible
late-time path, the GW170817/GRB 170817A timing analysis supplies the
conditional interval

```text
-3e-15 <= delta_T <= +7e-16,
```

subject to the source-emission-delay assumption in that analysis. Thus a
literal model that fixes a constant accessible `delta_T` outside that interval
is rejected by this gate; `delta_T=0` is the recovered common-cone control.
The gate neither constrains an inaccessible parent, an arbitrary scalar sound
speed, nor a redshift-, frequency-, or direction-dependent cone without a
separate propagation calculation.

Version 0.6 validates four compact GWOSC/Fermi/GCN source files and reconstructs
the published geocentric delay as `1.737774 s`, consistent with
`+1.74 +/- 0.05 s`. The rounded cone interval follows only after imposing the
paper's simultaneous-emission upper-bound assumption and ten-second source-lag
lower-bound assumption. That is a source-metadata and arithmetic reconstruction,
not a raw waveform/light-curve reanalysis or removal of the source-lag degeneracy.

## DSF-3: finite-domain causal calibration, not a bare variable `c`

The new nested-domain question is legitimate, but it needs one more category
boundary than the phrase “variable speed of light” usually carries:

> **DSF-3 — finite-domain causal calibration.** A finite transition may derive
> a different *relative* characteristic structure or clock/frequency transfer
> for its successor domain. The physical output must be dimensionless and
> operationally mapped. A different lapse, scale factor, coordinate speed, or
> numerical value assigned to dimensional `c` is not such an output.

DSF-3 is class **H/Q**. It currently has no separate numerical artifact. DSF-2
already supplies the smallest executable relative-cone control; a DSF-3 result
will be warranted only after an action and a transition map provide new inputs.

### A domain clock does not by itself vary local light speed

Write a domain-indexed homogeneous line element as

```text
ds_i^2 = -N_i(t)^2 c_ref^2 dt^2 + a_i(t)^2 dchi^2.
```

A radial null curve has coordinate slope

```text
dchi/dt = plus_or_minus N_i c_ref/a_i.
```

For a comoving local observer, however,

```text
dtau_i = N_i dt,
dell_i = a_i dchi,
dell_i/dtau_i = c_ref.
```

Thus an expanding or “thinning” geometry changes coordinate distances,
horizons, and travel times without automatically changing the locally measured
vacuum cone. For every smooth positive lapse, `tau_i=integral N_i dt` locally
removes `N_i`; it is clock parametrization until a physical clock sector says
otherwise. A finite proper-time endpoint, conformal boundary, geodesic
incompleteness, phase boundary, and finite handoff are also different global
claims. Each must be named rather than compressed into “time ends.”

Saying a created spacetime is “mapped end to end” can be a block-spacetime or
global-solution interpretation: one solution specifies the whole domain. It
does not imply that a creation event sends information instantaneously across
that solution or precomputes later events by a measurable superluminal signal.

### What a physical causal calibration would compare

At least two action-derived principal symbols are needed:

```text
G_a^(mu nu) k_mu k_nu = 0,
G_b^(mu nu) k_mu k_nu = 0.
```

In one operational frame the within-domain observable may be

```text
C_a|b(z,k,n) = c_a(z,k,n)/c_b(z,k,n),
delta_a|b = C_a|b - 1.
```

Across domains `i` and `j`, even this ratio requires a physical transition map
`M_(i->j)` that identifies observer worldlines, proper-time standards, phases,
sectors, energies, directions, and state normalization. Only then is a quantity
such as

```text
R_a|b^(i->j) = [(c_a/c_b)_j]/[(c_a/c_b)_i]
```

defined. Without that map, the causal settings are **uncompared**, not “faster,”
“slower,” or “infinite.” A corresponding spectral comparison must use the NGS
proper frequency `omega=-u^mu nabla_mu phi` and a dimensionless reference
ratio, not bare SI hertz.

This is where the finite-space intuition becomes a real programme: a
domain may begin and end in a declared proper-time or transition sense, and its
action may select its internal common cone. The scientifically new question is
whether the handoff derives `M_(i->j)` and a nontrivial `R_a|b^(i->j)` while
recovering the observed common cone inside our domain.

### VSL is a family of theories, not one density law

Albrecht and Magueijo explicitly studied a cosmological phase with faster early
light and modified evolution equations as an alternative route to the horizon,
flatness, and cosmological-constant problems. Their proposal does not derive
“more density or complexity means slower local light.” Magueijo’s review
separates hard Lorentz breaking, bimetric, locally Lorentz-invariant,
frequency-dependent, extra-dimensional, and vacuum-polarization VSL classes.
Those classes have different fields, symmetries, and observables; they cannot be
collapsed into a refractive-medium analogy. See the primary
[Albrecht--Magueijo model](https://doi.org/10.1103/PhysRevD.59.043516),
[VSL review](https://arxiv.org/abs/astro-ph/0305457), and distinct
[covariant construction](https://arxiv.org/abs/gr-qc/0007036).

For the fixed-coordinate metric `ds^2=-c^2dt^2+dx^2`, the formal
`c -> infinity` limit sends `g^(00) -> 0` and changes the causal structure to a
singular Galilean-type limit. It is not an ordinary finite Lorentzian domain
with instantaneous communication. FGC therefore gives a cleaner hypothesis:
a derived **large but finite dimensionless relative cone** may occur at a
coarser level, subject to hyperbolicity, stability, EFT validity, causal-loop,
local-recovery, and observational gates. A completed infinite speed is neither
needed nor admitted as a finite-domain result.

### QED vacuum effects are specific, and entanglement does not signal

The Lorentz-invariant Minkowski QED vacuum is not a universal material with a
refractive index determined by “density” or “complexity.” QED can alter photon
characteristics in specified backgrounds. Drummond and Hathrell derive
curvature- and polarization-dependent one-loop corrections in curved
spacetime; Barton derives a tiny direction-dependent parallel-plate result.
These are action-, background-, polarization-, and regime-dependent effects,
not a monotonic vacuum-thickness law or proof of usable faster-than-light
signalling. See [Drummond--Hathrell](https://doi.org/10.1103/PhysRevD.22.343)
and [Barton](https://doi.org/10.1016/0370-2693(90)91224-Y).

Quantum entanglement is likewise a correlation resource, not a controllable
instantaneous message channel. For a local trace-preserving operation
`E_A`, the remote reduced state obeys

```text
rho_B' = Tr_A[(E_A tensor I_B)(rho_AB)] = rho_B.
```

The parties need an ordinary causal channel to compare records. Entanglement
therefore cannot be used as evidence for the proposed coarse “top” domain.
Any modified top-layer dynamics must independently confront no-signalling,
Lorentz covariance or a declared preferred foliation, causal loops, and Bell
experiments. See the primary
[no-signalling analysis](https://doi.org/10.1103/PhysRevLett.87.170405).

### CCC is a neighbor, not an identity

Penrose’s Conformal Cyclic Cosmology joins the future conformal infinity of one
aeon to the conformally stretched beginning of the next. That makes CCC a
useful neighboring example of repeated domains and conformal handoff. It is not
the exact FGC/NGS nested-gradient proposal, and a finite conformal diagram does
not make every comoving observer’s proper time terminate at a finite value.
The crossover physics remains model-dependent; “time has been proved to stop”
is too strong. See the primary recent
[Meissner--Penrose formulation](https://arxiv.org/abs/2503.24263).

### DSF-3 failure gates

The branch fails or must be narrowed if:

- the effect disappears under `dtau=N dt`, a coordinate transformation, a
  unit change, or a common conformal rescaling;
- no healthy covariant action derives distinct characteristic metrics;
- a principal symbol loses Lorentzian signature, hyperbolicity, stability, or
  EFT control;
- no transition map defines proper time, clock, phase, observer, sector, and
  state comparisons across domains;
- the only prediction is a bare dimensional `c_i/c_j`, “thinness,” or
  complexity label;
- an accessible constant tensor/photon cone violates the stated GW170817 gate
  after its path and source-lag assumptions are applied; or
- entanglement correlations are treated as a selectable superluminal signal.

## DEB-1: boundary scaling is a conditional effective-fluid test

An alleged external or boundary scale must first be made operational in the
child domain. Let

```text
ell(a) = L(a)/L_ref = (a/a_ref)^s,
rho_B(a)/rho_B,ref = ell(a)^(-q).
```

Here `ell` is dimensionless. It does not license a comparison of an absolute
length in an unspecified parent with a child-domain length. For the explicitly
restricted case in which this component is separately conserved (`Q=0`), its
background identity is

```text
rho_B proportional to a^(-q s),
w_B = -1 + q s/3,
acceleration iff q s < 2.
```

These are conditional effective-fluid identities, not a derivation of a
boundary, an external energy source, or dark energy. They make several literal
subclasses fail-capable:

| Fixed subclass | Consequence | What the gate says |
| --- | --- | --- |
| `q=0` | Constant density and `w_B=-1` | In the v0.6 BAO-only profile it has `Delta chi2=1.2300` from the free constant-`w` minimum; it remains background-compatible, but this is no evidence that its origin is external or a boundary. |
| `q=1, s=1` | Constant-tension surface energy per volume, `w_B=-2/3` | It accelerates, but the real DESI DR2 BAO-only likelihood gives `Delta chi2=11.3325`; that restricted literal subclass is disfavored/rejected within the stated background model. |
| `q=2, s=1` | Inverse-area density, `w_B=-1/3` | It is nonaccelerating and the same profile gives `Delta chi2=69.8429`; it fails as a separately conserved explanation of the late accelerated background. |

The version-0.5 summary checks came from Table V of the primary
[DESI DR2 BAO cosmology paper, arXiv v3](https://arxiv.org/abs/2503.14738v3).
Version 0.6 instead evaluates the collaboration's public 13-observable BAO
distance vector and full covariance directly. Its free flat constant-`w`
profile minimum is `w=-0.911899`, `Omega_m=0.297693`, and `chi2=9.0410` for
ten nominal degrees of freedom. The calculation remains BAO-only: it does not
combine CMB or supernova likelihoods, use DR2 full shape, determine `H0`
without a sound-horizon model, or identify the physical origin of the fitted
background. See the [full model and likelihood scope](desi-dr2-bao-likelihood.md).

If `Q != 0`, the displayed `w_B` and dilution law are inapplicable until a
covariant stress tensor, interaction, and flux/junction law are derived. A
fitted transfer term cannot be used after the fact to rescue a failed
separately conserved boundary subclass.

## DST-1: one scalar sector can contain distinct effective components

### Controlled decomposition

DST-1 is not a claim that the existing CCT-1 stationary scalar already explains
the dark sector. It is a concrete next matter-sector ansatz to calculate. Let

```text
V(Phi) = V0 + m^2 varphi^2/2,
Phi = Phi0 + varphi,
```

in a homogeneous FLRW regime with a canonical scalar. Then

```text
rho_varphi = varphi_dot^2/2 + m^2 varphi^2/2,
p_varphi   = varphi_dot^2/2 - m^2 varphi^2/2,
rho_total  = rho_varphi + V0,
p_total    = p_varphi - V0.
```

For a nonzero constant term, `V0>0`, the component has

```text
p_V0 = -rho_V0 = -V0,       w_V0 = -1,
```

at the classical background level. At the allowed endpoint `V0=0`, the
component is absent and `p_V0/rho_V0` is undefined; the implementation returns
`None` rather than assigning it an equation of state. If the quadratic mode
oscillates rapidly,

```text
m >> H,
```

then its oscillation average satisfies

```text
<p_varphi> approximately 0,
<rho_varphi> proportional to a^(-3),
w_varphi approximately 0.
```

This is the standard coherent-oscillation result of M. S. Turner,
[“Coherent Scalar-Field Oscillations in an Expanding Universe”](https://doi.org/10.1103/PhysRevD.28.1243)
(*Physical Review D* 28, 1243–1247, 1983). It shows a background-level
possibility: one scalar potential can contain a vacuum-like term and a
matter-like oscillatory contribution.

That possibility is **not yet a unified dark-sector explanation**. The
background is degenerate: many distinct fluids, fields, and modified-gravity
models can produce similar `H(a)`. DST-1 must therefore demonstrate all of the
following before it can be called a dark-matter/dark-energy model:

1. **Dark-matter behavior:** a perturbation analysis showing viable clustering,
   lensing, sound speed/free-streaming behavior, halo phenomenology, abundance,
   and longevity—not merely `rho proportional to a^(-3)`.
2. **Dark-energy behavior:** the observed vacuum scale, radiative stability,
   technical naturalness or an alternative controlled explanation of its small
   value, and compatibility with late-time expansion constraints.
3. **Thermal history:** a reheating/production calculation, entropy transfer,
   and BBN/CMB consistency.
4. **Transition connection:** a derivation tying the scalar initial condition
   and parameters to a covariant collapse-to-child solution, rather than
   selecting them after observing the child cosmology.

### Dissipation is an interaction, not a name for balance

A phenomenological transfer term may be written

```text
Q = Gamma varphi_dot^2,

rho_varphi_dot + 3H(rho_varphi + p_varphi) = -Q,
rho_r_dot       + 4H rho_r                     = +Q.
```

Equivalently, in the simple homogeneous representation,

```text
varphi_ddot + (3H + Gamma) varphi_dot + m^2 varphi = 0.
```

This bookkeeping conserves the combined stress-energy at the background
level. It does **not** derive `Gamma`, particle production, thermalization, or
reheating. A viable model needs explicit couplings, e.g. to matter or radiation
fields, a quantum/nonequilibrium calculation, and constraints on every
interaction. Calling `Q` “energy balance” without that derivation hides the
physical work.

## Why black holes do not eat the universe

Black-hole gravity is not a rule that every ambient energy density must flow
inward. Growth requires a specified flux across a horizon. Schematically, for
a stationary-horizon energy current generated by `chi^mu`, the relevant input
is of the form

```text
dM/dlambda proportional to integral_H T_(mu nu) chi^mu dSigma^nu,
```

with the appropriate quasi-local/dynamical-horizon generalization in a real,
nonstationary spacetime. The sign and magnitude depend on the stress tensor,
trajectory, pressure, angular momentum, charge, boundary conditions, and the
global geometry—not merely on the statement that energy exists nearby.

In spherical symmetry the Misner–Sharp framework makes this dependence
explicit by relating quasi-local mass, fluid variables, work, and flux; it is
not a universal “everything falls in” rule. See Misner and Sharp,
[“Relativistic Equations for Adiabatic, Spherically Symmetric Gravitational Collapse”](https://doi.org/10.1103/PhysRev.136.B571).
Classical singularity results likewise have hypotheses about causal structure,
energy conditions, and trapped surfaces; they do not provide a cosmological
accretion law. See Penrose,
[“Gravitational Collapse and Space-Time Singularities”](https://doi.org/10.1103/PhysRevLett.14.57).

Finite black holes have finite mass, spin, and charge. Their interaction with
matter depends on capture cross sections and trajectories; pressure and angular
momentum can prevent direct radial infall, and cosmological expansion changes
the global question of which systems are causally connected. A cosmological
constant is described semiclassically by

```text
T_(mu nu)^Lambda = -rho_Lambda g_(mu nu),
```

not by an ordinary gas of particles streaming inward. It therefore does not
automatically supply ordinary particle accretion onto every black hole. Horizon
thermodynamics constrains accessible entropy, but it does not turn an ambient
vacuum energy into a universal inward mass flux; see the project’s existing
[horizon-thermodynamic discussion](../paper/recursive-horizons.md#5-horizon-thermodynamics-and-the-conditional-relation),
Bekenstein’s foundational [black-hole entropy paper](https://doi.org/10.1103/PhysRevD.7.2333),
and Bekenstein’s later [entropy-bound analysis](https://doi.org/10.1103/PhysRevD.49.1912).

There is nevertheless an exact GR anchor for the proposed language of
competing tendencies. In the weak-field limit of the spherical
Schwarzschild–de Sitter/Kottler solution with positive `Lambda`, radial
acceleration contains both terms

```text
a_r = -G M/r^2 + Lambda c^2 r/3.
```

They balance at

```text
r_balance = [3 G M/(Lambda c^2)]^(1/3).
```

Inside that radius the compact-mass term is stronger; outside it the
positive-`Lambda` term is stronger in this idealized geometry. Equivalently,
vacuum stress has `p=-rho`, so the FLRW acceleration source `rho+3p` is
negative. This is a real example in which attraction and accelerated
separation coexist and define a scale. The fixed-`M`, fixed-`Lambda` equality
is unstable because

```text
(d a_r/dr)_balance = Lambda c^2 > 0.
```

A small displacement is accelerated away from the balance radius. It is not a
stable equilibrium or a material wall, does not
prevent ordinary horizon-crossing accretion at small radius, and does not show
that the same balance creates spacetime or connects parent and child domains.
See Kottler’s original positive-`Lambda` spherical solution,
[doi:10.1002/andp.19183611402](https://doi.org/10.1002/andp.19183611402), and
the dynamical-horizon flux laws of Ashtekar and Krishnan,
[arXiv:gr-qc/0207080](https://arxiv.org/abs/gr-qc/0207080).

The correct question for this programme is not “why has every black hole not
eaten everything?” It is: **for a specified geometry and stress tensor, what
fluxes cross which horizons, and what global causal regions are connected?**

## Domain balance: the deeper idea in a form that can fail

The motivating intuition can be stated without pretending it is already a
theorem:

> **Domain-balance hypothesis (H/Q).** Stable causal domains may be dynamical
> equilibria or attractors of competing geometry, stress-energy, transport, and
> boundary terms. A domain transition may occur when that balance loses
> stability or crosses a threshold.

This branch uses the same **no unaccountable explanatory variable** rule as
NGS. A supposed balancing field, boundary pull, “dark” force, or relational
degree of freedom must enter a covariant equation, map to invariant
observables, change a predeclared prediction, and have a removal condition.
If the intended regime is a stable asymmetry, the equations must identify its
order parameter and show the relevant stability spectrum; if it is a driven
nonequilibrium regime, they must also close its finite work, flux, reservoir,
and entropy-production ledger. The words “balance” and “tension” do not supply
those missing quantities.

This is a promising bridge between the language of gradients and an actual
calculation. But it becomes physics only when the proposal supplies covariant
equations—such as Einstein/modified-gravity constraints, conservation laws,
matter equations, and explicit junction or flux conditions—and then solves
them. A diagram of forces “pulling from both sides” is not a stress tensor,
and a verbal equilibrium is not a junction condition.

GMF-1 makes that boundary concrete. A same-`Lambda` finite-`M` Kottler region
has a Misner–Sharp mass larger by `M` than the pure CCT-1 de Sitter core at the
same areal radius, so they cannot be directly Darmois matched. A shell may be
*prescribed* with Israel junction conditions, but its surface stress is then
distributional local data, not proof of a finite-invariant global transition,
topology change, stability, or an energy flux `Q`. The relevant audit is
[global-matching-flux-closure.md](global-matching-flux-closure.md), reproduced
by `python3 scripts/reproduce_gmf.py`; it does not identify an outside source
for dark energy or alter the local causal speed.

EC-1 sharpens rather than relaxes this gate. The action-motivated,
homogeneous Einstein--Cartan/Weyssenhoff reduction has a local finite turning
point, but its negative torsion-effective contribution scales as `a^-6` with
`w=+1`: it is neither a late-time vacuum component nor a derived outside
force. Its exact pressure-work boundary mass drift shows why a global
flux/taper/layer equation must be solved instead of assuming a fixed static
exterior. It derives no external `Q`, relative-cone signal, or variable locally
measured `c`; see [Popławski (2023)](https://arxiv.org/abs/2307.12190) and
[Ziaie et al. (2014)](https://doi.org/10.1140/epjc/s10052-014-3154-2).

For a spherical first pass, a completed model should state quantities such as

```text
nabla_mu T^(mu nu) = 0,
[extrinsic curvature and induced-metric matching data],
theta_+ and theta_-,
quasi-local mass and horizon flux,
delta S_gen or another declared stability functional.
```

It must determine whether a solution is stable, what perturbation triggers a
transition, whether energy conditions are modified or respected, and whether a
trapped region can connect to an expanding branch without a singularity or a
coordinate artifact. CCT-1 provides a deliberately limited closed-de Sitter
baseline for some of these questions; it does not yet supply the exterior,
junction, evolving order parameter, or flux law.

## Exact nonimplications

DSF-1/DSF-2/DSF-3/DEB-1/DST-1/domain balance do **not** imply any of the following:

```text
DSF-1  does not imply C10: a regular black-hole-to-child transition.
DSF-1  does not imply C13: a gravitationally derived information map.
DSF-1  does not imply C15: H-SAT or S_BH,parent = S_dS,child.
DST-1  does not imply observed dark matter, dark energy, or reheating.
NGS does not imply that unresolved gradient levels are dark matter or dark
energy, that every black hole adds one SI hertz, that friction universally
creates complexity or the golden ratio, or that unbounded continuation makes
observers logically necessary.
DSF-2 does not imply a variable local photon speed or a derived disformal action.
DSF-3 does not imply that every domain has a different cone, that space
expansion changes locally measured vacuum c, that a coarse domain has infinite
signal speed, that CCC is FGC, or that entanglement reveals a top layer.
DEB-1 does not imply an external origin merely because a background scaling fits.
Domain balance does not imply a literal Planck lattice or variable local c.
```

Conversely, a successful C10/C13/C15 construction could use a scalar sector or
domain-settings language, but it would still need its own derivation. The
finite IM-1 code map remains a mathematical demonstrator, and H-SAT remains
an optional coarse-entropy postulate rather than a consequence of unitarity or
scalar dynamics.

## Failure gates

DSF-1/DSF-2/DSF-3/DEB-1/DST-1 should be removed or narrowed if any of the following occurs:

- no generally covariant, hyperbolic, stable action derives the claimed
  effective setting/cone structure;
- the alleged effect is only a unit change, coordinate artifact, or conformal
  rescaling with no relative observable;
- a fixed accessible tensor/photon cone contradicts its declared
  multimessenger gate, or its source, redshift, frequency, and direction
  dependence has not been propagated;
- a claimed boundary scaling uses an unnormalized external length, assumes
  separate conservation without deriving it, or invokes `Q` without a covariant
  interaction and flux law;
- the recovered low-energy domain violates local Lorentz, multimessenger,
  laboratory, BBN, CMB, lensing, or structure-growth constraints;
- the scalar’s matter-like background fails clustering, lensing, abundance, or
  longevity tests;
- `V0` requires an uncontrolled radiative cancellation or fails expansion-data
  constraints;
- reheating or the radiation sector fails BBN/CMB entropy and thermal-history
  tests;
- no covariant flux/junction solution realizes the proposed domain balance; or
- the only agreement is obtained by fitting the same data used to select the
  settings.

The constructive standard is high but clear: derive the setting map, solve the
domain dynamics, recover known physics, and predict a dimensionless observable
that was not calibrated on the data used to test it.
