# Nested Gradient Spectrum

## Status

**Nested Gradient Spectrum (NGS)** is a provisional specialization of
[Finite Gradient Closure](finite-gradient-closure.md). It translates the idea
of “gradients within gradients within gradients” into a possible hierarchy of
finite causal domains, spectral bands, and complexity budgets. It is an H/Q
hypothesis, not an identification of dark matter, dark energy, black holes, or
our universe's ancestry.

FGC is the closure law being sought. NGS asks what such closure could *build*:
whether concentration at one level opens additional modes at the next level,
while coarser enclosing levels supply the boundary conditions and finite
capacity that prevent unlimited localization.

## The core picture

The proposed hierarchy is

```text
coarse, low-resolution domain D_0
  -> localized finite transition
  -> richer finite domain D_1
     -> localized finite transition
     -> richer finite domain D_2
        -> ...
```

Every realized `D_i` is finite in the operational sense required by FGC. The
index can continue without a final value, but no completed domain contains an
actually realized infinite energy, density, curvature, bandwidth, entropy, or
state description.

“Space upon space upon space” does not yet mean ordinary boxes embedded in one
larger Euclidean space. The scientifically conservative possibilities include:

- causally connected regions of one global spacetime;
- effective domains related by a transition surface or trapping horizon;
- quantum sectors related by an isometry or channel; or
- genuinely distinct child geometries joined only through a derived global or
  quantum completion.

The geometry and transition map must decide which meaning is physical.

## From “one extra hertz” to a spectral prediction

A hertz is defined relative to a clock. It cannot by itself be inherited
across domains with different geometry or redshift. For an observer with
four-velocity `u^mu` and a mode phase `varphi`, the local proper frequency is

```text
omega = -u^mu nabla_mu varphi.
```

Let a finite domain `D_i` have an accessible density of modes `g_i(omega)` over
a finite operational band. Its cumulative accessible mode count is

```text
N_i(Omega) = integral_0^Omega g_i(omega) d omega.
```

The precise version of “a localized gradient such as a black hole adds an
extra hertz” is therefore not necessarily

```text
omega -> omega + 2 pi (1 hertz).
```

It is the testable possibility that a finite transition changes the density of
states,

```text
g_(i+1)(omega) = T_i{g_i}(omega) + Delta g_i(omega),
```

with a finite, potentially quantized increase

```text
Delta N_i = integral Delta g_i(omega) d omega.
```

`T_i` includes redshift, mode conversion, boundary conditions, and any change
of the proper-time convention. A physical prediction must fix a dimensionless
ratio, spacing pattern, or transfer feature such as

```text
omega_(n+1)/omega_n,
Delta omega/omega_reference,
Delta N_i,
```

before looking at the data. A literal universal increment of one SI hertz is
not invariant and is not part of the active hypothesis.

## Causal calibration is separate from spectral richness

A domain with more accessible modes does not automatically have a different
locally measured photon speed. Mode count, dispersion, redshift, clock rate,
and characteristic cone are distinct objects. For a domain-indexed line
element

```text
ds_i^2 = -N_i(t)^2 c_ref^2 dt^2 + a_i(t)^2 dchi^2,
```

a null ray has `dchi/dt=plus_or_minus N_i c_ref/a_i`, but the local comoving
observer uses `dtau_i=N_i dt` and `dell_i=a_i dchi`, so

```text
dell_i/dtau_i = c_ref.
```

Expansion, dilution, or a different lapse can change global horizon size and
coordinate travel time while leaving the local common cone unchanged. A
physical causal-scale claim needs two principal symbols derived from an action,

```text
G_a^(mu nu) k_mu k_nu = 0,
G_b^(mu nu) k_mu k_nu = 0,
```

and a dimensionless within-domain ratio `c_a/c_b`. A comparison across levels
additionally needs the same transition map that NGS needs for proper
frequencies: observers, clocks, phases, sectors, and state normalization must
be identified. The mapped candidate is therefore

```text
R_a|b^(i->j) = [(c_a/c_b)_j]/[(c_a/c_b)_i],
```

not a bare `c_i/c_j`. Without the map, the levels are operationally
uncompared. The full requirements are the
[DSF-3 finite-domain causal-calibration gate](domain-settings.md#dsf-3-finite-domain-causal-calibration-not-a-bare-variable-c).

The proposed “top” level should consequently be described, if its equations
support it, by a **large but finite relative cone**. A literal `c -> infinity`
is a singular change of causal structure, not an ordinary finite Lorentzian
domain, and conflicts with FGC's refusal to turn a limiting process into a
completed physical infinity. Likewise, entanglement is compatible with
nonclassical correlations while obeying no-signalling; it does not establish
instantaneous transport through a coarser level.

## Why spectral language is more than metaphor

QFT represents fields through modes of a differential operator or background,
schematically

```text
Phi(x) = sum_n [a_n u_n(x) + a_n^dagger u_n^*(x)],
```

or through a continuum spectral measure. The geometry, boundary conditions,
state, and observer determine what counts as a mode and how its frequency is
measured. Wilsonian effective field theory then introduces a resolution scale:
declared high-energy/short-distance degrees of freedom above a cutoff are
integrated out, and their effects reappear in effective couplings and
operators. A separately specified infrared, open-system, or statistical
coarse-graining can trace over other declared variables. Renormalization-group
flow is therefore a concrete language for physics organized across nested
resolutions, but the eliminated and retained sectors must be named.

Thermodynamics and statistical field theory add occupation, free energy,
response, transport, and entropy. A spectrum is populated according to a state;
temperature enters through dimensionless combinations such as
`hbar omega/(k_B T)`; gradients of thermodynamic potentials drive flux; and a
spectral response function states which frequencies a system can absorb,
transmit, or dissipate.

This does not make every cosmological use of “frequency” correct. A complete
NGS model must specify:

- the operator whose eigenmodes define the spectrum;
- the state and observer defining positive frequency;
- the physical cutoff or operational resolution;
- the coarse-graining map between levels; and
- the invariant spectral quantities that survive a change of units or chart.

With those supplied, frequency, bandwidth, and resolution are model variables,
not decoration.

## Particle as inherited pole, resonance, or band

NGS records the following candidate interpretation at H/Q status:

> What we call a particle may be a daughter-domain pole, resonance, or
> spectral band produced by coupling to modes of an inherited or adjacent
> phase.

The familiar billiard-ball picture is an introductory approximation, not the
ontology of quantum field theory. QFT already describes particle content as
quantized field excitations and identifies it operationally through states,
correlators, scattering amplitudes, and localized detector events. That makes
wave-like propagation and discrete interaction outcomes compatible without a
literal hard marble or a classical material fluid. NGS asks a narrower,
additional question: can one derived finite-gradient substrate explain why
particular excitations localize, persist, and recur across scales?

This does not deny the empirical reality of Standard Model degrees of freedom.
A particle is already represented operationally through the poles, resonances,
continuum support, quantum numbers, and detector responses of correlation
functions. A deeper substrate could change their ontology while leaving that
effective description compulsory in its successful domain. A sharp isolated
spectral feature could behave as a long-lived daughter excitation; a broadened
feature could represent an unstable state coupled across a transition; a
continuum could encode modes that do not localize completely in the daughter
description.

To turn the interpretation into physics, a solved transition must provide an
operator or kernel such as

```text
K_(i->i+1)(omega, omega_prime; boundary data),
```

and derive the daughter spectral density

```text
rho_(i+1)(omega)
  = integral K_(i->i+1)(omega,omega_prime) rho_i(omega_prime) d omega_prime
    + rho_nucleated(omega).
```

Its discrete poles, widths, continuum thresholds, and residues must then
recover the observed masses, lifetimes, spin/statistics, gauge charges,
unitarity, local Lorentz behavior, and scattering/decay rates. Quarks need not
be free asymptotic objects for their effective field description to remain
physically required; short-lived resonances need not “stay” in one layer to be
measurable. The nested explanation earns distinct content only if it freezes a
new line-shape, coupling relation, missing-flux channel, threshold, or
scale-dependent residual that ordinary within-domain QFT does not already
predict. Without the transition kernel and that discriminator, “ancestral
particle” is an analogy rather than a result.

Nor may every excitation simply be called a soliton. A soliton-like particle
interpretation requires an actual nonlinear finite-energy solution, a
stability mechanism such as topology, conserved charge, gauge flux, or a
proved dynamical attractor, and a quantization map that recovers the observed
spin and statistics. Likewise, “surrounding pressure holds the particle” must
be replaced by a defined stress tensor or energy functional with a stable
finite stationary point. Nuclear fission is correctly described as a
reorganization of nuclear binding and rest energy into fragment kinetic
energy, radiation, and other products; it is compatible with field ontology
but does not by itself establish a deeper pressure medium.

## Tension, instability, nonlinearity, and dissipation

Structure often appears when opposed terms select and stabilize a finite band.
An illustrative free energy is

```text
F[Phi] = integral d^3x [
    (r/2) Phi^2
    - (alpha/2) (D Phi)^2
    + (beta/2) (D^2 Phi)^2
    + (u/4) Phi^4
],

alpha,beta,u > 0.
```

The negative `k^2` term favors a finite-wavenumber band, while the positive
`k^4` term suppresses arbitrarily fine structure. For the displayed quadratic
kernel, the minimum is `r-alpha^2/(4 beta)`: the homogeneous state becomes
linearly unstable only when that quantity is negative. In that regime the
positive `u` nonlinearity can saturate the ordered amplitude. A dissipative
relaxation law might be

```text
partial_t Phi = -Gamma delta F/delta Phi + xi,
```

where `Gamma` is a response coefficient and `xi` represents declared noise.
This is a pattern-selection control, not a gravitational action.

The word “tension” can correspond to a real gradient or surface-energy cost.
The word “friction” can correspond to dissipation and entropy production.
They should not be collapsed into one universal force. Conservative Hamiltonian
and unitary quantum systems can also generate correlations, interference, and
complex structure without friction. The more general ingredients are:

```text
difference/gradient
  + coupling or instability
  + nonlinear or boundary constraint
  -> selected finite structure.
```

Dissipation can then select an attractor, remove excess free energy, or make a
pattern long lived. It can also erase structure. Its sign and channel must be
derived.

## The golden ratio is a possible eigenvalue, not a premise

The golden ratio can arise naturally from a particular two-step recursion. If
derived transition scales obey

```text
L_(i+1) = L_i + L_(i-1),
```

then a fixed asymptotic ratio `q=L_(i+1)/L_i` satisfies

```text
q^2 = q + 1,
q = (1 + sqrt(5))/2.
```

That is an elegant, falsifiable route from nested inheritance to the golden
ratio. Friction alone does not derive this recurrence. A physical claim would
need the action, conservation law, transfer matrix, or extremization principle
to produce the two-step relation and to specify which observable scales follow
it. Otherwise `q` is a chosen pattern, not an explanation.

More generally,

```text
L_(i+1) = a L_i + b L_(i-1)
```

gives `q^2=a q+b`. Measuring or deriving `a` and `b` must come before naming
the resulting eigenvalue. NGS therefore treats the golden ratio as a candidate
spectral discriminator, never as automatic evidence of cosmic design or
friction.

## Cumulative gradients as the restoring response

NGS sharpens the suggestion that no single gradient need carry the whole
regularizing burden. In a localization variable `k`, a schematic stiffness
kernel at level `i` could be

```text
K_i(k)
  = -alpha_i k^2
    + sum_(n>=2) beta_(i,n) k^(2n)
    + B_(i-1)(k),
```

where:

- the first term represents the lower-order localization tendency;
- the positive leading higher-order response represents the cumulative
  restoring sector; and
- `B_(i-1)` represents boundary data inherited from the coarser domain.

This is a spectral control, not a covariant field equation. A viable action
must show that the combined response has the required physical sign, becomes
dominant at a finite scale, has a bounded Hamiltonian or proved degeneracy,
propagates constraints, and passes the Raychaudhuri/source gate. Calling one
outer gradient “weak” does not make it dynamically decisive; the derived sum,
couplings, and boundary conditions must make it so.

The conceptual role of the coarser, less dense level is then precise: it acts
as an infrared boundary condition or capacity constraint on the finer level,
while the finer level contains the wider ultraviolet spectrum. This is a
candidate ultraviolet/infrared relation, not an established fact about
gravity.

## Complexity inheritance

“Complexity” needs an operational measure. Candidate proxies include:

- the finite accessible mode count `N_i`;
- a coarse-grained entropy or generalized entropy;
- an effective Hilbert/code-subspace dimension;
- algorithmic description length relative to a declared encoding; or
- the number of long-lived distinguishable structures above a resolution
  threshold.

These quantities are not interchangeable. The working conjecture is only

```text
finite parent state and flux data
  -> finite inherited constraints
  -> concentration and mode conversion
  -> a richer but still finite accessible child spectrum.
```

A capacity inequality might eventually take a form such as

```text
C_(i+1) <= C_capacity(E_i, A_i, flux_i, couplings_i),
```

where `C` is one declared complexity measure. The optional horizon-entropy
saturation postulate is one possible capacity rule, not a consequence of NGS.
Energy concentration also does not universally imply greater complexity:
thermalization, horizon formation, decoherence, or a phase transition can
reduce a chosen measure. The action and state must determine the direction.

## A route to dark components—not an identification

At level `D_i`, unresolved finer or coarser gradients could appear after
coarse-graining as an effective source,

```text
G_mu_nu[g_i]
  = 8 pi G (
      T_mu_nu^visible
      + T_mu_nu^unresolved
      + T_mu_nu^boundary
    ).
```

This gives a respectable way to ask whether some phenomena called “dark” are
the gravitational residue of degrees of freedom outside the resolved
description. It does not establish that they are.

To qualify as a dark-matter account, the derived sector must reproduce more
than extra attraction: abundance, clustering, lensing, structure growth,
halo-scale behavior, and early-universe constraints must agree in one model.
To qualify as dark energy, it must reproduce the observed expansion and
perturbation behavior, explain or control its scale, and remain compatible with
local gravity and multimessenger propagation. A background equation of state
alone cannot distinguish a nested-gradient origin from an ordinary internal
field or modified gravity.

The word “dark” should therefore retain its observational meaning: a component
inferred primarily through gravitational or expansion effects, not proof that
it lives in another domain.

### Relational dark-sector posture: no unaccountable variable

The anti-epicycle rule is not “discard every variable that has not been
directly touched.” Physics legitimately infers latent fields and parameters.
The stronger rule is:

> A new degree of freedom may do explanatory work only if it has an
> operational map to observations, changes a predeclared prediction or
> likelihood, and has a stated falsifier. Otherwise it remains a nuisance,
> convention, representation, or interpretation—not evidence for a new cause.

For a proposed extension `M_1=(M_0,q)`, the minimum variable ledger is

```text
q = {status, units, invariant observable, measurement map,
     symmetry/gauge redundancy, prior domain, falsifier},

y = H(q,theta) + measurement error.
```

The predictive ledger then states, for every data block `D_r`, the likelihood
`p(D_r|theta,q,M)`, which inputs calibrated the model, which data are held out,
and what outcome would reject it. A physical branch must change at least one
predeclared operational distribution,

```text
p(D_*|M_1) != p(D_*|M_0).
```

If two models predict the same distributions for every declared experiment,
they are empirically equivalent over that experiment family no matter how
different their nouns sound. Fewer named substances is not automatically
simpler if the replacement imports unconstrained functions, initial data, or
transition maps.

This keeps both directions open. Astronomical and cosmological observations
strongly require additional gravitating behavior in standard GR-based fits,
and late-time acceleration is observed through a concordant data set. They do
not directly identify one dark-matter particle or prove that the fitted
`Lambda` is a literal vacuum substance. NGS may seek a geometric,
thermodynamic, or relational origin, but it must compare against the full
LambdaCDM-plus-GR baseline rather than one residual:

- CMB acoustic structure and lensing, BBN baryon constraints, galaxy and
  cluster dynamics/lensing, growth, nonlinear structure, and halo behavior for
  a dark-matter claim;
- supernova/BAO/CMB distances and calibration, growth and lensing, local
  gravity, stability, radiative control, and multimessenger propagation for a
  dark-energy or modified-expansion claim; and
- both packages with one action and parameter set if a unified origin is
  claimed.

The project uses “relational” in this operational sense: redshifts, detector
rates, curvature scalars, relative cones, clock ratios, and mapped transfer
features. That stance is compatible with Rovelli's Relational Quantum
Mechanics but does not assume it; RQM is a specific interpretation of quantum
states relative to physical systems, not a derived FGC spacetime mechanism.

## Planck length and absolute zero are different limits

Planck length is a length scale assembled from `G`, `hbar`, and `c`. Absolute
zero is a thermodynamic limit. They are not opposite endpoints of one
dimensionful axis. A chosen thermal mode can be compared through the
dimensionless ratio

```text
x = hbar omega/(k_B T),
```

and a length can be converted to a characteristic frequency only after a
propagation law and observer are specified, for example `omega ~ c/ell` in a
declared relativistic control. At `T=0`, thermal occupation can vanish while
ground-state energy and quantum fluctuations need not.

NGS can still ask whether nested closure produces a finite hierarchy spanning
extreme localization and extremely cold regimes. It must not claim that
Planck length and absolute zero are already the endpoints of one proven cosmic
frequency spectrum.

## Dynamic balance, not perfect equilibrium

Our observable universe contains expansion, gravitational collapse, radiation,
stellar burning, irreversible structure formation, and entropy production. It
is not in exact thermodynamic equilibrium. The useful version of “the pocket
is in perfect equilibrium” is a **metastable operating band** or dynamical
balance:

```text
inward/localizing flux
  + cumulative restoring response
  + enclosing boundary conditions
  -> long-lived finite regime.
```

A solution may have a slowly changing attractor, limit cycle, stationary flux,
or bounded excursion. It must define which observables are balanced and over
which proper-time interval. A planet's orbit is an analogy for competing
tendencies, not evidence for a cosmological gradient hierarchy.

The sharper candidate is **stable asymmetry**. Let `x` be a regulated state,
`A(x)` a declared asymmetry order parameter, and

```text
dot(x) = F(x; lambda, J),
```

where `J` denotes actual boundary or reservoir fluxes. A stable asymmetric
fixed state requires

```text
A(x_*) != 0,
F(x_*;lambda,J) = 0,
max Re eigenvalue[D_x F(x_*)] < 0.
```

A stable cycle instead needs all transverse Floquet multipliers inside the
unit circle. Symmetric equations may possess asymmetric attractors, so this is
not a claim that symmetry of the laws makes complexity impossible.

For the narrower *nonequilibrium generative* branch, declare fluxes `J_r`,
energy-valued generalized forces `X_r`, and reservoirs with temperatures
`T_r`. With that convention its entropy-production rate must satisfy

```text
sigma = sum_r J_r X_r/T_r >= 0,
at least one J_r X_r > 0,
0 < integral_I sigma dtau < infinity
```

on every declared finite operational interval `I`, with the energy/work ledger
closed. If a convention instead defines `X_r` as an affinity that already
contains inverse temperature, use `sigma=sum_r J_r X_r`; every channel must
state its units and sign orientation. This permits finite sustained throughput
without hardcoding mechanical friction. An equilibrium symmetry-broken phase
can have `A!=0` with `J=0` and `sigma=0`; it is a genuine structured state, but
it is not the throughput branch proposed for persistent life-like processing.

## Dissipative complexity and eternal opportunity

Life-like organization is not an equilibrium object. It persists by consuming
free energy, maintaining gradients, and exporting entropy to its environment.
That gives the intuition about friction a strong thermodynamic core:

```text
available free-energy gradient
  -> work and information processing
  -> internal organization
  + dissipated heat and exported entropy.
```

“Friction” here should mean the full irreversible response channel, not only a
mechanical drag force. The detailed chemistry, quantum dynamics, feedback, and
selection remain essential. Dissipation does not by itself manufacture life,
but sustained life-like complexity cannot be treated as cost-free equilibrium
order.

An unbounded genealogy can also turn a rare outcome into an almost-sure one
under explicit assumptions. In the simplest control, if each sufficiently
independent finite domain has a fixed probability `p>0` of producing a declared
complexity class, then after `N` trials

```text
P(at least one occurrence) = 1 - (1-p)^N,

limit_(N->infinity) P = 1.
```

This captures the intuition that “unlikely here” need not mean “unlikely in an
unbounded process.” It is not a logical consequence of eternity alone. A
deterministic process may never visit the relevant state; `p` may be zero;
trials may not be independent; the measure may change; or the genealogy may
not be ergodic. Probability one is also an almost-sure statement, not proof that
every possible history contains the event.

NGS must therefore derive the domain ensemble or deterministic recurrence,
measure, correlations, and nonzero accessibility of the complexity class
before saying that observers “must” arise. The conditional calculation is a
valuable bridge between the philosophical intuition and a falsifiable
population model.

## Neither perfect order nor complete chaos is the law

The philosophical intuition is that a perfectly featureless frozen state has
no usable gradient, while uncontrolled divergence destroys a finite
description. NGS explores the structured regimes between those limits. It
does not assert that exact ground states, highly ordered phases, maximally
mixed states, deterministic chaos, or mathematical idealizations are forbidden
by nature. Thermodynamic equilibrium removes net affinities for sustained work
extraction, but it can retain energy, fluctuations, correlations, entanglement,
or a broken-symmetry phase. It therefore does not imply zero complexity or
instantaneous destruction without a specified Hamiltonian, ensemble,
environment, and complexity measure.

The physical claim to test is narrower: a **generative causal domain** may need
a finite nonequilibrium gradient budget—enough difference to sustain dynamics
and structure, but constrained by restoring response, available free energy,
and boundary capacity.

## Falsifiable requirements

NGS fails as a physical explanation unless an implementation supplies:

1. a covariant definition of the levels `D_i` and their transition maps;
2. observer-defined spectral densities and dimensionless cross-level
   comparisons;
3. a finite mode/complexity budget at every realized level;
4. a derived coarse-to-fine inheritance and fine-to-coarse backreaction law;
5. the cumulative restoring sign and finite FGC activation scale;
6. a closed constraint, energy, quasi-local flux, and generalized-entropy
   ledger;
7. a demonstrated metastable or bounded regime, with a nonzero asymmetry order
   parameter, stability spectrum, and—if throughput is claimed—finite
   reservoirs, fluxes, affinities, energy balance, and entropy production;
8. if dark sectors are claimed, one joint background, perturbation, lensing,
   structure, and local-test calculation;
9. a dimensionless spectral spacing, transfer feature, or population law fixed
   independently of the observations used to test it; and
10. an observation capable of distinguishing nested-domain inheritance from a
    single-spacetime internal field or ordinary modified gravity; and
11. a variable/ontology and predictive ledger that identifies calibration
    inputs, priors, held-out targets, the LambdaCDM-plus-GR comparator, and the
    exact result that would make the added degree of freedom unnecessary.

If a golden-ratio relation is claimed, an additional requirement applies: the
cross-level dynamics must derive the two-step recurrence and its domain of
validity before the ratio is compared with observations.

## NGS-0 executable control

The repository now executes the three conditional pieces that can be tested
without inventing a cosmology:

- `finite_band_selection(alpha,beta)` verifies the stationary nonzero scale of
  the spatial `-alpha k^2+beta k^4` control;
- `two_level_recurrence_ratio(a,b)` reports the positive eigenvalue and treats
  `a=b=1` as one fixture rather than a privileged law; and
- `eventual_occurrence_probability(p,N)` evaluates `1-(1-p)^N` with the fixed-
  `p`, independent-trial assumptions and non-necessity labels exposed.

The machine-readable record is
[`results/nested-gradient-spectrum.json`](../results/nested-gradient-spectrum.json),
regenerated by
[`scripts/reproduce_nested_spectrum.py`](../scripts/reproduce_nested_spectrum.py).
NGS-0 is D-C for these declared equations only. It derives no physical
recurrence coefficients, nested domains, spectral transfer, dark component,
inevitable observer, or metaphysical conclusion.

## Immediate NGS programme

1. **NGS-1 spectral definition:** choose one solved or controlled FGC
   transition, define its proper-time mode basis, density of states, and finite
   bandwidth, and calculate `Delta N` without importing a universal SI hertz.
2. **NGS-2 complexity/capacity gate:** choose one complexity proxy and derive
   its bound from energy, area, constraints, or flux rather than analogy.
3. **NGS-3 coarse-graining gate:** integrate out declared degrees of freedom
   and calculate the resulting effective stress. Test whether it behaves like
   neither, either, or both dark sectors without fitting the conclusion into
   the premise.
4. **NGS-4 discriminator:** derive a dimensionless cross-level spectral or
   population signature that an ordinary single-domain model need not share.
5. **NGS-5 recurrence gate:** test whether the derived transfer operator has a
   stable two-level eigenvalue relation; report the resulting ratio even when
   it is not the golden ratio.

Until those gates pass, NGS is the project's disciplined translation of a
creative idea: concentration may open finite complexity inside finite
complexity, while inherited large-scale constraints could prevent a singular
endpoint or an unlimited realized spectrum if the FGC and NGS gates are met.
