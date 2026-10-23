# Inherited preparation, nonlinear concentration and a received clock

The [compensated reservoir](nsc-atlas-core-reservoir.md) supplied exact
initial data. This calculation evolves the field and its gravitational
geometry together, tests a diffuse preparation, and transports a photon
from a central proper clock to a specified receiver. It uses the same
canonical scalar/four-form action and frequency-2000 member throughout.

This is spherical, finite-time Einstein--Klein--Gordon evolution. The
Schrodinger--Newton approximation supplies the earlier coherent initial
profile; it is not the evolution equation here. The diffuse preparation
has no prescribed final core profile. The inherited parent, ambient field,
reservoir and radiation all remain in the calculation.

## What is reused and what is connected

The spherical Einstein--scalar evolution framework and oscillaton precedent
are attributed to [Balakrishna et al., Section II](https://arxiv.org/html/0710.4131).
The regular-centre issues and constrained radial approach are described by
[Alcubierre and Gonzalez](https://arxiv.org/html/gr-qc/0401113).
[Hayward's mass law](https://arxiv.org/html/gr-qc/9710089) provides the
independent flux/work identity. The numerical integrating-factor approach
is attributed to [Lawson](https://doi.org/10.1137/0704033).
The added work connects these ingredients to the actual inherited,
compensated parent data and an operational observer in its evolving metric.

The initial t=1 trial resolved the field and geometry for over thirteen
coherent-core binding times. A t=2 diffuse trial showed that initial
spreading reverses. The admitted continuation to t=5 uses a moving outer
grid whose boundary stays in the protected reference annulus. Six worker
processes evaluate independent preparations, resolutions and controls
within a 3600 CPU-second campaign budget. No P001 grids are reopened.

## Fixed action and initial preparation

Keep kappa=1, omega=2000, epsilon=1/4, q_parent=3/2, the existing membrane
charge and tension, and the same initial background. On this parent,

\[
U(\phi,q_p)=\rho_v+\frac{\omega^2\chi^2}{2},\qquad
\chi=\phi-\phi_{\min}(q_p),\qquad \rho_v=\frac{9}{32}.
\]

The initial ambient displacement is .00043301270189221924. The reference
has a0=6.62607707610 and H0=.442688940038. The prior coherent control uses
zeta=.08266217132541569 and Rc=.120974320414.

For the diffuse test, replace only the core displacement by a Gaussian
of width .03, with the same C2 compact taper and support Rc. Its amplitude
.00186251874314 is determined by matching the **total geometric core
excess** at Rc to the coherent preparation, approximately 9.57140113e-5.
The initial central excitation is about 10.54 against the ambient 3/8,
so this is already a nonlinear diffuse preparation, not an evolved small
P003 perturbation. This comparison fixes gravitational mass, not isolated
quadratic scalar norm; the ambient cross term differs between the two shapes. That is a
declared preparation choice, not an observational fit or a formed-core claim.
The same reservoir matching and pulse complete each initial slice.

The initial core, reservoir and pulse are real-field Cauchy data. The
outward characteristic condition on the pulse holds at the initial instant;
massive propagation can later carry energy in both directions. No absorber
removes that field from the evolving computational domain.

## A coordinate system tied to the parent expansion

Use areal radius R and

\[
ds^2=-\alpha^2dt^2+
\frac{(dR-\alpha H(t)R\,dt)^2}{B}+R^2d\Omega^2,
\qquad B=1+H^2R^2-\frac{2m}{R}.
\]

H(t) is the actual unaveraged homogeneous parent expansion. This is a
slicing choice for the interior, not a requirement that the local matter
remain homogeneous. The local metric, lapse and radial curvature respond
to that matter. Define v=HR, beta=sqrt(B), F=phi_R and Pi=n(phi), with
n=alpha^{-1} partial_t+v partial_R and s=beta partial_R.

The momentum constraint gives

\[
K_t=-H,\qquad K_r=-H-\frac{\kappa R}{2}\Pi F,
\qquad K=K_r+2K_t.
\]

The Hamilton constraint determines the mass radially:

\[
m_R=\frac{\kappa R^2}{2}(\rho-v\Pi F),
\qquad \rho=\frac{\Pi^2+BF^2}{2}+U.
\]

Preserving the angular slicing condition determines the lapse:

\[
B\alpha_R=\dot H R+
\alpha\left(H^2R+\frac{m}{R^2}+\frac{\kappa R}{2}p_r\right),
\qquad p_r=\frac{\Pi^2+BF^2}{2}-U.
\]

The lapse equals one on the untouched outer reference. The central lapse
is solved, so central proper time is the integral of alpha(t,0). Regularity
requires alpha_R(0)=0; it does not also require alpha(t,0)=1. In particular,
alpha_R/R has a finite central limit, representing acceleration divergence.
Dropping that limit would incorrectly impose a second clock normalization
or an additional central condition on the chosen H(t).

The scalar evolves by

\[
\phi_t=\alpha(\Pi-vF),
\]

\[
\Pi_t=\alpha\left[
BF_R+\left(\frac{2B}{R}+\frac{B_R}{2}
+\frac{B\alpha_R}{\alpha}\right)F
+K\Pi-\omega^2\chi-v\Pi_R\right].
\]

The reference variables obey the same parent equations already used in
the atlas: chi_dot=Pi, Pi_dot=-3H Pi-omega^2 chi, a_dot=aH and
H_dot=a^{-2}-kappa Pi^2/2. The zero-perturbation limit exactly recovers
that geometry in the discretized radial solve. The independently monitored
background Friedmann residual measures its finite time-integration error.

At the centre, the radial Laplacian is 3F_R(0). In the H=0 limit, the shift
vanishes and the equations recover the standard polar-areal system.
No horizon or singularity inference is made from coordinate failure;
the finite calculation requires positive B and lapse throughout.

## Full mass conservation instead of numerical damping

Subtract the evolving background density, not its initial value:

\[
M=m-\frac{\kappa\rho_{\rm bg}(t)R^3}{6},\qquad
B=1-\frac{R^2}{a^2}-\frac{2M}{R}.
\]

The equations are evaluated in perturbations of the homogeneous scalar
to retain exact cancellation when the source vanishes. Reconstruction of
M from the radial constraint is checked against the independent time law

\[
m_t=\frac{\alpha\kappa R^2}{2}
\left[(B+v^2)\Pi F-v(\Pi^2+BF^2)\right].
\]

The background derivative is subtracted as well. Six fixed-radius flux
integrals are advanced with the field, then compared with the independently
reconstructed mass changes. This includes pressure work and scalar flux.
It is not a statement that each enclosed mass stays constant.

Spatial derivatives use an eighth-order centred stencil with parity at
the centre. Radial constraints use fourth-order integration. A Lawson RK4
step integrates the frequency-omega rotation exactly and advances the
remaining gravitational, gradient and expansion terms. This avoids treating
ordinary explicit-step oscillator damping as physical radiation. It does
not average away the carrier oscillations or their gravitational stress.

## Follow the radiation while preserving core resolution

The radial grid uses

\[
R(t,x)=x+\left(\frac{a(t)}{a_0}-1\right)g(x).
\]

The smooth function g vanishes for x<=1 and equals x for x>=2. The core
therefore keeps fixed physical resolution, while the distant mesh follows
the expanding parent. With J=R_x and u=R_t, spatial derivatives transform
as F=phi_x/J and F_R=phi_xx/J^2-phi_x J_x/J^3; time derivatives acquire
the grid-advection term uF or uPi_R. J remains positive on the expanding
branch. The map changes coordinates, not the action or physical source.

The outer initial radius x=4 corresponds to theta_out=.648103368884.
The initial perturbed ball ends at theta_b=.073094296566. The established
reference conformal bound eta_*<=.515418396612 gives

\[
\theta_b+\eta_*<\theta_{\rm out}<\theta_w-\eta_*.
\]

Thus the outer boundary stays in the protected homogeneous domain. Its
reference scalar data can be supplied causally, with no reflected artificial
wall. The whole outgoing scalar disturbance stays inside the expanding
domain for the admitted run. A larger outer-domain control checks numerical
sensitivity. This supplies a coupled core/parent domain, not the still-open
global material-wall future.

## An operational signal after concentration

The emitter is the central geodesic pole. The receiver is actively held at
areal radius .6, with proper clock zero on the initial slice. The coherent
control emits at central proper age 0. The diffuse run emits at central
proper age 4, after its contraction has developed. The reception event and
receiver age follow from the null ray; coordinate time is not substituted
for either proper clock.

For an outward photon with normal-frame frequency E,

\[
\frac{dR_\gamma}{dt}=\alpha(v+\beta),\qquad
\frac{d\log E}{dt}=\alpha K_r-\beta\alpha_R.
\]

The receiver moves at peculiar velocity -v/beta relative to the normal,
so its proper-time rate and measured frequency are

\[
\frac{d\tau_o}{dt}=\alpha\sqrt{1-\frac{v^2}{B}},\qquad
\frac{\nu_o}{\nu_e}=\frac{E_o}{E_e}
\sqrt{\frac{\beta_o+v_o}{\beta_o-v_o}}.
\]

Its timelike condition B>v^2 is checked. An independent null Hamiltonian
calculation uses C=alpha(v+beta), p_R=E/beta and

\[
\frac{d\log p_R}{dt}=-\partial_R C.
\]

The two frequency equations agree only when the evolving geometry satisfies
the appropriate Einstein relation. Their symbolic equivalence and numerical
agreement supply a check independent of the radial mass reconstruction.

Emission is located by the central proper clock; the small partial emission
step interpolates metric data. Reception interpolates the arrival event.
These finite event errors are included in the numerical comparisons; the
clock calculation is not advertised as uniformly fourth order merely
because the bulk field stepper is fourth order.

The reported signal is the frequency departure from an independently
integrated homogeneous parent with the same source proper age and held
receiver radius. Its arrival proper-time difference is a second output.
The homogeneous null path also obeys R_gamma=a(t)sin(eta(t)-eta_em), giving
a separate direct comparison to the photon integrator. P003's linear
frequency-20 coefficients remain their own earlier conditional result.

## Numerical evidence and interpretation

Over parent times 4.5--5, the half excess-mass radius inside R=.12
averages .01849207, compared with .03930864 initially: a ratio .470433.
The late central excitation averages 148.3729 against 10.53893 initially,
a factor 14.0786. With the same field on the homogeneous metric it instead
averages .0477285. Thus the concentration depends on the computed local
gravitational response. The medium/fine spreads are 1.19e-7 in mean radius
and .00232 in mean central excitation.

| Resolution multiplier | Maximum integrated mass-flux gap, diffuse run |
| --- | ---: |
| 1 | 2.55e-7 |
| 2 | 1.27e-8 |
| 4 | 8.42e-10 |

The finest run keeps the local lapse above .98145. Its final exterior
mass offset is 2.04e-9, against the initial core excess 9.57e-5. These finite
residuals are numerical diagnostics, not exact-zero assertions. The
formation is followed by substantial breathing; no stationary endpoint
is imposed or inferred.

For the diffuse run emitted at central proper age 4, the fractional
frequency departure from the homogeneous comparison is

$$
\frac{(\nu_o/\nu_e)}{(\nu_o/\nu_e)_{\rm hom}}-1
=-0.0041375586\pm 1.1\times10^{-8}.
$$

This is approximately -0.413756 percent, with a finite numerical error
estimate of .0000011 percentage points. The receiver's additional proper
arrival time is .01703546, with numerical comparison error 2.2e-7 in atlas
time units. Both the emitter's accumulated proper clock and propagation
contribute to this delay. The independent canonical frequency differs by
less than 4.7e-11 on the finest run.

The coherent control, emitted at central age 0, has frequency departure
-.007456287 with numerical comparison error 2.1e-7 and proper arrival delay
.000476466. These use different declared emission ages; they are not a
same-epoch comparison of the two states. All errors here concern numerical
resolution and event calculation, not physical preparation uncertainty.


The [compact record](../results/summaries/nsc-atlas-core-evolution-v1.json)
stores the three resolutions, time and outer-domain controls, source hashes,
mass-flux comparisons and received signals. Its formation table distinguishes
initial spreading, later concentration and breathing. The same diffuse field
on a homogeneous fixed geometry is a diagnostic control: it intentionally
omits local gravitational response and is not a second full Einstein solution.

Finite concentration is evidence for a mechanism in this declared initial
family. It does not prove a stationary remnant, an attractor for arbitrary
preparations, nonspherical stability or a primordial population. Source
shape/phase and the frequency remain inputs. No empirical abundance or
physical SI length is fitted. Numerical error envelopes compare resolutions,
time steps, outer domains, independent photons and independent background
integration; they are not rigorous error bounds or model uncertainty.

The surrounding reference still has Lambda_vac*a0^2=12.3482523989285.
Rearrangement of its inherited excitation now develops into nonlinear local
structure in a finite sourced evolution. It does not select a new vacuum
coefficient or identify the observed cosmological constant.

## Reproduction

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_core_evolution.py --output results/summaries/nsc-atlas-core-evolution-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_core_evolution.py
```

The producer runs six independent worker processes. Full radial states are
hashed and regenerated when needed; transient plotting histories are not
committed. Next connections concern the inherited preparation/population
law, the longevity and nonspherical behaviour of the concentrations, and
additional proper-age signals from the same frozen state.

The completed campaign used 3431.34 process-CPU seconds including workers
and 1179.09 seconds wall time, with six independent worker processes. Grok
independently checked the expanding-slice equations and audited the moving
map, Lawson stages, regular centre and observer implementation. One proposed
central-lapse restriction was rejected by retaining the acceleration-divergence
term; the final implementation review confirmed the exterior clock normalization.
