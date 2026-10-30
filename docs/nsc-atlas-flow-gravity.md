# From unresolved circulation to geometry and horizon exchange

This connection follows Douglas's revised direction: define the offset by
what it does in the geometry, before assigning a particle interpretation.
Two explicit quantities now enter the atlas: the covariance of unresolved
motion and a field's frequency relative to a rotating horizon. They have
different equations and observables. Neither is renamed baryon asymmetry.

The [calculation](../results/summaries/nsc-atlas-flow-evolution-v1.json)
extends the retained vortex preparation through a finite evolution and
evaluates the associated gravitational boundary and curvature coefficients.
The separate rotating-horizon calculation gives the signed energy-transfer
law of the complex field derived in the [harmonic sector](nsc-atlas-inherited-charge.md).

## Choice, dependencies and finite question

We compared identifying the fluid velocity with gravity by name, inserting
the whole Navier--Stokes blowup into a black-hole interior, and using the
explicit fluid/gravity construction with a separately checked horizon flux.
The last route supplies actual equations at each connection. It introduces
a declared periodic domain and numerical resolution for the finite flow;
the gravity cutoff equals the existing viscosity. It does not introduce an
observational fit. The horizon comparisons use normalized mode parameters,
not an inferred source spectrum.

The periodic flow is a new finite-domain continuation of the existing
compact initial preparation, not a replacement of the older whole-space
initial-rate theorem. Its domain is [-2,2)^3, viscosity .001, mode m=4 and
0<=t<=.05. Grids 64,96,128 and one 128-grid half-step control were fixed
before running. Two processes use four FFT workers each, within 300 process
CPU seconds per case. Arrays are discarded; only summaries are retained.

The three main risks are a numerical rather than physical mean circulation,
a normalization error between fluid time and boundary proper time, and
extending a local gravity expansion into a new cosmic region without a
matching solution. Initial-rate/refinement checks, explicit time/energy
factors and separated spacetime sectors address those risks.

## 1. Sphere and compression

For positive isotropic surface energy gamma A at fixed enclosed volume,

$$
A^3\geq36\pi V^2,
$$

with equality at the sphere. Its local stationary pressure jump is
Delta p=gamma(k1+k2)=2 gamma/R. Inside and outside are the two regions;
the factor two here comes from the two equal principal curvatures. This
is the [isoperimetric input](https://arxiv.org/html/1503.06705v2) already used
in the book, not a theorem that any two forces generate a sphere.

At fixed angular momentum J, Cauchy--Schwarz gives E_rot>=J^2/(2I).
Homologous compression with fixed mass lowers I proportional to R^2 and
raises this rotational-energy lower bound. It does not generate angular
momentum from exactly spherical nonrotating data. The covariance source
below instead supplies an explicit way to redistribute existing motion
into oppositely rotating regions while preserving the total budget.

## 2. The offset is a retained correlation, with finite evolution

The [existing compact field](nsc-atlas-vortex-transfer.md) has zero rotational
mean but covariance R_ij=<u_i u_j>-U_i U_j. Its initial tangential mean obeys

$$
\left.\partial_t U_\theta\right|_0
=-\left(\partial_r R_{r\theta}+\frac{2R_{r\theta}}{r}
+\partial_zR_{z\theta}\right).
$$

The exact rates are -6.490606541844 and +2.336618355064 at r=.75,1.25.
Here R is a tensor with an energy cost, not a missing particle population.
This extracts a mechanism used in OpenAI's
[Navier--Stokes construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf).
That published forced singularity theorem remains an attributed input;
the present evolution is finite, unforced and nonsingular over the run.

The numerical continuation uses strictly 2/3-dealiased Fourier Galerkin
evolution in rotational form, orthogonal incompressibility projection and
SSPRK3. Circle readouts use periodic cubic interpolation. Pressure is
recovered from the Poisson constraint with zero spatial mean. An independent
Beltrami decay and Taylor--Green pressure solution test the implementation.

| Grid, step | Generated U_theta(.75) | Generated U_theta(1.25) | Final kinetic energy | Final viscous loss |
| --- | --- | --- | --- | --- |
| 64, .001 | -.34496010 | .10598034 | 5.37735345 | 1.60511344 |
| 96, .001 | -.31452503 | .11972764 | 5.39003952 | 1.64636435 |
| 128, .001 | -.32027250 | .11367327 | 5.39173769 | 1.65733166 |
| 128, .0005 | -.32027248 | .11367327 | 5.39173769 | 1.65733167 |

Generated means subtract the same grid's initial mean. On the fine grid
that initial projection/readout offset is -.00059380 and +.00341907.
The final spatial differences shrink by factors 5.30 and 2.27; the time
differences are 1.18e-8 and 6.55e-10. Twice the last spatial/time difference
gives empirical indicators .01149 and .01211. These are not rigorous
continuum error bars. Initial rate errors at the fine grid are .00475
and .03580; the latter is not monotone across all three grids. The signs
and finite circulation survive the tested refinements, without claiming
high precision for every local derivative.

The maximum energy-balance residual across all runs is 1.81e-9 for
E(t)-E(0)+integral D dt. The fine divergence RMS stays below 3e-15 at the
reported endpoint. The four cases used 310.19 aggregate process CPU seconds
and 206.73 wall seconds. No giant waveform JSON or additional scan is stored.

## 3. The same flow gives a gravitational boundary response

[Bredberg, Keeler, Lysov and Strominger](https://arxiv.org/html/1101.2451v2)
give a vacuum-gravity metric through cubic hydrodynamic order. We apply
their equations 14--25 with three fluid spatial coordinates, hence a
five-dimensional bulk, cutoff r_c=nu=.001, and convention 16pi G5=1.
The fluid time is tau; boundary proper time is sqrt(nu) tau. Introduce
slow coordinates y=epsilon x and slow time t=epsilon^2 tau, with
v=epsilon u(y,t), P=epsilon^2 p(y,t).

Their boundary tensor, in proper boundary time, gives

$$
\rho_{\rm BY}=\frac{\epsilon^2 |u|^2}{\nu^{3/2}},\qquad
\langle\rho_{\rm BY}\rangle
=\frac{\epsilon^2\operatorname{tr}R}{\nu^{3/2}}
\quad\hbox{when }U=0.
$$

Thus eliminating the fine field while retaining only its zero mean would
delete a nonzero geometric boundary response. The corresponding spatial
metric covariance at fixed mean pressure is

$$
\delta g_{ij}^{\rm cov}
=\left(1-\frac{r_{\rm bulk}}{\nu}\right)
\frac{\epsilon^2 R_{ij}}{\nu}.
$$

An averaged metric is not automatically another vacuum solution. The
calculated connection retains the covariance rather than making that
additional closure assumption. The momentum constraint on the cutoff
reduces to the same Navier--Stokes equation with viscosity nu; a different
viscosity would not satisfy the published regular-horizon construction.

At leading nonzero order, their Riemann-squared invariant gives, in three
fluid spatial dimensions,

$$
\mathcal K=\epsilon^4\mathcal K_4+\cdots,\qquad
\mathcal K_4=-\frac{3}{\nu^2}|\nabla_y\times u|^2.
$$

The contraction is Lorentzian and may be negative. It is not the Ricci
scalar or a vacuum-energy density. Periodic solenoidality then establishes
the connected energy/geometry identities

$$
\int\frac{\rho_{\rm BY}}{\epsilon^2}\,d^3y
=\frac{2E}{\nu^{3/2}},\qquad
\frac{d}{dt}\int\frac{\rho_{\rm BY}}{\epsilon^2}\,d^3y
=-\frac{2D}{\nu^{3/2}},\qquad
\int\mathcal K_4\,d^3y=-\frac{3D}{\nu^3}.
$$

For the continuum initial preparation the boundary-energy coefficient is
346366.913683 and the integrated curvature coefficient is -5.167803776e9.
The fine numerical endpoint gives 341003.433221 and -4.971995004e9.
Initial continuum and projected-grid values differ by the recorded initial
discretization error; the discrete energy identity uses its own initial state.
Epsilon powers, slow integration coordinates and the stated G5 convention
must accompany these coefficients. They are not SI or cosmological numbers.

The published metric solves the Einstein equation through cubic order;
this work does not bound all higher orders. Its local Rindler horizon is
not yet joined to the four-dimensional positive-vacuum parent. The
dimensional reduction, shared source and matching between those sectors
remain explicit edges of the atlas, rather than being silently assumed.

In particular, the time derivative of the boundary-energy coefficient
enters at fourth order in the hydrodynamic expansion. On a closed flat
cutoff, exact vacuum Codazzi conservation would conserve integrated
Brown--York energy. The kinetic-energy decrease cannot simply be declared
a completed gravitational energy loss. A work term or a further state/
geometry correction must close that next-order balance. Its calculation
is the immediate continuation of this finite result.

### The next-order averaged work relation

The [boundary-work calculation](../results/summaries/nsc-atlas-boundary-work-v1.json)
now supplies that average balance under a homogeneous boundary-strain ansatz.
Let its spatial metric be b(t)^2 delta_ij, with b=1+epsilon^2 b2 and b2(0)=0.
The mean density is 2epsilon^2 E/(Vnu^(3/2)), V=64; the background boundary
pressure is p0=1/sqrt(nu). Vacuum Codazzi conservation gives, at fourth order,

$$
\frac{d\bar\rho}{d\hat t}+3H_b p_0=0,\qquad
H_b=\frac{\epsilon^4}{\sqrt\nu}\frac{db_2}{dt}+\cdots.
$$

Terms H_b times the perturbed density/pressure enter at sixth order. Keeping
the fourth-order terms and the same energy history gives

$$
\frac{db_2}{dt}=\frac{2D}{3V\nu},\qquad
b_2(t)=\frac{2[E(0)-E(t)]}{3V\nu}.
$$

This is a calculated relation between energy depletion and increasing
boundary scale. At the fine endpoint it gives **b2=.877929728565**. The
96-to-128 difference is .00697875 and the half-step difference is 1.26e-8;
twice their maximum is an empirical .0139575 indicator. Direct integrated
dissipation agrees to 2.10e-9 in b2 for the half-step run.

The actual scale correction is epsilon^2 b2; epsilon remains the formal
hydrodynamic parameter, not a selected cosmological scale. This closes
the averaged Ward identity with the declared homogeneous strain and no
extra bulk source. It does not solve the local fourth-order Einstein
equations or select this strain ansatz over other possible source/work
completions. The positive work term comes from geometric boundary pressure;
the incompressible fluid pressure itself does no net work on its periodic
kinetic energy.

For zero-mean periodic flow, D>=2nu(2pi/L)^2 E and D=0 only at rest. Energy
therefore provides a monotone ordering of this dissipative branch, and b2
increases with the energy transferred. This gives an arrow along an already
specified evolution; it does not derive the existence of time. The run's
loss decreases slightly, so this is not a late-time accelerating-cosmology
prediction. It is a finite, source-linked energy/geometry relation.

## 4. The rotating horizon has a calculable transfer threshold

The complex harmonic field already derived from the internal geometry has
canonical scalar stress. On a stationary axisymmetric horizon let
k=xi+Omega_H eta be its future null generator and use a regular mode
psi=A exp[i(m phi-omega t)]. The external winding m is distinct from the
internal charge +/-1. Define

$$
\Delta_H=\omega-m\Omega_H,\qquad
\mathcal N=2\int_H|A|^2dA.
$$

Contracting the canonical stress with the same Killing generators gives

$$
\mathcal F_E=\mathcal N\omega\Delta_H,\qquad
\mathcal F_J=\mathcal N m\Delta_H,\qquad
\mathcal F_E-\Omega_H\mathcal F_J=\mathcal N\Delta_H^2\geq0.
$$

These are fluxes into the future horizon. The scalar potential/mass term
drops out of the contraction because the generator is null and orthogonal
to the horizon tangents. At m=1, normalized modes give:

| omega/Omega_H | Energy flux / (N Omega_H^2) | Angular flux / (N Omega_H) | Co-rotating energy flux / (N Omega_H^2) |
| --- | --- | --- | --- |
| .5 | -.25 | -.5 | .25 |
| 1 | 0 | 0 | 0 |
| 1.5 | .75 | .5 | .25 |

The zero is synchronization, the threshold also used by
[Herdeiro--Radu](https://arxiv.org/html/1403.2757v3) for their rotating scalar
clouds and asymptotically flat hairy black holes. Negative Killing-energy
flux is the superradiant extraction channel; see
[Brito--Cardoso--Pani](https://arxiv.org/abs/1501.06570). This is not escape
of matter already inside an event horizon. The table holds supplied horizon
amplitudes fixed; it is not a greybody calculation or a selected spectrum.
No asymptotically flat hairy solution is imported as a de Sitter child.

The dimensionless frequency ratio and flux signs survive a positive
rescaling of the time coordinate. Energy, omega and Omega_H must use the
same frame. Internal charge reversal at fixed positive omega and m leaves
these stress fluxes unchanged on a neutral rotating background. A spatial
circulation offset is therefore not automatically a baryon imbalance.

## 5. The atlas connection and its completion test

The calculated route is: specified fine motion -> finite circulation and
conserved energy budget -> gravitational boundary/curvature response in a
defined lift. The second route is: inherited field -> regular rotating
horizon mode -> signed energy and angular-momentum transfer. Both use
explicit tensors and conservation equations. The remaining edge joins
their geometry, source and state to an evolving parent and child.

Douglas's clarified one-coordinate reset is a reduced generating
description for another region with three spatial dimensions. It does not
require a change in the number of spatial coordinates. A new expanding
region still requires a regular interior solution and compatible transition
data; a horizon flux supplies part of that data, not the complete birth.

The later golden-ratio proposal concerns a statistical outcome of branching.
That is different from identifying every geodesic spiral or every radius
ratio with phi. A branching mechanism, energy allocation and relative clock
map must determine that average. The present finite work relation gives one
geometric response to depletion; it does not yet select a branching law or
show that the rate accelerates through a region's full history.

Independent source reviews by Grok check the boundary normalizations and
the horizon contractions. Numerical controls test known exact flows, direct
energy conservation and the fixed refinements. A later calculation that
fails either the boundary momentum constraint or the horizon flux identity
under these same assumptions would contradict the claimed connection.

## Local fourth-order continuation

The [local work entry](nsc-atlas-local-work.md) now identifies the earlier
average coefficient as the zero mode of the published higher-order scalar
constraint. It constructs its nonzero spatial modes on the retained vortex
and independently checks a regular fourth-order bulk witness. The full
numerical vortex metric and global attachment remain open. The affine-flow
description and the earlier shift-free expanding boundary are not silently
identified; their coordinate shift/lapse difference is explicit.
