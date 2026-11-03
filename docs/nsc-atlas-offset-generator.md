# One spatial seed and a repeated offset: a constructive geometric map

## Active start: parent compression, not an isolated offset

Douglas's correction, 2 November 2026: the child's proposed generating
seed arises through collapse inside an evolving parent, itself inside an
enclosing region. The offset must be connected to that inherited drive.
The isolated-loop transfer and fixed-alpha calculations below are retained
controls. Their open/closed comparisons are not a task to select a preferred
percentage before constructing the parent-child relationship.

### An additive measure from the enclosing geometry

Let F_p(u), F_c(u) describe corresponding local material-patch deformations
in a common comparison frame. Their orientation-preserving volume factors
are J_p=det F_p and J_c=det F_c. The label u indexes the supplied history;
it is not an independently postulated physical clock. For spatial dimension
d define normalized volume factors j_a(u)=J_a(u)/J_a(u_0), and

$$
\chi_{p\to c}(u)=\frac{1}{d}\log\frac{j_p(u)}{j_c(u)}.
$$

For the constructed three-dimensional patches d=3. This is logarithmic
relative LENGTH compression. If the child contracts in every direction
by half while its parent patch is unchanged, chi=ln 2. Its normalization
comes from taking the d-th root of a volume factor, not fitting a target.
It requires nonzero-volume comparison patches; it is not itself a birth
solution at a degenerate one-dimensional seed.

The logarithm is forced if a continuous scalar measure of a positive
relative scale must add under composition: C(ab)=C(a)+C(b) gives C(a)=k ln a.
Taking the logarithm of the relative length ratio fixes its geometric unit.
This does not establish a conversion to angular phase, seconds or energy.
Identifying inherited compression with the proposed generating offset
remains the explicit physical proposal being developed.

The increment delta chi is calculated from the deformation history.
It is not an independently chosen alpha. Where derivatives exist,

$$
\frac{d\chi}{du}=\frac{1}{d}
\left[\mathrm{tr}(F_p^{-1}F_p')
-\mathrm{tr}(F_c^{-1}F_c')\right].
$$

A common invertible frame change F_a -> B(u)F_a cancels in the relative
ratio, including its initial normalization. For three nested regions,

$$
\chi_{g\to c}=\chi_{g\to p}+\chi_{p\to c}.
$$

Thus zero ADDITIONAL child compression does not erase the inherited parent
contribution. This is an exact finite-chain log-ratio telescoping identity
(bookkeeping, not dynamical evidence for the mechanism); it assumes neither
a first universe nor a common increment for every level. It is not a
derivation of an infinite physical chain.

### Executed controls and the next connection

Analytic controls F_g=exp(-u)I, F_p=exp(-2u)I, F_c=exp(-4u)I give
chi_gp=u, chi_pc=2u and chi_gc=3u, with zero sampled composition residual.
These are test deformations, not selected collapse histories. Independent
controls check common-frame cancellation, reverse sign under exchanging
regions, ancestor inheritance under common compression, and log-Jacobian
evaluation at compression 400 without determinant underflow. Pure rotation
has unit volume factor: rotation and compression cannot be identified by
volume data alone.

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --inheritance-only
```

**Next:** obtain the actual relative deformation, including turning/shear,
from a shared parent-child compression construction and connect its inherited
increment to the drawing operation. A pressure-to-deformation law, angular
phase relation and receiving geometry have not been derived by this
kinematic identity. No free increment is used to stand in for those links.

## Coupled finite parent, receivers and characteristic continuation

The active dynamical completion uses the same angular phase and A0 stretch
cost, with explicit geometric inertia. It replaces prescribed enclosing
pressure and prescribed compression by a finite surrounding region whose own
energy does work. This is a candidate mathematical completion, not a proof
that the recurrence uniquely selects the inertia or the regional topology.

**Calculation map.** `parent_receiver_characteristics` evolves the connected
phase and geometry. `_parent_receiver_system` and `parent_receiver_evolution`
provide an independent finite-element Hamiltonian control. The existing CLI
`--evolution-only` runs three optical-step comparisons for a chain and a fork,
plus two short finite-element controls, in eight worker processes. It writes
no files. The readout remains in dimensionless offset units.

### One action, finite pressure and work counted once

For an ordered regional tree let q_i=ln R_i, S_i=R_i^3 and
$$
 V_i=S_i-\sum_{j\in\mathrm{children}(i)}S_j>0,\qquad
 \lambda_i=V_i^{1/3}.
$$
This is an available-volume reduction, not an embedded multi-hole metric.
The path-length/volume identification and isotropic log-boundary inertia
are declared premises. In particular the outer radius is free; no total
volume conservation or inexhaustible pressure reservoir is imposed. On one
generated material curve of length two per region, take
$$
 L={3\over2}\sum_i\dot q_i^2+
 {1\over2}\sum_i\int_0^2
 [\lambda_i\dot\psi_i^2-\lambda_i^{-1}\psi_{i,x}^2]dx.
$$
The kinetic term is the proposed dynamical addition to the static A0 cost.
The circle phase is the existing field, not newly added particle content.
At a contact its value is shared and outward phase tractions sum to zero;
there is no material membrane or point mass. Fixed phase endpoints close
the finite reference domain and can reflect waves. Initially each available
volume is one, all momenta vanish, the root carries the retained triangular
lift 2*pi*min(x,2-x), and receiving curves are quiet. The topology and
receiving domains are admitted data, not dynamically created spaces.

Canonical momentum pi_i=lambda_i*psi_dot gives regional phase energy
$$
 E_i={1\over2\lambda_i}\int_0^2(\pi_i^2+\psi_{i,x}^2)dx,
 \quad P_i=-\partial_{V_i}E_i={E_i\over3V_i},
 \quad\dot p_i=S_i(E_i/V_i-E_{p(i)}/V_{p(i)}).
$$
The outer enclosing pressure is zero. The initially excited parent drives
inward acceleration of its quiet inner boundary. The received phase energy
subsequently contributes an opposing pressure. No function P(t), compression
history, damping coefficient, or target abundance is supplied.

Writing F_i for net phase power into a region gives exact balances
$$
 \dot E_i=F_i-E_i\dot{\ln\lambda_i},\qquad
 \dot K_{\rm geometry}=\sum_iE_i\dot{\ln\lambda_i},\qquad
 \sum_iF_i=0.
$$
Thus sum(E_i)+K_geometry is conserved and parent work is counted once.
The phase part's conjugate pressure/density ratio is 1/3; it has not been
identified with the full cosmological stress. Also sum(p_i_dot)=sum(E_i):
local compression can coexist with increasing overall expansion momentum.

### The exponential grid break and the continuation

The initial finite-element calculation conserves its own Hamiltonian, but
the expanding parent stretches its last material cell enormously relative
to a receiver. Three meshes do not resolve the late fork transfer uniformly.
That endpoint is not a continuum claim. Instead change to each region's
optical offset s_i, with ds_i/dtau=1/lambda_i. The exact regional solution is
$$
 \psi_i=A_i(s_i-x)+B_i(s_i+x-2).
$$
Delayed outgoing histories and their integral/squared-integral primitives
carry the same phase into later regions without a spatial grid. They give
E_i directly as the sum of squared characteristic histories over the last
two optical units, divided by lambda_i. RK4 advances only the finite
geometry and those histories; optical steps and interpolation retain an
explicit numerical error. Before a front arrives, the next quiet region's
field is exactly zero in this representation. Finite-element precursors
are not interpreted as received energy.

At a degree-m junction the physical incoming/outgoing amplitudes satisfy
$$
 b=\left({2\over m}{\bf1}{\bf1}^{T}-I\right)a,\qquad
 S^TS=I.
$$
A two-curve junction transmits completely. One parent incident on two quiet
receivers reflects energy 1/9 and transmits 4/9 into each, hence 8/9 in total.
These are exact contact relations of the declared action. They are not the
earlier 88.4735% scalar-horizon result. Returning waves subsequently change
the instantaneous split. A continuous sign/orientation convention can be
carried on the graph; neither that sign nor this transparent contact proves
that the contact is the proposed one-way null horizon.

### Exact free-parent expansion and an optical future

For one room the full conformal phase budget Q0=4*pi^2 is constant. Let
a=sqrt(Q0/6). The same equations, from rest and unit initial size, have
$$
 q=2\ln\cosh(a\tau),\quad\lambda=\cosh^2(a\tau),\quad
 s={\tanh(a\tau)\over a},\quad
 E_\psi={Q_0\over\lambda},\quad K_g=Q_0\tanh^2(a\tau).
$$
This supplies expansion and decreasing local phase energy without an added
vacuum potential. In the declared curve metric -d(tau)^2+lambda^2 dx^2,
the future optical reach is finite, sqrt(6)/[2*pi]=0.3898484006. The proper
distance to the observer's causal-reach boundary is
$$
 \lambda(\tau)\int_\tau^\infty{dt\over\lambda(t)}
 ={1+e^{-2a\tau}\over2a}\longrightarrow0.1949242003.
$$
This finite expression continues the apparently troublesome product of an
unbounded scale and a vanishing optical tail. The late scale rate is
2a=5.1301993206. If that same scale is used for the earlier flat 3D metric
completion, its late geometric curvature has Lambda_curvature=3*(2a)^2
=8*pi^2 in these units. This is a conditional metric-curvature result,
not an SI vacuum prediction or a derived Einstein source equation.

### Zero-radius formation with finite radial inertia

For the logarithmic inertia and conserved budget H, each |d log R/dτ|
is at most sqrt(2H/3). Integrating bounds R below by its positive reference
radius times exp(-sqrt(2H/3)*|Δτ|). A leaf at R=0 is therefore at infinite
kinetic distance. An internal available volume can still vanish through
cancellation between positive enclosing volumes. Finite phase energy then
forces its end-to-end phase difference to vanish as O(V^(1/6)).

Consider the kinetic metric induced by mean squared radial displacement of
a spherical shell with normalized angular measure. This gives dR² and
admits zero radius at finite distance. In this candidate reduction keep two
affine parent segments of length y and one receiving segment of length r.
Use the continuous phase values (0,M,Θ,0). The root starts with y=1,
M=W=2π and zero velocities, preserving the triangular winding. The receiving
arc starts with r=Θ=0. Its admission is a premise of this continuation.

The kinetic and stretch terms are

$$
K=\frac{\dot y^2+\dot r^2}{2}
 +\frac{y}{6}(2\dot M^2+\dot M\dot\Theta+\dot\Theta^2)
 +\frac r6\dot\Theta^2,\qquad
Q=\frac{M^2+(\Theta-M)^2}{2y}+\frac{\Theta^2}{2r}.
$$

The phase mass matrix at birth has entries (2/3, 1/6, 1/3) in its upper
triangle. Eliminating the midpoint acceleration gives reduced mass
s=7/24 and s Θ''(0)=3W/2-k for limiting traction k=Θ/r. The radial
equation gives r''(0)=k²/2, hence k+s k³/2=3W/2. Its unique positive
root is 3.4477970166. With a=k²/4=2.9718260669, the opening has
r=aτ²+O(τ⁴), Θ=kaτ²+O(τ⁴) and receiving energy 4a²τ²+O(τ⁴).

For local existence divide r, Θ, y−1 and M−W by τ², and the two radial
velocities and two phase momenta by τ. The resulting system τX'=F(τ,X)
is analytic near the positive limiting coefficient a. At the limiting
state the variations w=δ(Θ/τ²)−k δ(r/τ²) and b=δ(Θ'/τ)−k δ(r'/τ)
satisfy w'=b−2w and b'=−b−(1/s+k²)w/a in logarithmic offset. The other
linear blocks have eigenvalues −1 and −2. All eigenvalues have negative
real part. Writing A=F_X at the limit, the integral operator for the
deviation is ∫₀¹ u^(−A) [F(uτ,X(uτ))−A(X(uτ)−X₀)] du/u. Its kernel is integrable
and its nonlinear part has arbitrarily small Lipschitz constant near the
limit. Contraction gives the local branch in the bounded rescaled class.
The forcing starts at order τ², giving the stated remainder.

Run `python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --birth-only`.
At offset 0.1 the radius is 0.028963565226, volume 2.429719051e-5 and
received energy 0.306711998348. Seeds 1e-5 and 1e-6 differ in radius by
4.5e-15 and the largest fine budget defect is 1.92e-13. These are numerical
checks of the affine continuation. The causal and spatial completion remains
open. If r is used as an isotropic volume radius, energy per volume grows
as τ⁻⁴ at birth even though energy per arc length stays finite.

### Continuum wave at the opening

Keep the radial-displacement inertia and the triangular parent profile.
Allow the phase to propagate on both material intervals. Write y for the
parent scale and r for the receiving scale. The parent incoming amplitude
at the contact is π/y until the first optical front arrives. A transparent
two-curve junction passes that characteristic into the receiver. In its
optical coordinate ds/dτ=1/r the outgoing derivative is A'(s)=πr/y. The
unit receiving interval has phase A(s−x)−A(s+x−2), with zero phase at x=1.
This branch uses zero incoming conformal wave energy in the optical past.

For optical lag z solve ∂_z t_z=−r(t_z), t_0=τ. Define
g_z=exp(−∫₀ᶻ r'(t_w)dw)=r(t_z)/r(τ). The exponential formula stays
regular at zero radius. The exact receiving force and the contact power are

$$
r''=\pi^2\int_0^2 \frac{g_z^2}{y(t_z)^2} dz,\qquad
J=\pi^2\left[\frac{1}{y(\tau)^2}-\frac{g_2^2}{y(t_2)^2}\right].
$$

The receiving phase energy is r r'' and its contact phase is
π∫₀² r(t_z)/y(t_z) dz. The parent satisfies
y² y''=4π²−∫₀ᵗ yJ du, from its characteristic flux ledger. These
equations give a local fixed point for the continuous accelerations.
To see this, reconstruct both radii by double integration on [0,T] in a
bounded box around the initial accelerations (2π²,4π²). Differences in
the radii and velocities are O(T²) and O(T) times the acceleration norm.
The lag flow stays in [0,τ] and depends Lipschitz continuously on r.
The exponential ratio then varies by O(T) times that same norm. The
parent integral and denominator also have Lipschitz bounds tending to
zero with T. Thus the update maps the box into itself and contracts for
sufficiently small T. The phase and radius work terms cancel in the total
budget by the same contact balance as for the positive-radius curves.

At birth g_z tends to one and y tends to one. Consequently

$$
r(\tau)=\pi^2\tau^2-\frac43\pi^4\tau^3+O(\tau^4),\qquad
\Theta/r\to2\pi,\qquad E_{\rm receiving}=4\pi^4\tau^2+O(\tau^3),
$$

where the last expression includes radial kinetic energy. The cubic
term follows from g_z=1−2π²zτ+O(τ²). The coefficient π² replaces
the affine phase approximation's 2.9718260669 in this continuum model.

`python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --continuum-birth-only`
evaluates the coupled equations through τ=.05 at 65, 129 and 257 samples.
The fine radius is .016736471416 and the receiving energy .299289538820.
The sampled budget discrepancy decreases from 2.33e-8 to 8.34e-10 to
2.79e-11. The last two radii differ by 2.5e-13. These are numerical
comparisons on the declared early interval.

The curve metric −dτ²+r² dx² has scalar curvature 2r''/r∼4/τ².
The flat isotropic three-dimensional extension would have Kretschmann
scalar ∼240/τ⁴. A change of coordinates cannot produce a C² extension
with these invariants. In the physical-distance chart X=rx the curve
metric is −dτ²+(dX−v dτ)², with v=r'X/r. Within the shrinking wedge
|v|≤|r'|→0. Extending v by clamping X/r to [0,1], and setting v=0
before birth, gives one continuous Lorentzian extension of the metric.
A physical phase and geometry matching law and selection of the opening
event remain open. The transparent contact used here is bidirectional.

### Shared mouth geometry and the returned signal

At positive offset, measure physical distance X from the mouth using
X=y(x−2) on the parent and X=rx on the receiver. Both propagation metrics
then have the form

$$
ds^2=-d\tau^2+(dX-H_i X d\tau)^2,\qquad H_i=\dot a_i/a_i,
\quad a_i=y\ \text{or}\ r.
$$

They agree at X=0, where the induced metric is −dτ² and the two null
slopes are ±1. The interface is a geodesic on each side. Its extrinsic
curvature vanishes on both sides. Sidewise coordinates
T=τ+H_i(τ)X²/2 cancel the linear cross term and yield matching first
metric derivatives there. Curvature can still jump. These statements
concern the positive-offset mouth, with the opening tip treated above.

With the fixed receiving endpoint, a phase packet emitted at τ_e returns
at τ_b satisfying ∫[τ_e,τ_b] dτ/r=2. Its amplitude reverses sign and
is multiplied by r(τ_e)/r(τ_b). Differentiation gives
dτ_b/dτ_e=r(τ_b)/r(τ_e), so the narrow-packet power and energy ratios are
respectively the square and first power of r(τ_e)/r(τ_b). Transparent
matching transmits the return into the parent. The existing continuum CLI
now includes a packet emitted at .02, returning to the mouth at
.029154034883 and reaching a parent detector at material distance .01
at .039385342912. Its mouth energy ratio is .502902410597. The final two
sample counts differ in detector offset by 2.5e-12. The detector and the
narrow probe on the calculated geometry are declared readout choices.

### Inertia criterion for birth and finite future optical reach

For a radial kinetic term M(r) rdot²/2 and a finite positive conserved
budget E with nonnegative remaining energy, |rdot|≤sqrt(2E/M). A monotone
opening therefore requires ∫₀^ε sqrt(M(r)) dr finite. An unbounded growing
branch with finite optical future requires ∫^∞ sqrt(M(r)) dr/r finite.
The latter criterion is also sufficient on a monotone expanding tail if
radial kinetic energy is bounded away from zero there. For unit radial
inertia, r(τ) grows at most linearly under the conserved budget, so its
optical future diverges on any forever-regular continuation. For
logarithmic inertia the birth integral diverges. The two tests concern
opposite ends of the scale range.

A candidate that passes both integral tests is M(r)=1/(1+r²/3). This is
the radial metric coefficient of a hyperbolic spatial geometry in its
areal-radius coordinate. Its small-radius normalization and large-radius
coefficient reproduce the two previous inertia limits. The choice of
this interpolation remains a modeling assumption. For the separate
free one-room solution starting from rest at r=1 with phase budget
Q₀=4π², energy conservation gives kinetic energy Q₀(1−1/r). Its entire
future optical reach, using u=1/r, is

$$
\frac1{\sqrt{2Q_0}}\int_0^1\frac{du}{\sqrt{(u^2+1/3)(1-u)}}
=0.2611589800633733.
$$

The integral is evaluated with u=1−z². This gives a concrete candidate
for a single inertia law spanning the two regimes. A coupled trajectory
from birth to an isolated late region remains to be established, with the
physical mouth and its transmission law part of that calculation.

### Last echo from a coupled hybrid opening

Apply M(a)=1/(1+a²/3) to both radii in the continuum opening, with the
same triangular initial phase and transparent mouth. In the canonical
radius u=√3 asinh(a/√3), kinetic energy is u'²/2 and
u''=E_phase cosh(u/√3)/a. The finite-inertia birth proof carries over,
with r/τ²→π² and parent y''(0)=16π²/3. Optical histories are carried
from the early solution into the retarded evolution.

The parent incoming amplitude π/y persists on this branch. Before its
optical distance from birth reaches one, its conformal phase budget
C=yE_parent obeys C'≥−π²/y, since the reflected power is nonnegative.
Thus C≥3π² on that interval. Integrating dK_parent/dy=C/y² gives
K_parent≥3π²(1−1/y). Substituting the hybrid inertia into the optical
integral bounds total parent optical distance by √2/π<1. This excludes
a first crossing of the triangular midpoint front. The receiver phase
energy is a positive history integral and the parent budget has the
lower bound just proved. Both canonical speeds increase from positive
values after birth and stay bounded by the conserved energy. Positive
radii and bounded canonical coordinates on finite intervals allow the
retarded equations to continue by steps for every finite positive offset.

At a checkpoint T let v_i=u_i'(T)>0, j be the other region and H=4π².
All subsequent speeds of region i lie between v_i and √(2H−v_j²).
Integrating the reciprocal radius along those two constant-speed curves
bounds its remaining optical distance between

$$
\frac{\mathrm{arsinh}(\sqrt3/a_i(T))}{\sqrt{2H-v_j^2}}
\quad\text{and}\quad
\frac{\mathrm{arsinh}(\sqrt3/a_i(T))}{v_i}.
$$

The inequalities apply to exact checkpoint states. The executable uses
floating-point states and an empirical energy allowance. Its refinement
differences quantify numerical sensitivity rather than certify an interval
for the exact checkpoint.

The receiver optical future decreases continuously from infinity at birth
to zero. The unique launch cutoff for a reflected echo is its crossing of
two optical units. The cutoff is a supremum of returning launch offsets.
At equality the return is at infinite offset. At the parent detector a
material distance .01 from the mouth, the echo needs another .01 units of
parent optical distance after it crosses the mouth. This gives an earlier
launch cutoff. The command `python scripts/run_python.py
scripts/check_nsc_atlas_offset_generator.py --last-echo-only` compares three
birth resolutions, tolerances and late checkpoints. Their cutoffs are
.22694238706 at the mouth and .11650293537 at the detector. The finest
two differ by less than 2.1e-11 in either cutoff. The fine budget
discrepancy is 3.94e-11 and the remaining receiver optical tail at T=24
is about 5.92e-10.

Both canonical velocities tend to positive limits. To obtain the asymptotic
rate, note that both radii grow at least exponentially after a positive
checkpoint. In the receiver phase energy integral, the factor r(u)/r(t)
decays exponentially with t−u and the incoming power falls as y(u)^−2.
The parent energy satisfies its work equation with source −J. Both terms
in J decay exponentially because the reflected emission tends to the finite
cutoff. This gives exponential phase energy decay, allowing polynomial
prefactors when rates coincide. Since cosh(u/√3)/a is bounded on each late
tail, the remaining acceleration and its first time moment are integrable.
Thus each radius is asymptotic to a constant times exp(κ_i τ), with
κ_i=v_i(∞)/√3. The computed rates are .90249212 and 5.05019337.

Let δ be the gap between a launch and its cutoff, and κ the rate of the
region in which the final leg ends. Its arrival offset is
−log(δ)/κ+O(1). Its energy ratio is asymptotic to κδ and its power ratio
to (κδ)². For the mouth this follows by substituting the optical tail
1/(κr) into the radius ratio. For the detector the derivative of the
mouth return time with respect to launch time is r_return/r_launch.
That factor cancels the receiver energy redshift when converting the
remaining parent optical distance into δ. Each decade closer to cutoff
adds about 2.5514 offset units at the mouth or .4560 at the detector.

There is also a null boundary inside the receiving region. In its material
coordinate x the edge of the mouth observer's causal past is
x_H(τ)=∫_τ^∞dt/r(t) wherever this lies in [0,1]. It enters at the far
endpoint near .635108 and its physical distance tends to 1/κ_receiver,
about 1.108043. The timelike mouth remains the phase contact. The atlas's
exterior (+), horizon (0), interior (−) convention concerns a physical
black-hole handoff. Identifying that handoff with this observer horizon
requires an additional geometric construction. The present results use
an admitted receiving branch, the chosen hybrid inertia and a fixed
reflecting endpoint.

### Scope and the next connection

The small [calculation record](../results/summaries/nsc-finite-parent-receiver-v1.json)
records the finite chain, branch, independent controls and empirical spreads.
The phase source is carried unchanged between regions. Child **formation**,
selection of the initial tree and geometric inertia, and identification of
the shared cut with a causal horizon remain distinct connections. A signed
radius alone does not supply them. Offset parametrizes this motion; its
ultimate origin and angular calibration are not established by the action.

The failed long fork comparison is retained as a diagnosis of interpolated
reflected fronts. The finite pre-reflection branch and exact scattering law
are separate results; no reflected cyclic history is claimed. No cosmological
target, eta rung or 95% floor enters these calculations.

The next work is the causal/contact geometry and its dynamically admissible
branch transition. Three foreseeable obstacles are an inertia not selected
by generated paths, loss of positive available volume, and discontinuous
reflected histories. They appear respectively as an assumption-dependent
trajectory, a volume guard, or noncontracting contact/energy errors. None is
a reason to launch additional shape or cosmic-number fits.

## Intrinsic deformation response and a conservative neck control

The next candidate uses the existing angular phase with a unit target
metric and the actual spatial metric:

$$
Q_g=\frac12\int g^{ij}\partial_i\psi\partial_j\psi dV_g.
$$

Choosing this functional and target metric is explicit. It is not the old
fixed-drawing-coordinate functional with a new name, and has no assigned
joule conversion. On a circle of radius R with phase m theta it gives
Q=pi m^2/R. Replacing theta by any smooth increasing periodic parameter
leaves this value unchanged when both the metric and volume factor are
transformed. The calculation checks this directly.

### How inherited compression changes the cost

For an affine deformation A of a reference patch, let
M=integral grad(psi) grad(psi)^T dV_reference. Then

$$
Q(A)=\tfrac12\det(A)\mathrm{tr}(A^{-T}MA^{-1}).
$$

The code represents M=C^T C, keeping cost positivity explicit. It does not
introduce another field. For uniform A=lambda I in dimension d,
Q(lambda)=lambda^(d-2) Q(1). Thus half-size costs relative to the original
are 2 in d=1, 1 in d=2, and 1/2 in d=3. These are scale controls, not
a selection of spatial dimension.

For the SAME carried reference phase pattern, this connects directly to
the enclosing-deformation result above:

$$
\frac{Q_c}{Q_p}=e^{-(d-2)\chi_{p\to c}}.
$$

This is a derived response to inherited compression, with no isolated alpha
or independently chosen damping rate. It is the phase-gradient portion of
a geometric budget, not a claim that total energy disappears.

### Work and directional response

Let Q_i be the contributions of the squared current phase gradient in
orthogonal directions. Changing a current principal length gives
partial Q/partial ln(lambda_i)=Q-2Q_i. These derivatives define the
generalized response; they are not inserted pressure labels. Under a
constant isotropic enclosing compressive load P>0, the load potential is
P V (equivalently, input work on compression is -P dV). In d=3,

$$
\frac{d}{d\lambda}[Q_0\lambda+P V_0\lambda^3]
=Q_0+3P V_0\lambda^2>0.
$$

Consequently this candidate alone has no finite isotropic holding radius
under that load. A one-dimensional circle's resistance cannot simply be
carried over to a three-dimensional region. The constant load is a symbolic
control; a physical enclosing load has not been obtained from this identity.

At FIXED volume, an isotropic M instead favours equal principal lengths
by the arithmetic/geometric-mean inequality. This is shape stability,
not stability of overall size. For A(s)=diag(exp(-s/2),exp(-s/2),exp(s)):

| Log axial stretch s | Q for normalized axial gradient | Q for M=I |
| --- | --- | --- |
| 0 | 0.5 | 1.5 |
| 1/2 | 0.1839397206 | 1.8326609913 |
| 1 | 0.0676676416 | 2.7859494701 |
| 2 | 0.0091578194 | 7.3982139184 |

An inherited directional phase pattern can therefore favour elongation and
transverse narrowing. These two reference patterns are analytic controls,
not selected preparations of a universe. A rigid rotation changes neither
cost; this functional alone supplies no torque that creates a vortex.

### Exact motion and its budget

For the normalized axial pattern Q(s)=exp(-2s)/2. Declare the geometric
inertial metric sum_i(d ln lambda_i)^2=(3/2)ds^2. From rest at s=0 its
conservative motion obeys

$$
\ddot s=\frac23e^{-2s},\qquad
s(\tau)=\log\cosh(\sqrt{2/3} \tau),\qquad
Q+\frac34\dot s^2=\frac12.
$$

This adds an explicit candidate inertia for deformation; it is not silently
deduced from the earlier phase kinetic rule. Offset tau is normalized and
has no calibration to the generating offset or seconds. Fixed volume makes
the isotropic load do zero work on this shape mode. The decreasing Q is
accounted for by increasing geometric motion:

| tau | Log axial stretch | Q | Geometric kinetic | Total |
| --- | --- | --- | --- | --- |
| 0 | 0 | 0.5 | 0 | 0.5 |
| 1 | 0.3017834121 | 0.2734288042 | 0.2265711958 | 0.5 |
| 2 | 0.9772952024 | 0.0708112364 | 0.4291887636 | 0.5 |
| 4 | 2.5742942146 | 0.0029037984 | 0.4970962016 | 0.5 |

### A0 least-stretch justification, enclosing compression minimization, and break continuation

Adopting the 3D volume-weighted functional $Q_g=\frac12\int|\nabla_g\psi|^2 dV_g$
without deriving it from the generated geometry dilutes the 1D phase gradient
by the 3D volume factor $\det(A)\propto\lambda^3$, so $Q_g(\lambda)=\lambda Q_0=e^{-\chi}Q_0$
decreases to zero under 3D compression and offers no holding resistance.

In the offset generator, the geometry is drawn by 1D phase curves (spokes/windings).
By Cauchy-Schwarz on a 1D curve of physical length $L_i=\lambda_i L_{i,0}$ carrying
phase change $\Delta\psi_i$, the least-stretch cost with respect to induced arc
length $ds_g$ is uniquely
$$
\frac12\int_0^{L_i}\left(\frac{d\psi}{ds_g}\right)^2 ds_g
\ge \frac{(\Delta\psi_i)^2}{2\lambda_i L_{i,0}}=\frac{m_i}{\lambda_i},
\qquad
Q_{\rm stretch}(A)=\tfrac12\mathrm{tr}\!\left[(A^TA)^{-1/2}C^TC\right].
$$
Under uniform compression $A=\lambda I=e^{-\chi}I$, $Q_{\rm stretch}(\lambda I)=\lambda^{-1}Q_{\rm stretch}(I)=e^{+\chi}Q_{\rm stretch}(I)$
in every dimension $d\ge1$, resisting compression rather than vanishing.

Minimizing total cost $E(\lambda_1,\lambda_2,\lambda_3)=\sum_{i=1}^3\frac{m_i}{\lambda_i}+PV_0\lambda_1\lambda_2\lambda_3$
under enclosing compression $P>0$ gives the exact unique equilibrium for $m_i>0$:
$$
\lambda_i^\star=\frac{m_i}{(m_1m_2m_3)^{1/4}(PV_0)^{1/4}},
\qquad
Q_{\rm stretch}^\star=3PV_0\lambda_1^\star\lambda_2^\star\lambda_3^\star=3W_{\rm ext}^\star.
$$

With $PV_0=1$, the three generator regimes give:

1. **Isotropic 3D winding ($m=(1/2,1/2,1/2)$) $\to$ uniform contraction:**
   $\lambda_1^\star=\lambda_2^\star=\lambda_3^\star=2^{-1/4}\approx0.8408964153$,
   $Q_{\rm stretch}^\star\approx1.7838106765$, $W_{\rm ext}^\star\approx0.5946035588$,
   $E_{\rm total}^\star\approx2.3784142353$.
2. **Raw two-accumulator moment ($m=(1/4,1/4,1/2)$) $\to$ localized prolate neck equilibrium:**
   $\lambda_\perp^\star=2^{-5/12}\approx0.7491535384$, $\lambda_3^\star=2^{7/12}\approx1.4983070769$,
   aspect ratio $\lambda_3^\star/\lambda_\perp^\star=2$, $s_\star=\tfrac13\log2\approx0.2310490602$.
3. **1D axial seed ($m=(0,0,1/2)$) $\to$ transverse neck-pinch break:**
   $\lambda_\perp^\star\to0$ ($s\to+\infty$), because zero transverse winding cannot
   hold the transverse radius against enclosing compression.

**Break continuation into the child frame.** In parent coordinates, the axial
neck motion $s(\tau)=\log\cosh(\sqrt{2/3} \tau)$ drives the transverse throat
radius $a(\tau)=e^{-s(\tau)/2}\to0$ while the normalized velocity fraction
$\dot s/\sqrt{2/3}=\tanh(\sqrt{2/3} \tau)$ hits the $0.999\dots$ wall
($0.9998888764$ at $\tau=6$, where $a(6)\approx0.1220964547$). Applying the
break rule via throat inversion $R_{\rm child}(\tau)=a(\tau)^{-2}=e^{s(\tau)}=\cosh(\sqrt{2/3} \tau)$
continues the trajectory into an expanding child frame obeying the exact
Friedmann-like balance
$$
\left(\frac{\dot R_{\rm child}}{R_{\rm child}}\right)^2+\frac{2}{3R_{\rm child}^2}=\frac23,
\qquad
H_{\rm child}(\tau)\to\sqrt{\frac23}\approx0.8164965809,
\qquad
\Lambda_{\rm child}=3H_\infty^2=2.
$$
At $\tau=6$, $R_{\rm child}(6)\approx67.0801295585$ and $99.9778$% of the initial
parent phase cost has transferred into child expansion kinetic budget.

**Turning and spiral comparison.** Both $Q_{\rm vol}(A)$ and $Q_{\rm stretch}(A)$
are $SO(3)$-invariant functions of symmetric stretch and produce zero rotational
torque on their own. When unit logarithmic compression $\lambda(\chi)=e^{-\chi}$
is paired with a supplied or inherited angular turn rate $\omega=2\pi\alpha$,
$z(\chi)=e^{(-1+2\pi i\alpha)\chi}$ has decay-per-radian $\gamma/\omega=1/(2\pi\alpha)$
and opposite-lobe asymmetry $\eta=(A_+-|A_-|)/(A_++|A_-|)=\tanh(1/(4\alpha))$.
For $\alpha_{\rm golden}=(3-\sqrt5)/2\approx0.3819660113$, $\gamma/\omega\approx0.4166730505$
and $\eta\approx0.5746972340$ ($|A_-|/|A_+|\approx0.2700854214$), whereas a golden
logarithmic spiral has $\gamma/\omega=2\ln\varphi/\pi\approx0.3063489625$ (which
would require $\alpha=1/(4\ln\varphi)\approx0.5195217303$).

- **Coupling $c$:** Still supplied along the 1D drawing coordinate; the 3D neck shape wave speed $\sqrt{2/3}$ follows only after declaring the log-strain kinetic metric.
- **Increment $\alpha$:** Still supplied (or inherited from external rotational/chiral parent history) because symmetric compression costs are $SO(3)$-invariant and generate zero torque.

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --intrinsic-only
```

**Next connection:** derive a chiral/rotational coupling from the two-sided
(`+`/`0`/`-`) contact so turning and $\alpha$ emerge from the compression
response rather than remaining supplied.

## Unified single-cost cycle, helical buckling, damped `+`/`0`/`-` asymmetry, and Goldilocks stabilization

To eliminate fluid-physics assumptions (parcel incompressibility $\lambda_1\lambda_2\lambda_3=1$,
surface tension, or an ad-hoc volume rule $V_-=V_0(1-e^{-3\chi})$), `unified_least_resistance_cycle`
minimizes **one single 5-variable cost functional** built purely from **A0**
(least resistance: $Q_{\rm stretch}+PV$), **A1** (quadratic winding cost $m\propto n^2$
from `winding_variation`, with parent convergence $\lambda_1=\lambda_2=\xi=e^{-\chi}$
in `+` and 3D offset distribution in `-`), and **A4** (turn budget split $n_++n_-=N_0$,
$n_+=(1-\theta)N_0, n_-=\theta N_0$ across `0`):

$$
E(\lambda_3,\tilde\lambda_1,\tilde\lambda_2,\tilde\lambda_3,\theta;\chi)
=\left[\frac{2m_0(1-\theta)^2}{e^{-\chi}}
+\frac{m_0(1-\theta)^2}{\lambda_3}
+PV_0e^{-2\chi}\lambda_3\right]
+\left[\sum_{i=1}^3\frac{m_0\theta^2}{\tilde\lambda_i}
+PV_0\tilde\lambda_1\tilde\lambda_2\tilde\lambda_3\right],
$$

where the first bracket is $E_+$ (infall region) and the second bracket is $E_-$ (across `0`).

Simultaneous stationarity $\partial E/\partial\lambda_3=0$, $\partial E/\partial\tilde\lambda_i=0$,
and $\partial E/\partial\theta=0$ (with strictly positive-definite $5\times5$ Hessian)
gives the unique closed-form global minimum (for normalized $m_0=PV_0=V_0=1$):

$$
\sqrt{\theta^\star(\chi)}=\frac{\sqrt{16e^{2\chi}+17}-3}{4e^\chi},\qquad
\lambda_3^\star(\chi)=(1-\theta^\star(\chi))e^\chi,\qquad
\tilde\lambda_1^\star=\tilde\lambda_2^\star=\tilde\lambda_3^\star=\sqrt{\theta^\star(\chi)}.
$$

All three stages emerge as **consequences** of this single minimization:
- **Parcel shrinkage AND 1D spaghettification on the way in (`+`):** Without
  assuming fluid incompressibility, the `+` parcel volume shrinks as
  $V_+^\star(\chi)=(1-\theta^\star)e^{-\chi}\sim\tfrac32e^{-2\chi}\to0$, while its
  aspect ratio stretches into a 1D filament for all $\chi>\tfrac12\ln2\approx0.3465735903$
  (where $\theta^\star=1/2$ and aspect ratio is exactly $1$):
  $\lambda_3^\star/\xi=(1-\theta^\star)e^{2\chi}\sim\tfrac32e^{+\chi}\to\infty$.
- **Rounding into a 3D sphere in `-` without surface tension:** By AM-GM on
  $\sum_{i=1}^3\theta^2/\tilde\lambda_i+\tilde\lambda_1\tilde\lambda_2\tilde\lambda_3$,
  the `-` scales satisfy $\tilde\lambda_1^\star=\tilde\lambda_2^\star=\tilde\lambda_3^\star=\sqrt{\theta^\star}$
  (aspect ratio $1$, a round 3D sphere).
- **Derived cooker volume takeover without assuming $V_++V_-=V_0$:** The `-`
  drop's share of the two-region volume, $f_{\rm drop}^\star(\chi)=V_-^\star/(V_+^\star+V_-^\star)$
  with $V_-^\star=(\theta^\star)^{3/2}$, rises from $37.90$% ($\chi=0$) to $65.05$%
  ($\chi=0.5$), $84.63$% ($\chi=1$), $97.73$% ($\chi=2$), and $99.95$% ($\chi=4$).

The calculation then connects the nested signed-zero crossing map, the
two-force `+ -> 0 -> -` zero-crossing, helical buckling, damped oscillation and
Goldilocks stabilization through least resistance and budget conservation:

- **Nested signed-zero crossing map and the $\sim95$% exponential horizon (`nested_zero_crossing_map`).**
  Iterating the openly declared signed-zero crossing axiom ($R$ passes through
  `0` into `-` with $\mathrm{sign}(\det F_k)=(-1)^k$ at each one-way `0`
  wall) places every local room $k$ between an **outer birth `0` valve** (inflow
  from parent $k-1$) and **inner black-hole `0` valves** (outflow to children
  $k+1$). Because the room is continuously fed from outside while draining
  inward, multiple compressions cannot empty the room of its visible matter.
  At one 3D volumetric e-fold ($3\chi_1=3$, i.e. $\chi=3$):
  - The 3D exponential horizon is $1-e^{-3}\approx0.9502129316$ ($95.02$% dark
    sector) with visible floor $e^{-3}\approx0.0497870684$ ($4.98$%), split by
    the 3D Virial lock ($3:1$) into $\tfrac34(1-e^{-3})\approx71.27$% Dark
    Energy-like expansion stretch and $\tfrac14(1-e^{-3})\approx23.76$% Dark
    Matter-like boundary load.
  - Evaluating the 5-variable single-cost minimum $E^\star(\chi)=3(\theta^\star)^{3/2}+[(\theta^\star)^{3/2}+2(1-\theta^\star)e^{-\chi}]+2(1-\theta^\star)^2e^\chi$
    directly at $\chi=2.5$ and $\chi=3.0$ gives:
    - $\chi=2.5$: $68.21$% Dark Energy-like, $23.24$% Dark Matter-like ($91.46$% combined dark sector), $8.54$% Visible filamentary.
    - $\chi=3.0$: $70.72$% Dark Energy-like, $23.79$% Dark Matter-like ($94.51$% combined dark sector), $5.49$% Visible filamentary.
    - $\chi\approx3.1045$: $71.05$% Dark Energy-like, $23.95$% Dark Matter-like ($95.00$% combined dark sector), $5.00$% Visible filamentary.
  - **Antimatter check and spectral band shift (`spectral_band_and_baryon_harmonic`).**
    The one-way signed-zero crossing selects a single bulk orientation $(-1)^k$
    inside each child room. Furthermore, each room $k$ at relative compression
    $\chi$ has a finite wavelength band $[\lambda_{\min}(\chi),\lambda_{\max}(\chi)]=e^{-\chi}[\lambda_{\min,0},\lambda_{\max,0}]$
    stretching from its shortest wavelength before `0` transitions to `-`
    (the room's Planck cutoff $\lambda_{\min}$) to its longest wavelength as
    wide as the space itself (the room's Absolute-Zero horizon $\lambda_{\max}$).
    Both endpoints shift down by $e^{-\chi}$ inside more compressed children
    and up by $e^{+\chi}$ in the parent, while keeping the room's dynamic range
    $\mathcal{B}=\lambda_{\max}/\lambda_{\min}$ invariant. Every mode of
    wavelength $\lambda\in[\lambda_{\min},\lambda_{\max}]$ riding the room's
    horizon decay envelope $e^{-\gamma\tau}$ has half-cycle damping
    $e^{-\pi(\lambda/\lambda_{\max})}$ under the spectral identification
    $\gamma/\omega=\lambda/\lambda_{\max}$, giving the amplitude-area asymmetry curve
    $$
    \eta(\lambda)=\tanh\!\left(\frac{\pi}{2}\frac{\lambda}{\lambda_{\max}}\right).
    $$
    At $\lambda=\lambda_{\max}$ (Absolute-Zero horizon), $\eta=\tanh(\pi/2)\approx0.917152$;
    at the harmonic ratio $\lambda_B/\lambda_{\max}=(2/\pi)\mathrm{artanh}(\eta_B)\approx3.88338\times10^{-10}$
    ($\approx31.26$ octaves below the horizon, well above the Planck cutoff),
    the same curve yields $\eta(\lambda_B)=6.1\times10^{-10}$.
  - **Bidirectional branching chain and cyclic vortex-diffusion ring (`branching_cosmos_chain_evolution`, `cyclic_vortex_diffusion_ring`).**
    When an expanded parent (`O`) births $b$ child black-hole seeds (`o`),
    transferring fraction $\theta^\star(\chi)$ across the `0` valves drains the
    parent from `O` into a compressed remnant `o` while ballooning the children
    into expanded `O` rooms (`Oooo -> oOOOooooooooo -> ...`) with zero budget
    loss. On a closed ring of coupled rooms (`oo / OO oO / oO Oo / oo`),
    vortex concentration (`O`, shedding into nested mini-hurricanes toward
    $\lambda_{\min}$) and outward diffusion into still air (`o`, relaxing toward
    $\lambda_{\max}$) circulate around the ring so diffused air re-accumulates
    and spins back up into an enormous hurricane `O` every full cycle.
  - **Non-local to local Atlas readout, two-valve fixed point $\chi_\star$, and $K=7$ baryon harmonic (`nonlocal_to_local_atlas_readout`).**
    Solving the zero-parameter two-valve fixed point where the 5-variable
    single-cost visible filament share $\Omega_{\rm vis}(\chi_\star)$ equals the 3D
    volumetric floor $e^{-3}\approx0.0497870684$ gives $\chi_\star\approx3.0669683518$
    ($\Omega_{\rm DE}(\chi_\star)\approx71.1441$%, $\Omega_{\rm DM}(\chi_\star)\approx23.8772$%,
    combined dark sector $1-e^{-3}\approx95.0213$%, visible filament share
    $e^{-3}\approx4.9787$%; at integer $\chi=3$, $\Omega_{\rm DE}(3)\approx70.8827$%,
    $\Omega_{\rm DM}(3)\approx23.8134$%, $\Omega_{\rm vis}(3)\approx5.3039$%).
    Along the discrete 3D volumetric e-fold harmonic ladder
    $\eta_K=\tanh\!\left(\frac{\pi}{2}e^{-2\chi_{\rm spaghetti}}e^{-3K}\right)=\tanh\!\left(\frac{\pi}{4}e^{-3K}\right)$
    (where $e^{-2\chi_{\rm spaghetti}}=1/2$ is the exact spaghettification onset
    from `unified_least_resistance_cycle`), the $K=7$ harmonic ($3K=21$ e-folds)
    gives $\eta_7=\tanh\!\left(\frac{\pi}{4}e^{-21}\right)\approx5.955329\times10^{-10}$,
    inside the observed BBN window $[5.8\times10^{-10},6.5\times10^{-10}]$.
  - **Circle-pendulum duality, why `-` exists, why `(-) × (-) = (+)`, and why the radius of space must cross `0` (`signed_circle_pendulum_duality`).**
    Seen edge-on, the Layer-1 closed circle $z(\tau)=R_0e^{i\omega\tau}$ is the
    harmonic pendulum $R(\tau)=\mathrm{Re}z(\tau)=R_0\cos(\omega\tau)$,
    which swings $+R_0\to0\to-R_0\to0\to+R_0$ with invariant energy
    $\tfrac12\dot R^2+\tfrac12\omega^2R^2=\tfrac12R_0^2\omega^2$. Without
    negative coordinates ($R<0$), $50$% of the circle ($\pi R_0$ out of
    $2\pi R_0$) is missing and the loop cannot close. Multiplying by $-1=e^{i\pi}$
    is a half-turn (in 3D, spatial inversion $\mathcal{P}=-I_3$ with
    $\det\mathcal{P}=(-1)^3=-1$), so two nested crossings satisfy
    $\mathcal{P}^2=(-I_3)(-I_3)=+I_3$ and $\det(\mathcal{P}^2)=(-1)(-1)=+1$.
    When standard physics restricts the radius of space to an unsigned size
    $r(\tau)=|R(\tau)|\ge0$, folding $R<0$ onto $r>0$ at $\tau_0=\frac{\pi}{2\omega}$
    creates an artificial velocity cusp $\Delta\dot r=+2R_0\omega$ and a
    Dirac-delta acceleration singularity $\int_{\tau_0^-}^{\tau_0^+}\ddot r d\tau=2R_0\omega$,
    whereas the signed coordinate $R(\tau)$ crosses `0` into `-` with zero
    velocity jump ($\Delta\dot R=0$) and zero acceleration ($\ddot R(\tau_0)=0$).
  - **Vacuum bubble crystallization, surrounding-medium flow energy, $180^\circ$ ripple annihilation, shared-flow resonance, and discrete bursts (`vacuum_bubble_crystallization_and_mass_spectrum`).**
    Integrating the radial phase-flow gradient in the surrounding medium $r\ge R$
    outside a bubble of radius $R$ with integer winding $n\in\mathbb{Z}$ gives
    $E_{\rm surrounding}(R)=m_0n^2/R$ ($50$% in $[R,2R]$, $90$% in $[R,10R]$,
    $0$% in a point singularity). Balancing against cavity pressure work
    $W_{\rm pressure}(R)=PV_0R^3$ crystallizes the wave at
    $R_n^\star=(m_0n^2/(3PV_0))^{1/4}\propto|n|^{1/2}$ and $E_n^\star\propto|n|^{3/2}$,
    with $75$% of the mass-energy in the surrounding medium flow and $25$% in
    cavity pressure work. Two opposite bubbles ($n_+=+1, n_-=-1$) have net
    winding $0$, unlocking $R^\star\to 0$ and releasing $2E_1^\star$ into the open wave
    equation as two equal back-to-back ripples ($E_{\rm left}=E_{\rm right}=E_1^\star$
    at $180^\circ$). Superposing two bubbles at separation $d$ gives interaction
    energy $\Delta E_{\rm int}=-(2m_0n_1n_2/d)\cos((\omega_1-\omega_2)\tau+\Delta\phi_0)$,
    which time-averages to $0$ off-resonance and yields an attractive $-1/d^2$
    shared-flow force when resonance-matched ($\omega_1=\omega_2$). Because
    loop closure $z(1)=z(0)$ forces $n\in\mathbb{Z}$, bursting transitions shed
    sharp quantized line energies $\Delta E=E_{n_i}^\star-E_{n_f}^\star$ rather than a
    smear.
  - **Derived 3-accumulator volumetric bubble ladder vs. quarantined post-hoc coincidences (`derived_three_accumulator_bubble_ladder`).**
    Iterated **A1** offset accumulation yields binomial windings $W_k(n)=\binom{n}{k}$
    and 3D equilibrium energy $E_{\rm 3D}^\star(n)=4(PV_0)^{1/4}m_0^{3/4}\sqrt{\binom{n}{1}\binom{n}{2}\binom{n}{3}}$.
    This proves from **A0+A1** that $n=1$ ($W_2=W_3=0$) is a 1D loop, $n=2$ ($W_3=0$)
    is a 2D sheet, and **$n=3$ ($W_1W_2W_3=9>0$) is the minimum integer supporting
    a 3D volume**, with untuned ratio $E_{\rm 3D}^\star(3)/E_{\rm 1D}^\star(1)=6$ and
    3D ladder $1\text{ (1D)}:6:19.60:44.72:84.85$ (which misses $1836.15$ by
    $\approx 306.03\times$ without a continuum/nesting scale factor). Post-hoc
    numerical coincidences (Lenz 1951 $6\pi^5+e^{-3}/\sqrt{2}$, $n\text{--}p$ gap
    $2+\tfrac12+e^{-3}/\sqrt{2}$, rung selection $K=7$ on $\eta_K$, CMB projections
    using measured $H_0/G$ and $-4/5$ in `thermodynamic_gravity_and_cmb_temperature`,
    and the $71.14$% / $23.88$% / $4.98$% split vs. $\Lambda\text{CDM}$ $68.5$%)
    are strictly quarantined in `numerical_coincidences_not_derived`.
  - **Finite-Facet Sphere, Planck/Absolute-Zero Resolution Band, and Stellar Collapse into a Black Hole (`stellar_resolution_and_neutron_star_collapse`, Proposition 9.2b).**
    Under **Axiom A3**, a sphere sits between the Absolute-Zero floor ($\lambda_{\max}=2R$) and the Planck ceiling ($\ell_P=\sqrt{\hbar G/c^3}$) and has finite facet count $F<\infty$ ($A_F^3/(36\pi V_F^2)=1+\pi\sqrt{3}/F+O(F^{-2})>1$). A self-gravitating star saturates its Planck/Compton facet resolution at $\mathcal N_{1\mathrm{D}}=\lambda_n/\ell_P=M_P/m_n=1.29942\times 10^{19}$, $F_\star=(\lambda_n/\ell_P)^2=1.688\times 10^{38}$ facets, $N_\star=(\lambda_n/\ell_P)^3=2.194\times 10^{57}$ baryons, and $\sqrt{R_\star\ell_P}=\lambda_n$, fixing the base resolution mass $M_{\star,0}=M_P^3/m_n^2=1.8481\ M_\odot$, the neutron-star radius $R_{\mathrm{NS}}=4\lambda_n^2/\ell_P=10.916\ \mathrm{km}$ ($R_S/R_{\mathrm{NS}}=1/2=\theta^\star(\chi_{\mathrm{spaghetti}})$), and the TOV black-hole collapse threshold $M_{\mathrm{TOV}}\in[2.0414,2.1454]\ M_\odot$.
  - **Universal Exponential-Wall Bypasses, the `(1, 1, 1)` Orthogonal Triad Winner, and Hard-Core Nuclear Repulsion (`exponential_wall_bypasses_and_triad_winner`).**
    Whenever compression or winding hits an exponential or superlinear resistance
    wall, **A0 (Path of Least Resistance)** bypasses the wall by opening a new
    orthogonal degree of freedom or crossing `0` into `-`:
    1. *1D Compression Wall $\to$ Helical Buckling Plateau:* Straight 1D compression
       costs $E_{\rm 1D}(\chi)=e^{+\chi}\to\infty$; buckling sideways into a helix
       bypasses the wall and caps the cost at the finite plateau
       $E_{\rm helix}^\star(\chi)=\tfrac32-\tfrac12e^{-2\chi}\to 1.5$.
    2. *3D `+` Compression Wall $\to$ Signed `0`-Crossing Plateau:* Staying in `+`
       ($\theta=0$) costs $E_{\rm stay,+}(\chi)=2e^{+\chi}+2e^{-\chi}\to\infty$;
       opening the `0` valve into a `-` child sphere caps the 5-variable cost at
       the finite plateau $E_{\rm 5\text{-var}}^\star(\chi)=4-\tfrac92e^{-\chi}+O(e^{-2\chi})\to 4.0$.
    3. *Single-Axis Winding Wall $N^{3/2}>N$ $\to$ Orthogonal `(1, 1, 1)` Triad Ground State & Nuclear Hard Core:*
       Piling $N\ge 2$ windings onto a single axis costs $N^{3/2}>N$, so windings
       split across the 3 orthogonal axes $(n_1,n_2,n_3)\in\mathbb Z_{\ge 0}^3$
       with $E^\star(n_1,n_2,n_3)=4(PV_0)^{1/4}m_0^{3/4}\sqrt{n_1n_2n_3}$. States
       `(1,0,0)` ($N=1$) and `(1,1,0)` ($N=2$) have zero 3D volume ($V^\star=0$);
       **`(1,1,1)` ($N=3$) is the unique minimum-cost 3D ground state** ($E^\star=4.0$,
       cost per winding $4/3\approx 1.3333$, round sphere $1:1:1$); adding a 4th
       winding forces prolate **`(2,1,1)`** ($4:1:1$, $E^\star=4\sqrt{2}\approx 5.6569$,
       cost per winding $\sqrt{2}\approx 1.4142>4/3$). Two `(1,1,1)` triads attract
       at medium range ($d>2R^\star$) via shared surrounding flow ($-2m_0/d$), reaching
       a bound minimum $E_{\rm contact}=7.0$ at $d=2R^\star$, but **strongly repel at
       short range** against merging into `(2,2,2)` because $E^\star(2,2,2)=8\sqrt{2}\approx 11.3137>2E^\star(1,1,1)=8.0$
       (hard-core barrier $+3.3137$).

0. **Two-force zero-crossing across `0` (`two_force_zero_crossing`).**
   When inward compression $F_{\rm in}$ wins against outward resistance
   $F_{\rm out}$ on a sphere of signed coordinate $R(\tau)$ with $R(0)=R_0>0$,
   the net force $F_{\rm net}=F_{\rm in}-F_{\rm out}>0$ drives $\ddot R=-F_{\rm net}<0$,
   so $R(\tau)=R_0+v_0\tau-\tfrac12F_{\rm net}\tau^2$. Comparing the three
   candidate paths at the horizon $R=0$:
   - Busting back into the parent `+` sphere ($\Delta R>0, R>0$) is blocked
     uphill against the winning force $F_{\rm net}>0$.
   - Stopping at $R=0$ requires an infinite deceleration impulse to destroy
     the horizon kinetic budget $K(\tau_0)=\tfrac12\dot R(\tau_0)^2>0$ in zero distance.
   - Crossing `0` into $R<0$ (`-` region) has zero impulse cost, flips the 3D
     Jacobian determinant $\det[(R/R_0)I_3]=(R/R_0)^3$ from `+` to `-` through `0`,
     and converts the inward implosion velocity $\dot R<0$ directly into the
     outward expansion rate $\dot r_-=-\dot R>0$ of a new sphere $r_-=|R|=-R>0$
     growing from scratch without reversing the driving force.

1. **Helical buckling under axial compression (`helical_buckling_under_compression`).**
   A 1D seed compressed axially to $\lambda=e^{-\chi}$ that turns laterally at
   azimuthal speed $v_\perp=r_\perp\omega$ has Euclidean arc-length factor
   $\ell=\sqrt{\lambda^2+v_\perp^2}$. Minimizing the Cauchy-Schwarz stretch
   cost plus lateral enclosing confinement,
   $$
   E(v_\perp;\lambda)=\frac{Q_0}{\sqrt{\lambda^2+v_\perp^2}}+\frac12P_\perp v_\perp^2,
   $$
   shows that at $\lambda_{\rm crit}=(Q_0/P_\perp)^{1/3}$ ($\chi_{\rm crit}=-\tfrac13\log(Q_0/P_\perp)$),
   the straight seed ($\omega=0$) loses transverse stability
   ($\partial^2E/\partial v_\perp^2\vert_{v_\perp=0}=P_\perp-Q_0/\lambda^3<0$) and
   buckles into a right- or left-handed helix with exact turning rate
   $$
   r_\perp\omega^\star(\chi)=\sqrt{\lambda_{\rm crit}^2-e^{-2\chi}},
   \qquad
   \alpha^\star(\chi)=\frac{\omega^\star(\chi)}{2\pi},
   $$
   lowering total cost from $Q_0/\lambda$ to $\tfrac32P_\perp\lambda_{\rm crit}^2-\tfrac12P_\perp\lambda^2$.

2. **Coupled `+`/`0`/`-` damped oscillation and asymmetric remainder (`damped_horizon_oscillation`).**
   Coupling the two-sided horizon displacement $q(\tau)=A_0e^{-\gamma\tau}\sin(\omega\tau)$
   to the receiving geometry budget via $\dot E_{\rm geom}=2\gamma\dot q^2$ gives
   $$
   \ddot q+2\gamma\dot q+(\gamma^2+\omega^2)q=0,
   \qquad
   E_{\rm wave}(\tau)+E_{\rm geom}(\tau)=\frac12A_0^2\omega^2=\text{const.}
   $$
   Because each `+` lobe ($q>0$) precedes its compensating `-` lobe ($q<0$) by
   $\pi/\omega$, the peak-amplitude ratio, single-cycle area ratio, and infinite
   cumulative area ratio all share the exact invariant asymmetry
   $$
   \eta=\frac{I_+-I_-}{I_++I_-}
   =\frac{I_0^+-|I_1^-|}{I_0^++|I_1^-|}
   =\tanh\!\left(\frac{\pi\gamma}{2\omega}\right)>0,
   \qquad
   \Delta I=I_+-I_-=\frac{A_0\omega}{\gamma^2+\omega^2}>0.
   $$

   The same oscillator fixes a second signed readout. Write
   $\rho=\exp(-\pi\gamma/\omega)$ for the amplitude multiplier over a half
   cycle. At every zero crossing the oscillator energy is purely kinetic
   and multiplies by $\rho^2$. Consequently the first two lobes deposit
   $W_+=E_0(1-\rho^2)$ and $W_-=\rho^2 W_+$ in the receiving geometry.
   For positive damping, each adjacent pair has the same normalized work
   imbalance, as does the cumulative sum starting with the positive lobe:

   $$
   \eta_{\rm work}=\frac{W_+-W_-}{W_++W_-}
   =\tanh(\pi\gamma/\omega)
   =\frac{2\eta_{\rm area}}{1+\eta_{\rm area}^2}.
   $$

   Thus the recorded $K=7$ amplitude-area residue $5.955329\times10^{-10}$
   corresponds to a work residue $1.191066\times10^{-9}$ for the same
   oscillator. A spectral half-cycle exponent $-(\pi/2)\lambda/\lambda_{\max}$
   would instead give the smaller amplitude-area expression
   $\tanh[(\pi/4)\lambda/\lambda_{\max}]$. The exponent used above preserves
   the existing numerical amplitude-area curve. The oscillator supplies
   a signed amplitude residue and a signed geometric work residue. The
   missing-antimatter interpretation still requires a rule relating one
   of these quantities to a matter count and its normalization.

   If a two-daughter split deposits its entire available release into
   this oscillator, and the oscillator fully relaxes, the signed difference
   of the deposited work is $E_N^\star(1-1/\sqrt{2})\eta_{\rm work}$.
   This states the budget for that proposed coupling. The split calculation
   leaves the excitation amplitude, damping ratio and phase at release to
   the transition dynamics.

3. **Goldilocks stabilization (`buckling_to_goldilocks_calculation`).**
   As $E_{\rm geom}(\tau)$ populates the child's isotropic 3D stretch budget
   against enclosing load $PV_0$, its least-resistance holding scale expands as
   $$
   \lambda_{\rm child}(\tau)=\left(\frac{E_{\rm geom}(\tau)}{3PV_0}\right)^{1/4}
   \xrightarrow{\tau\to\infty}
   \lambda_{\rm Goldilocks}=\left(\frac{E_{\rm total}}{3PV_0}\right)^{1/4}.
   $$
   For $\omega=2\pi$, $\gamma=\tfrac12\ln2$ and $A_0=1$, the child reaches
   $98.40$% of its Goldilocks equilibrium scale $\lambda_{\rm Goldilocks}\approx1.6015928510$
   by offset $\tau=4$, with zero conservation residual throughout.

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --buckling-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --atlas-only
```

### Exact relation between the three budget shares

The normalized five-variable minimum gives an exact relation among the
stretch, load and filament shares across its entire compression family.
Write these shares as $D,M,V$ respectively, and put $B=M-D/3$ for the
radial infall share. Let $x=e^{-\chi}$, $u=\sqrt{\theta^\star}$ and
$a=1-u^2$. The minimum obeys $3u=x+2a/x$. With common total cost $E$,

$$
D=3u^3/E,\qquad B=2ax/E,\qquad V=2a^2/(xE).
$$

Multiplying the stationarity equation by $a/E$ gives
$(V+B/2)/D=a/u^2$. Multiplying the radial and filament shares gives
$BV/D^2=(4/9)(a/u^2)^3$. Eliminating the ratio proves

$$
\boxed{9DVB=4(V+B/2)^3,\qquad B=M-D/3.}
$$

The identity uses the same cost and stationarity equation as the existing
split and holds at any compression in the declared family. It can also be
read in reverse at any finite point with positive shares:

$$
u^2=\frac{D}{D+V+B/2},\qquad e^{-\chi}=\frac{3uB}{2V+B}.
$$

In deep compression, $a=3x/2+O(x^2)$ and $E=4+O(x)$. Consequently
$V=9x/8+O(x^2)$ and $B=3x^2/4+O(x^3)$, giving the intrinsic scaling

$$
M-D/3=\frac{16}{27}V^2+O(V^3).
$$

Thus the excess above one third of the stretch share falls quadratically
with the filament share. The exact identity applies before taking that
limit. The three budget shares encode the compression that generated them.

The existing `nested_zero_crossing_map(1, chi)` calculator gives a maximum
absolute cubic residual below $7\times10^{-17}$ at the sampled compressions
0, 1, 2, 3.0669683518, 4, 6, 8 and 10. At the existing $V=e^{-3}$ point,
$B=.00162558745$ and the inverse recovers $\chi=3.06696835185$. At $\chi=10$,
$B/V^2=.59264976$ approaches the exact limit $16/27=.59259259$. The
equation concerns the atlas's geometric budget shares and applies equally
to points chosen without the filament-floor condition.

### Spherical fragmentation preserves area and releases budget

Use the existing bubble cost $Q_n(R)=m_0n^2/R+PV_0R^3$. Compare relaxed,
separated spheres at the same pressure and couplings so their costs add.
The stationary radius gives $(R_n^\star)^2=\sqrt{m_0/(3PV_0)} |n|$ and
$E_n^\star=4PV_0(R_n^\star)^3$. For a positive winding $N$ split into positive
integers $n_i$ with $\sum_i n_i=N$, these identities imply

$$
\sum_i4\pi (R_{n_i}^\star)^2=4\pi (R_N^\star)^2,\qquad
\frac{E_{\rm daughters}}{E_N^\star}
=\frac{\sum_i (R_{n_i}^\star)^3}{(R_N^\star)^3}
=\sum_i\left(\frac{n_i}{N}\right)^{3/2}<1
$$

for any nontrivial split. Thus the equilibrium surface area is additive
in same-sign winding. Fragmentation reduces the total enclosed volume
and minimized energy by the same fraction. The energy difference is
available to other degrees of freedom in a conserving evolution.

For a fixed number $k$ of integer daughters the least cost occurs when
their windings differ by at most one. Indeed the forward differences of
$n^{3/2}$ strictly increase, so moving a unit from a winding at least two
larger than another lowers the sum. When $k$ divides $N$, each radius is
$R_N^\star/\sqrt{k}$ and the available release fraction is $1-1/\sqrt{k}$.
For $N=8$, successive balanced divisions into 2, 4 and 8 daughters give
cumulative releases of 29.2893%, 50%, and 64.6447%. The smallest positive
winding in this sector is one. The fraction is independent of pressure
in this comparison, while the absolute release scales as $P^{1/4}$ at
fixed couplings and winding.

This is an equilibrium budget law for admitted same-sign partitions.
The onset, separation dynamics and spectrum of released waves require
a transition calculation. The opposite-winding annihilation construction
below addresses cancellation of signed winding. The Navier–Stokes
[concentrating vortex example](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf)
provides a separate mathematical example of increasing local intensity
with bounded total kinetic energy. The fragmentation law above follows
directly from the atlas bubble functional.

### A downhill winding-two neck and split

Extend the bubble's exterior radial phase $n/r$ to two unit centers
separated by $2c$. Use the cavity $\Phi\geq 2/R$ with
$\Phi=1/|\mathbf{x}-c\mathbf e_z|+1/|\mathbf{x}+c\mathbf e_z|$.
This specifies a candidate phase and shape family. The unit labels are
carried as fluxes of the exterior phase gradient. For $c=0$ the boundary
is the winding-two sphere. The dimensionless separation is $s=c/R$.
The midplane neck has radius $R\sqrt{1-s^2}$ for $s<1$. The components
touch at $s=1$ and separate for $s>1$. Each separated component encloses
one center. The constant boundary phase keeps the exterior field
compatible with both components.

Green's identity on the exterior gives
$\int |\nabla\Phi|^2 d^3x=(2/R)(8\pi)$, since the total normal
flux is $8\pi$ and the contribution at infinity vanishes. With the
normalization $m_0/(4\pi)$ this is the original stretch cost $4m_0/R$.
Only the cavity volume changes at fixed $R$.

For the volume, use prolate coordinates $\xi=(d_1+d_2)/(2c)$ and
$\eta=(d_1-d_2)/(2c)$, where $d_i$ are distances to the centers. The
volume element is $c^3(\xi^2-\eta^2)d\xi d\eta d\varphi$.
Writing $q=(1+\sqrt{1+4s^2\eta^2})/2$, the boundary is $s\xi=q$.
Integration gives the volume factor relative to a sphere of radius $R$,

$$
v(s)=3\int_{\eta_0}^1
\left(q^2-\frac23 q^3-s^3/3+s^3\eta^2\right)d\eta,
$$

where $\eta_0=0$ for $s\leq1$ and $\eta_0=\sqrt{1-1/s}$ otherwise.
At $s=0$ use the continuous limit $v=1$. On the connected part the
derivative of the last two terms integrates to zero. The first two have
derivative $-2s^2\eta^2\partial_s q<0$ away from the axis of integration
$\eta=0$. Thus volume decreases strictly before the pinch. After the
pinch, translate coordinates with one center and consider its half-space.
Increasing the separation lowers the other center's potential in the
common half-space. In the newly added strip between the old and new
midplanes, both axial distances are at least $R$ so that the potential
cannot exceed $2/R$. The translated lobes therefore shrink as separation
increases, converging to two spheres of radius $R/2$ and $v=1/4$.

At fixed pressure, minimizing $4m_0/R+PV_0R^3v(s)$ over the overall
scale gives $R^\star(s)/R^\star(0)=v(s)^{-1/4}$ and
$E^\star(s)/E^\star(0)=v(s)^{1/4}$. This relaxed path is downhill throughout.
At the pinch $v=.6331341144$ and $E^\star/E^\star(0)=.8920186380$, leaving
10.7981% available. In the separated limit the available fraction is
$1-1/\sqrt2$, and each daughter has radius $R^\star(0)/\sqrt2$.
Near the original sphere, $v(s)=1-3s^4/5+O(s^6)$, so the first energy
decrease is quartic in separation. Pressure sets the scale and absolute
budget. The downhill sign holds for every positive pressure in this
family.

Run `python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --neck-split-only`
for the volume integral and budget along the path. The path gives the
energy available to outgoing degrees of freedom as the shapes relax.
A time evolution, outgoing spectrum and pressure-triggered onset remain
to be derived. The phase and geometry in this split are symmetric under
interchange of daughters. The damped-oscillation residue interpreted as
missing antimatter is a separate NSC hypothesis. Its identification with
a matter count needs a signed observable and a conversion law.

### The split launches its own phase waveform

The downhill neck path now has a coupled dynamical calculation. Use the
same exterior phase action with kinetic and gradient terms
$\int(\Phi_t^2-|\nabla\Phi|^2)d^3x/(4\pi)$, in units $m_0=1$ and
initial sphere radius one. The existing pressure normalization is then
$PV_0=4/3$. The near field follows the two-center equipotential family
and its relaxed radius. This collective-coordinate reduction is a stated
approximation. The outer field satisfies its wave equation with outgoing
conditions and acts back on the changing shape.

Use $q=c^2$ as the shape coordinate. The first phase departure about the
sphere is $2qP_2(\cos\theta)/r^3$. Direct integration of its squared
time derivative gives total quasistatic inertia $8/15$. The potential
is $16/3-(4/5)q^2+O(q^3)$. Thus a small admitted deformation has an
unstable direction with a definite restoring-force sign. The initial
reference takes $s=c/R=.05$ and zero velocity. The initial outer field
is the static field of that deformation. The seed is explicit preparation
data and sets the timing of departure from the exact sphere.

At a matching sphere $r=B$ the even angular coefficients are
$g_l(q)=2q^{l/2}/B^{l+1}$. The total flux fixes a constant monopole,
and reflection symmetry removes odd angular modes. The retained dynamical
modes satisfy the radial wave equation outside B. For Laplace variable p,
their exact normal response is

$$
D_l(p)=p+\frac{l+1}{B}
-p\frac{P_l'(pB)}{P_l(pB)},\qquad
P_l(z)=\sum_{k=0}^{l}\frac{(l+k)!}{k!(l-k)!2^k}z^{l-k}.
$$

Here P_l is the outgoing radial polynomial, distinct from the Legendre
angular polynomial. If its roots are rho_j, states h_j satisfy
$\dot h_j=(\rho_j/B)h_j+g_l$. The normal response is
$D_lg_l=\dot g_l+g_l/B-\sum_j\rho_jh_j/B^2$.
Initial h_j=−Bg_l(0)/rho_j recovers the prepared static exterior.
This construction supplies causal memory and energy transfer without
selecting a damping coefficient or a sinusoidal output.

The projected near-field inertia M(q) is the integral of
$(\partial_q\Phi)^2/(2\pi)$ over the near exterior. Untransported higher
harmonics retain their quasistatic inertia. The retained exterior static
energy is subtracted once from the full potential U(q). Variation gives

$$
M\ddot q+\tfrac12M'\dot q^2+U'
+\sum_l\frac{2B^2}{2l+1}g_l'
\left[D_lg_l-\frac{l+1}{B}g_l\right]=0.
$$

The exterior work is $\sum_l 2B^2(D_lg_l)\dot g_l/(2l+1)$. Including
that field energy closes the combined budget. Independently reconstructing
the wave at finite radii and integrating its kinetic and gradient energies
checks this account. The source of the phase pulse is the evolving shape.

At the reference B=2.5, l<=8 calculation, the neck pinches at offset
3.960294 and reaches the declared separated state s=1.2 at 4.264066.
The potential falls from 5.333328 to 4.513152. The projected kinetic
budget is .726247. A receiver at radius 10 records positive polar phase
departure .0506724 and negative equatorial departure −.0171683 at the
last arrival in this window. Its integrated outward energy is .0553464.
All these quantities use the source's normalized length and offset units.
The signed angular mean of the emitted departure is zero because every
varying mode has positive angular degree. The two antipodal poles have
equal signals by the source symmetry.
This angular integral and the earlier temporal lobe-area residue use
different readout operations. This calculation assigns no value to the
missing-antimatter identification of that temporal residue.

Six bounded controls compare shape quadrature, angular cutoffs and the
size of the quasistatic near field. At fixed B=2.5 the last two shape
discretizations change the l<=6 receiver energy by about 2.13e-8, while
adding l=8 changes it by 5.99e-6. The independent exterior energy integral
differs from boundary work by less than 8e-11 across the cases. Changing
B from 2.3 to 3 changes receiver energy from .0717211 to .0314597.
This latter variation measures dependence on the near-field reduction
and dominates the numerical differences. It must accompany any quoted
wave amplitude. A fully evolving near phase field remains to be calculated.

The level set itself continues through the pinch with a smooth local
phase normal form. Writing rho for cylindrical distance and z for axial
distance, near the origin at s=1 its equation is
$2z^2-\rho^2=2R^2(s-1)$ to leading order. Neck radius has a square-root
coordinate behavior, while the underlying phase and its shape derivative
remain finite. The q coordinate carries the evolution across that event.

The [compact result](../results/summaries/nsc-neck-waveform-v1.json) and
`python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --split-waveform-only`
give the finite waveforms and all six controls. The calculation reaches
two separated phase lobes. Relaxation to the infinitely separated energy
limit, an inherited seed mechanism and a signed matter-count readout
remain open. In the calculated finite window the pole records a pulse,
and its opposite angular signs cancel in the full-sphere phase integral.

### Whole-exterior causal response at the start of the split

The matching-radius sensitivity can be removed exactly in the quadratic
small-deformation problem. Use the same phase kinetic term over the entire
exterior r>=1. Let the quadrupole be f(r,t)P2(cos(theta)), with boundary
shape coordinate q=f(1,t)/2. The previous static path gives potential
16/3−4q²/5+O(q³). A static exterior quadrupole f=2q/r³ has gradient
energy 12q²/5, fixing the remaining quadratic shape term to −16q²/5.
The resulting quadratic action is

$$
L_2=\frac15\int_1^\infty
[r^2(f_t^2-f_r^2)-6f^2]dr+\frac45 f(1,t)^2.
$$

Variation gives the exterior wave equation and the free shape condition

$$
f_{tt}-f_{rr}-2f_r/r+6f/r^2=0,\qquad f_r(1,t)+4f(1,t)=0.
$$

This is the B->1 limit of the preceding quadratic model. The matching
sphere and projected near-field inertia have disappeared from this onset
calculation. The same canonical phase inertia and static shape cost remain
its premises. Nonlinear boundary motion requires its further variation.

The outgoing radial solution at Laplace rate p is proportional to
$\exp[-p(r-1)](p^2r^2+3pr+3)/r^3$. Substituting it into the shape condition
gives p³−3p−3=0. Its positive root is 2.1038034027, compared with
1.7630684000 in the B=2.5 near-field approximation. Its other roots are
−1.0519017014 +/- .5652358517 i. The same causal boundary response
therefore contains growing separation and a decaying oscillatory transient.
Both rates follow from the field and geometric cost.

Retain the initial static perturbation f(r,0)=2q0/r³ and zero initial
bulk velocity. Releasing the shape produces a finite-energy characteristic
corner. Its boundary velocity immediately after release is q'(0+)=q0,
as fixed by the solution. This pointwise boundary value coexists with the
zero-velocity initial data in the open bulk. The earlier finite near-field
inertia had instead imposed zero initial shape velocity.

The exact Laplace boundary response is

$$
\widehat q(p)=q_0\frac{p(p+1)}{p^3-3p-3}.
$$

If lambda_j are the three roots and D'(lambda)=3lambda²−3, the complete
field after the front t=r−1 is

$$
f(r,t)=\frac{2q_0}{r^3}\sum_j
\frac{\lambda_j^2r^2+3\lambda_jr+3}{\lambda_jD'(\lambda_j)}
e^{\lambda_j(t-r+1)}.
$$

Ahead of that front the field retains its initial static value. The
quadratic conserved energy is the spatial field energy minus 4f(1,t)²/5.
Direct integration gives the constant −4q0²/5. Its negative sign is the
second-order departure from the sphere's positive total budget. This
balance checks the kinetic, wave and shape terms independently.

For the existing s=.05 seed, q0=.0025000046764, the shape coordinate
reaches .0133275776068 at offset1. The largest sampled energy discrepancy
is 4.7e-18 and the boundary residual is below 1.4e-17. These are numerical
checks of the exact quadratic solution. Run the existing CLI with
`--split-onset-only` to reproduce them. The nonlinear waveform through
pinch still retains the preceding near-field approximation. This exact
onset response supplies a reference for extending that calculation.

### Inherited daughters and the direction of the next split

The first split's endpoint s=1.2 now supplies the next geometry, its
velocity and its outgoing memory. The identical phase can be written as
two sources of winding two or as four coincident unit sources. Let the
centers of the two groups lie at z=+s and z=−s in units of the current
equipotential scale R. Splitting each group in the transverse x direction
places unit sources at (x,z)=(+b,+s),(-b,+s),(+b,-s),(-b,-s). The
boundary is the same level Phi=4/R. At b=0 this reproduces the first
field exactly, including its derivatives with respect to the original
shape coordinate. The new splitting velocity is initially zero, so its
addition preserves the existing field's time derivative.

The calculator retains the first split's shape velocity and every outgoing
memory coefficient. The existing winding-two reference converts to winding
four at the same pressure by multiplying lengths and offsets by sqrt(2),
phase by sqrt(2), and total energy by 2sqrt(2). The squared separation and
outgoing memory each multiply by two. The record includes that full scaled
handoff explicitly. Shape values below use the parent's own reference units.
In a moving local chart centered on one daughter,
the shared potential has a self term and a sibling term. If their center
distance is d, the latter has the exact local expansion

$$
\frac{n_s}{|d\mathbf e_z+\mathbf x|}
=n_s\sum_{l=0}^{\infty}\frac{(-1)^l\rho^l}{d^{l+1}}P_l(\cos\theta),
\qquad \rho<d.
$$

Its time derivative follows from the inherited d and d_dot. The same
parent state therefore supplies the daughter's angular preparation. Direct
projection of the complete boundary at s=1.2 gives a normalized dipole
−.07340185 and quadrupole +.02750552, with higher modes retained. The
quadrupole radius coefficient is .02060825 in parent reference units,
and its inherited offset derivative is −.05818555. These are evaluated
source data from the first split, with the previously declared near-field
approximation carried along.

The next least-resistance direction also follows from this field. Write
zeta=b². The volume derivative for an axial second split is the surface
integral of its phase perturbation divided by the gradient magnitude.
The phase perturbation is one half of the axial second spatial derivative
of the original two-center phase. Harmonicity makes the trace of this
directional tensor zero. Axial symmetry therefore gives

$$
\partial_\zeta v(\beta)
=A(s) \frac{3\cos^2\beta-1}{2},
$$

where beta is the angle to the first splitting axis. At the inherited
endpoint A(1.2)=.152756298060, independently checked by a finite variation
of the four-source volume. The fixed-pressure relaxed energy is proportional
to v^(1/4), so its normalized directional derivatives are +.0630224551
along the old axis and −.0315112276 in the transverse plane. That plane
minimizes the directional cost. Its azimuth remains degenerate under the
original rotational symmetry. The calculation fixes a plane, and an x-axis
representative is used for the subsequent four-center path.

The second pair of necks reaches the level-set critical points at
b/R=.627460420792 and z/R=+/-1.15970809670. At these points Phi=4/R
and the phase gradient vanishes. Finite source distances keep the phase
regular. The resulting four lobes each enclose one unit source. The full
shared phase includes the cross terms throughout the path.

| Stage in the same four-source family | Available equilibrium energy / original winding-four budget |
| --- | ---: |
| First split's s=1.2 endpoint, b=0 | .1537839232 |
| Second necks reach zero | .2554477258 |
| Four separated lobes at b/R=1 | .3386902228 |
| Fully separated four-unit reference | .5 |

The finite volumes use two angular quadratures and the coincident-pair
limit agrees with the earlier analytic two-center integral. The listed
volume comparisons are numerical and the local direction identity is
analytic. A monotonicity proof between all finite four-source cuts has
not been supplied.

At the same fixed pressure, the equilibrium scale map for a binary child
has a=1/sqrt(2). It acts on the local phase and its coordinates as

$$
\mathbf x'=a\mathbf x,\qquad t'=at,\qquad
\Phi'(\mathbf x',t')=a\Phi(\mathbf x,t).
$$

First spatial and temporal derivatives are invariant. The source flux
scales as a², the energy as a³, and the pressure is invariant. The
normalized wave operator R²(partial_t²−Delta) is inherited. Frequencies
in a common unit increase by sqrt(2). Two binary generations give four
units with total reference energy one half of the original. These are
same-pressure similarity relations. The finite sibling fields and their
time-dependent corrections remain in the inherited state above.

`--cascade-inheritance-only` on the existing offset CLI reproduces the
[handoff and direction calculation](../results/summaries/nsc-inherited-daughter-v1.json).
It extends the first wave calculation to the second generation's initial
field, shape and static continuation, and evaluates both the **three-axis
orthogonal cascade theorem ($\mathbf e_z\to\mathbf e_x\to\mathbf e_y$)**
and the **joint $2\times 2$ kinetic handoff**:

1. **Three-axis orthogonal cascade from $\Delta\Phi=0$ (`three_axis_orthogonal_cascade`).**
   On any separated lobe of the four-source level set $\Phi=4/R$ at
   $(\pm b,0,\pm s)$, splitting each unit source by $\pm\sqrt{\zeta_i}$
   along axis $i\in\{x,y,z\}$ perturbs the phase by
   $\delta\Phi=\frac12\zeta_i\partial_{ii}\Phi+O(\zeta_i^2)$. Because $\Phi$
   is harmonic away from the sources ($\Delta\Phi=\partial_{xx}\Phi+\partial_{yy}\Phi+\partial_{zz}\Phi=0$),
   the directional volume-derivative trace vanishes identically at every
   separation:
   $$
   \partial_{\zeta_x}v+\partial_{\zeta_y}v+\partial_{\zeta_z}v=0.
   $$
   At $s=1.2$, for all $b\ge 0.7$ after the second neck pinches ($b_{\rm pinch}=0.62746042$),
   the third orthogonal axis $\mathbf e_y$ is the unique global minimum
   of the directional volume and energy derivatives. By $b=1.0$, **both**
   already-split axes in the $xz$-plane are uphill
   ($\partial_{\zeta_x}v=+0.15537634>0$, $\partial_{\zeta_z}v=+0.03404474>0$),
   leaving $\partial_{\zeta_y}v=-0.18942108<0$ as the **only downhill axis in 3D**.
   At the symmetric four-lobe square $b=s=1.2$, four-fold in-plane symmetry
   ($\partial_{\zeta_x}v=\partial_{\zeta_z}v=+0.05787268$) and $\Delta\Phi=0$
   fix the exact harmonic ratio $-2$:
   $$
   \partial_{\zeta_y}v(s,s)=-2 \partial_{\zeta_x}v(s,s)=-2 \partial_{\zeta_z}v(s,s)=-0.11574536<0.
   $$
   Together with the first-to-second split ratio $-1/2$, harmonicity and
   least resistance force the cascade through all three mutually orthogonal
   axes $\mathbf e_z\to\mathbf e_x\to\mathbf e_y$ ($1\to 2\to 4\to 8$).
   Continuing along $\mathbf e_y$ on $\Phi=8/R$ pinches the four lobes into
   eight separated 3D octant lobes at third-split half-separation
   $c_y/R=0.360449324968$ (for $b=1.0, s=1.2$) and $c_y/R=0.347227979949$
   (for $b=s=1.2$).
2. **Joint $2\times 2$ kinetic matrix and outgoing-memory reaction at handoff (`joint_two_coordinate_kinetic_handoff`).**
   With $q_1=c^2$ and $q_2=\tilde b^2$, the transverse phase derivative at
   $q_2=0$ is $\partial_{q_2}\Phi\vert_{q_2=0}=\frac12\partial_{xx}\Phi_0=-\frac14\partial_{zz}\Phi_0+\frac14(\partial_{xx}-\partial_{yy})\Phi_0$.
   Because $\partial_{zz}\Phi_0\equiv\partial_{cc}\Phi_0$ for the two-center
   potential, the axisymmetric ($m=0$) exterior multipole coefficients
   satisfy the exact identity
   $$
   \left.\frac{\partial g_{l,0}}{\partial q_2}\right\vert_{q_2=0}
   =-\frac{l-1}{2}\frac{\partial g_{l,0}}{\partial q_1}
   =-\frac12 l(l-1)q_1^{l/2-1}B^{-(l+1)},
   $$
   while the quadrupolar ($m=2$) component expands as
   $\sum_{\text{even }l\ge 2}\frac12 q_1^{l/2-1}r^{-(l+1)}P_l^2(\cos\theta)\cos(2\varphi)$.
   Summing the inner prolate quadrature and outer solid-harmonic series gives
   the $B$-independent full quasistatic kinetic matrix at $s=1.2, b=0$
   (matching to $10$ digits across $B\in\{2.5,3.0,4.0\}$):
   $$
   M_{\rm full}=\begin{pmatrix}0.2074628302 & -0.0054817211 \\ -0.0054817211 & 0.6751131619\end{pmatrix},
   \qquad \det M_{\rm full}=0.1400308435>0,
   $$
   and at the $B=2.5, l_{\max}=8$ transported reduction gives
   $M_0=\begin{pmatrix}0.1695831644 & 0.0176210767 \\ 0.0176210767 & 0.6564127651\end{pmatrix}$
   ($\det M_0=0.1110060515>0$). Even with zero initial transverse velocity
   $\dot q_2(0)=0$, the inherited axial velocity $\dot q_1=2.92661678$ gives
   nonzero transverse conjugate momentum $p_2(0)=M_{0,21}\dot q_1=+0.05157014>0$,
   and the first split's outgoing wave memory exerts a direct transverse
   reaction $R_2=-0.06917181<0$ that reinforces the downhill potential
   slope $-\partial_{q_2}U=+0.12034460>0$ for a combined positive transverse
   drive $-\partial_{q_2}U-R_2=+0.18951641>0$.

## Retained foundation control: opposite windings, 2 November 2026

The geometry-only calculation is:

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --winding-only
```

It does not evaluate the supplied parent action, select an observational
target or write another result file. The older parent comparison below is
retained as a conditional construction, outside these premises.

### Closed orientations and their contact

Take a circle of radius R and a normalized drawing coordinate x in [0,1].
The two orientations are z_+(x)=R exp(2 pi i m x) and
z_-(x)=R exp(-2 pi i m x). Closure requires integer m. Their degrees are
+m and -m. Define the geometric squared variation, explicitly in this
drawing coordinate, by Q[z]=(1/2) integral |dz/dx|^2 dx. Each loop has
Q=2 pi^2 R^2 m^2: reversing orientation preserves its positive value.
This functional has no supplied conversion to joules or particle mass.

Join the loops at one contact and use a combined coordinate x in [0,2].
Let f(x)=x for x<=1 and f(x)=2-x for x>=1. The joined path

$$
Z_h(x)=R\exp[2\pi i m(1-h)f(x)],\qquad 0\leq h\leq1,
$$

is closed, continuous and never passes through zero. It is piecewise
smooth with finite squared variation at the join. At h=0 it follows the
two loops; at h=1 it is constant. During contraction the two halves need
not remain separately closed: allowing their common contact to reconnect
is an explicit operation, not a consequence of the offset update alone.

Its variation is exactly

$$
Q[Z_h]=4\pi^2R^2m^2(1-h)^2.
$$

Cauchy-Schwarz proves this is the least variation at each h for a fixed
circle radius, fixed endpoint lift and midpoint lift: on each length-one
half, integral (phase derivative)^2 >= (phase change)^2, with equality
for the constant phase slopes used above. This is a sharp constrained
minimum, not a claim that nature selects the constraints.

For R=m=1:

| Contraction h | Q |
| --- | --- |
| 0 | 39.4784176044 |
| 1/4 | 22.2066099025 |
| 1/2 | 9.8696044011 |
| 3/4 | 2.4674011003 |
| 1 | 0 |

This supplies a regular geometric continuation from opposite winding to
zero variation. The parameter h is not yet the offset's selected history.
No received energy, new region or outgoing signal is inferred from the
decrease without a transfer law.

### What forces the 1:1 and 180-degree pair?

Suppose exactly two outputs have positive scalar amounts a,b, unit
directions u,v, total a+b=Q_out, and balanced weighted directions
a u+b v=0, with no remaining receiver. Taking norms gives a=b, hence
u=-v. Therefore each amount is Q_out/2: **ratio 1, angle 180 degrees**.
This is an exact geometric consequence of those explicit output rules.

Balance alone does not select two outputs. Three equal outputs separated
by 120 degrees, or four separated by 90 degrees in a plane, also balance.
The calculation constructs all three cases. Nor does opposite orientation
select a contact rule: multiplying conjugate unit-circle maps gives a
constant, adding them gives 2 cos(phase), and retaining both keeps two
loops. The offset recurrence alone specifies none of these interactions.

**Next missing connection:** the contact composition/evolution law that
selects output channels and transfers the lost geometric variation.
Deriving that law, its coupling to the offset and an instrument readout is
what can distinguish this construction from a change of names. No 511 keV
scale, electron charge or matter excess has been derived here. Integer
winding is conditional on the circle target; it is not a charge assignment.

## Conservative contact evolution, 2 November 2026

The preceding contraction supplies a path but not its motion or a receiver.
Three distinct choices were considered: prescribing h directly, descending
Q by a geometric gradient, and giving the same circle phase conservative
inertial motion. The first prescribes the history; the second dissipates Q
without specifying a receiver. The third permits a closed transfer ledger
and an explicit open/closed comparison, so it is the worked candidate here.

### The added rule and its exact solution

Write Z=R exp(i psi). Its squared spatial variation is
Q=(R^2/2) integral psi_x^2 dx. **Add** geometric inertia using the same
circle metric: K=(R^2/2) integral psi_tau^2 dx. Stationarity of the
integral of K-Q, with fixed variation at the endpoints, gives

$$
\psi_{\tau\tau}=\psi_{xx}.
$$

This equation is derived from that newly declared local quadratic
inertial rule, not from the offset recurrence or from boundary minimization
alone. A ratio of inertia and spatial coefficients can be absorbed into
the dimensionless tau unit; no physical propagation speed or duration is
thereby predicted. No new field species, potential or observed scale is used.

Let g(x)=2 pi m f(x), with the triangular f above. Start with psi(x,0)=g(x)
and psi_tau(x,0)=0. Compare two drawing-coordinate continuations: extend g
by zero outside [0,2], or identify that interval periodically. These are
different domain assumptions for the SAME circle-valued target.
For either continuation the exact weak solution is

$$
\psi(x,\tau)=\tfrac12[g(x-\tau)+g(x+\tau)].
$$

The moving corners have finite energy and need no singularity cutoff.
Derivatives exist almost everywhere; integrals are evaluated separately
on their exact intervals. The two characteristics are derived from the
second-order equation. Equal amounts follow from zero initial velocity;
they are not selected by opposite winding alone.

### Where the variation goes

The local conserved density and flux are

$$
e=\frac{R^2}{2}(\psi_\tau^2+\psi_x^2),\qquad
j=-R^2\psi_\tau\psi_x,\qquad
\partial_\tau e+\partial_x j=0.
$$

Equivalently, e=e_right+e_left and j=e_right-e_left, with
e_right=(R^2/4)g'(x-tau)^2 and e_left=(R^2/4)g'(x+tau)^2.
Each integrates to 2 pi^2 R^2 m^2: exactly half the initial budget.
These are geometric conserved amounts, not calibrated joules.

For R=m=1, the initial total is 39.4784176044:

| Domain | Offset | Spatial Q | Inertial K | Total |
| --- | --- | --- | --- | --- |
| Open | 0 | 39.4784176044 | 0 | 39.4784176044 |
| Open | 1/2 | 9.8696044011 | 29.6088132033 | 39.4784176044 |
| Open | 1 | 19.7392088022 | 19.7392088022 | 39.4784176044 |
| Periodic | 1/2 | 0 | 39.4784176044 | 39.4784176044 |
| Periodic | 2 | 39.4784176044 | 0 | 39.4784176044 |

On the open line, the supports [tau,2+tau] and [-tau,2-tau] separate
at tau=1. Each outgoing packet carries 19.7392088022, with j=+e or j=-e.
Across each end of the original drawing interval [0,2], the outward
flux is pi^2 R^2 m^2 for 0<tau<2 and zero afterward, apart from corner
representatives that do not affect integrals. Thus the explicitly receiving
left and right intervals each acquire

$$
B_{\rm left}=B_{\rm right}=\pi^2R^2m^2\min(\tau,2),\qquad
B_{[0,2]}=4\pi^2R^2m^2-2B_{\rm left}.
$$

At tau=2 the original interval is empty of this geometric budget and each
receiving side holds exactly half. Direct boundary-flux and spatial
integration independently check these identities. These drawing-domain
boundaries are not being identified with the black-hole horizon (0).
Every fixed bounded interval eventually clears. The opposite propagation
directions have angle 180 degrees if read along a straight spatial axis;
the solution alone does not identify that axis or a three-dimensional
detector embedding.

On the periodic coordinate, psi is constant at tau=1/2 but its velocity
is not. Zero spatial variation is not zero total budget. The full initial
state recurs at tau=2: no irreversible drain or child has been produced.
This applies the continuation question concretely: changing the domain
changes where the conserved amount can go, rather than erasing it at a
flat drawing. Neither continuation has been derived to be an event horizon.

Reproduce without files or the conditional parent:

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --transport-only
```

**Next connection:** select the inertial contact rule, initial velocity and
receiving geometry from the offset construction. The present result is an
exact conditional transfer mechanism. Its 1:1 ratio and opposite directions
do not establish electron/positron identification or the 511 keV scale.

## Retained offset-history construction

Douglas's starting proposal is one spatial seed plus an offset operation,
applied repeatedly to generate a circle, a sphere and continued structure.
Here the offset is the primitive update; no independent proper-time
coordinate is used to generate the spatial map. The construction below
realizes part of that proposal with an explicit metric and a proof.

The [calculation](../results/summaries/nsc-atlas-offset-generator-v1.json)
uses one continuous seed coordinate and one irrational offset. A history
rule supplies the two phase values rather than choosing two independent
initial phases. Retaining the earlier draws, choosing a cross-pass distance
law and taking its metric completion are explicit assumptions.

## Preparation and choice

We compared supplying three independent spatial coordinates, claiming a
single smooth instantaneous trajectory fills a volume, and retaining an
offset history whose distance completion has higher dimension. The third
formulates the proposed drawing without making the false instantaneous-rank
claim. The earlier [suspension and golden-map controls](nsc-atlas-golden-vortex.md)
remain valid in their scope; this entry replaces independently supplied
phase coordinates with a single deterministic history.

The most likely mistakes are hiding a second freely supplied phase,
confusing a dense history with its completed space, and using an embedding's
component count as the intrinsic spatial dimension. The formula, completion
proof and induced-metric calculation check those separately. Finite coverage
uses fixed counts 100,200,1000,4000 and at most 120 CPU seconds; no physical
target selects a phase, count or distance rule.

## 1. The same accumulation operation is applied twice

Let alpha be irrational, theta_0=phi_0=0, and set, modulo one,

$$
\theta_{n+1}=\theta_n+\alpha,\qquad
\phi_{n+1}=\phi_n+\theta_n.
$$

Equivalently, Delta theta=alpha and Delta phi=theta. One fixed offset feeds
the first accumulator, whose history feeds the second. Thus

$$
\theta_n=n\alpha\quad\bmod1,\qquad
\phi_n=\frac{n(n-1)}{2}\alpha\quad\bmod1.
$$

The index n counts applications. It is not a separately supplied clock in
seconds. Both phases are determined by that count and alpha; there is no
second independent phase input. The inverse update exists, so the map
alone does not select a thermodynamic arrow or explain why updates occur.

**Density proof.** For any nonzero integer pair (k,l), the corresponding
Fourier character is the polynomial phase

$$
\alpha\left[kn+\frac{l n(n-1)}{2}\right].
$$

If l=0, its average tends to zero by the geometric-sum formula. If l is
nonzero, subtracting the value at n from the value at n+h gives a linear
phase with irrational slope l h alpha for every positive integer h.
Its average again tends to zero. The van der Corput bound then gives an
upper limit of order 1/H for the squared original average for every H;
letting H grow makes it zero. Weyl's criterion proves joint equidistribution
on the two-phase torus. The general polynomial theorem is attributed to
Weyl; [Tao's notes](https://terrytao.wordpress.com/2010/03/28/254b-notes-1-equidistribution-of-polynomial-sequences-in-torii/)
give the character/differencing framework, and
[Yifrach](https://arxiv.org/abs/2201.07138) records the classical statement
before developing a further generalization. No generalization beyond this
quadratic special case is needed here.

A rational offset need not work: alpha=1/4 returns to the initial state
after eight updates and supplies only finitely many directions. A zero
offset leaves the original line. Nonresonant accumulation is an essential
condition, not a word added after the calculation.

## 2. Line, circular boundary, spherical boundary

Put Theta_n=2pi theta_n and Phi_n=2pi phi_n. A first rotation of the retained
seed interval r in [-R,R] gives lines (r cos Theta_n,r sin Theta_n). Their
completion is a disk with a circular boundary.

For the second stage define the direction correlation

$$
C_{nm}=\cos\Theta_n\cos\Theta_m
+\sin\Theta_n\sin\Theta_m\cos(\Phi_n-\Phi_m),
$$

and the distance rule

$$
d^2[(r,n),(s,m)]=r^2+s^2-2rsC_{nm}.
$$

This is a declared metric kernel. Its Gram factorization is

$$
v_n=(\sin\Theta_n\cos\Phi_n,\sin\Theta_n\sin\Phi_n,
\cos\Theta_n),\qquad C_{nm}=v_n\cdot v_m.
$$

The distance is therefore Euclidean distance between r v_n and s v_m;
coincident labels, including every representation of the origin, are
identified. Density of the phase pairs implies density of directions on
the sphere. Every point of the ball is a limit of retained seed points.
Conversely every seed point lies in the ball, which is complete and closed.
**The metric completion is exactly the three-dimensional closed ball.**

This is a conditional constructive theorem: one spatial interval and one
offset history suffice for the stated completion. Finite histories and
the countable union of lines have zero three-dimensional volume before
completion. If motion were restricted to walking along the spokes through
their shared origin, a different distance would result. The distance law,
memory and completion must not be hidden inside the phrase "one plus time."

## 3. The same history supplies the closed parent geometry

For a radial seed r in [0,pi R], choose the positive-curvature distance law

$$
d=R\arccos\left[
\cos(r/R)\cos(s/R)+\sin(r/R)\sin(s/R)C_{nm}\right].
$$

Its completion is the round closed three-dimensional space used by the
parent model. A four-component representation is
(R cos(r/R), R sin(r/R) v_n), with one fixed-norm constraint. The intrinsic
metric, calculated by differentiation, is

$$
d\ell^2=dr^2+R^2\sin^2(r/R)
\left(d\Theta^2+\sin^2\Theta d\Phi^2\right).
$$

There are three intrinsic spatial directions. Four embedding components
do not make four independent spatial dimensions. Likewise a separate
five-coordinate auxiliary fluid/gravity model is not the count of primitive
variables in this generator.

The radius is still a parameter of this geometric construction. The
existing constant-field parent at its retained minimum supplies a compatible
choice through its Hamiltonian constraint. For zero extrinsic curvature,

$$
{}^{(3)}R=\frac{6}{R^2}=2\Lambda_\star,\qquad
R=\sqrt{\frac{3}{\Lambda_\star}}=3.684522144189
$$

in the current model units, with Lambda_star=.220983024213. Momentum
constraints vanish for those constant fields. This joins the generated
spatial geometry to existing admissible parent initial data. It is not a
primordial formation calculation or an SI scale selection.
The earlier action's phase ratio and radius are compatibility inputs;
this construction has not derived that action or its auxiliary
internal-dimensional interpretation from the offset alone.

The scalar Laplacian then has eigenvalues l(l+2)/R^2. The first two positive
values are .220983024213 and .589288064567. Under radius rescaling,
R^2(-Delta) preserves its dimensionless spectrum. This is an explicit
normalized operator inheritance law across separately specified regional
metrics; it does not yet attach those regions physically.

## 4. Finite resolution and measure

For N complete radial lines, let h be the spherical covering radius of
their two endpoint directions. The largest distance from the unit ball
to that line family is exactly sin(h). For N oriented north-to-south
meridians in the closed parent, it is min(h_oriented,pi/2) times R.
Spherical Voronoi vertices and a nearest-neighbor check evaluate these
finite errors; all phases are calculated at 80-digit precision before
converting the direction coordinates.

| Offset | Passes | Ball distance / R | Closed-parent distance / R |
| --- | --- | --- | --- |
| Retained action ratio | 100 | .35147195 | .50145491 |
| Retained action ratio | 1000 | .16691160 | .19978947 |
| Retained action ratio | 4000 | .08015006 | .12788791 |
| Golden-angle fraction | 100 | .39247396 | .57017792 |
| Golden-angle fraction | 1000 | .15218010 | .22441514 |
| Golden-angle fraction | 4000 | .08859494 | .13519849 |

The retained action offset is alpha=.0573586470504441, from the previously
derived ratio of linear mode frequencies. The golden comparison uses
alpha=(3-sqrt(5))/2=.381966011250105. Both are irrational; neither was fitted
to these errors. Their physical identification as the fundamental offset
is still a hypothesis. This construction does not select a cosmic pass
count, the 10^60 hierarchy, or a golden branching average.

Uniform phase passes are also not uniform physical volume. Their direction
second moments tend to (1/4,1/4,1/2). The metric area weight |sin Theta|
instead gives (1/3,1/3,1/3). The closed-parent volume element is
R^2 sin^2(r/R)|sin Theta| dr dTheta dPhi; the ball has r^2 in place of
R^2 sin^2(r/R). This is the metric measure used by the covariant action,
not an identification of update frequency with energy or particle density.

## What this adds to the atlas

The generator is now explicit and does not require an independent second
starting phase. Its completed spatial geometry and a compatible parent
slice are proved under the stated construction rules. Existing spacetime
calculations can be read as tests on a generated geometry, rather than
being asserted as the primitive origin of it.

What remains is to derive the physical selection of this distance kernel,
the offset law, two accumulation levels followed by nesting, and a duration
or irreversible arrow. Further accumulation levels can be constructed, so
exactly three is not selected by irrationality alone. A coherent physical
theory must also supply regional stress, formation and observable predictions.

## A conditional dimension-selection lemma

A more restrictive candidate can select three local directions. Let a
nontrivial real inner-product space V, dim(V)>=2, have a bilinear
antisymmetric bracket satisfying Jacobi and ad-invariance, with one c>0:

$$
\langle[x,u],v\rangle=-\langle u,[x,v]\rangle,
\qquad |[x,y]|^2=c^2(|x|^2|y|^2-\langle x,y\rangle^2).
$$

Then dim(V)=3. Ad-invariance and polarization give
ad_x²=-c²(|x|² Id-x tensor x). For orthonormal x,y define z=[x,y]/c.
The three vectors are orthonormal and satisfy [x,y]=cz,[y,z]=cx,[z,x]=cy.
Their orthogonal complement W is invariant under these adjoint maps.
On W put I=ad_x/c,J=ad_y/c,K=ad_z/c. The squared identity implies
I²=J²=K²=-Id and IJ+JI=0. Jacobi gives [I,J]=K, so K=2IJ and
K²=-4Id, contradicting K²=-Id unless W=0. This proves the claim.
Dimensions0 and1 are excluded because their zero brackets satisfy the
norm equation vacuously.

The quaternionic reference frame realizes this algebra with c=2r.
[Baez, Section4.1](https://arxiv.org/html/math/0105155v4) supplies the
related quaternion/octonion cross-product context; the seven-dimensional
norm-compatible cross product fails Jacobi. The short argument above
does not require a classification theorem.

These are candidate composition/isotropy postulates. The existing single
inheritance identity and offset generator do not derive them, and the
physical sourced metric need not be isotropic. The lemma selects neither
global topology nor radius, internal d=12 or an SI scale. The physical
dimension-selection gate therefore remains open.
