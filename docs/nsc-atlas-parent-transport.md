# Real parent sources: causal access, arrival and received energy

The [regular preparation](nsc-atlas-scalar-preparation.md) determines a
nearly empty homogeneous child. Its tiny Euclidean pole amplitude is
not a transmission coefficient for real energy already in the parent.
This calculation closes a prerequisite to that different preparation:
which parent events can affect the central child, and how radial
freely falling test bodies carry energy between those local frames.

The result uses the solved scalar/gravity cap profiles and charged
membrane, with unchanged couplings. For a source .3 L_ref into the
parent on the time-symmetric slice, the central-arrival threshold is
gamma_s=1.290485734. A supplied gamma_s=2 gives

$$
\frac{c\tau_r}{L_{\rm ref}}=4.285837228,\qquad
\frac{E_r}{mc^2}=1.546392954,\qquad
\frac{E_r}{E_s}=.773196477.
\tag{T1}
$$

This is a geometric response to supplied source data. It selects no
population density, adds no identified dark-matter species, and does
not solve scalar/metric/membrane perturbations. A wave interpretation
requires its own controlled geometrical-optics and interface analysis.
There is no direct nongravitational interaction between the test body
and the membrane in this calculation.

## Scope, reusable mathematics and choice of calculation

The programme needs a substantial evolving component from a physical
parent preparation. Three routes were inspected: continuing to refine
the regular state's small tail; exciting the field through colliding
charged walls; and transporting real, specified parent data. The
first is already decided by the regular branch. The simple additive
collision control below does not supply a homogeneous excitation.
The third requires a causal map before selecting or evolving a packet.
That map is the bounded calculation here; it is not the packet solution.

[Dong and Harlow](https://arxiv.org/html/1109.0011), Section 2.2,
supply the wall/open-chart continuation used in the preparation.
The radial geodesic equations, a de Sitter embedding and the inherited
cap solution supply independent checks. No new field-grid solver,
coupling, observational fit or P001/P002 trajectory refinement is used.
Four supplied source speeds, one threshold root, two existing background
tolerances and independent trajectory integration fit a 120 CPU-second
budget. Evidence is a compact summary; curves are regenerated for the
figure. Failure of the comparison checks stops the calculation.

## The geometry fixes the source domain

In the wall chart, with the global radial coordinate increasing from
the child pole toward the parent,

$$
ds^2=d\xi^2+a(\xi)^2
\left[-d\eta^2+\cosh^2\eta\,d\Omega_2^2\right].
\tag{T2}
$$

The child pole is xi=0, the wall is xi=L_i, and the parent pole is
xi=L_i+L_o. The stored parent profile instead uses its own distance
r from the parent pole: xi=L_i+L_o-r. Thus parent integrals on an
inward ray run from r=L_o-d to r=L_o, then down the child cap.
The regular pole normalizes the conformal coordinate without a free
additive clock constant:

$$
x(\xi)=\ln\xi+
\int_0^\xi\left[\frac{1}{a(s)}-\frac{1}{s}\right]ds,
\qquad e^x\sim\xi.
\tag{T3}
$$

Lengths in the logarithm are expressed in L_ref. Its apparent unit
dependence cancels in the arrival relations. The implementation splits
the integral at the membrane and uses the regular pole expansion.
The scalar layer and its gravitational corrections stay in a(xi).

An inward null ray has constant N=exp(eta+x). Across the child cone,
N matches exp(zeta+chi), where d zeta=d tau/b and zeta-ln(tau) tends
to zero at the regular origin. In the vacuum comparison child,

$$
b(\tau)=\frac{\sinh(h\tau)}{h},\quad
h^2=\frac{\kappa\rho_i}{3},\quad
e^\zeta=\frac{2}{h}\tanh\left(\frac{h\tau}{2}\right).
\tag{T4}
$$

Consequently a source reaches the central worldline at finite time
only if N<2/h. A finite label alone is insufficient: the child has
finite future conformal time. At eta_s=0 this bounds the accessible
parent distance from the wall by

$$
d<3.700717257\,L_{\rm ref}.
\tag{T5}
$$

The parent cap length is 7.070494445 L_ref. The accessible portion is
.775076272 of the proper parent volume on this particular eta=0
slice, obtained from the integral of a^2 dxi. This fraction is not a
cosmic abundance or a slicing-independent fraction of the universe.
The open cosmology beyond the other pole is a different causal region;
it cannot simply be used as an incoming source for this central child.

For d=.3, the source has a_s=2.862664737 and x_s=1.190196801.
A light signal emitted at eta_s=0 arrives at tau_r=3.420193817,
with 1+z=b_r/a_s=1.294206955. At that same source location the
latest central null emission time is eta_s=1.091977295; equality
means asymptotic reception, not a finite observation.

## Massive transport and the wall

For a radial timelike geodesic, eta is cyclic in the reduced radial
metric. The conserved charge per unit rest mass is

$$
K=a^2\frac{d\eta}{ds}=a_s\gamma_s,\qquad
\left(\frac{d\xi}{ds}\right)^2=\frac{K^2}{a^2}-1.
\tag{T6}
$$

This is a radial symmetry statement; partial_eta is not a Killing
field for arbitrary angular motion in T2. Inward passage requires
K greater than the maximum a on the route. Otherwise there is a
turning point or limiting orbit before the child cone.

The metric is continuous at the membrane and has a finite derivative
jump. There is no delta impulse in this minimally coupled test-body
geodesic. Its continuous tangent gives the same wall-frame gamma=K/R_*
on both sides. This is not a claim that the field transmission
coefficient is one: field stresses, wall bending and wave reflection
are absent from the test-body approximation.

The delay relative to a null ray is

$$
D(K)=\int_{\rm source}^{\rm child\ pole}
\frac{|d\xi|}{a}
\left[\frac{1}{\sqrt{1-a^2/K^2}}-1\right],\qquad
L=\exp(\eta_s+x_s+D).
\tag{T7}
$$

The integrand is regular at the child pole. The calculation evaluates
its rationalized form a/[K^2 sqrt(1-a^2/K^2)(1+sqrt(1-a^2/K^2))].
The common regular-cone limit of the two charts is
xi exp(eta)=tau exp(chi). Radial continuation then gives

$$
b^2\left|\frac{d\chi}{ds}\right|=K,\qquad
\frac{d\chi}{d\tau}=-\frac{K}{b\sqrt{b^2+K^2}}.
\tag{T8}
$$

Write chi=ln L-I(tau), with I-ln(tau) tending to zero at the
origin. For T4, exact integration gives

$$
e^{I(\tau)}=
\frac{2/h}{\coth(h\tau)+
\sqrt{\coth^2(h\tau)+(hK)^{-2}-1}}.
\tag{T9}
$$

Since I is increasing and exp I tends to 2K/(1+hK), central reception
is equivalent to

$$
L<\frac{2K}{1+hK}.
\tag{T10}
$$

A trajectory failing T10 can still enter the child. It approaches a
nonzero comoving distance and misses the central observer. This is
distinct from the pre-cone turning condition in T6.

At central reception the local energy follows from the mass shell:

$$
\gamma_r=\sqrt{1+\frac{K^2}{b_r^2}},\qquad
\frac{E_r}{E_s}=\frac{\gamma_r}{\gamma_s}.
\tag{T11}
$$

These are energies in specified local frames, not a globally conserved
energy. For massive bodies an envelope's arrival-time derivative is
not a photon frequency-redshift law and is not used as T11.

| Source gamma | Central arrival c tau/L_ref | E_r/(mc^2) | E_r/E_s |
| --- | --- | --- | --- |
| 1.1 | No finite central arrival | — | — |
| 1.5 | 5.816966363 | 1.160849530 | .773899687 |
| 2 | 4.285837228 | 1.546392954 | .773196477 |
| 5 | 3.523508127 | 3.863701164 | .772740233 |

The primary received speed is .762773146 c. This transport therefore
does not immediately supply cold matter. A distribution's subsequent
redshifting, angular stress and density must be calculated from its
preparation, not inferred from one radial trajectory.

### What this nesting changes against a specified alternative

Remove the membrane and parent vacuum change, replacing the complete
geometry with a single de Sitter space at the same child h. Hold the
source's initial areal radius a_s and its inward gamma_s=2 fixed on
the time-symmetric slice; choose the hemisphere nearest the receiving
pole. This is an explicit comparison with equivalent local source data,
not a second fit. T12 gives

| Quantity | Prepared parent/wall/child | Unwalled de Sitter |
| --- | --- | --- |
| Central threshold gamma_s | 1.290485734 | 1.232271209 |
| Reception c tau/L_ref | 4.285837228 | 4.014164210 |
| E_r/E_s | .773196477 | .811509668 |

The nested reference delays reception by .271673018 L_ref/c and reduces
the received energy by 4.721224 percent relative to this control. These
are consequences of the different intervening geometries with the local
source held fixed. They do not establish nesting observationally: unknown
source radius, velocity or preparation could mimic a changed response,
and this one comparison does not exhaust alternative geometries.

![Accessible parent and central trajectories](figures/nsc-atlas-parent-transport.png)

The first panel shows the causal portion of the parent on eta=0.
The second shows massive trajectories after entering the child. The
gamma_s=1.1 trajectory enters that chart but never reaches its center.

## Independent checks and approximation bounds

For a single de Sitter space with no physical wall, let theta=h xi_s.
Its embedding in five-dimensional Minkowski space gives a geodesic
X(s)=X_s cosh(hs)+U_s sinh(hs)/h. Setting its radial coordinate to
zero yields, independently of T7--T9,

$$
\tanh(hs_r)=\frac{\tan\theta}{\sqrt{\gamma_s^2-1}},\qquad
\sinh(h\tau_r)=\gamma_s\sinh(hs_r),\qquad
\gamma_r=\gamma_s\cos\theta.
\tag{T12}
$$

It requires theta<pi/2 and gamma_s>sec(theta). The implementation
recovers this embedding arrival, energy and threshold across an
arbitrary chart split. Direct integration of T8 on the actual
transmitted labels agrees with T9 inversion within 7e-14 L_ref/c
for the three arriving examples. Changing background tolerance from
1e-7 to 1e-9 changes the primary arrival by 1.06e-11, the energy ratio
by 8.74e-13 and the threshold by 1.35e-12. A conservative numerical
comparison tolerance of 2e-8 is retained; these are numerical checks,
not interval-certified or observational error bars. A common length
rescaling preserves all dimensionless energy ratios and thresholds.

The child contains a residual excitation E0=9.6245e-47. The already
proved comparison bounds are b_hminus<=b_actual<=b_hplus, where
hminus^2=kappa rho_i/3 and hplus^2=hminus^2+kappa E0/3. The T8
integrand decreases with b. Thus, for fixed transmitted K and L,

$$
\tau_{h_-}\le\tau_{\rm actual}\le\tau_{h_+}
\tag{T13}
$$

whenever both comparison trajectories reach the center. T9 is used
at hminus numerically. These analytic brackets isolate the original
regular-state correction without resolving another tiny tail. They
do not bound backreaction from a new finite-density population, and
the threshold's asymptotic limit is not a finite-time error estimate.

## A collision control: what free passage does not supply

[Giblin, Hui, Lim and Yang](https://arxiv.org/abs/1005.3493) develop
the high-speed free-passage construction for scalar bubble collisions.
For the current quadratic potential, the simpler fixed-metric scalar
equation is linear:

$$
(\Box_g-\omega^2)\phi=\mu q,\qquad
\phi_*(q)=-\frac{\mu q}{\omega^2}.
$$

On one prescribed metric, continuing prescribed membranes with
additive swept flux permits superposition. If the two interiors
reach their minima,

$$
\phi_*(q_p-e_L)+\phi_*(q_p-e_R)-\phi_*(q_p)
=\phi_*(q_p-e_L-e_R).
\tag{T14}
$$

For q_p=1.5 and equal charges .5, the overlap has q=.5, vacuum
density .03125 and zero homogeneous displacement energy at the
minimum. This identity cannot provide the earlier quench energy
.09375. It does not prove that an interacting gravitational collision
radiates nothing: dynamical membrane deflection, scalar tails and
metric perturbations were not included. Terminating the membranes
would additionally require an annihilation and energy-transfer law.

## Next physical connection

Prepare a finite parent scalar packet inside the accessible region,
state its supplied energy and solve its gravitational constraints.
Derive the scalar, metric and membrane-displacement response together;
then compute the transmitted stress and its angular/time dependence.
The source must be fixed before comparing its output density with data.
The smooth-wall perturbation framework of
[Garriga, Montes, Sasaki and Tanaka](https://arxiv.org/html/astro-ph/9811257)
is relevant bulk mathematics, but its interface conditions cannot be
imported unchanged to an independently charged thin membrane.

The most likely obstructions are constraint-incompatible source data,
extra interface dynamics omitted by a scalar probe, and a transmitted
stress that stays anisotropic or too fast to cluster as required.
Each has a definite next equation to test. None is resolved by raising
the precision of the regular pole amplitude or choosing a child density.
The vacuum value Lambda L_ref^2=.125 remains the supplied branch's
value; transport has not newly selected it.

Evidence: [compact record](../results/summaries/nsc-atlas-parent-transport-v1.json).

```
python scripts/run_python.py scripts/check_nsc_atlas_parent_transport.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_parent_transport.py -q
python scripts/run_python.py scripts/plot_nsc_atlas_parent_transport.py
```
