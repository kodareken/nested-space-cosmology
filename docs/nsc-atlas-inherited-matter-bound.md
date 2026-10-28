# One inherited scalar supplies a controlled vacuum-and-matter background

This continues the [shared scalar parent](nsc-atlas-shared-scalar-slice.md)
and the [high-frequency hierarchy](nsc-atlas-hierarchy-evolution.md).
The inherited potential already has a vacuum floor and an excited real
field. The missing connection is a controlled pressureless-background
approximation with its normalization fixed by that same state.

## Preparation and decision

Three possible continuations were considered: solve the full interacting
charged-wall future, infer the scalar parameters from geometric time jets,
or first bound the background approximation already used by the structure
calculation. The third directly joins the existing source, expanding parent
and observable-distance calculations without a new state or dust component.
The protected parent is an actual smooth domain of the shared initial atlas;
its calculation does not require prescribing the distant walls' future.

The classical oscillating-scalar matter mechanism is borrowed, not a claim
of priority: [Turner (1983)](https://doi.org/10.1103/PhysRevD.28.1243) gives
the rapid-oscillation energy law; [Hwang and Noh, equations (2)--(6)](https://arxiv.org/html/0902.4738)
state the background averaging and its H/m expansion. Here the target is an
explicit, uniform bound in the actual closed, self-gravitating parent,
followed by its light-cone consequence. Clustering perturbations and the
full halo approximation retain their separate equations and error questions.

Keep kappa=1, mixing=1/4, fluxes (1/2,1,3/2), charge=1/2, tension=1/5,
the uniform grandchild-minimum initial field and zero scalar momentum.
The same initial geometry is independent of frequency. Evaluate only the
already used frequency members 20, 2000 and 10^9. Their parent vacuum is
rho_v=9/32 and excitation E0=3/8. No SI ruler or observation chooses them.

The finite calculation will: derive the corrected-energy identity;
prove exact rational comparison factors from independently checked input
bounds; evaluate density/equality and luminosity-distance predictions;
then compare those enclosures with the original unaveraged Einstein-scalar
equations in proper and expansion clocks. The light-cone protocol uses
A_observer=2 and the primary redshift z=1/4: the
emitter at z=1/4 then has A_emitter=8/5, so the comparison integral admits
an exact rational polynomial/tail enclosure. This choice is fixed before
evaluating any target and stays within the previously solved one-e-fold
domain. It changes no source coupling or state.

Two frequency-20 evolutions, one short frequency-2000 control, quadrature
and symbolic checks fit 180 CPU seconds and 180000 total RHS calls,
with up to four useful workers. Store summaries and hashes, not histories.

Failure conditions are an incorrect pressure/scale factor, loss of positive
comparison energy, or an unaveraged source/observable escaping its stated
enclosure. A failure is diagnosed before another calculation is designed.
Success closes this background approximation and its specified light-cone
map; it does not select the ancestral flux, produce primordial seeds, or
establish perturbation/halo averaging errors.

## Exact connection to calculate

In the shared parent, A=a/a0 and

$$
E=\frac{\Pi^2+\omega^2\chi^2}{2},\qquad
p=\frac{\Pi^2-\omega^2\chi^2}{2},\qquad
\dot E=-3H\Pi^2.
$$

Define D=A^3 E, S=A^3 chi Pi and I=D+(3/2)H S. The candidate identity is

$$
\dot S=2A^3p,\qquad
\dot I=\frac{3}{2}\dot H S,\qquad |S|\leq D/\omega.
$$

The comparison uses the actual turning-point initial state S(0)=0.
The initial inherited field is at the grandchild minimum, so

$$
E_0=\frac{1-\epsilon}{2}(q_p-q_g)^2=\frac38,\qquad
\rho_v=\frac{\epsilon q_p^2}{2}=\frac9{32},\qquad
\frac{E_0}{\rho_v}=\frac43.
$$

This is the normalization of the already specified state. No independent
dust density has been added. It does not select that ancestral preparation
from other possible states.

## A uniform comparison on the actual parent

The initial wall gives H0 and curvature C=1/a0^2. The existing positive
vacuum bound and decreasing E imply H_minus<=H<=H_plus for all A>=1.
Independent exact rational evaluation of the unsquared initial junction
verifies the conservative choices

$$
H_-=\frac{133}{500},\quad H_+=\frac{117}{250},\quad
H_0<\frac{443}{1000},\quad C<\frac{57}{2500}.
$$

Since dot H=C/A^2-kappa*Pi^2/2 and A>=exp(H_minus*t),

$$
\int_0^t|\dot H|\,dt
\leq H_0-H(t)+2\int_0^t\frac{C}{A^2}\,dt
\leq V_H\equiv\frac{443}{1000}-\frac{133}{500}
+\frac{57/2500}{133/500}.
$$

For delta=3H_plus/(2omega)<1, the corrected energy obeys
(1-delta)D<=I<=(1+delta)D. Its initial value is I0=E0 because
the inherited momentum is zero. Thus

$$
|\log(I/I_0)|\leq\Gamma,
\qquad\Gamma=\frac{3V_H}{2\omega(1-\delta)}.
$$

**Uniform source theorem.** If Gamma<1, the exact rational factors

$$
L=\frac{1-\Gamma}{1+\delta},\qquad
U=\frac{1}{(1-\Gamma)(1-\delta)}
$$

satisfy, on every future time in the protected parent,

$$
L\leq\frac{A^3E}{E_0}\leq U,\qquad
\frac43LA^{-3}\leq\frac{E}{\rho_v}\leq\frac43UA^{-3}.
$$

**Proof.** The scalar equation gives dot S=2A^3 p. Energy conservation
then cancels the pressure term in dot I, giving the stated identity.
The inequality omega*|chi*Pi|<=E bounds S. Positivity of I bounds its
logarithmic derivative by 3|dot H|/[2omega(1-delta)]. Integrate and use
exp(-Gamma)>=1-Gamma and exp(Gamma)<=1/(1-Gamma). Every step is uniform
in the future interval. Zero excitation remains identically zero.
For nonzero initial momentum the additional factor I0/E0 must be retained;
it cannot silently be replaced by one.

The excitation decreases strictly across any nonzero time interval of this
nontrivial solution and tends to zero. Its unique energy-equality event
E=rho_v therefore has

$$
\left(\frac43L\right)^{1/3}\leq A_{\rm eq}
\leq\left(\frac43U\right)^{1/3}.
$$

This equality concerns densities, not the phase-dependent instantaneous
acceleration. Instantaneous scalar pressure can have amplitude comparable
to E. Neither pressureless perturbations nor a small error in dot H is
asserted.

## The same bound reaches an interior light cone

Keep the exact parent curvature and define the comparison

$$
H_d^2=\frac{\kappa\rho_v}{3}
+\frac{\kappa E_0}{3A^3}-\frac{C}{A^2}.
$$

Let epsilon_E=max(1-L,U-1), d=kappa*E0/3, and
r=(d/H0^2)*epsilon_E. For A>=1, the denominator of
d/(A^3 H_d^2) is increasing because 3*(kappa*rho_v/3)*A^2-C>0.
Consequently

$$
1-r\leq\frac{H^2}{H_d^2}\leq1+r.
$$

For comoving endpoints, z=A_o/A_e-1 and

$$
\eta=\frac1{a_0}\int_{A_e}^{A_o}\frac{dA}{A^2H},\qquad
D_L=a_0A_o(1+z)\sin\eta.
$$

Integrating the H bounds, then using the increasing, concave sine on the
declared angular domain, proves

$$
\frac1{\sqrt{1+r}}\leq\frac{D_L}{D_{L,d}}
\leq\frac1{\sqrt{1-r}},\qquad
\frac1{1+r}\leq\frac{\Lambda_vD_L^2}{\Lambda_vD_{L,d}^2}
\leq\frac1{1-r}.
$$

The protocol A_o=2, A_e=8/5 has eta_d/sqrt(1-r)<1<pi/2;
the original causal parent margin contains the complete ray. For u=1/A,
H_d^2 is a positive cubic polynomial P(u). Expand P^(-1/2) about the
midpoint of its endpoint range. Twelve binomial terms integrate exactly
as rational polynomials, with the entire tail bounded by
Delta_u*q^13/(1-q). An alternating sinc series then bounds the closed-space
distance. Thus the saved decimal interval is rounded outward from exact
rationals, rather than estimated from floating quadrature agreement.

The [frozen analytic calculation](../results/summaries/nsc-atlas-inherited-matter-forecast-v1.json)
was constructed without an unaveraged target solve:

| Frequency | Uniform relative density ceiling | Relative ceiling for Lambda_v D_L^2 | Enclosed Lambda_v D_L^2 at A_o=2, z=1/4 |
| --- | --- | --- | --- |
| 20 | .05798115 | .03840303 | [.242823328986305, .261473628319700] |
| 2000 | .00054834 | .00034988 | [.251715570867755, .251891708276777] |
| 10^9 | 1.09608e-9 | 6.99120e-10 | [.251803608594062, .251803608946144] |

For the high-frequency member, the approximate printed energy-equality
range is A_eq in [1.10064241589608,1.10064241670034]; the cube-root formula
above is the bound itself. The vacuum length product remains the inherited
Lambda_v*a0^2=12.3482523989285. This work controls the evolving component
and the observer reduction without replacing that supplied vacuum source.

The unaveraged proper-time, independent expansion-clock and short
high-frequency checks are the next execution step. A successful comparison
will be recorded separately, leaving the analytic forecast unchanged.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_inherited_matter.py
python scripts/run_python.py scripts/check_nsc_atlas_inherited_matter.py --validate
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_inherited_matter.py
```
