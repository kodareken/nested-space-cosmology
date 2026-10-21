# A shared scalar atlas with an evolving protected parent

The [geometric observer](nsc-atlas-geometric-observer.md) established a
two-sector ambiguity: different inherited laws can reproduce one interface.
This entry supplies the missing third region on a compact initial slice,
rules out the alternative law on that geometry for every canonical scalar
state, and evolves a causally protected part of the actual parent.

The whole initial atlas satisfies the bulk and material constraints and
the instantaneous temporal Israel condition. A smooth parent domain evolves
for all future proper time. The coupled future of both walls and their weak
scalar corners remains a separate existence and transmission problem.

## Rethinking, inputs and calculation boundary

Requiring a second time-symmetric O(4) preparation or a completed wall PDE
first would postpone both composition and cosmology. The alternatives were
a double preparation saddle, full one-wall evolution, and a moving common
Cauchy slice. The third route permits the same scalar state to occupy three
charge sectors, with one shared middle mass. The parent domain of dependence
then supplies genuine evolution without prescribing the walls' future.

The existing [flux-transition equations V5--V9](nsc-atlas-flux-transition.md)
supply the signed thin-wall identities. The smooth bulk constraint equations
are standard Einstein-scalar mathematics; see
[Choquet-Bruhat, Isenberg and Pollack, Section 2](https://arxiv.org/html/gr-qc/0610045).
Their smooth constraint results are not invoked as a thin-shell evolution
theorem. The spherical constraint calculation, wall joins, sector obstruction
and causal bound below are evaluated explicitly for this action.

The new data are two radii, one middle mass and a specified uniform scalar
preparation. They are initial conditions, not new couplings. A quintic blend
owns the middle slicing choice, and a positive lapse owns the initial gauge
extension. Neither is a hidden force or fitted waveform. One finite connector
quadrature, an independent proper-distance ODE, symbolic identities and two
parent clocks fit a 120 CPU-second budget. No spacetime grid scan is needed.

## One action and three regions

Keep the scalar/four-form action and constant charged tension fixed:

\[
U(\phi,q)=\rho_b+\frac{m^2\phi^2}{2}
              +\frac{(q+\mu\phi)^2}{2},\qquad
\omega^2=m^2+\mu^2,\qquad \epsilon=\frac{m^2}{\omega^2}.
\]

Use kappa=1, epsilon=1/4, omega=20, mu=10 sqrt(3), charge e=1/2,
tension sigma=1/5 and physical bare density zero. Inside outward, label the
grandchild, middle child and parent by g,m,p, with fluxes (1/2,1,3/2).
On the initial slice set the field everywhere to the grandchild minimum,
phi0=-sqrt(3)/80, with zero momentum and spatial gradient. Then

\[
(U_g,U_m,U_p)=\frac{(1,7,21)}{32},\qquad
(\rho_{v,g},\rho_{v,m},\rho_{v,p})=
\left(\frac1{32},\frac18,\frac9{32}\right).
\]

The outer regions are excited: the initial field forces are (0,S,2S),
where S=mu e=5 sqrt(3). They cannot be evolved as vacuum Kottler regions.
At this single instant, however, their stress is isotropic with p=-U,
so the constraint representation is

\[
f_j(R)=1-\frac{2M_j}{R}-\frac{\kappa U_jR^2}{3},\qquad
(M_g,M_m,M_p)=(0,3/25,0).
\]

The algebraic reference flux is q+mu phi0; the reference constant term
m^2 phi0^2/2=3/128 is part of this field's existing potential, not a change
to the physical bare density. Set wall radii R1=1 and R2=19/10.

## Two moving walls and one middle geometry

The common normal points inside to outside. Write k=kappa sigma/2=1/10.
The ordinary branch of the existing Israel equation gives

\[
\beta_i=\frac{M_o-M_i}{kR^2}
 +\frac{\kappa(U_o-U_i)+3k^2}{6k}R,\qquad
\beta_o=\beta_i-kR,\qquad v^2=\beta_i^2-f_i.
\]

Both walls use the positive square root v. The two middle entries are
literally the same region object and the same mass M_m.

| Initial quantity | Inner wall | Outer wall |
| --- | ---: | ---: |
| Radius R | 1 | 1.9 |
| Areal speed v | 1.204916145077 | .841108986073 |
| Inside beta | 1.5625 | 1.148006694367 |
| Outside beta | 1.4625 | .958006694367 |
| Pressure difference P | .1875 | .4375 |
| Areal acceleration | -3.173177083333 | 1.401482769823 |

The first areal speed exceeds one while the wall remains timelike; it is
not a local Lorentz speed . Both walls initially expand, and the inner one
initially decelerates. Persistence and stability do not follow.

For wall extrinsic eigenvalues (A,B,B), B=beta/R, the pressure condition is
sigma times mean K=P. It determines

\[
\overline{A}=P/\sigma-(B_i+B_o),\qquad
A_i=\overline{A}+k/2,\quad A_o=\overline{A}-k/2.
\]

The inner A pair is (-2.0375,-2.1375); the outer pair is
(1.129071900666,1.029071900666). Both temporal jumps are -k. An independent
acceleration calculation uses the full moving Misner--Sharp mass:

\[
m_j=\frac R2(1+v^2-\beta_j^2)
 =M_j+\frac{\kappa U_jR^3}{6},\qquad
\ddot R=-\frac{m_j}{R^2}-\frac{\kappa R p_{n,j}}2+A_j\beta_j.
\]

Here p_n=-U initially. Omitting v^2 in the mass produces an incorrect
acceleration. The tests explicitly reject that resting-wall substitution.

## Explicit compact spatial geometry

Write the spatial metric as d ell^2+R(ell)^2 d Omega^2, and use L for its
extrinsic curvature, distinct from the wall's K. The walls initially comove
with the slice normal. Put H1=v1/R1 and H2=v2/R2.

The grandchild is a regular hyperbolic cap with sectional curvature
c_g=U_g/3-H1^2=-369/256 and L=-H1 times the identity. Its proper radius is
.846543496782. The parent is a large spherical cap with
c_p=U_p/3-H2^2=.022776502368 and L=-H2 times the identity. Measured from its
regular pole, its length is 18.889384373634 and

\[
a_0=c_p^{-1/2}=6.626077076100,\qquad
R(r)=a_0\sin(r/a_0),\qquad
\theta_w=\pi-\arcsin(R_2/a_0)=2.850764359770.
\]

The global outward normal is opposite to increasing parent-pole r.
For the middle annulus R1<=R<=R2, set z=(R-R1)/(R2-R1) and

\[
H(R)=H_1+(H_2-H_1)(6z^5-15z^4+10z^3),\quad v_s=RH(R),
\]

\[
\frac{dR}{d\ell}=\beta_s=\sqrt{f_m+v_s^2},\qquad
L_r=-\frac{dv_s}{dR},\quad L_\theta=-\frac{v_s}{R}.
\]

The endpoint first derivatives of H vanish, so both extrinsic eigenvalues
match the adjacent caps. The endpoint minimum of f_m is a global positive
lower bound .610455043860 because f_m is concave. The connector is therefore
regular without a sampled positivity assumption. Its length is
.657495129391, and total pole-to-pole length is 20.393422999807. R is positive
between regular poles; the compact topology is S3, with the intended two
material walls and no extra junction or arbitrary outer cutoff.

Direct substitution proves the bulk constraints:

\[
{}^{(3)}R=-4R_{\ell\ell}/R+2(1-\beta_s^2)/R^2,
\quad R_{\ell\ell}=f_m'/2+v_sv_s',
\]

\[
{}^{(3)}R+(\operatorname{tr}L)^2-L_{ij}L^{ij}=2\kappa U_m,
\qquad -2\partial_\ell L_\theta
 +2\beta_s(L_r-L_\theta)/R=0.
\]

At each material wall, [beta]=-kR gives the scalar-curvature delta
coefficient -4[beta]/R=4k=2 kappa sigma. The continuous L_theta gives no
surface momentum, as required in the comoving slice. The scalar value and
first normal derivative are continuous, and the outward flux jump is e.

A globally unit lapse would not implement these temporal junctions at
fixed wall labels. Use disjoint proper collars of width delta=.1 and put

\[
\log N=A_{\rm side}\delta z(1-z^2)^4,
\quad z=(\ell-\ell_w)/\delta,\quad |z|<1;
\qquad N=1\ \text{outside the collars}.
\]

N=1 on each wall, and its one-sided logarithmic derivatives are exactly
A_i,A_o. Its global bounds are .956493750123<=N<=1.045485137642. It is
smooth through each collar edge and one near both poles. Together with
metric time derivative -2NL and zero shift, it realizes the continuous
initial metric and the required wall first jets. This fixes an initial
gauge, not a solved future lapse or a new physical source.

## Local observations now reach an actual third region

The existing observer inverse at both walls recovers the same phi0 and
zero momentum; squared-gradient residuals are at roundoff. Subtracting
the inferred excitation gives ground densities 1/32 and1/8. At the parent
side, subtracting excitation 3/8 gives 9/32. Initial Ricci density is not
automatically the ground vacuum. In particular, the raw density second
difference is 1/4, while the ground-vacuum second difference is 1/16.

At the common initial field, P_outer-P_inner=e^2=1/4. At later differing
wall fields the exact relation is instead e^2+S(phi_outer-phi_inner).

Candidate B from the preceding observer entry has epsilon=25/28,
e=sqrt(7)/2 and parent flux 3sqrt(7)/10. Every integer sector satisfies

\[
q_n=e(3/5-n),\qquad
U_{\min,B}(n)=\frac{25}{32}(3/5-n)^2\ge\frac18.
\]

The minimum is at n=1, including all positive and negative integers.
Canonical kinetic and gradient energy are nonnegative. But the new
grandchild's geometry independently requires

\[
\rho_{\rm geometry}
=\frac{{}^{(3)}R+(\operatorname{tr}L)^2-L_{ij}L^{ij}}{2\kappa}
=\frac1{32}.
\]

Thus every B state misses this geometry by density at least 3/32, or
Hamilton residual at most -3/16. Its prescribed next sector has the larger
ground value 49/32. No scalar preparation or choice of integer sector can
remove the obstruction under B's fixed action. This rejects B for the
constructed geometry; no astronomical observation has been supplied.

Conditional on a bounded-geometry weak development, the two corners have
inward flux 18.75 tau1^2 and 56.25 tau2^2 to leading order. These use separate
wall proper clocks. Mass change additionally includes pressure work, and
the middle's initial vacuum-subtracted integration constant is .12, not
zero. Subtracting these local flux coefficients is not a middle energy
budget. The record preserves their distinct pressure and mass coefficients.

## A parent cosmology inside the actual domain of dependence

The parent cap data are homogeneous closed-FLRW data for the same scalar.
Use its local proper time, chi=phi-phi_min(q_p), and
rho_v=9/32. The initial excitation is E0=3/8. The smooth equations are

\[
\dot\chi=\Pi,\quad \dot\Pi=-3H\Pi-\omega^2\chi,
\quad \dot a=aH,\quad \dot H=a^{-2}-\kappa\Pi^2/2,
\]

\[
E=(\Pi^2+\omega^2\chi^2)/2,\quad
3(H^2+a^{-2})=\kappa(\rho_v+E),\quad \dot E=-3H\Pi^2.
\]

Construct the homogeneous smooth solution and restrict it to the domain
of dependence of the actual parent cap. Smooth causal uniqueness fixes
this region independently of a material-wall continuation. The nonconstant
initial lapse collar is a coordinate choice and does not change that
development. This argument neither assumes nor proves a full weak
development of the two shells.

Let h_min^2=kappa rho_v/3. Since a0 h_min>1, the initially positive H cannot
reach zero: while expanding, E>=0 and a>=a0 imply
H^2>=h_min^2-a0^-2>0. Also E<=E0, so H is bounded above. The scalar and
expansion equations remain regular at every finite future time. Moreover

\[
\int_0^\infty\frac{dt}{a(t)}
\le\int_{a_0}^\infty\frac{da}{a\sqrt{h_{\min}^2a^2-1}}
=\arcsin\frac1{a_0h_{\min}}=.515418396612.
\]

Every comoving worldline with parent-pole angle theta<2.335345963158 stays
in this domain for all future proper time. Its initial volume fraction is
at least .907035570245 of the parent cap. This is a sufficient lower bound,
not the exact event horizon. It protects against all causal wall signals;
it does not supply the interacting exterior of this domain.

![Compact shared scalar geometry and protected parent](figures/nsc-atlas-shared-scalar-slice.png)

The chosen reference length is the initial parent spatial curvature radius
L_ref=a0, determined by the declared initial wall geometry independently
of any target observation. The inherited ground contribution is

\[
\Lambda_{\rm parent}L_{\rm ref}^2
=\kappa\rho_v a_0^2=12.348252398928.
\]

This is a nonzero conditional dimensionless value. It depends on the action
and initial data; it is not a state-selection theorem, an unexplained extra
vacuum generated by averaging, or the observed cosmic Lambda in SI units.

Two independent integrations, in proper time and in N=log(a/a0), reach
one e-fold with:

| Quantity | Parent result |
| --- | ---: |
| Proper time | 2.888845513139 |
| Scale factor | 18.011544909930 |
| H | .311446897816 |
| Scalar excitation E | .018994903671 |
| Integrated scalar expansion work | .356005096325 |
| E a^3 / (E0 a0^3) | 1.017394237405 |
| Vacuum density fraction | .966503113854 |
| Scalar excitation fraction | .065275141493 |
| Spatial-curvature fraction | -.031778255360 |

The scalar is within 1.74% of pressureless dilution over this interval.
That comparison does not establish clustering or a dark-matter explanation.
The fractions use the actual local H and sum to one within numerical error;
they are outputs of this benchmark, not fits to measured cosmic fractions.
The two clocks agree to below 1.5e-11 on the compared endpoints. Maximum
Friedmann and work residuals are below 3.3e-11; the declared numerical
comparison tolerance is 1e-9 in these normalized variables. This is a
numerical consistency tolerance, not an interval proof or model uncertainty.

## Reproduction and next connection

Run the existing isolated runtime:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_shared_scalar_slice.py
python scripts/run_python.py scripts/plot_nsc_atlas_shared_scalar_slice.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_shared_scalar_slice.py
```

The [compact record](../results/summaries/nsc-atlas-shared-scalar-slice-v1.json)
contains parameters, wall jets, bounds, endpoint comparisons and producer
hashes; it stores no trajectory arrays. Eleven symbolic identities and
eleven targeted tests include an independent derivative check, proper-radius
integration, exact lattice floor, lapse derivatives and both evolution clocks.
Grok independently checked the constraints, lapse and energy-floor argument;
the parent causal/evolution review is recorded when complete.

Next: evolve scalar and metric perturbations on this protected parent and
calculate whether the same inherited excitation can cluster, at which
scales, and with what observer transfer. That advances the cosmology and
dark-matter questions using the constructed source. The wall transmission,
formation and full atlas future remain important distinct tasks; none is
silently promoted from these initial constraints or the protected domain.
