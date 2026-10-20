# The shared parent on a transported clock

This connects the [shared two-wall history](nsc-atlas-connection-slice.md#a-shared-two-wall-history-and-its-characteristic-origins)
to the [geometric average and its conservation law](schur-memory.md#from-nested-geometry-to-a-conserved-mean-response).
The full retained radial parent now evolves on a single specified clock,
including both material interfaces and both bridge minima. The result is a
finite evolving effective scalar source, with an evaluated trace coefficient
and energy exchange. It also decides the large-domain limit of this route.
It uses the same quiet piecewise Einstein geometry; the regional vacuum
values are still inputs.

## Decision, domain and calculation

The previous obstruction was the absence of a transported averaging leaf.
Reassigning the old tapers at new wall radii would hide that choice.
Three routes were considered: a congruence of initially normal geodesics,
radar simultaneity, and a four-volume curvature/action projection. The first
continues the existing normal-frame average without a new material law.
Radar coordinates do not cover the full bridge-connected domain from the
cap's single observer. A four-volume projection would answer a different
question from the already defined spatial effective density and pressure.

The new clock is the proper time of freely falling test observers initially
normal to the original slice. Their initial proper-distance label is x.
This is a declared observation/averaging prescription, not a selected cosmic
rest frame. No field state, coupling, wall tension or target Lambda is fitted.
The outer boundary follows the original R=4 ell normal ray; its subsequent
areal radius is not artificially fixed.

The finite domain is 0<=t/ell<=.01, across the entire centre-to-cut slice.
The wall events stay inside their previously established proper-time
intervals [0,.04 ell]. Equal normal times are now simultaneous by this
specified foliation. Their two wall proper times are different.
Take L_ref=2/k=ell, fixed by the existing membrane scale before comparison
with any cosmological data.

The calculator permits 3,000,000 charged scalar quadrature/RHS evaluations
and 120 process-CPU seconds. It uses two fixed spatial quadrature rules,
32 and 64, at three normal times. A two-step time difference checks the
independent evolution identities. It is not a P001 field-grid campaign.
Only scalar summaries and source hashes are saved. The first exploratory
calculation preceded this completed record; these are synthetic calculated
consequences, not preregistered observations.

The main failure checks are: a lost signed junction or crossing branch,
a focal/singular normal map, and a finite trace being misidentified as
vacuum. They are exposed respectively by independent neighbouring-ray
derivatives and wall matching, positive radius/Jacobi bounds, and the
separate pressure and drift calculations below.

## A horizon-regular normal evolution

In each smooth region write

$$
ds^2=-dt^2+A(t,x)^2dx^2+R(t,x)^2d\Omega^2,\quad
V=R_t,\quad B=A_t,\quad E=R_x/A.
$$

The signed E is the normal geodesic's Killing energy on a static overlap.
It is the original Rprime, not the original V. Its sign distinguishes the
two sheets. For regional mass M and supplied curvature H_i,

$$
R_t=V,\qquad V_t=-\frac{M}{R^2}+H_i^2R,\qquad
A_t=B,\qquad B_t=\left(\frac{2M}{R^3}+H_i^2\right)A.
\tag{N1}
$$

Initially A=1, B=V_R on each smooth original segment. The first integral
is V^2+f=E^2, f=1-2M/R-H_i^2R^2. Differentiating the radial equation with
respect to x proves R_x=EA throughout each such congruence: both R_x and
EA solve the same linear Jacobi equation with the same initial data.
Consequently E_t=0 and

$$
{}^{(3)}R=-\frac{4E_x}{AR}+\frac{2(1-E^2)}{R^2},\quad
\theta=\frac BA+\frac{2V}{R},\quad
\sigma^2=\frac13\left(\frac BA-\frac VR\right)^2.
\tag{N2}
$$

This computes intrinsic curvature independently of the Hamiltonian
constraint. The latter, R3+4(B/A)(V/R)+2(V/R)^2=6H_i^2,
is then checked. Differentiating that first integral in x supplies the
same constraint; the spherical momentum constraint follows from R_x=EA.
Thus (N1) transports the actual Kottler vacuum, rather than fitting a new
stress to a prescribed radial profile.

No equation divides by f or E. At a bridge minimum E=0, while A remains
positive: a spatial minimum is not a caustic. At a cosmological horizon
f=0 the same equations stay regular. Static Killing clocks are used only
on the exterior wall's static collar, where f is bounded away from zero.
The signed mass/embedding conventions are inherited from the
[original spherical construction](nsc-atlas-connection-slice.md) and
[Kuchar's canonical treatment](https://arxiv.org/html/gr-qc/9403003v1);
(N1)--(N2) are derived directly here for the retained regional metrics.

## Crossing the two moving walls

Let U,n be the wall's unit future tangent and outward unit normal.
With N,s the normal-clock frame, set
U=gamma N+b s, n=b N+gamma s, b=gamma v.
The signed Israel law is [beta]=-kR. It implies, at the same wall event,

$$
[V]=kRb,\qquad [E]=-kR\gamma,\qquad
\frac{dt}{d\tau}=\gamma,\qquad \frac{dx_W}{dt}=\frac vA.
\tag{N3}
$$

Brackets in (N3) are outward spatial jumps. A normal geodesic crosses
in the direction opposite the sign of b. Along its future trajectory,

$$
V_{\rm after}=V_{\rm before}-kR|b|,\qquad
A_{\rm after}=A_{\rm before},\qquad
B_{\rm after}=B_{\rm before}+\frac{kA}{|b|}.
\tag{N4}
$$

The last term is necessary. To derive it, use wall Gaussian distance eta.
The radial curvature impulse of the pure-tension Israel jump is
R_{N s N s}^{sing}=-k delta(eta). This two-dimensional sectional
curvature is unchanged by the radial boost. The Jacobi equation gives
A_tt^{sing}=k A delta(eta); N(eta)=-b converts the impulse to kA/|b|.
The angular Jacobi impulse similarly gives the radial-velocity jump in
(N4). A coordinate or sampling prescription cannot omit these impulses.
The existing [surface-source conventions](nsc-atlas-connection-slice.md#distributional-shell-constraints)
and [shell-frame reference](https://arxiv.org/html/gr-qc/9708008v2)
fix the signs; no dust-shell dynamics is substituted for membrane tension.

For this state the inner wall has v<0: incoming normals originate in the
cap and enter the host. The outer wall has v>0: incoming normals originate
in the exterior and enter the host. Neither statement supplies a causal
signal between the two walls. They inherit one initial slice and its
already matched host chart origins.

For inner crossings the incoming solution is exact:
R(t,x)=x exp(t/ell), A=exp(t/ell). Integrating dt/dtau=gamma on the wall
locates its original label. For outer crossings, integrate the geodesic
from its original exterior-taper radius and match both R and the signed
static Killing time to the wall. The initial time offset is fixed by
dT/dR=V/(fE), and the geodesic has T_t=E/f. This gives its parent time,
initial label and wall proper time without resetting a host clock.

The outgoing signed energy is constant after crossing. Its original-label
derivative is also needed in (N2), rather than being chosen to satisfy the
constraint. If s denotes crossing time and d=sign(b), then

$$
E_{\rm after}=E_{\rm before}+d\,kR\gamma,\quad
\gamma_s=b\left(\frac{K^\tau{}_{\tau,\rm in}}{\gamma}
                   -v\frac BA\right),
$$
$$
(E_x)_{\rm after}=(E_x)_{\rm before}
 +d\,\frac{kA}{v}\left(\dot R_W+R\gamma_s\right).
\tag{N5}
$$

Here Rdot_W is the derivative with wall proper time.
Equation (N5) follows by differentiating E_after along the wall and
dividing by dx_W/ds=v/A. Independent differences of neighbouring
transmitted rays verify both R_x=EA and this E_x.

The original label of an initial wall also separates two smooth congruence
pieces after t=0. Their one-sided label Jacobians can differ: the original
normal has a spatial derivative kink. The same geometric tangent and
signed E agree there, so no additional material atom is inserted.
Piecewise proper charts and their proper-volume Jacobians are used.

## The finite normal domain

All original non-central pieces have B(0)>=0. With M>=0, H_i^2>=0,
(N1) and the positive impulses (N4) preserve A>0 and B>=0.
The cap centre is handled by its exact regular expanding solution.
The smallest non-central initial radius is the host throat
rb/ell=.020008009615388216.

As long as R>=rb/2, R_tt>=-4M/rb^2 for an original host ray,
whose initial V is nonnegative. Thus

$$
R(t)\geq r_b-\frac{2Mt^2}{r_b^2}
 \geq .01501201201859\,\ell>r_b/2
\quad (t\leq.01\ell).
\tag{N6}
$$

A first-exit argument closes that bound. The outer throat has the much
larger bound .799875 ell. Transmitted rays enter at R>=.8 ell with
positive V on this interval and remain well outside these minima.
Ordered initial labels and positive one-sided Jacobi scales preserve the
radial normal map on the retained sheets.

The crossing signs are protected by the wall rapidity equation in (N5).
For the outer incoming collar the original label radius lies in
[1.497,1.5] ell, its taper derivative is exponentially small and
0<=B/A<.003/ell. The wall acceleration factor stays between
-1.37/ell and -1.35/ell. Starting at atanh(.13348702063824), a
time of .01 ell cannot reduce the rapidity to zero. For the inner
cap B/A=1/ell, Ktau_minus=1/ell-.01 ell^2/R^3, and the initial
rapidity atanh(-.28088655126719) cannot reach zero in .01 ell.
The resulting radius bounds R_inner<.806 ell, R_outer<1.517 ell
and label bounds close the collar estimates. These estimates specify
the finite construction, rather than extrapolating its sampled minima.
No infinite-time foliation, stability or regular formation is claimed.

## Integrating the actual moving geometry

Original rays are integrated only over the labels still in their original
region. Each crossed strip is added once, with crossing time s as coordinate:

$$
dV(t;s)=4\pi R(t;s)^2 A(t;s)
             \frac{|v(s)|}{A(s;s)}\,ds.
\tag{N7}
$$

This is ordinary change of variables dx=|v|ds/A at crossing. It includes
the swept volumes that a collection of independent fixed-region averages
would miss. For each material wall the singular integrals are

$$
E_W=4\pi R^2\sigma\gamma,\quad
S_W=4\pi R^2\sigma(\gamma-3/\gamma),\quad
C_W=2\kappa E_W,\quad \kappa=8\pi.
\tag{N8}
$$

They use the actual incoming normal's boost at that event, not a newly
chosen flat radial slice. Bulk E_i=3H_i^2/kappa and S_i=-3E_i retain
their original interpretation. Combining (N2), (N7), and (N8) gives
the same seven moments as the composition theorem, now at common future
times. The Hamiltonian and Raychaudhuri averages have P=0 in this
geodesic clock, including their explicit distributional sources.

Define Q=2 Var(theta)/3-2<sigma^2>, R3=<three-curvature>,
H_D=<theta>/3, and the geometric effective trace

$$
\Lambda_{g,\rm tr}=\frac{Q-R3}{4},\qquad
{\cal A}=3Q+R3=-3\kappa(\rho_g+p_g).
\tag{N9}
$$

The trace including the explicit regional and wall sources is instead
Lambda_total,tr=Lambda_g,tr+kappa(<E>-<S>)/4. These are different
quantities. Neither one is automatically a constant vacuum component.

| t / L_ref | H_D L_ref | Lambda_g,tr L_ref^2 | A L_ref^2 | Lambda_total,tr L_ref^2 |
| ---: | ---: | ---: | ---: | ---: |
| 0 | .0926090487390 | -.0153506330754 | 2.11796240047 | .443183309166 |
| .005 | .0940093784337 | -.0143471868939 | 2.13982145265 | .449242562255 |
| .01 | .0954271417243 | -.0134681818918 | 2.16143093935 | .455231112462 |

At the endpoint the wall radii are .805273647888 ell and
1.516505403648 ell. Their proper times are .009613643764 ell and
.009919162878 ell. This is an actual simultaneous normal-clock cut
of the shared history, rather than equal proper times on the two walls.

At t=.005 ell, direct differentiation of the constructed moments gives

$$
\Lambda_{g,\rm tr,t}\ell^3\simeq .18865034,\qquad
J\ell^3\simeq .0437522915,\quad
J=\langle E\rangle_t+3H_D\langle E\rangle+H_D\langle S\rangle .
\tag{N10}
$$

Inserting the independently differentiated Q,R3,E into the previously
derived conservation identity gives

$$
Q_t+R3_t+2H_D{\cal A}=2\kappa J,\qquad
\Lambda_{g,\rm tr,t}={\cal A}_t/4+H_D{\cal A}-\kappa J .
\tag{N11}
$$

The first residual decreases from 7.22e-9 to 1.80e-9 ell^-3 when the
time-difference step halves from 1e-4 to 5e-5 ell. The independent
H_D derivative/Raychaudhuri difference decreases from 4.05e-9 to
1.02e-9 ell^-2. Thus the previously uncomputed J and trace drift now
come from the same evolving geometry and moving sources.

The 32/64 change in Lambda_g,tr L_ref^2 is at most 1.28e-10.
These are convergence indicators, not observational uncertainty or
rigorous rounding enclosures. The finite prediction remains conditional
on the supplied geometry and the declared normal-clock averaging rule.

### The same evolving wall supplies an observer signal

The inner normal clock is the cap's flat-slicing time, equal to central
proper time. An inward negligible flash emitted at its wall event t_e
therefore has

$$
t_{\rm rec}=t_e-H^{-1}\log(1-HR_W(t_e)),\qquad
Z=\frac{\beta_-+\dot R_W}{1-H^2R_W^2}.
\tag{N-observer}
$$

This reuses the established cap light-cone protocol with its calibrated
wall proper frequency. Its reception domain is supplied by the quiet cap's
already established continuation; the entire normal foliation need not
extend to that later reception time.

| Emission t_e / ell | Reception t_rec / ell | 1+z |
| ---: | ---: | ---: |
| 0 | 1.60943791243410 | 3.74639337610612 |
| .005 | 1.62760502067277 | 3.81927466927820 |
| .01 | 1.64616002849800 | 3.89500911178850 |

Directly differentiating t_rec with wall proper time recovers Z. Thus
the same physical history supplies both the effective-source table and
a delayed clock observable. This is a conditional joint forward result;
the redshift alone does not measure the whole parent average or select
its foliation, and no sky data have yet tested it.

## What happens on the large parent domain?

There is an analytic answer without a cutoff scan. On the original
outer time-symmetric Schwarzschild end, let r denote the initial areal
radius. Its geodesic starts at rest. At any fixed bounded regular time,
the radial equation gives

$$
R=r-\frac{Mt^2}{2r^2}-\frac{M^2t^4}{12r^5}+O(r^{-8}),\quad
A=\partial R/\partial r
 =1+\frac{Mt^2}{r^3}+\frac{5M^2t^4}{12r^6}+O(r^{-9}).
\tag{N12}
$$

Although x is initial proper distance, R_x=E R_r, so A=R_r here.
It follows that theta=O(r^-6), sigma^2=O(r^-6), and R3=O(r^-6).
Their proper integrals over the end converge; the source integrals are
confined to the retained interior and walls. Total volume grows as the
cube of the initial outer cut. Therefore, for this finite evolved family,

$$
\lim_{R_{\rm cut}\to\infty}\Lambda_{g,\rm tr}L_{\rm ref}^2=0,
\qquad
\lim_{R_{\rm cut}\to\infty}\Lambda_{\rm total,tr}L_{\rm ref}^2=0.
\tag{N13}
$$

The regional vacuum value remains nonzero inside its supplied regions.
Equation (N13) concerns the effective spatial average on this isolated
asymptotically Euclidean parent, not a pointwise change of that value.
The same outer Schwarzschild spacetime also retains its exactly zero
local nesting-induced bulk Lambda difference.

This is the connected conclusion for Step 7's present route: nesting
produces a calculable finite geometric response, with nonzero evolution,
but this isolated parent's averaging does not generate a nonzero bulk
cosmological constant. A cosmological parent requires a different large-scale
geometric specification, such as a maintained nonzero regional volume
fraction and a compatible non-AE connector/history, or a genuinely new
source mechanism. A volume fraction chosen to reproduce the observed
Lambda would be a fit, not its derivation.

## Reproduction and remaining connection

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_parent_evolution.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_parent_evolution.py -q
```

The [calculator](../src/recursive_horizons/nsc_atlas_parent_evolution.py),
[application](../scripts/check_nsc_atlas_parent_evolution.py),
[independent controls](../tests/test_nsc_atlas_parent_evolution.py), and
[compact evidence](../results/summaries/nsc-atlas-parent-evolution-v1.json)
include exact de Sitter recovery, parametric Schwarzschild free fall on
both sheets, the E=0 throat, neighbour-ray spatial derivatives across
both walls, the historical initial moments, evolving conservation and
explicit domain/compute stops.

The full plan still needs the cosmological parent, its source selection
and its observer light cone. P001's sourced finite-amplitude convergence,
the quantum response and the infinite operator's physical identification
keep their separate domains. This quiet normal evolution supplies none
of those by relabelling. Its completed connection is shared geometry to
compatible future slices to actual evolving effective sources; the next
outward connection is the cosmological domain and its physical volume law.

The completed validation run passes 68 connected scientific controls and
14 manuscript-renderer controls. The finite application uses 739,848 charged
scalar evaluations and 10.0 process-CPU seconds in the recorded run; its
compact JSON is 11,770 bytes. Runtime varies by host.
