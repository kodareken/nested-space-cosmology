# Valve formation: compression, closed surfaces and a receiving region

This owner covers only roadmap step4: compression sources, circulation,
future trapping and a compatible receiving geometry. Eta is parked at
`dfda6e5a`; horizon/GH/boson refinement and downstream cosmology are frozen
except for the finite evolution needed to test formation. The compact-parent
action is retained. The earlier two-mode preparation and the new one-mode
occupation below remain explicitly distinct states.

## Current result: a same-action trapped neck

The old two-mode finite searches below found no admitted trap; their
occupation was supplied rather than selected. Continuing that optimizer,
adding an external pressure term, and changing the occupation are materially
different routes. The last route retains the action, introduces no coupling,
and supplies a directly calculable formation threshold. It is a declared
new state, not a spherical average of the nonspherical field.

Set chi=epsilon A X1 with kappa A²=24/(z+3), all normal field momenta and
Kij zero, and s=sstar. Its round-reference excitation energy equals the
old two-mode normalization. Its solved proper energy is recalculated.
The [initial calculation](../results/summaries/nsc-atlas-valve-neck.json),
[finite continuation](../results/summaries/nsc-atlas-valve-spherical-v1.json)
and [independent full-GH control](../results/summaries/nsc-atlas-valve-neck-gh-v1.json) supply numerical evidence of a
closed future-trapped region. No regular receiving child is admitted.

### Geometric threshold and local-development statement

With u=X1=cos(chi), the reflection-even conformal factor solves

$$
-8\Delta_0\Psi+6\Psi
-\kappa\epsilon^2A^2(1-u^2)\Psi
-[6+\kappa\epsilon^2A^2zu^2]\Psi^5=0.
$$

The fixed set u=0 is a totally geodesic closed S2. Both initial null
expansions vanish. Write its proper areal radius R=H^(-1)Psi² sin(chi).

Existence is not inferred from Newton convergence alone. Reuse the weighted
equivariant compactness proof in [finite parent data](nsc-atlas-finite-parent-data.md)
with antipodal symmetry, whose minimum orbit size is2. Put
a=kappa epsilon² A²=24epsilon²/(z+3). The elliptic weights are
A_w=6-a(1-u²), B_w=6+azu². For z>1 and epsilon<=1, A_w>0, while

$$
\frac{B_{\max}}{\langle B\rangle}
=\frac{6+az}{6+az/4}<\frac52,
\qquad
\left(\frac{\mathcal Q(1)}{T}\right)^3
=\left(\frac{\langle A_w\rangle}{6}\right)^3
\frac{B_{\max}}{4\langle B\rangle}<\frac58<1.
$$

The inherited Hebey--Vaugon concentration argument therefore supplies at
least one smooth positive antipodally even solution for every0<=epsilon<=1
in this one-mode family. It does not prove uniqueness or enclose the
numerically continued solution. The SO(3)-invariant calculation and the
strict negative slope below retain that distinction.

Independent ADM evolution and radial warped curvature give

$$
\dot\theta_\pm={\rm Ric}_h(\nu,\nu)-\kappa(\nu\chi)^2,
\quad
\frac{{\rm Ric}_h(\nu,\nu)}{H^2}
=\Psi^{-4}(2-4\Psi_{uu}/\Psi),
$$

and the constraint at u=0 reduces this to

$$
\frac{\dot\theta_\pm}{H^2}
=3-\left(1+\frac{\kappa\epsilon^2A^2}{2}\right)\Psi(0)^{-4}.
$$

For an exact nearby smooth constraint solution retaining the strict negative
sign, local Einstein--scalar development and reflection symmetry imply
both expansions are negative on the same sphere for sufficiently small
positive proper time. The residual/refinement calculations below support
that nearby solution; they are not an interval-certified continuum enclosure.

| Initial degree | Numerical epsilon_crit | Independent full conformal RMS |
| --- | --- | --- |
| 16 | .8339066729 | 8.44e-5 |
| 24 | .8339066144 | 1.57e-7 |
| 32 | .8339066143 | 2.52e-10 |

At epsilon=1 the finest neck has Psi=.69449811, versus the sign threshold
.76010912, and dot(theta_pm)/H²=-1.30467886. The overall trace initially
grows at3H²: radial stretching accompanies tangential compression. This
is an internal prepared-field/geometry threshold, not an external minus-side
pressure or a selected offset count n. No physical membrane is introduced.

### Finite continuation and independent full-metric check

The exact SO(3) reduction uses proper synchronous time,

$$
ds^2=-dt^2+q^2d\chi^2+b^2\sin^2\chi\,d\Omega^2,
\qquad q=b+(1-u^2)c,\quad k_a=k_b+(1-u^2)d.
$$

These variables enforce regular poles. Both scalars and their backreaction
are evolved; no homogeneous average replaces the physical fields. The
expansions are theta_pm=-2kb plus/minus2R_chi/(qR). At proper .02:

| Spherical modes | Both theta_pm/H | Independent Hamiltonian RMS/H² |
| --- | --- | --- |
| 24 | -.00707953272 | 8.21e-5 |
| 32 | -.00707954029 | 3.76e-7 |
| 48 | -.00707954032 | 8.86e-11 |

Independent quadrature evaluates constraints outside the retained spatial
nodes. At proper .3 both expansions/H=-.10672195 and the finest Hamiltonian
RMS/H²=7.85e-11. The radius is1.77330393 and radial scale1.80256812.
The source begins on a bi-marginal slice; an inherited earlier preparation
or collapse history has not been derived.

The existing full nonspherical GH engine, supplied the same one-mode
initial data, gives theta_pm/H approximately-.0074003, -.0071164 and
-.0070798434 at degrees12,16,24 and matched neck proper .02. Its finest
discrepancy from the spherical result is3.03e-7, or .0043% of the magnitude;
initial Hamiltonian RMS/H² falls .2644 ->.02763 ->.0002085. All sampled
points of the full reflection S2 have negative expansions. This checks
the symmetry reduction and lapse/proper-time handling without changing
the older saved-parent engine or its source pins.

### The other side: interior time and an actual regional ledger

The areal invariant

$$
1-\frac{2m}{R}=\nabla R\cdot\nabla R
=-\frac{R^2\theta_+\theta_-}{4}
$$

makes grad(R) timelike inside the strict trap. Decreasing R is a local
future time there, while the synchronous metric remains regular through
the marginal boundary. At proper .3 its neck norm is-.0006595557 and
m_geo=.88723676. This is the break-question continuation; it does not
create a new spatial coordinate or prove a new child. At the initial
bi-marginal neck grad(R)=0, so areal time alone is degenerate there.

Reuse [Hayward, Section II, (2.2) and (2.7)--(2.9)](https://arxiv.org/html/gr-qc/9710089v2)
for spherical geometric mass and supply/work. In the present conventions
j=-sum(Pi phi_l) and k_a,k_b are mixed extrinsic-curvature eigenvalues.
The full proper-energy density satisfies

$$
\partial_t(4\pi qR^2\rho)+\partial_\chi(4\pi R^2j)
=4\pi qR^2(k_ap_r+2k_bp_t).
$$

For moving marginal boundaries, retain the boundary term
[4piR²(rho q chi_dot-j)]_left^right. At proper .2 the connected trapped
band has volume16.77650134, energy3.77629301 and mean pressure-.21738048,
so w=-.96573118. Subtracting the supplied reference vacuum gives
excitation energy.06897101 and excitation w=.87628261. Near-vacuum total
stress here does not derive a new Lambda.

| Energy-rate contribution | Action units |
| --- | --- |
| Pressure work | .1584497891 |
| Physical material influx across both band boundaries | .0856197876 |
| Boundary sweeping existing energy | 24.1826846267 |
| Direct derivative | 24.4267542203 |
| Absolute ledger residual | 1.69e-8 |

The marginal boundary moves at1.30618 relative to normal observers; it is
a spacelike marginal tube, not matter moving faster than light. Reflection
makes the total material flux through the central timelike neck exactly0.
The two hemispheres' energy changes close with pressure work, not a supply
from one hemisphere into the other. A band receiving energy is not yet a
regular expanding child. Neither its proper energy nor geometric mass is
identified with Parke E=3/8 or the different88.47% reference transfer.

The [independent horizon classification](../results/summaries/nsc-atlas-valve-horizon-classification-v1.json)
checks that the surrounding marginal tubes are future outer, not inner
cosmological boundaries. With reciprocal null normals of inner product-2,
Hayward surface gravity is kappa_H=m/R²-kappa R(rho-pr)/4 and
L_in(theta_out)=-4kappa_H/R at marginality. At proper .02,.2,.3 the
computed kappa_H is+.0849121,+.0843238,+.0804675, and cross expansion is
negative at both boundaries. Independent directional and finite-difference
calculations agree within1.6e-12 and9.2e-10. The material-current component
changes sign over the carrier oscillation; full horizon crossing includes
its motion through the matter. An outward current relative to a spatial
slice is not proof of causal escape from the trapped region.

The [bounded continuation to proper2](../results/summaries/nsc-atlas-valve-neck-continuation-v1.json)
tests what the pinch becomes. Modes48/64/80 reach the endpoint in under
three seconds each. Finest Hamiltonian RMS/H²=8.80e-11. Neck radius
falls10.71% to1.58681277 while its radial scale grows73.79% to3.08848102;
both expansions/H=-.98374109. Its areal radius remains a spatial minimum,
and neither a bounce nor metric degeneration is observed in this window.
The hemisphere volume grows .57435% and its energy .22352%, with exact
zero central transfer. This is a calculated instance of radial stretching
coexisting with angular compression, not an isotropically inflating child.
Zero central flux does not exclude symmetric accretion through the two
surrounding marginal boundaries. The positive material influx above is
real. The next physical question is whether this receiving trapped band
admits a regular expanding-child continuation.

### Same-state continuation in a local null frame

The [conformal continuation](../results/summaries/nsc-atlas-valve-conformal-v1.json)
changes the numerical frame at the later global loss of accuracy. It
inherits the actual80-mode proper2 Cauchy state. Proper radial distance X
on that slice gives a=0, a_T=-ka, R_T=-kbR and phi_T=Pi for phi=(S,chi).
There is no resetting of fields, energy, occupation or source normalization.
Use ds²=exp(2a)(-dT²+dX²)+R²dOmega². The same Einstein--scalar equations are

$$
R_{TT}-R_{XX}=-\frac{R_T^2-R_X^2}{R}
-\frac{e^{2a}}R+\kappa WRe^{2a},
$$
$$
a_{TT}-a_{XX}=\frac{R_T^2-R_X^2+e^{2a}}{R^2}
-\frac\kappa2\sum_A(\phi_{A,T}^2-\phi_{A,X}^2),
$$
$$
\phi_{A,TT}-\phi_{A,XX}
=-\frac{2(R_T\phi_{A,T}-R_X\phi_{A,X})}{R}-e^{2a}W_{,A}.
$$

A shrinking characteristic stencil stays inside the Cauchy domain of
dependence, discarding one point at each side per step. Its initial
proper-radial half-width is3.45839; the common conformal endpoint is1.14819.
No exterior values or new incoming state are inserted at the patch edge.
All three calculations use the same inherited state and domain.

| Initial points | Center proper time | Radius | Both theta_pm/H | Hamiltonian RMS/H² |
| --- | --- | --- | --- | --- |
| 2001 | 3.945529 | .65495433 | -11.85692 | .0379729 |
| 4001 | 3.945509 | .65496885 | -11.85045 | .0094763 |
| 8001 | 3.945504 | .65497252 | -11.84883 | .0023727 |

Successive radius differences shrink by3.97; full fine-minus-middle spread
is3.66e-6. This is an empirical refinement indicator, not a certified error.
The finest constraint residual is about2e-5 of the local curvature scale;
the reference-H normalization alone obscures the growing physical curvature.
Exact de Sitter identities check the PDE factors to1.56e-15, and its field
errors converge at second order. Tiny finite-difference vacuum constraint
noise is nonmonotone and is not sold as an enclosure. Proper-time-aligned
overlap at2.2 differs from the previous frame by1.22e-8 in radius and
1.79e-7 in theta/H; full-vs-decimated actual-history clock indicators are
reported separately. The entire retained local patch is future trapped,
but it does not include the whole receiver's moving outer boundary.
No global child volume, cosmological epoch or new vacuum follows from it.

### What an actual zero radius would mean in both frames

At the reflection center, R_l=s_l=0. The sphere's spacetime sectional
curvature is K_Omega=(1+R_tau²)/R²>=1/R². Reflection eliminates the
odd-time mixed components in the spherical Riemann contraction, giving

$$
\mathrm{Kretsch}_4\geq4K_\Omega^2\geq4R^{-4}.
$$

Thus a genuine zero of this reflection sphere cannot be a regular center
or a removable coordinate zero. This statement does not assume zero is
reached by the finite computed history.

The [declared internal reconstruction](nsc-atlas-inherited-charge.md)
is ds16²=exp(-12s)g4+exp(2s)gamma12, with R[gamma]=15. Put x=exp(2s),
so the external areal radius is x^(-3)R. Its gradient is timelike or zero
at reflection, hence K_ext,16>=x^6/R². In12 dimensions positive scalar
curvature15 guarantees at least one internal sectional plane with
K_gamma>=15/(12*11)=k=5/44. The internal radius gradient is likewise
timelike or zero there, so K_int,16>=k/x. Consequently, for every volume roll,

$$
\max(K_{\rm ext,16},K_{\rm int,16})
\geq k^{6/7}R^{-2/7},\qquad
\mathrm{Kretsch}_{16}\geq4k^{12/7}R^{-4/7}.
$$

The bound follows by minimizing max(x^6/R²,k/x) over x>0. Odd-time mixed
curvature components also vanish in this reflection-symmetric warped
ansatz, so their Lorentz signs cannot cancel the positive square. Roundness
of gamma is unnecessary for this obstruction; its fixed positive scalar
curvature suffices. A volume/conformal redistribution cannot regularize
both factors if R actually tends to zero within this ansatz.

If the spatial minimum, radial stretching and canonical null-energy signs
persist, R_tautau<=0 and the endpoint values conditionally imply no more
than .622 further proper units before zero. Persistence has not been proved.
The next physical test is a neck-type/receiving-geometry transition before
that limit; neither coordinate failure nor a curvature blow-up is renamed
as a child. This is a scoped constraint on the current candidate, not on
all nested geometries or other published regular-interior mechanisms.

### Scope and next connection

The conditional source has a resolved future-trapped development and a
checked local receiving-band budget. Its inherited occupation, nonzero
directional supply history, regular child continuation and rotational formation
mechanism remain open. The completed threshold and this finite window
are frozen. The next calculation must follow the actual receiving band
and its horizons, retaining the parent-supplied state. A directional asymmetry,
if required, must be inherited rather than inserted to force a result.

A [later bounded probe](../results/summaries/nsc-atlas-valve-neck-late-probe-v1.json)
loses global synchronous constraint accuracy near proper2.9 before any
trustworthy bounce or neck-type change. Subsequent numerical metric
degeneration is not a physical endpoint. The next vantage point is a local
conformal/double-null characteristic development from trusted proper2
Cauchy data, with no new incoming boundary history. Its causal patch must
remain inside the inherited data domain. No further global resolution scan
is admitted by this diagnostic.

At a reflection neck the exact spherical Einstein equations also give
R_tt+R_ll+ka R_t=-kappa R(rho+pr)/2. If R_t=0 and the neck remains a
strict spatial minimum R_ll>0, canonical rho+pr>=0 implies R_tt<0.
A temporal bounce retaining that neck type is excluded under this action.
Changing coordinate names cannot remove this invariant sign condition;
a continuation must change neck type or its source assumptions.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_valve_neck.py
python scripts/run_python.py scripts/check_nsc_atlas_valve_spherical.py
python scripts/run_python.py scripts/check_nsc_atlas_valve_neck_gh.py --degree 24 --threads 4 --budget 180
python scripts/run_python.py scripts/check_nsc_atlas_valve_neck_continuation.py 80
python scripts/run_python.py scripts/check_nsc_atlas_valve_conformal.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_valve_neck.py tests/test_nsc_atlas_valve_spherical.py tests/test_nsc_atlas_valve_conformal.py
```

## Calculation and stopping rules

1. Derive a dimensionless geodesic focusing ratio from the actual stress;
   distinguish it from the sufficient conditions for a trapped surface.
   Test whether external pressure and accumulated offset winding have
   separate stress terms in the current action before treating them as
   compression parameters.
2. Evaluate closed two-surfaces with the full spatial metric and extrinsic
   curvature. Solve the deformed theta_plus=0 equation and independently
   evaluate theta_minus. A cosmological marginal surface is not admitted
   as a future black-hole boundary.
3. First numerical window: epsilon=1, zero initial normal momenta and Kij=0,
   degrees12,16,20, proper interval0 to.3 in the existing model units.
   This reaches beyond one massive carrier period. It is an initial
   formation window, not a claim about the entire future. Use four threads
   per evolution, two concurrent workers, 900 wall seconds per case,
   rtol3e-9, atol3e-11 and maximum coordinate step.01. Passive matter-work
   accumulation is omitted; the coupled Einstein/scalar equations are
   unchanged. No observational quantity selects these numerical settings.
4. Compare constraints, both expansions and the full surface residual at
   independent denser angular nodes. Surface residual convergence alone
   does not establish evolution convergence. Search multiple initial
   radii and centers; a search miss is not a proof that all closed surfaces
   or all admissible states are untrapped.
5. A resolved future trap triggers the interior continuation and received
   energy/volume/pressure calculation. Otherwise report the precise
   tested domain, the missing source or unresolved dynamics, and the next
   necessary change. Neither a solver divergence nor a coordinate horizon
   supplies an invented child. The break question requires a regular
   continuation or a diagnosed invariant obstruction.

The likely obstacles are confusing local focusing with trapping,
mistaking a cosmological root for a local future trap, and trusting an
under-resolved metric because a surface solve converges. The independent
null signs, known-geometry controls and three metric resolutions test
those separately. A later finer evolution or a surface outside the tested
shape/chart family can overturn a finite search conclusion.

## Fixed-action pressure diagnostic

For the two canonical real fields and full potential W,

$$
R_{nn}=\kappa(\Pi_S^2+\Pi_\chi^2-W).
$$

For a hypersurface-orthogonal geodesic normal congruence,

$$
n(\theta)=\kappa W-\kappa(\Pi_S^2+\Pi_\chi^2)
-\sigma_{ij}\sigma^{ij}-\theta^2/3.
$$

When W>0, the normalized right-side competition is

$$
\mathcal F=\frac{\Pi_S^2+\Pi_\chi^2}{W}
+\frac{\sigma_{ij}\sigma^{ij}+\theta^2/3}{\kappa W}.
$$

Expansion decreases when F>1. This is a geodesic focusing criterion;
it is not a universal TOV, Jeans or black-hole mass threshold. Accelerated
slice normals also require their acceleration-divergence term. At the
time-symmetric reference event F=0 at every amplitude in the existence
theorem. The potential and gradients do not supply an initial inward
normal velocity by themselves.

The horizon is smooth with no surface stress. An independently adjustable
external pressure has not been specified. The outside's already present
bulk gravity is included in the common parent; assigning a new stress
to the minus label would change the model. The phase accumulator likewise
does not appear as an independently winding material field or forcing in
the action. Its index is an event label, not an accumulating energy source.
The actual oscillating scalar kinetic stress remains a possible dynamical
focusing source and is included in the finite evolution.

## An internal pressure-release threshold already in the action

The retained volume potential and its existing wave coupling supply a
specific candidate release mechanism without an added force or fit:

$$
V=7e^{-12s}-\frac{15}{2}e^{-14s}+\frac{25}{32}e^{-24s},
\qquad W=V+200\chi^2e^{-12s}.
$$

There are two different, explicitly conditional tests. At fixed
instantaneous chi, the stationary-volume polynomial has its unique
positive minimum at x=(14/15)^(1/5), where x=exp(-2s). Both extrema
disappear when

$$
\chi^2>\chi_c^2
=\frac{(175/24)(14/15)^{1/5}-7}{200}
=.000958715247985449.
$$

At the prepared scalar maximum, this occurs for epsilon>.571125282882.
The wave oscillates, so an instantaneous loss of the well is not evidence
of permanent escape.

For slow volume evolution, the correct other test holds the fast
oscillator action fixed rather than holding chi fixed. Neglecting local
spatial transport and holding external volume fixed, its density scales
as rho_chi=M exp(-6s). The effective volume potential is
V_eff=V+M exp(-6s). Its stationary points satisfy

$$
\mathcal M(s)=-14e^{-6s}+17.5e^{-8s}-3.125e^{-18s}.
$$

The minimum and barrier merge at M'(s)=0. Thus

$$
225x^6-560x+336=0,\qquad
x_{\rm fold}=.998727268724,\qquad
s_{\rm fold}=.000636770943,
$$

$$
\mathcal M_{\rm crit}=.375079333827.
$$

The unchanged initial two-mode source has exact maximum
f_max=(a²+b²)/(2b)=.0542142471364. At epsilon=1, its local oscillator
loading is M_local=(mu_star² f_max²/2) exp(6s_star)=.728021263069.
The dimensionless release test is therefore

$$
\mathcal B=\frac{\mathcal M_{\rm local}}{\mathcal M_{\rm crit}},
\qquad \mathcal B>1,\qquad
\epsilon_{\rm crit}=.717776996933.
$$

The epsilon=1 preparation has B=1.94097940732; the old epsilon=.25
control has B=.121311212957. This is a derived local adiabatic
stabilization threshold in the fixed action. Expansion, transport,
oscillation and gradients determine whether it actually releases the
volume. It is separate from the F>1 focusing test and theta_plus=0
surface test. No one of these substitutes for the other two.

The older Parke reference E0=3/8 is .0211512% below M_crit. These are
compared reference coefficients, not a calculated energy received by a
new child. The finite parent's normalized total excitation is another
quantity; it is not identified with E0 merely because numbers are close.

The unexcited bare maximum is s_top=.060000535042 with
V_top=.354536677477. The existing reconstruction is
ds_16²=exp(-12s) ds_4,E²+exp(2s) gamma_internal. Its radius and volume
factors are exp(s-s_star) and exp[12(s-s_star)]. This supplies a precise
internal-volume channel to track. It does not construct a new
three-dimensional regional child. The full evolution, regular higher
geometry and omitted internal modes retain their own requirements.

## Closed-surface definition

With the inherited convention Kij=-(Lie_n h)ij/2 and outward unit spatial
normal v, evaluate

$$
\theta_+=H_S+K(v,v)-\operatorname{tr}K,\qquad
\theta_-=-H_S+K(v,v)-\operatorname{tr}K.
$$

The finder uses a finite harmonic graph on S2 in geodesic polar charts of
S3 and its actual physical metric. Its independent validation uses denser
angular nodes. Strictly negative theta_minus is required at a candidate
theta_plus root, followed by checks on nearby surfaces and metric
resolution. A finite search does not certify the global outermost surface
on a compact parent.

An independent necessary sign test is useful on each leaf. If lambda1,
lambda2,lambda3 are the ordered physical eigenvalues of K, then for
every surface normal v,

$$
\operatorname{tr}_S K=\operatorname{tr}K-K(v,v)
\leq\lambda_2+\lambda_3.
$$

Where that upper bound is nonpositive, theta_plus+theta_minus cannot
be negative for any orientation. A sampled spatial bound is a numerical
diagnostic, not a certified global no-trapping theorem. This test never
replaces the actual closed-surface equation.

The expansion equation and elliptic-surface approach are published inputs:
[Thornburg's review](https://doi.org/10.12942/lrr-2007-3) and
[AHFinderDirect](https://einsteintoolkit.org/thornguide/EinsteinAnalysis/AHFinderDirect/documentation.html).
The regular-frame S3 implementation, physical-metric adapter and this
specific parent comparison are the atlas's calculations.

## Computed initial window

The [source-bound three-grid calculation](../results/summaries/nsc-atlas-valve-formation-early-v1.json)
reaches proper.3 with epsilon=1. These are new finite-window runs, not
the earlier under-resolved long histories. At proper.075 the fixed scalar
maximum has F=5.66663596 on the finest metric: intrinsic scalar kinetic
stress develops enough local focusing. Negative slice-normal expansion
also develops, so the initially expanding state does not stay everywhere
expanding. This does not yet establish a trapped closed surface.

| Metric degree | Final Hamiltonian RMS/H_reference² |
| --- | --- |
| 12 | .0338456813 |
| 16 | .00262947587 |
| 20 | .000175834904 |

The surface search at proper.225 uses both the scalar maximum and a
scalar node as centers, three round starts and deformed shape degrees4
and8. No sampled future-trapped or future-MOTS candidate is found in
this search. The finest maximum-centered candidate has full outgoing
RMS 7.85824e-5 and ingoing expansion range [-.0891827,+.0786602] in model
units. The node-centered candidate likewise has mixed ingoing signs.
Neither reaches the prescribed max-residual tolerance or qualifies as a
future horizon. No additional precision is pursued on these nonfuture
candidates.

The same maximum-centered graph is evaluated on all three metrics. Its
adjacent outgoing-expansion RMS gaps/H_reference shrink
.000295633 -> .0000227065, and ingoing gaps shrink
.000293449 -> .0000241497. The corresponding node-centered gaps contract
as well. This checks physical metric/field resolution separately from
surface fitting; it is not a continuum no-trapping proof.

The [passive dense-history onset calculation](../results/summaries/nsc-atlas-valve-onset-events-v1.json)
locates the first two local events without another evolution. The
geodesic-relaunch focusing score reaches 1 at proper duration .02822386;
the slice-normal trace at the tracked observer subsequently reverses
from expansion to contraction. Its three-resolution values are:

| Metric degree | First normal contraction, proper duration |
| --- | --- |
| 16 | .0538292690 |
| 20 | .0538328333 |
| 24 | .0538334449 |

The adjacent contraction-time differences decrease 3.56434e-6 -> 6.11599e-7.
Differences in the earlier F=1 event lie below the declared 1e-8 root
tolerance and are not treated as resolved extra digits. The early
observer-phase readout gives 1.94076 and 3.70068 protocol passes at the
two finest events, respectively. These are synchronization counts, not
an inherited nesting depth or a formed-child threshold. The measured
trace is that of slice normals at the geodesic observer's position, not
the expansion of a separately evolved geodesic bundle.

Independent new histories extend to proper1.2, with dense polynomial
outputs, at degrees12,16,20,24. The degree12/16 late constraints become
large; the finer runs are used to identify the interval where the geometry
is reliable. Repeating the proper.3 checkpoint with maximum step.025
instead of.01 supplies the time-step comparison without another physical
preparation. A later unresolved leaf is not admitted merely because its
integrator reached the endpoint.

## Circulation and the inward-channel question

The existing current obeys dot J+3hJ=5q1q2/a² and sources a regional
vector metric response. Its prepared coexact angular source is orthogonal
to the six reference rotation Killing fields: local circulation is
compensated globally. This permits local counter-rotation but does not
select a unique inward axial valve. A future collapsing channel must be
tested with its actual inward flux and null expansions, not inferred from
the presence of a nonzero circulation coefficient.

The [same-state vortex calculation](../results/summaries/nsc-atlas-valve-vortex-v1.json)
now makes a more specific connection. Use the canonical quadrature map
Z=chi+i Pi_chi/mu_reference, with positive constant mu_reference. In the
two-mode comparison, D=q1p2-q2p1 is its transverse orientation determinant.
For D nonzero, its zeros are the two S3 circles X1=X2=0 and X1=X4=0.
Their fixed-point source symmetries force chi=Pi_chi=0 there in the full
nonlinear evolution too, though the nonlinear winding must be checked.
The geometric eigenvalue gap generates
dot D+3hD=-5q1q2/a² from initially zero D.

This phase geometry is tied to the actual momentum density:

$$
j_\chi=-\Pi_\chi\,d\chi,\qquad
dj_\chi=d\chi\wedge d\Pi_\chi,
\qquad \oint j_\chi=\frac12\oint(\chi\,d\Pi_\chi-\Pi_\chi\,d\chi).
$$

Both sampled cores and both transverse radii have winding -1 at metric
degrees 16,20,24 at proper .9. Finest 256/512-point loops confirm the sampled
index with a nonzero boundary margin; this is not a certified continuum
nonvanishing bound. The phase winding is not itself a quantized material
circulation. For example chi=x,Pi=y gives a phase vortex while j=-y dx
is a shear flow. Circular fluid streamlines and an inward drain require
additional measurements.

At the two selected cores, the chi current vanishes but the volume field
has opposite axial momentum densities of magnitude 1.92824e-4 in model
units. Thus the quadrature zero is not an endpoint of the full stress:
the other field carries a finite current there. The cores' normal
expansion is about +.171428. A direct three-resolution round-surface search
about these source-derived centers finds no sampled future trap among
the thirteen radii .2 to 2.6. This does not rule out other deformed surfaces,
other centers or later events. It also does not identify a child merely
because the phase coordinate is undefined at its core.

[The 2026 forced-fluid construction](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
Theorem1.1 and its early construction, supplies radial inflow with axial
stretching/outflow in an incompressible constant-density fluid and smooth
forcing. Those inspected passages do not construct the present source,
a relativistic horizon or a child beyond the singular limit. The broader
proof is not used as an NSC formation theorem.
[Bredberg and collaborators](https://arxiv.org/abs/1101.2451) supply a
fluid/gravity map with specified boundary and expansion assumptions;
those assumptions do not automatically hold in the strong-field limit
tested here. The useful shared question is concentration and transfer;
their mathematical source and regime remain explicit.

Published [massive-scalar collapse](https://arxiv.org/abs/gr-qc/9709014)
and [rotating scalar collapse](https://arxiv.org/abs/gr-qc/0405101)
show how preparation, pulse scale and angular data affect thresholds.
Their isolated/spherical or complex-field assumptions are not silently
substituted for this nonspherical compact two-real-field state.

## The break question in the declared internal frame

The breathing-mode reconstruction provides a concrete second vantage
point for the same fields. For a deformed surface Sigma2 wrapped over the
fixed internal compact geometry M12, its full area density is

$$
dA_{14}=e^{-12s}dA_E\,e^{12s}dV_{\rm int,0}
=dA_E\,dV_{\rm int,0}.
$$

With the null normals normalized as n plus/minus the unit surface normal,
k_16=exp(6s) k_E. Consequently

$$
\theta_{16,\pm}=e^{6s}\theta_{E,\pm}.
$$

The wrapped trapping signs are exactly preserved, even with nonuniform s
on the surface. Looking only at external Jordan-frame two-area would
instead give

$$
\theta_{J,\pm}=e^{6s}[\theta_{E,\pm}-12k_{E,\pm}(s)].
$$

Its timelike three-volume trace is
vartheta_J=exp(6s)[vartheta_E-18n_E(s)]. Internal expansion and external
three-volume contraction can therefore coexist. The full area identity
prevents this reallocation of volume from being misidentified as a new
wrapped trapped surface. These are exact kinematic identities of the
block-diagonal ansatz; localization, omitted internal modes and a new
regional child still require their own dynamical solution.

At the degree24 F=1 event, the same state gives
vartheta_E/H_reference=.07507030 and
vartheta_J/H_reference=-.10515381, with n(s)=.00309498365 in model units.
At the subsequent normal-contraction onset, vartheta_E=0 while
vartheta_J/H_reference=-.21686585. These are a calculated example of
internal-volume growth coexisting with external three-volume contraction
in the declared reconstruction. The wrapped null signs remain governed
by the full area identity, so neither number establishes a hidden child.

No child has been admitted in the completed window. Received child
energy, volume, pressure and w remain undefined. The already computed
88.47% is a regional Killing-energy crossing in the earlier prepared
background, and E0=3/8 is a separate reference coefficient. Neither is
reported as the energy of an unconstructed child. This finite search
also establishes no universal prohibition for all amplitudes, surfaces
or future times, so there is no justified one-line universal blocker.

## Extended search, current limit and next connection

The [extended finite diagnostic](../results/summaries/nsc-atlas-valve-formation-extended-v1.json)
uses metric degrees 16,20,24 at matched proper .9. Hamiltonian residuals
there decrease .00680 ->.00169 ->.000221 in reference H² units. Its
maximum-centered candidate has outgoing expansion range [-.0277647,.0233364]
and ingoing range [-.215176,.293958]. The node-centered candidate also has
mixed signs. No prescribed seed, finite deformed fit or sampled core
surface supplies an admitted future trap. The shape4 fits remain
nonconverged; neither optimizer success nor a small average is a horizon.
The proper1.2 leaf has larger residuals and remains a diagnostic rather
than a new formation claim.

The results are: a local source-driven focusing threshold, actual local
contraction, an internal stabilization-loss threshold, and phase vortices
with finite longitudinal volume-field current. A child receiving energy
has not been constructed. The current sources are a broad standing
preparation; no independently specified external compression or inward
characteristic history is selected by the offset map. The finite search
cannot establish that no admissible state ever forms a trap, so it cannot
justify a universal blocking assumption or a necessary action change.
A localized inward preparation within the same action is a possible next
candidate, provided its constraints and parent inheritance are derived.

Before another global resolution campaign, the next formation connection
is the null-generator history on the existing resolved region. With
l_plus/minus=n plus/minus the unit surface normal,

$$
R_{ll}=\kappa\sum_A(\Pi_A\mathbin{\pm}v^iD_i\phi_A)^2.
$$

The potential cancels from this null projection, while still affecting
the field and opposite-expansion evolution. Integrated Raychaudhuri
requires the actual null shear, affine or boost weight and a common
closed cut with negative opposite expansion. Timelike K shear and the
old 88.47% Killing-energy fraction do not supply these data.
[Massless spherical positive-Lambda formation](https://arxiv.org/abs/2005.03434),
[anisotropic vacuum short pulses](https://arxiv.org/abs/1302.5951) and
[nonsymmetric massless-scalar formation](https://arxiv.org/abs/2304.01695)
provide criteria to reuse when their source, characteristic and norm
hypotheses hold. Those inspected theorem/setup sections do not certify
this massive compact standing state. In particular its old modes have
mu/k about 52.69 and 32.27, outside a massless short-pulse scale limit.

Reproduction uses the source-bound evolution command already admitted
above, saving sparse and dense caches under .build; then:

```powershell
foreach ($valveDegree in 12,16,20) {
  python scripts/run_python.py scripts/evolve_nsc_atlas_parent_gh.py --degree $valveDegree --epsilon 1 --proper-endpoint 0.3 --threads 4 --budget-seconds 900 --rtol 3e-9 --atol 3e-11 --max-step 0.01 --no-work --state-cache ".build/valve-e1-L$valveDegree-early-v1.npz" --output ".build/valve-e1-L$valveDegree-early-v1.json"
}
foreach ($valveDegree in 16,20,24) {
  python scripts/run_python.py scripts/evolve_nsc_atlas_parent_gh.py --degree $valveDegree --epsilon 1 --proper-endpoint 1.2 --threads 4 --budget-seconds 2400 --rtol 3e-9 --atol 3e-11 --max-step 0.025 --no-work --state-cache ".build/valve-e1-L$valveDegree-extended-v1.npz" --dense-cache ".build/valve-e1-L$valveDegree-extended-v1-history.npz" --output ".build/valve-e1-L$valveDegree-extended-v1.json"
}
python scripts/run_python.py scripts/analyze_nsc_atlas_valve_onset.py --output .build/reproduced-valve-onset.json
python scripts/run_python.py scripts/check_nsc_atlas_valve_vortex.py --core-surfaces --output .build/reproduced-valve-vortex.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_valve_formation.py tests/test_nsc_atlas_trapped_surfaces.py tests/test_nsc_atlas_trapped_parent.py
```

The surface producer takes explicit record paths and independently
admitted raw-cache SHA256 values. The compact records retain the original
paths, source commits, parameters, checkpoints, numerical tolerances and
cache hashes; full fields regenerate locally. The goal/frontier and all
other open programme items remain in [the progress record](atlas-progress.json).

## Completed null-sheet searches of the older occupation

The [compact causal record](../results/summaries/nsc-atlas-valve-null-search-v1.json)
preserves actual geometric-null transport on the source-bound dense parent
histories, both expansions of each common cut, and the time-tilted search.
Its curvature comes from metric jets; it is not replaced by the canonical
source to force focusing. Direct Ricci/source disagreement is reported
independently. Contracting outgoing rays coexist with positive or mixed
opposite expansion. The best degree4 time-profile search still has maximum
expansion/H=+.156706, not a future trap. This is a finite local search,
not a nonexistence theorem for the parent or all preparations. Raw rays and
dense arrays regenerate locally; the old searches are frozen.

A useful conditional geometric obstruction avoids further blind fitting.
Let a smooth caustic-free null sheet have cuts with reciprocal null normals
L,Lbar, g(L,Lbar)=-2 and angular metric q. A graph cut lambda=F(x) has, at
a critical point of F,

$$
\bar\theta_F=\bar\theta_{\rm base}(F,x)+2\Delta_qF.
$$

At the global minimum of F, Delta_q F>=0. Uniform nonnegative base incoming
expansion across a sheet window therefore rules out a strictly future-trapped
closed graph cut in that window. A connected closed spacelike cut of an
embedded S2 times interval sheet is a graph: its angular projection is a
local diffeomorphism, hence a covering of simply connected S2. In coordinate
time with L=N(n+v), the normalized shift is2N Delta_q F, with N>0.
The flat-null-plane and static-S3 controls check the sign independently.
Finite sampled positive expansions do not establish the uniform hypothesis;
this lemma is not promoted to an all-sheet physical no-go theorem.
