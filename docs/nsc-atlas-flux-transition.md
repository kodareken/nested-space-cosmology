# An inherited vacuum source, two moving boundaries and an expanding parent

The [programme reassessment](nsc-cosmological-atlas.md#programme-reassessment-the-generative-and-cosmological-connection)
asks for a physical parent-to-child source and its observer consequence.
The retained conformal-field support and trace identities already exclude
the proposed additional classical bulk vacuum in that candidate. This
calculation therefore takes the specified alternative: an Einstein/four-form
action with charged membranes. It does not rename the former Schur response.

**Result.** One action fixes the vacuum changes across two nested boundaries,
their normal forces and their common intermediate spacetime. That spacetime
has a regular core and a finite joint evolution, and its exterior joins an
expanding dust-plus-vacuum parent. A local free-fall curvature measurement
isolates the inherited vacuum from the shells' mass. The absolute vacuum,
charge and initial state remain inputs; the observed cosmological constant
has not been selected.

![One flux law, a shared evolution and a curvature measurement](figures/nsc-atlas-flux-transition.png)

## Borrowed construction and finite calculation

Brown--Teitelboim's membrane mechanism is the input, not an NSC discovery.
[Bousso--Polchinski](https://arxiv.org/html/hep-th/0004134), Sections 2.1--2.3,
specifies the four-form action, fixed-flux boundary term and the distinction
between a quantized jump and an absolutely quantized flux.
[Hirano](https://arxiv.org/html/1804.09985), equations (2.1)--(2.4) and
(2.28)--(2.31), gives the coupled flux, junction and membrane equations.
Our extrinsic curvature has the opposite sign to Hirano's convention.
The [existing parent owner](nsc-atlas-parent-observables.md) supplies the
Einstein--Straus matching method; here the same vacuum also continues into
the dust region.

The missing connection is simultaneous compatibility across two physical
transitions and the exterior cosmological join, followed by an observer
relation. Inputs are fixed below before calculating. Budget: two nine-variable
ODE solves, at most 12000 RHS calls each, six exact symbolic identities,
one independent 40-digit separated-quadrature calculation, and 120 process
CPU seconds. No field grid, data fit, parameter scan or nucleation simulation.
Signed junction or conservation failure rejects this assembly. Success
completes this finite connection; it does not authorize further wall
refinement without a new question.

## The action determines the source and the inheritance law

Use signature (-+++), c=hbar=1 and kappa=8 pi G. With F=dA a four-form,

$$
S=\int\sqrt{-g}\left[
\frac{\mathcal R-2\Lambda_b}{2\kappa}
-\frac{F_{\mu\nu\rho\sigma}F^{\mu\nu\rho\sigma}}{2\,4!}
\right]d^4x
-\sigma\sum_j\int_{\Sigma_j}\sqrt{-\gamma}\,d^3y
+e\sum_j\int_{\Sigma_j}A+S_{\rm dust}+S_{\rm boundary}.
\tag{V1}
$$

S_boundary contains the gravitational boundary terms and the fixed-flux
Legendre term. Vary before inserting the on-shell four-form: omitting that
boundary term and varying a substituted on-shell Lagrangian gives the
wrong vacuum sign. The charge coupling is metric-independent and adds
no separate surface tension. The fixed surface stress is S_ab=-sigma gamma_ab.

Away from a membrane, d(*F)=0. Write F=f epsilon. Direct metric variation gives

$$
T^F_{\mu\nu}=\frac1{3!}F_{\mu\alpha\beta\gamma}
F_\nu{}^{\alpha\beta\gamma}
-\frac{g_{\mu\nu}}{2\,4!}F^2
=-\frac{f^2}{2}g_{\mu\nu},\qquad
\Lambda(f)=\Lambda_b+\frac{\kappa f^2}{2}.
\tag{V2}
$$

Orient each membrane so its inside flux is f_in=f_out-e. Across consecutive
equal-charge transitions,

$$
f_n=f_0-ne,\qquad
\Lambda_n=\Lambda_b+\frac{\kappa}{2}(f_0-ne)^2,
\tag{V3}
$$
$$
\Lambda_{n+1}-\Lambda_n=-\kappa e(f_n-e/2),\qquad
\boxed{\Lambda_{n+2}-2\Lambda_{n+1}+\Lambda_n=\kappa e^2.}
\tag{V4}
$$

The composition f_(n+m)=f_n-me is exact for any supplied finite sequence.
It does not establish which transitions occur, unlimited spacetime
continuation, or absolute quantization of f_0. V4 is a corollary of the
borrowed source law. Its second difference removes the bare vacuum and
the parent flux; its first difference does not.

In these units [kappa]=L^2, [f]=[e]=L^-2, [sigma]=L^-3 and [Lambda]=L^-2.
No change of units selects a coupling. The same action holds in every
region, while the state f changes. This is a physical transition map;
identification with the earlier inherited Schur operator remains open.

## Both junction equations and distributional conservation

Each vacuum region is Kottler,

$$
ds_n^2=-F_n(r)dT_n^2+\frac{dr^2}{F_n(r)}+r^2d\Omega^2,\qquad
F_n=1-\frac{2M_n}{r}-\frac{\Lambda_n r^2}{3}.
\tag{V5}
$$

M_n is a geometric mass, not a freely chosen second copy at the next
interface. Use a normal pointing from inside to outside and
K_ab=+h_a^c nabla_c n_b. On the ordinary orientation,

$$
\beta_-=\sqrt{\dot R^2+F_{\rm in}},\quad
\beta_+=\sqrt{\dot R^2+F_{\rm out}},\quad
\beta_--\beta_+=kR,\quad k=\kappa\sigma/2.
\tag{V6}
$$

Dots denote that membrane's proper time. Let delta_M=M_out-M_in and
delta_Lambda=Lambda_out-Lambda_in. Then

$$
A=\frac{\delta_M}{k},\quad B=\frac{\delta_\Lambda}{6k}+\frac k2,\quad
\beta_-=\frac A{R^2}+BR,\quad
\beta_+=\frac A{R^2}+(B-k)R,
$$
$$
V(R)=1-\frac{A^2}{R^4}
-\frac{2(M_{\rm in}+AB)}R
-\left(\frac{\Lambda_{\rm in}}3+B^2\right)R^2,\qquad
\dot R^2+V=0,\quad \ddot R=-V'/2.
\tag{V7}
$$

The temporal equation is checked independently using
K^tau_tau=(ddot_R+F'/2)/beta. Together with V6 it gives
[K^a_b]=-k delta^a_b. The normal force is

$$
\sigma\overline K=e\,\frac{f_{\rm out}+f_{\rm in}}2
=\frac{\Lambda_{\rm out}-\Lambda_{\rm in}}{\kappa}.
\tag{V8}
$$

This is also the distributional normal conservation equation:
[T_nn]=-e(f_out+f_in)/2=-sigma K_bar. Tangential exchange is zero because
the bulk source is proportional to the metric. Hence constant tension
satisfies d(sigma Area)+p d(Area)=0 with p=-sigma. The four-form pressure,
wall force and Einstein jumps are one conserving system.

Given inner mass, radius and proper speed, the unsquared equation fixes

$$
M_{\rm out}=M_{\rm in}
+\frac{(\Lambda_{\rm in}-\Lambda_{\rm out})R^3}{6}
+kR^2\sqrt{\dot R^2+F_{\rm in}(R)}
-\frac{k^2R^3}{2}.
\tag{V9}
$$

This is applied from the regular center outward. It is not two independently
chosen exterior masses. Positive beta_+ is retained; squaring alone would
admit unwanted orientations.

## A composed finite geometry with calculated vacuum changes

Set an independent reference length L_ref and declare

$$
\Lambda_bL_{\rm ref}^2=0,\quad
\sqrt\kappa f_0 L_{\rm ref}=3/2,\quad
\sqrt\kappa eL_{\rm ref}=1/2,\quad kL_{\rm ref}=1/10.
$$

Both initial proper speeds vanish. Inner radius is 0.4 L_ref, outer radius
0.8 L_ref, and the core mass is zero. All membrane couplings are shared.
No observation selected these benchmark inputs.

The reference calculation uses normalized variables
r_hat=r/L_ref, f_hat=sqrt(kappa) L_ref f and
sigma_hat=kappa L_ref sigma. The classical action has an overall factor
L_ref^2/kappa, which drops out of these field equations. The record's
kappa=1 is the coefficient in those normalized equations; it does not
set a physical Planck-length ruler or remove that factor from a future
quantum nucleation exponent. Choosing an SI L_ref still requires physical
input. The general evaluator also supports dimensionful unit rescalings.

| Region, outside to inside | sqrt(kappa) f L_ref | Lambda L_ref^2 | M/L_ref |
| --- | ---: | ---: | ---: |
| Parent vacuum collar | 1.5 | 1.125 | .01523145571612779 |
| Shared child region | 1 | .5 | .01162657748024113 |
| Regular grandchild core | .5 | .125 | 0 |

Thus the inward vacuum changes are -0.625 and -0.375 in L_ref^-2.
The next curvature is not a separately adjustable cap input. The scalar
vacuum stress is nonzero and has exactly vacuum form throughout each bulk.
Nevertheless, the outer vacuum is inherited from the supplied f_0; adding
an inner membrane does not change the outside flux or create acceleration
in an otherwise unchanged parent.

One initial trial with outer radius L_ref instead of 0.8 L_ref has outer
geometric mass -.00753563494260208 L_ref. It cannot be mass-compensated
by positive dust. Reducing the initial radius gives the admissible positive
mass above. This is a matching restriction, not observational calibration.

**Common clock and finite existence.** Both boundaries use the same middle
Kottler T. Their individual proper-time rates are

$$
\frac{d\tau_{\rm outer}}{dT}
=\frac{F_1(R_{\rm outer})}{\beta_{-,{\rm outer}}},
\qquad
\frac{d\tau_{\rm inner}}{dT}
=\frac{F_1(R_{\rm inner})}{\beta_{+,{\rm inner}}}.
\tag{V10}
$$

For this benchmark A, D=M_in+AB and C=Lambda_in/3+B^2 are positive.
The acceleration derivative is
10 A^2/R^6+2D/R^3+C>0. At the initial radii the outer acceleration
is +.9991346271/L_ref and the inner acceleration is -2.931761288/L_ref.
Their respective outward and inward motions therefore persist until
the declared endpoint.

Moreover |dR/dT|=F_1 |dot_R|/sqrt(dot_R^2+F_1)<F_1<1.
For 0<=T/L_ref<=.25, consequently
0.15<R_inner/L_ref<=.4 and .8<=R_outer/L_ref<1.05.
All adjacent lapses have positive lower bounds on these closed radius
boxes (the smallest simple bound exceeds .54). The algebraic beta signs
stay positive since A>0 and B-k>0. Velocities and the smooth ODE coefficients
are bounded there. The ordinary ODE continuation theorem therefore gives
a joint solution through the full interval. The shells separate, and
the regular core excludes the singular centers of the two massive charts.
This finite continuation argument is independent of sampled residuals.

At T=.25 L_ref:

| Quantity | Outer wall | Inner wall |
| --- | ---: | ---: |
| R/L_ref | .82656089965304 | .31785194531276 |
| Proper time / L_ref | .22977119877507 | .21669213569815 |
| dR/dtau | .23281518506626 | -.93188354329553 |

The proper radial derivative can exceed one without violating timelikeness;
the coordinate/light-cone comparison is V10. Equal wall proper times
would select different middle-time slices.

## The same exterior joins a cosmological parent

Add pressureless dust outside R_b(0)=1.2 L_ref, with the *same* outer
four-form flux. No new membrane is placed at that join. The Darmois
conditions give

$$
M_0=\frac{\kappa\rho_m R_b^3}{6},\qquad
H_p^2=\frac{\kappa\rho_m+\Lambda_0}{3},\qquad
R_b=R_b(0)a_p,\qquad \rho_m=\rho_m(0)a_p^{-3}.
\tag{V11}
$$

Writing h=sqrt(Lambda_0/3) and a_p(0)=1, the exact expanding history is

$$
a_p(\tau)=\left[\cosh(3h\tau/2)
+\frac{H_p(0)}h\sinh(3h\tau/2)\right]^{2/3}.
\tag{V12}
$$

Both extrinsic forms match: beta_b=1 and
ddot_R_b+F_0'(R_b)/2=0. The common exterior clock satisfies
dT_0/dtau_b=1/F_0(R_b), and is synchronized with the outer membrane's
T_0. At the displayed endpoint tau_b/L_ref=.10872542374877,
a_p=1.07034528578242 and R_b/L_ref=1.28441434293890.
The join stays outside both membranes. Its static cosmological horizon
requires infinite T_0, so a finite smooth T_0 cannot reach it.

Initially H_p L_ref=.62660114879518 and kappa rho_m L_ref^2=.05288699901433.
These follow from the flux and compensated geometry, not from measured
cosmic density fractions. The parent accelerates because V2 supplies
vacuum stress. The selected parameters are not an explanation for the
observed dark-energy fraction.

## A curvature measurement and a cosmological cross-check

In each vacuum collar, release an infinitesimal orthogonal array of freely
falling probes momentarily at rest in the static orthonormal frame. The
relative-acceleration eigenvalues per unit separation are

$$
\mathcal A_r=-F''/2=\frac{2M}{r^3}+\frac{\Lambda}{3},\qquad
\mathcal A_t=-F'/(2r)=-\frac M{r^3}+\frac{\Lambda}{3}.
$$
$$
\boxed{\mathcal A_r+2\mathcal A_t=\Lambda},\qquad
\frac{\mathcal A_r-\mathcal A_t}{3}=\frac M{r^3}.
\tag{V13}
$$

The first combination removes the Weyl mass, while the second measures it.
In SI, relative acceleration divided by separation equals c^2 A.
This is a free-fall measurement, not the acceleration needed to hold a
static detector. Finite instruments need separation/time corrections;
the exact statement concerns the geodesic-deviation limit away from walls.

Laboratories at r/L_ref=1.1,.6,.1 remain in their respective vacuum
regions throughout this interval. Their traces are exactly 1.125,.5,.125
in L_ref^-2 units. With charge set to zero while the parent flux, tension,
gravity and initial radii are fixed, all three traces instead equal 1.125;
the shell masses change but cannot counterfeit this trace difference.
The control changes one stated coupling; it is not a claim that charged
and neutral membranes are the same action.

A comoving dust observer measures trace(A)=Lambda_0-kappa rho_m/2, not
Lambda_0 alone. In the flat dust parent the combined expansion measurement
gives a second relation:

$$
\Lambda_0=H_p^2+2\ddot a_p/a_p
=3\ddot a_p/a_p+\kappa\rho_m/2.
\tag{V14}
$$

It must agree with the adjacent vacuum laboratory. V13--V14 are standard
Einstein curvature identities now evaluated in one joined construction;
they are not uniquely diagnostic of NSC. The equal-charge three-region
relation V4 is an additional conditional test. Comparing spatially separated
laboratories still requires data transport; this finite interval does not
claim that one observer has already received all three measurements.

## Evidence, limits and next connection

The [compact record](../results/summaries/nsc-atlas-flux-transition-v1.json)
contains the inputs, outputs, source hashes and error diagnostics. Six
symbolic identities vanish exactly. Independent 40-digit first-integral
quadratures reproduce the endpoint radii and clocks within 1.0e-15 L_ref.
Two integration tolerances differ by less than 3.0e-12 in the state vector
on 65 common cuts; sampled junction, force and energy residuals are below
4.2e-13 in the normalized calculation. These are numerical checks, not
interval-certified endpoint errors or observational uncertainty.

Reproduce with
```
python scripts/run_python.py scripts/check_nsc_atlas_flux_transition.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_flux_transition.py -q
python scripts/run_python.py scripts/plot_nsc_atlas_flux_transition.py
```

The action, bare vacuum, parent flux, membrane charge and tension, initial
radii, dust state and physical ruler are explicit inputs. The calculation
does not derive those from empty premises. It establishes the conditional
source/geometry/observer connection which the former candidate lacked.
It does not select the observed Lambda, give dark-matter or baryon
fractions, establish general stability, or produce a reheating history.
The supplied two-wall Cauchy data are not claimed to arise from two
independent nucleation instantons. Classical finite evolution and
quantum formation probabilities are different problems.

The [21 October selection and conversion calculation](nsc-atlas-flux-conversion.md)
now decides what the fixed flux sector and isolated bounce formulas can
supply, and studies an explicitly separate scalar response. It retains
this completed geometry and its limits.

The next useful question is whether a specified transition/preparation
history selects flux sectors and a physical scale while leaving a
predictive observer relation. The literature's membrane formation
conditions and empty-universe/spacing restrictions must be applied to
that history. Another refinement of these two trajectories cannot
answer it.
