# Valve formation: compression, closed surfaces and a receiving region

This owner covers only roadmap step4: compression sources, circulation,
future trapping and a compatible receiving geometry. Eta is parked at
`dfda6e5a`; horizon/GH/boson refinement and downstream cosmology are frozen
except for the finite evolution needed to test formation. The existing
compact-parent action and constrained nonspherical source are retained.

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
