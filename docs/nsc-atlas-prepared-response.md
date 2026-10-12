# The interior response to a finite exterior preparation

The [preparation certificate](nsc-atlas-field-preparation.md) supplies independent
null data in the retained Schwarzschild white-hole sector IV, a finite
exterior-to-wall operator, and its connection to the regular interior observer.
This successor evaluates that operator. Its question is whether the prescribed
compact exterior data produce a nonzero interior signal after the input ends,
and whether that signal can be distinguished from numerical error.

The completed evaluation finds a bipolar central waveform followed by a
positive field tail within the admitted window. The tail agrees with the
independent analytic bracket below and is resolved under the prescribed
spatial and geometric refinements. A first unsegmented attempt stopped at
its call limit; the event-aligned run passed without changing the theory,
input or accuracy thresholds.

The physical source is the declared finite IV preparation. The ordinary right
parent remains causally excluded in this isolated reference. This calculation
does not supply a global source, formation history, quantum state, energy flux,
permanent detector memory, or a unique reconstruction of the exterior.

## Fixed input, domain and borrowed methods

All physical inputs remain those of the preparation certificate: H ell=1,
kappa ell=2, M/ell=0.4, a real classical massless scalar with xi=1/6 and angular
mode zero, linear probe evolution, and no extra surface scalar interaction.
The wall expands through 0.45<=R/ell<=0.55. The retained exterior strip has
0.4<=R/ell<=0.55. Zero cap data are supplied at t_a=t_W(0.45); the cap has
no centre-to-wall return during the whole window. Consistency with an earlier
zero cap additionally requires no earlier incoming drive, as the owner explains.

Write x=R/ell, m=0.4, U_a=U_W(V_a), U_b=U_W(V_b), and

$$
F(x)=(1-x/(2m))e^{x/(2m)},\qquad UV=F(x),
\qquad U_0=F(0.4)/V_a.
$$

The domain is D={V_a<=V<=V_b, U_0<=U<=U_W(V)}. Both wall null coordinates
increase. The independent input is

$$
g(V)=\epsilon B\left(\frac{4(V-V_{\rm mid})}{V_b-V_a}\right),
\qquad V_{\rm mid}=(V_a+V_b)/2.
$$

Here B(z)=(1-z^2)^4 for |z|<1 and zero otherwise. It is C3, not infinitely
smooth. Epsilon is a positive formal small amplitude. Set DeltaV=V_b-V_a,
V_s=V_a+DeltaV/4 and V_e=V_a+3DeltaV/4. Thus the input and its derivative
vanish for V>=V_e. This fixes the source independently of any desired response.
The bottom characteristic data on V=V_a are zero.

The scalar equation and matching are borrowed unchanged. Numerical evaluation
uses [SciPy's DOP853 implementation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.DOP853.html)
of the established explicit Runge-Kutta method and a clipped composite
trapezoidal spatial integral. Geometry uses
[cubic Hermite interpolation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicHermiteSpline.html)
for U_W(V) and T(V), using their known derivatives, with monotonicity checked;
[PCHIP interpolation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html)
supplies the positive Vdot coefficient.
The reproducible evidence records the implementation and its environment.
These methods are integration tools, not new field dynamics or a general
PDE framework.

## An exact extension for a finite numerical representation

The exterior equation and inherited wall condition are

$$
\chi_{UV}+K\chi=0,\qquad
K=\frac{16m^4}{x^4}e^{-x/(2m)}>0,\qquad
\chi_U|_W=0,\qquad \chi=R\phi.
$$

Changing the integration order in the admitted wedge operator gives

$$
(Af)(U,V)=\int_{V_a}^{V}dq
\int_{U_0}^{\min(U,U_W(q))}dp\,K(p,q)f(p,q),
\qquad \chi=g-A\chi.
$$

For numerical storage only, extend chi constantly in U from U_W(V) to U_b.
If U>=U_W(V), then U>=U_W(q) for every q<=V. Both the integral above and
g(V) are consequently independent of U there. Thus the integral equation
itself defines this extension. It agrees with the retained solution and its
zero U derivative at the wall; it does not apply the vacuum equation on the
excised side or identify that side with another exterior region.

Differentiation at fixed U gives

$$
\chi_V(U,V)=g'(V)-\int_{U_0}^{\min(U,U_W(V))}
K(p,V)\chi(p,V)dp.
$$

At U=U_b the stored value equals the moving wall value. The total wall
derivative has no additional term because chi_U|_W=0. Therefore

$$
J(V)=\frac{d\chi_W}{dV}=\chi_V|_W
=g'(V)-\int_{U_0}^{U_W(V)}K(p,V)\chi(p,V)dp.
$$

Numerically evolve z=chi-g, initially zero, so that

$$
z_V(U,V)=-\int_{U_0}^{\min(U,U_W(V))}
K(p,V)[z(p,V)+g(V)]dp.
$$

This removes cancellation between a large prescribed source derivative and
its numerical time integration. Evaluate J directly from this right-hand
side plus g'; do not obtain a small tail by differencing stored field values.
The last spatial cell ends at the actual moving wall. Whole cells beyond it
contribute nothing. Its endpoint value uses the constant extension, rather
than setting the field to zero or propagating vacuum dynamics across the cap.

## An independent sign and magnitude check after the input

The existing analytic norm bound is

$$
\|A\|_\infty\le q_*
=16e^{-1/2}[\sqrt{F(0.4)}-\sqrt{F(0.55)}]^2
=0.13882099234357653<1.
$$

The nonnegative source and positive kernel give a stronger pointwise
comparison than the norm statement alone. Reorder the positive integrations
in A^2g. For every fixed inner point (p,q), the remaining integration is
over its future inside the original wedge, and its total kernel integral
is bounded by q_*. Consequently

$$
0\le A^2g\le q_*Ag.
$$

Positivity of A then yields A^(n+1)g<=q_* A^n g for n>=1. The uniformly
convergent Neumann series alternates with pointwise decreasing terms from
Ag onward. Its first two bounds are

$$
g-Ag\le\chi\le g-Ag+A^2g.
$$

After the source ends, g=g'=0. It follows that

$$
-Ag\le\chi\le-(1-q_*)Ag<0
\quad\hbox{for }U>U_0,\ V_e\le V\le V_b.
$$

This negative reduced exterior field produces a positive derivative trace
J=-integral(K chi). The trace's leading post-input contribution is second
order in the potential: one K creates -Ag and another converts it into J.
A first-order transmitted-pulse calculation would miss it.

Define the positive reference integral

$$
\mathcal B(V)=\int_{V_s}^{V_e}g(q)dq
\int_{U_0}^{U_W(q)}dp\,K(p,q)
\int_p^{U_W(V)}du\,K(u,V).
$$

For V>=V_e, monotonicity ensures U_W(q)<=U_W(V), so every factor and
integration length has the stated sign. Reordering integral(K Ag) gives
exactly this expression, without a factor of two. The exact derivative
trace therefore obeys

$$
(1-q_*)\mathcal B(V)\le J(V)\le\mathcal B(V).
$$

This bracket requires neither a computed wave solution nor numerical
differentiation of it. At fixed negative V, implicit differentiation of
UV=F(x) gives

$$
K\,dU=-\frac{4m^2}{Vx^3}dx,\qquad
\int_a^b K(U,V)dU=\frac{2m^2}{V}
[x(b,V)^{-2}-x(a,V)^{-2}].
$$

The latter is positive: x increases with U, while V is negative. Using it
for the innermost integral leaves an independent positive double quadrature
for the reference bracket.

There is also a conservative nonzero lower bound using only endpoints.
Put K_min=K(0.55), and

$$
G=\int g(V)dV=\frac{64\epsilon\Delta V}{315}.
$$

For U<=U_a, the entire source history retains [U_0,U], giving
Ag(U,V)>=K_min(U-U_0)G. Restricting the current J integral to
[U_0,U_a] proves

$$
J(V)\ge\frac{(1-q_*)K_{\min}^2(U_a-U_0)^2G}{2}>0,
\qquad V_e\le V\le V_b.
$$

The recorded geometric endpoints give approximately 6.53246e-5 epsilon.
The formula establishes the positive bound; the quoted decimal inherits
the geometric quadrature diagnostics and is not an interval-certified
decimal enclosure. It supplies a useful absolute scale below which a
numerical tail error must fall.

## Interior observer and unit scaling

The inherited conformal matching convention is I=2 Vdot J. The regular
centre and the no-return condition give

$$
\phi_C[T(\tau)]=\frac{I(\tau)}{T'(\tau)},\qquad
T=t_{\rm flat}-\ell\log(1-R/\ell).
$$

Both Vdot and T' are positive. Thus the positive post-input J gives a
positive central field over the corresponding finite arrival interval.
Report the dimensionless observable

$$
Y(T/\ell)=\frac{\ell\phi_C}{\epsilon}
=\frac{2\ell\dot V}{T'}\frac{J}{\epsilon}.
$$

Changing ell rescales proper times and phi at fixed reduced amplitude;
the null geometry, z/epsilon, J/epsilon, and Y as a function of T/ell are
unchanged. Epsilon scales the linear field, not the background geometry.
The admitted clock window remains
-0.581758951595344<=T/ell<=-0.250230700560700. No conclusion extends it to
T=0, to infinite time, or to a permanent remnant. The old compact
wall-input calibration and these compact exterior-input data are different
preparations, so their differing support is not a contradiction.

## Bounded numerical evaluation and acceptance

The separate [admission](prior-art-and-open-claim.md#prepared-scalar-response-admission-12-october-2026)
authorizes the following finite evaluation. The theory, domain and source
are fixed before execution. The complete recorded experiment has at most
60 CPU seconds and 10 million kernel evaluations. Each evolution permits
at most 6,000 right-hand-side evaluations. Direct geometric evaluation
uses at most 257 knots per grid and at most 300 points per adapter call.
Use spatial grids of 129, 257 and 513 points and 257 common output points.

Compare 129 and 257 direct geometric knots, retain 257 for the reported
geometry, and disable extrapolation. Geometry interpolation, spatial
integration and temporal tolerance are separate error sources. DOP853's
local error control is not a total-error certificate; these piecewise cubic
interpolants are C1, the source is C3, and clipped cell crossings also limit
smoothness. No eighth
order convergence claim follows merely from choosing that integrator.

The final spatial refinement must meet
|Y_513-Y_257|<=1e-6+1e-4|Y_513| at every output point; geometry refinement
must meet the same mixed threshold. Temporal tolerances must be tighter
than that reported threshold. Report the coarser differences as evidence
of convergence, not just the final pass flag. After the input, additionally
require derivative-trace refinement errors below ten per cent of the
elementary positive lower bound above, rather than accepting a tail solely
because of relative agreement near zero.

Independently evaluate the positive reference integral with Gauss orders
32 and 64 at V=V_e+f(V_b-V_e), for f=1/4, 1/2 and 1. Their difference must
also be below ten per cent of that lower bound. The numerical J must lie
between (1-q_*) times the reference integral and the reference integral,
allowing the separately recorded refinement discrepancies. This is a
bounded numerical check supported by an analytic sign proof, not a
rigorous interval enclosure of every output value.

Two borrowed or elementary controls check different parts of the method:

- Setting K=0 must give chi=g, z=0, J=g', and no post-input tail. This is a
  numerical control, not an alternative Schwarzschild theory.
- On a fixed rectangle with constant K, zero bottom data and g(V)=V-V_a,
  put t=V-V_a and w=K(U-U_0)t. The exact solution is
  chi=t J_1(2 sqrt(w))/sqrt(w), with its continuous value t at w=0;
  its V derivative is J_0(2 sqrt(w)). The standard
  [Bessel power series](https://dlmf.nist.gov/10.2#E2) verifies both the
  integral equation and normalization independently. This manufactured
  rectangle tests the numerical integration, not the moving-wall geometry.

Also check amplitude linearity, unit rescaling, retained-side clipping,
source/bottom compatibility, no-return, and rejection of clock-domain exits.
The owner of the numerical implementation is
[nsc_atlas_prepared_response.py](../src/recursive_horizons/nsc_atlas_prepared_response.py).
The record must distinguish analytic bounds, numerical measurements and
unresolved extensions. No evaluated waveform is claimed by this admission.

## First execution stop and algorithm clarification

The first admitted prototype stopped at the 6,000-right-hand-side-call limit
on its 129-point spatial grid, near normalized source coordinate
s=(V-V_a)/(V_b-V_a)=0.575. It did not reach the end of the window and did not
produce an accepted waveform. The reported attempt took about 6.6 seconds;
the call cap, rather than a changed physical result, stopped it.

The moving wall crosses fixed spatial nodes. At each crossing the discrete
clipped integral changes which cell is partial. Its coefficients are
continuous but have reduced derivative smoothness, causing the adaptive
ODE driver to spend steps locating known events. This is a numerical
representation issue; it does not justify adjusting the field, background,
source, requested tolerances, or scientific domain.

Before the next attempt, split the normalized V interval at every actual
crossing U_W(V)=U_i, and at the two source support edges. Obtain crossing
coordinates by a bracketed scalar root solve on the existing monotone
Hermite wall interpolant. Carry the same correction state continuously
between segments. Apply DOP853 separately on each segment, starting with
its segment length as the initial step and retaining the same tolerances.
The method still integrates exactly the same spatially discretized wedge
operator; no extra physical boundary or interpolation across the cap is
introduced. Source queries and all recorded output remain inside the
original finite window.

The redirected run retains all original per-record caps, including 6,000
total right-hand-side calls **summed across segments for each solve**, not
6,000 per segment. Its total kernel and CPU limits remain 10 million and
60 seconds. Record the unsuccessful attempt separately. Another exhausted
or unresolved attempt stops this branch of numerical evaluation without
weakening the acceptance thresholds.

## Measured finite response and validation

The event-aligned run completed the full declared window. The canonical
[prepared-response record](../results/development/nsc-cosmological-atlas-v8/prepared-response.json)
stores all 257 output samples, analytic checks, refinement diagnostics and
budget counts. The result remains conditional on the stated classical
preparation and fixed background.

| Quantity | Recorded value |
| --- | ---: |
| Positive sampled extremum of Y=ell phi_C/epsilon | 47.29833694 |
| Negative sampled extremum of Y | -45.69595941 |
| Y at the final admitted arrival time | 0.00129049806 |
| Minimum sampled J/epsilon strictly after the input | 0.000716524633 |
| Analytic uniform lower bound for post-input J/epsilon | 0.0000653246206 |
| Maximum sampled change in Y from the zero-kernel control | 0.26595212110 |

These extrema are sampled field amplitudes, not energy fractions or bounds
on unsampled continuum extrema. The zero-kernel comparison omits the
potential as a numerical diagnostic; it is not a second physical spacetime.
The terminal value occurs at T/ell=-0.250230700560700. No terminal zero,
later continuation, permanent memory or asymptotic decay law is inferred.
The end of the exterior input maps to T/ell=-0.340772594886934; the
post-input claim concerns its finite arrival interval through that terminal
event. A field tail here does not assert a particular long-time decay law.

The independent positive-integral comparisons are

| Source coordinate s | Analytic lower value (1-q_*) B/epsilon | Measured J/epsilon | Analytic upper value B/epsilon |
| --- | ---: | ---: | ---: |
| 0.8125 | 0.0006371321074 | 0.0007361136575 | 0.0007398370162 |
| 0.875 | 0.0006563785613 | 0.0007575130409 | 0.0007621859747 |
| 1 | 0.0006965878198 | 0.0008020390045 | 0.0008088769160 |

The largest Gauss 32/64 difference in B is 7.18e-17. The mathematical
inequality is exact for the admitted continuum operator; these evaluations
of its reference integral and the field are numerical. The agreement is
independent evidence for the predicted sign and size, not a rigorous
interval enclosure of the numerical waveform.

| Refinement or consistency diagnostic | Maximum difference |
| --- | ---: |
| Y, spatial grid 129 to 257 | 1.11704166e-6 |
| Y, spatial grid 257 to 513 | 2.79286258e-7 |
| Post-input J/epsilon, spatial grid 257 to 513 | 3.70648451e-9 |
| Y, geometry knots 129 to 257 | 3.40622250e-8 |
| Post-input J/epsilon, geometry knots 129 to 257 | 1.84677577e-14 |
| Signed observer area versus twice the wall endpoint difference | 2.50474214e-10 |

The maximum mixed-threshold ratios are 0.01005 for the spatial comparison
and 8.91e-6 for the geometry comparison, both below the admitted limit one.
The post-input spatial discrepancy is also well below ten per cent of the
analytic lower bound. The roughly fourfold spatial reduction supports the
expected trapezoidal refinement behavior in this case; it does not prove
a general convergence order. The observer-area diagnostic uses
integral(phi_C dT)=2 Delta chi_W over this finite window. It checks the
clock and normalization independently of peak sampling and does not
measure energy.

The four evolution call counts were 1,465, 2,192, 4,352 and 4,352 for the
three field grids and the separate geometry refinement. The complete
record used 5,597,264 kernel values, below 10 million, and passed the
60-CPU-second guard. The 6,000-call cap was applied to each whole evolution,
including every event segment. The failed prototype is retained as a
separate stopped attempt, rather than included as accepted evidence.

The [15 focused tests](../tests/test_nsc_atlas_prepared_response.py) passed.
They include the exact transverse primitive checked by independent
quadrature, the zero-potential control, a constant-potential Bessel
benchmark for both chi and its independently known J_0 derivative trace
in the same solve, fixed-theory amplitude linearity, regional unit rescaling,
positive-tail brackets, budget enforcement and interpolation/clock-domain
rejection. The inherited preparation tests retain corner compatibility,
full retained-side membership and no-return checks. The completed result
closes the finite preparation-to-wall-to-observer calculation for these
data; it does not remove the preparation and inference questions below.

## Stop conditions and remaining connections

A calculation stops if its budget expires, refinement cannot separate the
tail from error, the independent positive bracket disagrees, or the input
or geometry fails an inherited compatibility check. It reports that failure
without adjusting the coupling, source profile, domain or background to
preserve a desired result. A successful finite response closes the
preparation-to-wall-to-observer connection for this candidate only.

The likely obstacles are unresolved spatial or interpolation error in a
small tail; accidental integration across the excised cap; and an apparent
late response inferred outside the certified clock or no-return window.
The corresponding checks are the absolute tail bound, exact moving-boundary
clipping, and explicit domain rejection. A verified signal still leaves its
remote preparation and observational distinguishability open. Other atlas
geometries and communicating-mouth routes remain independent candidates.
