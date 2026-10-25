# Formation needs an inherited growth window

## Decision and admission

The [spatial correlation](nsc-atlas-compensation-response.md) is complete
within its declared preparation. Formation cannot be inferred by evolving
the same already dense core again. The older omega2000 diffuse-cloud
experiment also began with a nonlinear overdensity. Compare three routes:
another nonlinear cloud, a new boundary-driven preparation, and the growth
available to a small ambient disturbance in the same omega1e9 parent.
The third tests the missing prerequisite without choosing a final core.

Complete: the all-wavelength finite-future bound and compact earlier
growth window are established below. Their
[nonlinear continuation](nsc-atlas-seed-collapse.md) now tests two explicit
finite-amplitude preparations, preserving the linear calculation's domain.

The only new dynamical state is the linear density/phase perturbation.
The parent action, vacuum, excitation and curvature remain unchanged.
An earlier initial scale is explicitly a comparison on the averaged
background continuation, not a constructed earlier shared-wall spacetime.
Use the existing envelope, direct density/phase equations, an independent
dust quadrature and conserved transfer determinant. At most120 CPU seconds
for a small set of physically named modes and epoch controls, four workers,
no PDE grid, source fitting, empirical amplitude or primordial spectrum.

The likely failures are a nonregular infinite-time endpoint, treating the
pressureless comparison as an unproved finite-k bound, and confusing an
earlier homogeneous continuation with the completed atlas past. Reciprocal
scale, a Volterra comparison proof and explicit epoch/state scope address
these separately. A linear threshold identifies when nonlinear equations
become necessary; it is not a collapse theorem.

The bounded calculation also includes a compact compensated seed, built
from a smooth bump's autocorrelation rather than the stationary core.
Five spectral/tolerance/epoch controls evaluate it on the same background.
No nonlinear evolution is started before the growth window is known.

## Density, phase and the actual inherited background

Keep omega=1e9, Lambda=9/32, D=kappa Eb=3/8 and C=.0227765023676067.
Use the already derived averaged parent,

\[
H^2=\frac{\Lambda}{3}+\frac{D}{3A^3}-\frac{C}{A^2}.
\]

In comoving atlas distance y=R/A the same total envelope can be written
psi=b(1+delta/2+i theta) to first order. Its continuity and phase equations
for a local comoving wave number k are

\[
\dot\delta=\frac{k^2\theta}{\omega A^2},\qquad
\dot\theta=-\frac{k^2\delta}{4\omega A^2}
+\frac{\omega D\delta}{2Ak^2}.
\]

Eliminating phase gives

\[
\ddot\delta+2H\dot\delta+
\left(\frac{k^4}{4\omega^2A^4}-\frac{D}{2A^3}\right)\delta=0.
\]

The derivative term follows from the existing classical real-scalar
envelope, without a new quantum interpretation. The related averaged
massive-scalar equation and its Jeans scale are attributed to
[Hwang--Noh,0902.4738, equations (19)--(28)](https://arxiv.org/html/0902.4738).
Their flat-background derivation is not substituted for our curved
inherited expansion; this calculation linearizes our already specified
local envelope directly. Physical use needs the original averaging,
weak-field and local/subhorizon domain. No global flat parent is inferred.

An initially real envelope has theta=0 and zero peculiar current. A
growing-mode velocity is additional state information, not an automatic
consequence of specifying density. We retain the independent input
v_i=delta_dot_i/H_i in the transfer matrix.

## Infinite expansion has a regular endpoint in reciprocal scale

Set u=1/A, nu=k^2/(2 omega) and p=A^2 delta_dot. Then

\[
\frac{d\delta}{du}=-\frac{u p}{H(u)},\qquad
\frac{dp}{du}=\frac{\nu^2u-D/2}{H(u)}\delta.
\]

For this background H never vanishes and tends to sqrt(Lambda/3).
The entire averaged expanding future is therefore the finite regular
interval u_i to zero. No stopping age or guessed future growth tail enters
the result. The matrix in (delta,p) has zero trace, so its fundamental
determinant is exactly one. With columns (delta_i,v_i), the determinant
instead equals A_i^2 H_i; the code checks the correct normalization.

This is an infinite-future property of the stated reduced equations.
It does not establish an infinite-time error bound relative to full GR,
nor the unresolved complete future of both material walls.

## A proved comparison bound for every wavelength

Introduce kinetic time d eta=dt/A^2, starting at zero. The equation becomes

\[
\delta_{\eta\eta}+\nu^2\delta=\frac{DA}{2}\delta.
\]

Its Volterra integral equation has free part
delta_i cos(nu eta)+p_i sin(nu eta)/nu and source kernel
sin[nu(eta-s)]/nu. Since the absolute kernel is at most eta-s and the
source coefficient DA/2 is nonnegative, iteration of the positive Volterra
kernel bounds every finite-k solution by the zero-frequency solutions:

\[
|\delta_k(\eta)|\leq
D_0(\eta)|\delta_i|+E_0(\eta)|p_i|,
\]

where D0 has initial data (1,0) and E0 has (0,1) in (delta,p).
Both are nonnegative and increasing. The integrable endpoint singularity
DA/2 in eta causes no problem: the u system is regular, and
integral (DA/2) d eta = (D/2) integral du/H is finite.
This proves the comparison; it does not assume that a positive pressure
term always lowers an arbitrary signed oscillatory solution.

An elementary bound needs no perturbation integration. For H>=Hmin>0,
the integral of the positive endpoint kernel is at most
D u_i^3/(12 Hmin^2). Gronwall therefore gives, for p_i=0,

\[
|\delta_k|\leq
\exp\left(\frac{D u_i^3}{12H_{\min}^2}\right)|\delta_i|.
\]

For u_i=1, the minimum is
Hmin^2=Lambda/3-4C^3/(3D^2)=.0936379692352. The elementary factor is
1.39616911774, safely below1.397. The sharper positive comparison is

\[
|\delta_k(t)|\leq1.26934144123|\delta_i|
+0.68728524662|v_i|\qquad(A_i=1).
\]

The defining comparison coefficients are exact solutions; the printed
decimals are numerical evaluations, not interval-certified constants.
Independent quadrature and the ODE agree within8e-15; tighter integration
changes the two coefficients by less than3.4e-12 using twice the difference.
The elementary bound remains available without those evaluated coefficients.

The independent pressureless basis is H(A) and
H(A) integral_0^A da/[a^3 H(a)^3]. This established result is attributed to
[Heath (1977)](https://doi.org/10.1093/mnras/179.3.351), whose abstract and
bibliographic record were read. We verify H's equation symbolically and
derive the second solution by reduction of order; the numerical integral
does not call the perturbation solver.

Parseval gives the same bound on the spatial L2 norm for zero initial
current. A pointwise bound requires a specified Fourier norm; it does not
follow from an arbitrary pointwise density bound alone. The compact seed
below has a nonnegative Fourier weight and admits exactly that stronger
pointwise conclusion.

Consequently a small zero-current perturbation introduced at A=1 cannot
reach order-one contrast within this linear evolution. The necessary
linear-extrapolation amplitude is at least .78781, already outside the
small-perturbation regime. This rules out that proposed linear route,
not all nonlinear configurations, externally driven seeds or earlier epochs.

## The inherited supply and core scales have different growth histories

The modes below use the previous geometry, not newly fitted wavelengths.
The core wave number is one divided by its initial half-mass radius;
it is not a claim that an entire core equals one Fourier mode.

| Mode starting at A=1 | Final zero-current density gain | Jeans crossing A |
| --- | ---: | ---: |
| Positive pressureless comparison | 1.26934144 | comparison only |
| k=1/Rd | 1.26933642 | .0000247930 |
| k=khalf from the spatial turnover | 1.26798790 | .00668511 |
| Initial Jeans wave number | 1.07277044 | 1 |
| k=1/Rhalf,core | .44810375 | 52670.8 |

The crossing is A_J=k^4/(2 omega^2 D). It only marks a sign change in the
linear restoring coefficient. It is neither sufficient growth nor a
collapse time; even the late crossing of the small core scale does not
evade the finite-future bound. Nonlinear concentration must collect from
larger scales or arrive with already nonlinear data.

## A compact seed without an installed core

Choose, in comoving local distance,

\[
g(r)=(1-r^2/L^2)^4_+,\qquad
\delta_i(r)=\epsilon\frac{-\Delta(g*g)(r)}{\int|\nabla g|^2d^3x}.
\]

The convolution has support2L. Integration by parts gives delta_i(0)=epsilon;
its integral vanishes. Its Fourier transform is proportional to
k^2 |g_hat(k)|^2, hence nonnegative. Therefore |delta_i(r)|<=epsilon, and
the proved all-mode bound yields a uniform pointwise future bound for this
entire compensated disturbance. Its linear Poisson potential is compact
as well. No stationary-core solver is called to select this seed shape.

Take L equal to the previous comoving reservoir radius .000481562379847,
epsilon=.05 and zero phase. These are explicit benchmark preparation
inputs, not a derived primordial state or empirical fit. At A_i=1 its
central future gain is1.26488829 and the general pointwise bound remains
1.26934144. A 5% seed stays below about6.35% in the linear model.

Now use A_i=.01 on the same averaged homogeneous branch. This changes the
preparation epoch, not the action or vacuum. For q=kL the normalized
transform of g is10395 j5(q)/q^5. The radial spectral weight for delta is
q^4 times its square. Parseval fixes its integral to

\[
I=\frac{64\pi B(5/2,7)}{B(3/2,5)^2}=295.53715505942.
\]

No truncated spectrum is renormalized. The explicit sine/cosine polynomial
for j5 bounds the omitted weight beyond Q by

\[
\frac{10395^2}{I}\int_Q^\infty q^{-8}
\left(1+\frac{15}{q}+\frac{105}{q^2}+\frac{420}{q^3}
+\frac{945}{q^4}+\frac{945}{q^5}\right)^2dq.
\]

Multiplying by the positive dust gain bounds the omitted **pointwise**
response at every future time. At Q=60 this bound is3.06712e-6 for the
unit-normalized seed. Quadrature and ODE controls add a measured gain
comparison of1.02e-8. The compact seed's final central linear gain is
35.60903223; a conservative presentation is35.609032 +/- .000004, where
the spectral part is bounded and the integration part is a comparison
estimate, not an interval enclosure of the whole calculation.

The declared 5% seed reaches delta=.1 near A=.0957400. Its **formal linear
extrapolation** reaches delta=1 near A=1.0190483 and tends to1.78045.
The unit-contrast crossing is outside linear validity: it identifies a
nonlinear continuation to perform, not a formed core or a completed
collapse proof. The seed at the earlier epoch is still supplied, and its
shape and age have not been selected by an inherited probability law.

The largest retained k/(omega A_i) is .01246 and H_i/omega=3.53e-7.
These identify the reduction hierarchy; they are not a total GR error bar.
The earlier expanding background is a conditional continuation, not the
completed past of the existing shared material walls or a thermal history.
A local use also needs a parent domain large enough for the intended
causal calculation. Radiation, other matter and a different early state
could change this window and would define a different comparison.

## Verification and next calculation

The [compact evidence](../results/summaries/nsc-atlas-formation-window-v1.json)
contains twelve named mode/epoch/tolerance cases and five compact-packet
controls, computed in four workers in about4 CPU seconds. Symbolic
elimination, the independent dust basis, direct phase transport, conserved
transfer determinant, an exact no-gravity oscillator, real-space Parseval
normalization and the analytic spectral tail check the result. The exact
homogeneous scalar remains homogeneous by uniqueness; no numerical noise
is advertised as an inherited seed mechanism.

A native Grok review independently recovered the density equation and the
late zero-current gain1.2693414412, while distinguishing an independent
initial velocity from density data. The all-wavelength comparison and
compact-packet error bound above are established by their stated proofs
and executable checks, not by the reviewer's agreement.

The next admitted nonlinear calculation should start from this explicit
small compact density/current preparation on a declared earlier parent,
solve its compatible initial source/constraint dictionary and follow the
emerging concentration. It should not install the old stationary core
or choose a new spectrum to match its mass. Formation, retained mass and
the spatial turnover can then be tested as outputs of the same history.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_formation_window.py --output results/summaries/nsc-atlas-formation-window-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_formation_window.py
```
