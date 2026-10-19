# What an interior signal determines about its preparation

The [prepared response](nsc-atlas-prepared-response.md) connects independent
exterior data to a measured finite interior waveform. The next connection is
an inverse: **if this geometry is given, what exterior data are required by a
complete central field history?** That question must precede interpreting a
waveform as evidence for the geometry itself.

On the same retained wedge, compatible signed field histories determine a
unique preparation. The inverse is explicit as a convergent integral
resolvent, and it preserves an initially quiet interval. It does not preserve
the original pulse's late cutoff, positivity or polynomial shape. Those
restrictions give a genuine test of whether a proposed history belongs to
the allowed response family.

This is a conditional application of established contraction and resolvent
mathematics. No new waveform is evolved and no second geometry is shown to
have the same observations. The [admission](prior-art-and-open-claim.md#observer-preparation-inverse-admission-12-october-2026)
limits this cycle to proof, source comparison and evidence packaging.

## Fixed model and data class

Retain H ell=1, kappa ell=2, M/ell=0.4, the classical massless conformal
spherical scalar, no extra surface interaction, zero initial cap, and the
no-return wall interval 0.45<=R/ell<=0.55. Use the exact smooth wall and
clock, not a numerical interpolant as the definition of the geometry.
The [preparation certificate](nsc-atlas-field-preparation.md) supplies

$$
D=\{(U,V):V_a\le V\le V_b,\ U_0\le U\le w(V)\},
\qquad w(V)=U_W(V),\qquad U_a=w(V_a).
$$

Both wall null coordinates increase, R/ell>=0.4 everywhere in D, and

$$
K(U,V)=16M^4R^{-4}e^{-R/(2M)}>0,\qquad f_{UV}+Kf=0.
$$

Here f=chi/epsilon, g(V)=f(U_0,V) and h(V)=f(w(V),V). All three are
signed classical reduced fields. Set L(p)=V_a for p<=U_a and
L(p)=w^{-1}(p) otherwise. The inherited forward equation is

$$
f=Eg-Af,\qquad
(Af)(U,V)=\int_{U_0}^{U}dp\int_{L(p)}^Vdq\,K(p,q)f(p,q),
$$

where E extends a function of V constantly in U. The inherited analytic
bound is ||A||<=q_*<1 in the uniform norm, with

$$
q_*=16e^{-1/2}[\sqrt{F(0.4)}-\sqrt{F(0.55)}]^2,
\qquad F(x)=(1-x/0.8)e^{x/0.8}.
$$

The recorded q_* is approximately 0.138821; this cycle reuses that bound
without a new geometric or kernel evaluation.

Fix V_s=V_a+(V_b-V_a)/4, the original preparation's onset. Define the
enlarged preparation and wall-trace class

$$
\mathcal S=\{a\in C^2([V_a,V_b]):a=0\text{ on }[V_a,V_s]\}.
$$

Members need not stop at V_e=V_a+3(V_b-V_a)/4, remain positive, or keep the
original profile. The original C3 pulse belongs to this class, but the
class is larger. The formal small physical amplitude epsilon must still
be compatible with the linear-probe approximation.

## Inverse from a wall trace

Subtract the forward equation at w(V) from its value at U. This gives

$$
f=Eh+\mathsf Bf,\qquad
(\mathsf Bf)(U,V)=\int_U^{w(V)}dp\int_{L(p)}^Vdq\,K(p,q)f(p,q).
$$

Every integration point lies in D. The same rectangular majorant used for
A bounds ||B||<=q_*<1. The borrowed Neumann/contraction argument yields

$$
f=(I-\mathsf B)^{-1}Eh,\qquad
g=S(I-\mathsf B)^{-1}Eh,\qquad
\|g\|_\infty\le\frac{\|h\|_\infty}{1-q_*},
$$

where S takes the source trace at U_0. This is a convergent analytic
representation; no inverse iteration has been run.

The following checks establish that it is the inverse of the physical
finite problem, not just an abstract integral equation. Direct differentiation
gives

$$
f_U=-\int_{L(U)}^V K(U,q)f(U,q)dq,
\qquad
f_V=h'(V)+\int_U^{w(V)}K(p,V)f(p,V)dp.
$$

The moving upper-limit term in the second identity is zero because
L(w(V))=V. Thus f_UV=-Kf, f_U=0 on the wall, and f(w(V),V)=h(V).
At V_a the integral vanishes and h(V_a)=0 supplies the zero bottom.
Integrating f_U from U_0 recovers f=Eg-Af with g=Sf. Conversely every
forward solution gives the B equation by subtraction. Uniqueness proves
that the constructions are mutual inverses. Source and wall data have
not been independently imposed in the same forward problem.

Both operators depend only on V-prefixes. Uniqueness on each prefix
implies that an initially zero h gives f=g=0 on that entire initial
strip; the analogous statement holds in the forward direction.

### Regularity at the initial characteristic and wall corner

The lower limit L has a derivative switch at U_a. Its smoothness must not
be silently assumed. Put H=Kf. Away from the switch,

$$
f_{UU}=-\int_{L(U)}^V H_U(U,q)dq+L'(U)H(U,L(U)).
$$

The extra term has zero one-sided limit at U_a because
f(U_a,V_a)=h(V_a)=0. In the chosen initial-zero class, f is identically
zero near this corner, so the term is zero in a neighborhood of the
switch. The integral of H_U is continuous across it. The first-derivative
identities and smooth K give f_UV=-Kf and

$$
f_{VV}=h''+w'K(w(V),V)h(V)
+\int_U^{w(V)}(K_Vf+Kf_V)dp.
$$

Hence f has continuous second derivatives up to the retained boundary.
The same argument with the forward identities proves C2 regularity from
g in S. Its wall restriction is C2. This proves a bijection S<->S between
source and wall traces. Smooth joining to the zero cap is compatible
because the full initial open interval and its derivatives vanish.
No theorem for arbitrary incompatible corner data is claimed.
This regularity statement concerns the reduced exterior field and the C1
central trace below. It does not assert C2 regularity of the full
four-dimensional scalar at the centre for every C2 input: recovering that
scalar from the odd reduced radial field can lose a derivative. The
original C3 pulse and the smooth exclusion example have greater regularity.

## Inverse from the central observer

Use the known strictly increasing clock theta(V)=T(V)/ell and the signed
observable Y(theta)=ell phi_C/epsilon. The inherited central response is

$$
Y(\theta(V))=\frac{2h'(V)}{\theta'(V)},\qquad
h(V)=\frac12\int_{\theta_a}^{\theta(V)}Y(s)ds.
$$

The zero initial cap fixes the integration constant. Consequently the
preparation map is a bijection from S to

$$
\mathcal Y=\{Y\in C^1([\theta_a,\theta_b]):
Y=0\text{ on }[\theta_a,\theta(V_s)]\}.
$$

For two such histories within this same fixed model, linearity gives

$$
\|g_1-g_2\|_\infty\le
\frac{\theta_b-\theta_a}{2(1-q_*)}\|Y_1-Y_2\|_\infty.
$$

This bounds reconstructed source amplitude in response to uniform errors
in complete compatible field histories. It does not bound source derivatives,
frequency content, energy, backreaction or uncertainty in the geometry.
The forward map differentiates h, so this is not a claim of uniform-norm
stability in both directions. Intensity-only detections, finitely sampled
traces, missing early history and unknown clocks are different inverse
problems. No inference beyond the recorded finite clock interval follows.

## The preparation restriction is observable

Let P be the source class allowed by a physical preparation mechanism.
A proposed complete trace is allowed by this fixed model precisely when

$$
S(I-\mathsf B)^{-1}E
\left[\frac12\int_{\theta_a}^{\theta(V)}Y(s)ds\right]\in P.
$$

For example, keeping the original source cutoff requires this reconstructed
function to vanish throughout [V_e,V_b]. A single zero-area condition does
not enforce that functional restriction. Positivity and a fixed profile
are further restrictions. The original single normalized polynomial
preparation admits only its own response (or its amplitude multiples if
that amplitude is left free).

There is also a direct exclusion for the larger class of nonzero
nonnegative sources supported in [V_s,V_e]. The positive-kernel argument
already given in the [response owner](nsc-atlas-prepared-response.md#an-independent-sign-and-magnitude-check-after-the-input)
uses positivity and support, not the particular polynomial. It therefore
gives, after the source ends,

$$
h(V)<0,\qquad h'(V)>0,\qquad
Y(\theta(V))>0,\qquad
\int_{\theta_a}^{\theta(V)}Y(s)ds=2h(V)<0.
$$

These are joint constraints on the finite signed trace, not on energy.
The strict signs follow because the continuous nonzero source is positive
on an interval and the kernel and retained transverse interval are
positive. They do not extend past V_b.

To see that the enlarged response class is genuinely larger, choose a
smooth, nonzero, nonnegative Y supported strictly between theta(V_s) and
theta(V_e). Its integrated h is nonnegative and becomes a positive
constant. Positivity of B makes f=(I-B)^(-1)Eh>=Eh and therefore g>=h>0
after V_e. This trace has a unique initially quiet preparation in S,
but that preparation cannot switch off at the prescribed cutoff. This is
an exact within-model counterexample to unrestricted fitting with a
compact source, not a calculated alternative universe or a new wave run.

## What this unlocks for the atlas

Three questions remain distinct:

1. **Fixed geometry, unknown preparation:** the complete compatible signed
   field trace has the unique preparation just derived.
2. **Different geometries, different permitted preparations:** a trace is
   degenerate only if it lies in both restricted response ranges. No second
   compatible geometry is supplied by this result.
3. **Different geometries, the same independently identified preparation:**
   comparing response operators is a stronger experiment. Identifying the
   physical input across regions is part of that experiment.

The next useful connection is an independently justified preparation
restriction and its response-range test, or an additional field channel
whose preparation is physically tied to the first. A second detector in
the unchanged cap need not help: a fixed full wall history and zero cap
already determine the interior solution. An extra observable must have
an independent dependence before it earns another calculation.

This inverse is a reconstruction by the analyst; it supplies no remote
source mechanism and creates no causal route from the excluded right
parent. Formation, quantum-state selection, descendants, connected mouths
and exterior-response alternatives retain their separate open connections.

## Borrowing, independent review and evidence

[Brunner's *Volterra Integral Equations*, Section 1.2.1](https://assets.cambridge.org/97811070/98725/excerpt/9781107098725_excerpt.pdf)
supplies the standard Neumann/resolvent method (equations 1.2.4-1.2.10,
Theorem 1.2.3 and Remark 1.2.4). Its one-variable triangular-kernel theorem
is not imported as a ready-made thin-wall theorem: the retained-domain
norm and differentiated identities above supply this application.

An independent analytic review recovered the inverse, PDE and wall
conditions, source/bottom compatibility, the corner correction and central
normalization. The forward-map cross-check is h=(I-Q)g with
Q=RA(I+A)^(-1)E, where R is wall restriction. Here ||Q||<=q_*/(1-q_*)<1,
so its resolvent gives the same inverse by uniqueness. Its sufficient
q_*<1/2 condition is stronger than the B proof's q_*<1 and is not a
necessary physical bound. The review checks the C2 statement used here;
no computation is substituted for the proof.

The [literature library](arxiv-studies.md#observer-inverse-library)
contrasts arbitrary bulk-forcing nullspaces, controlled nonlinear
source-response maps, and recorded stochastic boundary inputs with this
one passive finite trace. Those results preclude a universal claim that
geometry can never be inferred. No full proof audit of those external
inverse theorems or exhaustive priority search has been performed.

The creation-only [evidence record](../results/development/nsc-atlas-observer-inference-v1/inference.json)
binds this note to the original v7/v8 inputs and explicitly distinguishes
an analytic result from numerical reconstruction. Its
[manifest](../results/development/nsc-atlas-observer-inference-v1/manifest.json)
and [checker](../scripts/record_nsc_atlas_observer_inference.py) authenticate
the declarations and prerequisite bytes. Passing that checker authenticates
the record; the mathematical argument above supports the theorem.

## A physical radar and wall-clock inverse

The preceding theorem inverts a field preparation at fixed geometry. This
separate connection instead measures geometry using the existing quiet
cap/wall states. It joins the atlas, its physical membrane law and an
observer protocol without identifying the abstract inherited tree with an
Einstein operator. The scalar is still zero. An ideal negligible radar
transponder and a calibrated proper clock ride on the membrane. They are
test instruments, not a new reflective law for the transparent scalar.

Use c=1 and central proper time T. Inside the regular static cap,

$$
ds^2=-f\,dT^2+f^{-1}dR^2+R^2d\Omega^2,\qquad
f=1-H^2R^2>0,\qquad r_*(R)=H^{-1}\operatorname{artanh}(HR).
$$

The flat limit is r_*=R. For one reflection event on a future timelike
wall, put v=dR/dtau, beta=f*dT/dtau>0 and beta^2-v^2=f. A central radar
ray is sent at t_-=T-r_* and received at t_+=T+r_*. Its entire cap-side
path exists while this quiet wall history exists: in cap conformal
coordinates the outward null ray overtakes a timelike wall only at the
reflection; the inward ray falls strictly behind it. The original
monotone wall has its established quiet history. The inner large-radius
branch extends through its regular radial turning event; in the static
cap its two approaches to R=1/H occur at infinite static time. Thus the
finite central send and receive events used here stay inside the wall.
No path through the host bridge is assumed.

Define four readouts, all in the same central clock:

$$
D=\frac{t_+-t_-}{2},\quad Z=\frac{dt_+}{d\tau},\quad
K=\frac{dt_+}{dt_-},\quad s=\frac{d\log Z}{dt_+}.
$$

Z is the ratio of emitted proper frequency to received frequency. K is
the stretch of a sequence of central radar pulses. The transponder has
known zero delay; a drifting or uncalibrated proper standard is not part
of this inverse. The null derivatives give

$$
\frac{dt_+}{d\tau}=\frac{\beta+v}{f}=\frac1{\beta-v},\qquad
\frac{dt_-}{d\tau}=\frac{\beta-v}{f}=\frac1{\beta+v},\qquad
\boxed{f=K/Z^2}.
$$

This cancellation of the Doppler motion is more general than the cap,
and does not require a static geometry. In any spherical radial patch
with unique outward and returning null paths, send/receive times define
null coordinates u=t_- and v=t_+. Put T=(v+u)/2 and D=(v-u)/2. Then

$$
ds_2^2=-C(u,v)\,du\,dv=C(T,D)(-dT^2+dD^2),\qquad
1=C\frac{du}{d\tau}\frac{dv}{d\tau},\qquad C=K/Z^2.
$$

The normalizing central clock fixes the null-coordinate scale. This
measures C **on the observed transponder worldline**. It also gives
dD/dT=(K-1)/(K+1) and the two proper-clock derivatives directly.
The `radar_clock_metric` component accepts any C>0. A separate test on
the explicitly evolving metric C=exp(2*a*T) succeeds even where the
quiet-cap inverse rejects C>1. Thus failure of that physical model
does not invalidate the radial metric readout.

Such clocks alone do **not** determine the angular-radius function,
the metric away from that worldline or the field equations. Unique
radial paths and calibrated transponder timing are required; caustics
and multiple echoes require separately identified branches. The
following inverse is conditional on the declared vacuum-cap family.

This also identifies a readout of the existing connected conformal
candidate, without declaring it an Einstein cap. Its owned radial metric
is N(t,x)^2*(-dt^2+dx^2), N=r*Q. For fixed-coordinate reference and test
stations separated by d, reflection at t gives coordinate send/return
times t-d and t+d. Write N_s,N_r for the reference lapse at those
events and N_b for the beacon lapse at reflection. Direct proper-clock
conversion gives Z=N_r/N_b and K=N_r/N_s, hence
C_radar=N_b^2/(N_s*N_r). Thus the clock readout measures the evolving
metric variable already used by the [finite causal observer calculation](schur-memory.md#causal-radial-distances-and-the-observers-clock),
with both clock Jacobians retained. It is not generally equal to N_b^2,
and this coordinate transformation adds no vacuum source or new state.

For 0<f<=1 and D>0, put x=sqrt(1-f). Then

$$
H=\frac{\operatorname{artanh}x}{D},\qquad
R=D\frac{x}{\operatorname{artanh}x},\qquad
\beta=\frac{Zf+Z^{-1}}2,\qquad v=\frac{Zf-Z^{-1}}2.
$$

At x=0 the limits are H=0 and R=D. H^2 is well-conditioned at zero
curvature; H itself has square-root sensitivity there. Expanding motion
is K>1; contraction uses the same formulas. The result reconstructs the
supplied cap Lambda=3H^2, rather than deriving its physical origin.

### The signed junction becomes a clock relation

For a zero-mass cap and the existing constant-tension membrane law,

$$
\beta(R)=a/R^2+cR,\qquad a=M_+/k,\qquad
c=\frac{k^2+H_+^2-H^2}{2k},\qquad k=4\pi\sigma.
$$

The acceleration is vdot=beta*beta_R+H^2R. Differentiating
Z=1/(beta-v) with wall proper time, then dividing by dt_+/dtau=Z,
proves s=beta_R/Z+H^2R. Consequently the same readouts give

$$
\beta_R=Z(s-H^2R),\qquad
c=\frac{\beta_R+2\beta/R}{3},\qquad
a=\frac{R^2(\beta-R\beta_R)}3.
$$

These are two physical junction combinations. With independently specified
host H_+, retain every positive solution of
k^2-2ck+(H_+^2-H^2)=0 and set M_+=ka. Schwarzschild host H_+=0
has one positive root k=c+sqrt(c^2+H^2). Equal vacua give k=2c.
If the host curvature is unknown there is a continuous ambiguity. If it
is larger than the cap curvature, two positive tensions can produce the
same a,c and **all** cap-side clock histories. Later cap measurements
cannot remove that branch degeneracy.

Conversely, if k is fixed independently by the common membrane law, the
same interior clock data infer the host itself:

$$
H_+^2=H^2+2kc-k^2,\qquad M_+=ka.
$$

At the common k*ell=2, the two examples below return H_host^2*ell^2=0
and 1. This is a conditional inference of surrounding geometry through
its effect on the material boundary. It requires neither a new incoming
waveform nor a signal through the host bridge. Unknown k or a different
wall material restores the stated ambiguity; the finite clock prediction
tests continuation of the supplied material law.

### Finite predictions in the established states

The [calculator](../scripts/check_nsc_atlas_radar_inverse.py) uses the
original one-level quiet wall at R/ell=.5, M/ell=.4, H_host=0 and the
nested inner wall at R/ell=.8, M/ell=.01, H_host*ell=1. Both have
H*ell=1 and k*ell=2. These are two established states under the same
membrane law, not two observations linked through a constructed parent
history. Inversion recovers those inputs to roundoff.

| Central readout | One-level wall | Nested inner wall |
| --- | ---: | ---: |
| D/ell | .549306144334055 | 1.09861228866811 |
| Z | 2.62549201161791 | 3.74639337610612 |
| K | 5.16990622730208 | 5.05276679827145 |
| s*ell | -.433158428652097 | 1.06171003724629 |
| Predicted Z after .04 ell of receiver time | 2.58466102249367 | 3.90852634473171 |

Infer H,a,c at the first event, then integrate only
dt_+/dR=Z/sqrt((a/R^2+cR)^2-1+H^2R^2). This predicts the later Z
without further calibration. An independent second-order proper-time
integration of the original signed junction and its null clock agrees
within 3.56e-15 in Z; 32/64 Gauss values agree at displayed precision.
An independent 80-digit nested calculation also agrees within 3e-14.
Those are numerical controls, not experimental or total model errors.

The finite consequence is not an inversion identity. Replacing beta(R)
by beta(R)+.5*(R-R0)^2/ell^2 preserves beta, its first derivative and
all first-event readouts, but changes the later Z by .0002413953 and
.0001163056 respectively. Fixed bulks and the Israel jumps give this
alternative a varying material law and positive surface density on the
short interval. Thus the later clock reading tests the proposed constant
membrane law. It does not discriminate candidates sharing the same a,c.

At the declared ell=10000m, .04 ell/c is **1.33425638079 microseconds**.
The inner drift is 31829.2661749 per second; its one-way radar time is
36.6457614043 microseconds. Neither scale is an inferred astronomical size.

### Explicit input uncertainty and the remaining connection

For input bounds [D0,D1], [Z0,Z1], [K0,K1], let
f0=K0/Z1^2 and f1=K1/Z0^2, with 0<f0<=f1<=1.
The function y(f)=artanh(sqrt(1-f)) decreases, while
r(f)=sqrt(1-f)/y(f) increases: artanh(x)/x is increasing for 0<=x<1.
Therefore H^2 lies in [y(f1)^2/D1^2,y(f0)^2/D0^2], and R lies
in [D0*r(f0),D1*r(f1)]. Correlated errors can narrow this enclosure.
The code evaluates the analytic endpoint bounds in float64; it is not
a directed interval proof of transcendental rounding.

For an explicitly **assumed** independent 100ppm bound on D,Z,K,
the nested example gives Lambda_cap*ell^2 in
[2.99837630274001,3.00162434110127] and R/ell in
[.799850976738802,.800149048509163]. These are synthetic instrument
bounds, not observed uncertainties or bounds on model/preparation error.
The [compact evidence](../results/summaries/nsc-atlas-radar-inverse-v1.json)
contains 1480 scalar evaluations, no field grids, source hashes and the
reproduction command. The physical source/origin of H, the common parent,
and a sky light cone remain unresolved. This proof supplies a conditional
observer and shared-law test in the actual atlas rather than another
vacuum coefficient selected by averaging convention.

The likely failures are an uncalibrated/delayed beacon (which corrupts
F), leaving HR<1 or selecting one of two positive tension roots (which
breaks identifiability), and a changing metric or material law (which
fails the finite-time clock prediction). Each is visible in its own
readout/domain/consistency check; none requires reopening P001.
