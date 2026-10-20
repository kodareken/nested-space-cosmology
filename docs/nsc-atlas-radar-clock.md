# A physical radar and wall-clock inverse

The [field-inverse theorem](nsc-atlas-observer-inference.md) inverts a field preparation at fixed geometry. This
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
