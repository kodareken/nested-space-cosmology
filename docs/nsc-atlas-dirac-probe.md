# Dirac probe on the actual quiet cap and material wall

The flat ratio `Psi/Phi=.9989169187321566` in the
[spectral-response owner](schur-memory.md#covariant-curvature-channels-and-a-static-gravitational-observable)
is an imported fixed-profile Dirac calculation about flat space. It is not
the response of this atlas. The [classical scalar action](nsc-atlas-cap-wall-response.md) does not contain a
Dirac field. Here a **massless, minimally coupled probe with no direct
surface spinor interaction** is an explicit extension, not a consequence of
scalar transmission. Work first on the known zero-scalar reference; the
unconverged P001 endpoint cannot supply a controlled corrected metric.

For a smooth spherical metric
`ds^2=C(-dt^2+dz^2)+r^2 dOmega^2`, the two-dimensional spin connection
contributes `d log(C)/4` and the sphere warp contributes `d log(r)`.
Consequently the spinor rescaling

$$
\psi=C^{-1/4}r^{-1}\,u
$$

removes both derivatives from the radial first-order operator. In an
angular spinor sector `k=+/-1,+/-2,...`, its angular coefficient is
`W=k*sqrt(C)/r`. There is no angular zero mode. This is an operator on the
curved metric, not the flat profile evaluated at a replacement momentum.
The [component](../src/recursive_horizons/nsc_atlas_dirac_probe.py) and
[controls](../tests/test_nsc_atlas_dirac_probe.py) keep these statements
independent of any quantum-state choice.

In the static cap, `A=1-H^2*r^2`, `C=A` and
`z=atanh(H*r)/H`. The radial Hamiltonian and its two squared channels are

$$
\mathsf{H}_k=-i\sigma_2\partial_z+\sigma_1 W,
\qquad W=\frac{k\sqrt{A}}{r}=\frac{kH}{\sinh(Hz)},
\qquad V_\pm=W^2\pm\partial_z W,
\qquad \partial_z W=-\frac{k\sqrt{A}}{r^2}.
$$

The original spinor norm `r^2*dr/sqrt(A)` becomes `dz` after the rescaling.
At `R/ell=.5`, `H*ell=1`, `k=1`, the
[one reference evaluation](../results/summaries/nsc-atlas-dirac-probe-reference-v1.json)
gives `W*ell=sqrt(3)`, `V_plus*ell^2=-.4641016151377548` and
`V_minus*ell^2=6.464101615137754`. Flat space at the same areal radius
has `W*ell=2`, `V_plus=0` and `V_minus*ell^2=8`. These are **bulk curvature
differences**, not a nesting-induced gravitational slip. The control which
continues the same de Sitter cap smoothly has precisely the same cap bulk
operator. A negative partner-potential value alone is not an unstable Dirac
mode: stability depends on the full squared operator and its domain.

The material interface is a distinct change. Use Gaussian normal distance
`n` increasing from cap to exterior, `K_ab=partial_n(gamma_ab)/2`, and
`[X]=X_plus-X_minus`. In a common spin frame parallel along the normal,
the local Dirac operator has normal part
`gamma^n*(partial_n+K/2)`. Transparent first-order transmission makes `psi`
continuous. Requiring the same domain for `D psi` gives

$$
[\partial_n\psi]=-\frac{[K]}{2}\psi.
$$

This follows directly from the first-order spin connection. It is a domain
condition for the squared operator, not an independently imposed reflecting
wall. The quiet reference's Israel jump is
`[K_ab]=-k_wall*gamma_ab`, where `k_wall*ell=2`; hence

$$
[K]=-\frac6\ell,\qquad
[\partial_n\psi]=\frac3\ell\psi.
$$

Here `k_wall` is the material parameter, distinct from the spinor angular
integer. Also `[n(r)]=beta_plus-beta_minus=-k_wall*r`, so the normal
derivative of `k/r` jumps by `k*k_wall/r`; it is `4/ell^2` for the same
lowest channel and radius. Both jumps vanish for the smooth-cap control.
The signed `beta` and its derivative recover the `K` jump without division
by `beta_plus` at its zero. Ten controls independently check the tetrad
rescaling, norm, squared Hamiltonian, flat/angular limits and interface
domain. The reference's Schwarzschild `f_plus=-.6` at this cut: its usual
static coordinate is spacelike. Extending the cap's static Hamiltonian into
that region is not a valid joined evolution.

The full loop calculation must propagate spinors through the actual regular
exterior and this moving domain. For a discretized matched problem its
cap resolvent has the Schur form
`G_cap=(D_cap-Sigma_wall)^(-1)`, with `Sigma_wall` obtained from the exterior
operator and interface. The smooth-cap comparison keeps `D_cap` and the
common source/readout fixed and changes that exterior/domain response.
Its difference isolates nesting; comparing directly with a flat cap also
changes local curvature. Geometric coefficients alone do not specify the
exterior state or the loop's finite interface terms. Relevant imported
methods are the Dirac and transmittal calculations in
[Vassilevich, Sections 3.3 and 6.3](https://arxiv.org/html/hep-th/0306138).

For a Lorentzian calculation the missing two-point state is a mathematical
input. Even in a finite free-fermion representation, the retarded response
of a bilinear vertex is `-i*Tr(C*[V(t),V(0)])`, with initial covariance `C`.
Changing `C` while keeping the geometry changes that expression whenever
the traced commutator is nonzero. The
[existing causal implementation](../src/recursive_horizons/nsc_causal_common.py)
already evaluates this dependence; no second response framework is needed.
The physical metric response also includes local contact terms, wall
displacement and the variation of the domain, not only a bulk commutator.
The [closed-time-path derivation](https://pmc.ncbi.nlm.nih.gov/articles/PMC5660882/)
provides the relevant state-dependent influence method.

Thus **no nested `Psi/Phi` or nesting-induced slip difference has yet been
computed**. Required inputs are a joined Dirac two-point state, bulk/interface
renormalization conditions, and a conserved perturbing source with a
gauge-invariant motion/clock/lensing protocol on this evolving background.
The flat fixed-profile action supplies none of these choices for the joined
cap--wall spacetime. This calculation identifies and evaluates the actual
geometry/domain changes; it does not label them as a completed loop or
replace the scalar candidate with the separate periodic Dirac model.

## Joined initial surface at the source anchor

The earlier [complete sourced slice](nsc-atlas-connection-slice.md) uses a
flat expanding cap at `Rwall/ell=1.5`. It is not initial geometry at the
current source anchor `.45`. The
[companion construction](../src/recursive_horizons/nsc_atlas_dirac_slice.py)
instead fixes `Tcap=0` in the quiet static cap at that anchor, keeping
`H*ell=1`, `k_wall*ell=2`, `M/ell=.4` and zero scalar field. The cap data are

$$
r=H^{-1}\sin(H\lambda),\qquad
\lambda_w=H^{-1}\arcsin(HR_w),\qquad K^i{}_j=0.
$$

Write `v=Rdot`, `Aminus=1-H^2*Rw^2`, and use the existing signed wall
coefficients `beta_minus, beta_plus`. In the common wall frame the cap's
slice normal and spatial tangent are the boost

$$
u_\Sigma=\gamma u_w-\gamma_v n_w,\quad
e_\Sigma=\gamma n_w-\gamma_v u_w,\quad
\gamma=\frac{\beta_-}{\sqrt{A_-}},\quad
\gamma_v=\frac{v}{\sqrt{A_-}}.
$$

The continued slice therefore has exterior wall data

$$
E=u_\Sigma(r)=\gamma v-\gamma_v\beta_+>0,\qquad
J=e_\Sigma(r)=\gamma\beta_+-\gamma_v v<0,\qquad
J^2-E^2=1-\frac{2M}{R_w}.
$$

Hold `V=u_Sigma(r)=E` initially along the Schwarzschild portion. With
`a=sqrt(1+E^2)` and `rmin=2M/a^2`, an exact proper-distance parametrization is

$$
r=r_{\min}\cosh^2q,\quad r'=a\tanh q,\quad
\lambda=\lambda_{\min}+\frac{r_{\min}}a
(q+\sinh q\cosh q).
$$

Start at `qwall=-acosh(sqrt(Rw/rmin))`, aligning the proper distance with
`lambda_w`. The minimum at `q=0` is regular and is a property of this slice;
the spacetime event horizon remains `r=2M`. On the increasing-radius branch
use the existing smooth step to taper `V` from `E` to zero between
`r/ell=1` and `1.2`. For either a constant or varying `V(r)` the vacuum data

$$
r'^2=1-\frac{2M}{r}+V^2,\quad
r''=\frac{M}{r^2}+VV_r,\quad
K_{\lambda\lambda}=-V_r,\quad K^\theta{}_\theta=-\frac Vr
$$

satisfy the spherical Hamiltonian and momentum constraints. After the
taper, `K=0` and the same Schwarzschild data extend to one AE end. The cap
centre closes the other direction. This supplies sourced initial data,
without asserting a regular formation history or global past completeness.
The taper is a slicing choice; using it later to select a quantum state
would make it an explicit preparation choice.

No division by the Schwarzschild factor is needed at the horizon. In
outgoing Eddington--Finkelstein coordinates
`ds2=-f*du2-2*du*dr+r2*dOmega2`, the spatial tangent and unit normal have

$$
e_\Sigma^u=-\frac1{V+r'},\quad e_\Sigma^r=r',\qquad
u_\Sigma^u=\frac1{V+r'},\quad u_\Sigma^r=V.
$$

Their norms are `+1,-1` and their mutual contraction vanishes, including
at `f=0` on the outgoing segment. At the material interface,

$$
[r']=-k_{\rm wall}R_w\gamma,\quad
[K^\theta{}_\theta]=-k_{\rm wall}\gamma_v,\quad
{}^{(3)}R_{\rm singular}=-4[r']/R_w.
$$

With `k_wall=4*pi*sigma`, these are precisely the surface Hamiltonian and
momentum sources for the independently projected energy `sigma*gamma` and
momentum `sigma*gamma_v`. The
[controls](../tests/test_nsc_atlas_dirac_slice.py) also differentiate the
embedding and fields independently, check the horizon and minimum directly,
and verify length-scale covariance.

## Intrinsic spatial operator versus quantum preparation

On this initial three-metric the intrinsic radial operator in angular
sector `k` is unitarily represented by

$$
D_{\Sigma,k}=-i\sigma_2\partial_\lambda+\sigma_1 W_\Sigma,\quad
W_\Sigma=k/r,\quad \partial_\lambda W_\Sigma=-kr'/r^2.
$$

Here `u=r*psi` and the norm is `integral |u|^2 d(lambda)`. This is a
different operator from the static-cap optical Hamiltonian above. At the
wall `W_Sigma` is continuous while its first derivative has the finite jump
`k*k_wall*gamma/Rw`. The squared **intrinsic spatial** operator has partner
potentials `W_Sigma^2 +/- d_lambda W_Sigma`, with a step and no delta.
Its transparent square has continuous `u` and `d_lambda u`. For the
first-order spatial operator only continuity of `u` is required; its domain
allows square-integrable piecewise derivatives.

These spatial domain statements are not matching conditions for the
four-dimensional squared operator. Its normal is the moving material
wall's normal, and its spin connection, time derivatives, lapse and frame
transport must be retained when pulling it back to a slice. Selecting the
negative spectral projector of `D_Sigma` would be an additional state
ansatz, not a derivation of the physical covariance. The smooth-spacetime
[Dirac Hadamard construction](https://arxiv.org/html/2108.11630v3), Sections
3.2, 6 and 7, supplies methods to inspect; its smoothness and evolution
hypotheses have not been established at this material wall. The joined
surface removes a geometric obstruction to specifying the state, without
silently supplying that missing state or a value of `Psi/Phi`.

The [source-anchor evaluation](../results/summaries/nsc-atlas-dirac-source-slice-v1.json)
uses 231 geometric points and `.046875` process-CPU seconds. The largest
algebraic constraint/embedding discrepancy is `8.89e-15`; the fixed 32/64
outgoing-distance comparison differs by `2.03e-13*ell`. These are numerical
consistency comparisons, not a physical uncertainty enclosure. It gives
`rmin/ell=.40535925102607046`, while the event horizon is at `.8`, and
`E=.9866904367219993`. At the same wall and lowest angular sector:

| Intrinsic spatial coefficient | Cap side | Exterior side |
| --- | ---: | ---: |
| `W_Sigma*ell` | 2.2222222222222223 | 2.2222222222222223 |
| `d_lambda W_Sigma*ell^2` | -4.410017555430062 | 2.1850400461857395 |
| `V_plus*ell^2` | .5282540495082104 | 7.123311651124012 |
| `V_minus*ell^2` | 9.348289160368335 | 2.7532315587525327 |

The derivative jump is `6.595057601615801/ell^2`, agreeing with
`k*k_wall*gamma/Rw`. A smooth continuation of the same cap has no such jump.
This is a quantified spatial interface change on the joined initial data,
not a loop polarization or a nesting-induced `Psi/Phi` value.

## Physical time generator at the moving interface

The intrinsic spatial operator is insufficient even to supply the initial
time generator at the wall. On the slice above take zero shift, radial
metric length `a=1` and `K_lambda_lambda=0` on each side. With
`ds2=-N^2*dt^2+a^2*dlambda^2+r^2*dOmega2`, write
`a_Sigma=partial_lambda log(N)` and `s=gamma_v/gamma`. The wall-frame
mixed curvature is `Ktau^tau=n_w dot nabla_(u_w) u_w`. Its acceleration is

$$
K^\tau{}_{\tau,w}=\partial_\tau\theta+
\gamma a_\Sigma-\gamma_v K_{\lambda\lambda},\qquad
\sinh\theta=\gamma_v.
$$

The boost and its derivative along the wall are common. Therefore its
temporal Israel condition requires

$$
[a_\Sigma]=-\frac{k_{\rm wall}}\gamma.
$$

This is a necessary lapse jet for this initial slice and frame. A globally
constant lapse would lose this part of the material-wall geometry. Normalize
the cap's static clock to `Nwall=1` at the initial instant:
`Nminus=sqrt(1-H^2*r^2)/sqrt(1-H^2*Rw^2)`. Its wall acceleration is
`a_Sigma_minus=-H^2*Rw/sqrt(1-H^2*Rw^2)`, and the jump fixes the exterior
first jet. It does not select the exterior lapse farther away or its higher
time derivatives.

The half-density `chi=sqrt(a)*r*psi_Sigma` has norm
`integral |chi|^2 dlambda`. A direct diagonal-tetrad contraction gives the
radial spin connection

$$
\frac{\gamma^0}{N}\left(\frac{a_t}{2a}+\frac{r_t}r\right)
+\frac{\gamma^1}a\left(\frac{N_\lambda}{2N}+\frac{r_\lambda}r\right).
$$

The time derivative of the half-density removes the trace-of-slice-curvature
term. With `c=N/a` and `alpha=sigma_2`, the one-sided massless generator is

$$
H=-i\alpha\left(c\partial_\lambda+
\tfrac12\partial_\lambda c\right)+\sigma_1\frac{Nk}r.
$$

The convention is `alpha=-gamma_hat_t*gamma_hat_lambda` with spacetime
Clifford signature `(-,+,+,+)`. At the initial instant `c=N`; its derivative
companion contains the step in
the lapse gradient. The [component](../src/recursive_horizons/nsc_atlas_dirac_evolution.py)
retains it. The [controls](../tests/test_nsc_atlas_dirac_evolution.py) compute
the connection from the full four-dimensional diagonal tetrad, verify the
half-density and Green identity, and recover the earlier cap optical
Hamiltonian by its unitary coordinate/density change. The independent
[3+1 derivation](https://arxiv.org/html/2503.03918), Eq. IX.29, also retains
the lapse-gradient and slice-curvature terms before rescaling.

Transparency is continuity along the **moving** wall. Differentiating it
and using the one-sided equation gives

$$
[\chi_t]=-N_ws[\chi_\lambda],\qquad
[\chi_\lambda]=-
\frac{[a_\Sigma]}2(I-s\alpha)^{-1}\chi
=\frac{k_{\rm wall}\gamma}2(I+s\alpha)\chi.
$$

Here `|s|<1`, so the inverse exists. Requiring both the time generator and
spatial derivative of `chi` to be continuous would be the wrong moving
domain. Restoring the derivatives of `r` using the checked geometric jumps
gives, in the slice spin frame,

$$
[\partial_\tau\psi_\Sigma]=0,\qquad
[\partial_n\psi_\Sigma]=
\left(\frac{3k_{\rm wall}}2 I+
\frac{k_{\rm wall}s}2\alpha\right)\psi.
$$

The normal component of this frame's two-dimensional connection has jump
`[omega(n)]=gamma_v*[a_Sigma]=-k_wall*s`. Its spinor term
`alpha*[omega(n)]/2` cancels the matrix contribution. In a frame parallel
along the material normal the result is exactly the earlier condition
`[partial_n psi]=3*k_wall*psi/2=-[K_wall]*psi/2`. No surface spinor
interaction has been added. Setting the lapse-gradient step to zero would
instead give `k_wall*psi`, missing the temporal part of the connection.

These are compatible initial one-sided equations and wall jets. They do not
yet supply a full foliation, propagated covariance, renormalized polarization
or gravitational slip.

The [one-point evaluation](../results/summaries/nsc-atlas-dirac-wall-time-v1.json)
uses `.015625` process-CPU seconds. At the source anchor `s=.7388177118` and
the required logarithmic lapse-gradient jump is `-1.3478106524/ell`.
The two eigenvalues of the density-derivative jump matrix are
`.3875652529/ell` and `2.5802106678/ell`; these are derivative coefficients,
not quantum probabilities or reflection amplitudes. The temporal Israel,
moving continuity and normal-frame checks differ by at most `2.67e-15` in
the reference units. In the parallel wall-normal frame the recovered
coefficient is `3/ell`, exactly the earlier four-dimensional result.

## Finite joined time chart

The next connection supplies a time-dependent background for propagation.
It is a coordinate pullback of the same quiet de Sitter and Schwarzschild
metrics, with the same moving material wall. It does not use the failed
finite-amplitude P001 endpoint as a corrected metric.

Let `v(R), beta_minus(R), beta_plus(R)` be the existing signed quiet wall
data and `A_w=1-H^2*Rw^2`. Define the constant initial clock normalization

$$
Q_0=\frac{\beta_{-,0}}{v_0\sqrt{A_{w,0}}},\quad
R_w(t)=R_{w,0}+t/Q_0,\qquad R_{w,0}/\ell=.45.
$$

The original proper wall velocity and static-cap clock then give

$$
\frac{d\tau_w}{dt}=\frac1{Q_0v},\qquad
\frac{dT_{\rm cap}}{dt}=\frac{\beta_-}{Q_0 A_wv},\qquad
\lambda_w=H^{-1}\arcsin(HR_w),\quad
\frac{d\lambda_w}{dt}=\frac1{Q_0\sqrt{A_w}}.
$$

Use the coordinate `y=lambda-lambda_w(t)` pinned to the wall, where `lambda`
is cap proper distance. The exterior extends the initial spatial label;
it need not remain proper distance there. On `-lambda_w<=y<=0` the cap is
explicitly

$$
r=H^{-1}\sin(H(y+\lambda_w)),\quad
a_-=1,\quad N_-=\sqrt{1-H^2r^2}\,T_{{\rm cap},t},\quad
b_-=\lambda_{w,t}.
$$

Here the shift `b` is distinct from the signed normal wall quantities
`beta_plus, beta_minus`. At the wall its value is `Nw*s`, with
`Nw=gamma/(Q0*v)` and `s=gamma_v/gamma`. Initially `Nw=1`.

On the exterior use `X=(U,r)` in outgoing Eddington--Finkelstein coordinates
with `g2=-f*dU^2-2*dU*dr`, `f=1-2M/r`. The initial slice already fixes
`X0(y)`, `r0_y`, `V0=u_Sigma(r)` and `U0_y=-1/(V0+r0_y)`. Select the
declared initial coordinate lapse and shift

$$
S(y)=1-\operatorname{step}(y/\ell),\quad
N_0(y)=\exp(a_{\Sigma+,0}yS(y)),\quad
b_0(y)=s_0S(y).
$$

This uses the same smooth step as the initial-slice taper. It is flat at
the endpoints, and `S=0` for `y>=ell`. It fixes gauge data, not a quantum
state. The initial coordinate velocity in the exterior is

$$
\dot U_0=\frac{N_0-b_0}{V_0+r_{0,y}},\qquad
\dot r_0=N_0V_0+b_0r_{0,y}.
$$

The exact material trajectory in the same exterior coordinates has

$$
\frac{dU_w}{dR_w}=\frac1{v(\beta_++v)},\qquad
e_w(t)=\left(-\frac1{E+J},J\right),\quad
E=\frac{k_{\rm wall}R_wv}{\sqrt{A_w}},\quad
J=\sqrt{A_w}-\frac{k_{\rm wall}R_w\beta_-}{\sqrt{A_w}}.
$$

No division by `beta_plus` or `f` is made. With `F0=S` and `F1=y*S`, set

$$
\begin{split}
X(t,y)={}&X_0(y)+t\dot X_0(y)\\
&+\bigl(X_w(t)-X_0(0)-t\dot X_0(0)\bigr)F_0(y)\\
&+\bigl(e_w(t)-X_{0,y}(0)-t\dot X_{0,y}(0)\bigr)F_1(y).
\end{split}
$$

The values and first derivatives of `F0,F1` enforce the wall embedding and
common spatial tangent exactly: `X(t,0)=Xw(t)`, `X_y(t,0)=ew(t)`. The
initial slice is `X0`. Recovering the prescribed initial coordinate velocity
also requires `e_w,t(0)=dot(X0)_y(0)`; this is checked using the temporal
Israel lapse jet above. Future slices need not retain constant `V` or
vanishing `K_lambda_lambda` at the wall.

Differentiate the embedding to extract the full radial metric:

$$
h=-fU_y^2-2U_yr_y,\quad
B=-fU_tU_y-U_tr_y-r_tU_y,\quad
D=U_tr_y-r_tU_y,\qquad
a=\sqrt h,\quad N=D/a,\quad b=B/h.
$$

The chart is admitted only where `h>0,D>0`; taking the absolute value of
`D` would conceal a reversal of its time orientation. The signed lapse
makes the time gradient timelike. The coordinate vector `partial_t` need
not itself be timelike for a valid foliation. Since each exterior point is
a pullback of Schwarzschild, no extra bulk matter is introduced. The metric
and its derivatives still have to be computed consistently; the spherical
constraints and geometric mass provide independent checks.

For this metric the half-density generator acquires shift transport:

$$
H=-i\left[\alpha\left(c\partial_y+\frac{c_y}2\right)
-I\left(b\partial_y+\frac{b_y}2\right)\right]
+\sigma_1\frac{Nk}r,\qquad c=N/a.
$$

The interface is fixed at `y=0`. Put `B_D=c*alpha-b*I`. Its transparency
condition implies

$$
[\chi_y]=-\tfrac12 B_D^{-1}[B_{D,y}]\chi.
$$

At the wall `a=1`, `b=N*s` and `|s|<1`, so this inverse exists. Restoring
the spatial half-density derivatives and the normal spin connection must
recover `3*k_wall/2` in the parallel material-normal frame. The shift term
and its half-derivative are required in that check.

The [admitted reference calculation](prior-art-and-open-claim.md#a-finite-joined-time-chart-for-the-actual-quiet-dirac-background)
uses fixed wall-radius cuts and geometric nodes. Sampled positivity is a
numerical check, not a uniform interval enclosure. No propagated quantum
covariance or loop polarization is implied by the coordinate construction.
The [component](../src/recursive_horizons/nsc_atlas_dirac_time_chart.py),
[nine controls](../tests/test_nsc_atlas_dirac_time_chart.py) and
[bounded recorder](../scripts/measure_nsc_atlas_dirac_time_chart.py) preserve
the signed coordinate-domain test; a failed cut does not trigger a gauge
scan or a claim of full time coverage.

The [frozen reference](../results/summaries/nsc-atlas-dirac-time-chart-v1.json)
records an **unsuccessful full-domain chart**. Five predeclared cuts through
`Rw/ell=.5` pass the sampled metric checks. The next cut,
`Rw/ell=.5125`, has signed Jacobian `-.2130723253` at initial spatial label
`y/ell=.5214966935`, with areal radius `r/ell=.5977107183`. The sampled
spatial norm remains positive at that cut; the failed determinant is a
time-orientation/coordinate fold. This is the first failed predeclared cut,
not a located exact onset radius. The remaining bulk cuts were not run.

On the successful metric calls the smallest lapse is `.01403670785`,
smallest spatial norm `.06659155747`, and smallest Jacobian `.01219117082`.
Mass, Hamiltonian and momentum consistency discrepancies are below
`6.79e-13`; the wall and spinor checks differ by less than `6.22e-15`.
The largest 32/64 coordinate comparison is `1.08e-12*ell`. These comparisons
are not a physical-error enclosure. The application uses 9024 quadrature
points and `.0625` process-CPU seconds, including setup and source binding.
No gauge profile was changed and no field evolution or polarization was run.

The fixed quiet spacetime and its material trajectory still supply the wall
data. This particular wide time chart does not cover the entire strip. A
subsequent response calculation must establish the chart coverage actually
needed for its source and observer, using compatible patches where required.
The failed chart cannot be used as a full-strip Dirac evolution or as an
explanation of P001's separate projection error.

## Conformal calculation variables remove the singular wall jet

For this massless probe a useful exact simplification follows from
four-dimensional conformal covariance. On `r>0`, use

$$
\Omega=\ell/r,\qquad \widetilde g=\Omega^2g,\qquad
\widetilde\psi=\Omega^{-3/2}\psi,\qquad
\widetilde D\widetilde\psi=
\Omega^{-5/2}D\psi.
$$

This is a change of calculation variables for the probe. The gravitational
action, physical metric and observer clocks remain those of `g`. The
conformal weights are checked directly in the
[controls](../tests/test_nsc_atlas_dirac_conformal.py); the smooth one-sided
operator identity also appears in
[Gerard--Stoskopf, 2.7.2](https://arxiv.org/html/2108.11630v3). Its use here
is local on each side and in the explicit interface limits; no smooth
global-state theorem is applied to the material wall.

The original pure-tension jump obeys
`[K^a_b]=-k_wall*delta^a_b` and `[n(log r)]=-k_wall`. Consequently

$$
[\widetilde K^a{}_b]=\Omega^{-1}
\bigl([K^a{}_b]+\delta^a{}_b[n\log\Omega]\bigr)=0.
$$

Both the temporal and angular components cancel independently. The normal
spinor jump transforms in the same way:

$$
[\widetilde n(\widetilde\psi)]
=\Omega^{-5/2}\bigl([n(\psi)]
-\tfrac32[n\log\Omega]\psi\bigr)=0.
$$

Here common spin frames are identified at the wall; tangential derivatives
of the conformal factor are shared. Thus the coefficient `3/ell` in the
original variables is not a reflection amplitude or an independent surface
spinor interaction. In the conformal variables the corresponding first
normal jets are continuous. This does not remove curvature variation away
from the wall or establish a vanishing scattering or gravitational response.

The transformed metric is a product of the radial two-metric and a sphere
of constant radius `ell`. For a smooth one-sided radial factor
`g2=-f*dt2+dr2/f`, its radial curvature is

$$
\widetilde R_{(2)}=\ell^{-2}
\bigl(-r^2f''+2rf'-2f\bigr).
$$

This expression has no division by `f` and continues through a coordinate
horizon. In particular,

$$
\widetilde R_{(2),-}=-2/\ell^2,\qquad
\widetilde R_{(2),+}=(-2+12M/r)/\ell^2.
$$

A direct Ricci contraction of the transformed two-dimensional metric checks
these results independently. The cap operator in these variables lives on
the constant-curvature radial geometry, with sphere eigenvalues `k/ell`;
the nested exterior retains a finite radial curvature change `12M/(r*ell^2)`.
A smooth continuation of the same cap has no such change. The transformed
Weyl invariant is zero in the cap and `48*M^2/(r^2*ell^4)` outside. No
curvature distribution is squared in these calculations.

The conformal map is singular at the physical centre and asymptotic end.
Regularity and incoming-state data at those ends still have to be imposed;
they cannot be erased by the local operator simplification. At the quantum
level the measure and renormalization also generate a conformal anomaly;
[Vassilevich, 7.1](https://arxiv.org/html/hep-th/0306138) gives its effective-action
relation. A two-dimensional mode calculation must retain a consistent
four-dimensional regulator and renormalization prescription. None of these
classical cancellations establishes a quantum Weyl symmetry or a nested
value of `Psi/Phi`.

The [one source-anchor evaluation](../results/summaries/nsc-atlas-dirac-conformal-interface-v1.json)
uses `.03125` process-CPU seconds. All transformed wall first-jet jumps are
zero in the reference arithmetic. At `Rw/ell=.45` its radial curvature is
`-2/ell^2` in the cap and `8.6666666667/ell^2` outside, a finite difference
`10.6666666667/ell^2`. The transformed exterior Weyl square is
`37.9259259259/ell^4`. The independent curvature identity discrepancy is
`1.78e-15` in these reference units. These are operator/geometry inputs,
without an inferred quantum state or a reflection/slip amplitude.

## Local one-loop Dirac trace on the original physical geometry

The [local trace component](../src/recursive_horizons/nsc_atlas_dirac_trace.py)
uses the original `g`, with one free massless four-component Dirac field.
Free Dirac quantization and the general Laplace-type heat coefficient are
imported. The spin and statistics contractions are calculated explicitly:
`E=-R/4`, `tr I=4` and `tr(Omega_mu_nu*Omega^mu_nu)=-Riemann^2/2`.
The fermionic determinant reverses the sign of the traced heat coefficient.
For the two quiet smooth regions, `Box R=0`, giving

$$
\frac{\langle T^\mu{}_\mu\rangle}{\hbar c}
=\frac{7R_{\alpha\beta\gamma\delta}R^{\alpha\beta\gamma\delta}
+8R_{\alpha\beta}R^{\alpha\beta}-5R^2}{360(4\pi)^2}.
$$

The [controls](../tests/test_nsc_atlas_dirac_trace.py) contract explicit
four-component Euclidean Clifford matrices and reproduce the corresponding
Weyl/Euler coefficients `c=1/20`, `a=11/360`. The imported determinant and
heat-coefficient framework is described in
[Vassilevich, 3.3, 4.1 and 7](https://arxiv.org/html/hep-th/0306138).
The one-sided Lorentzian curvature invariants are already checked directly
in the [probe controls](../tests/test_nsc_atlas_dirac_probe.py). Substitution
therefore yields

$$
\frac{\langle T^\mu{}_\mu\rangle_-}{\hbar c}
=-\frac{11H^4}{240\pi^2},\qquad
\frac{\langle T^\mu{}_\mu\rangle_+}{\hbar c}
=\frac{7M^2}{120\pi^2r^6}.
$$

A smooth continuation of the same cap retains the cap value. A flat
background has zero curvature trace. This actual geometric local quantum
term is distinct from the old flat nonlocal polarization.

For the trace projection only, let `lP^2=hbar*G/c^3`, `h=H*ell`,
`m=M/ell`, `x=r/ell`. The corresponding contribution has

$$
\Delta\Lambda_{\rm trace}\ell^2
=\frac{l_P^2}{\ell^2}
\begin{cases}
11h^4/(120\pi),&\text{cap},\\
-7m^2/(60\pi x^6),&\text{exterior}.
\end{cases}
$$

This projects `T` onto `-T*g/4`; it does not prove that the full stress has
vacuum form. The exterior coefficient varies with radius. The independently
declared `ell=10 km` fixes an illustrative physical scale, as in P001; it
is not calculated from a target Lambda. The cap's `3/ell^2` remains an
input of the quiet geometry. State-dependent traceless stress, interface
renormalization, causal polarization and backreaction are not supplied by
this local trace. Material-wall curvature distributions are not squared.
The scheme-dependent total derivative cannot be ignored in a wall or
nonconstant-curvature calculation merely because it vanishes here.

The [actual bulk-trace evaluation](../results/summaries/nsc-atlas-dirac-bulk-trace-v1.json)
at the same source anchor and fiducial `ell=10000 m` gives:

| Smooth one-sided region | `trace*ell^4/(hbar*c)` | `DeltaLambda_trace*ell^2` |
| --- | ---: | ---: |
| Cap | -.004643887583607148 | 7.622217590518666e-80 |
| Exterior at `r/ell=.45` | .1138837996360018 | -1.8692250517106162e-78 |

The shape factors multiplying `lP^2/ell^2` are `.029178406233514146`
and `-.7155530165987104`. Identity comparisons differ by `1.12e-16` or
less. The measured process time rounds to zero at this host's timer
resolution; no field evolution or nonlocal polarization was evaluated.
The constants, source binding and interpretation limits are retained in
the compact record. These numbers quantify the curvature trace of the
declared probe, not a prediction of the full observed cosmological Lambda.

## Euler-curvature contribution at the material interface

The same imported Dirac type-A coefficient also acts on the actual wall's
Euler curvature. This calculation specifies a thin limit: smooth the kinked
conformal factor while the conformal product curvature remains bounded.
It does not assume equivalence with arbitrary microscopic walls carrying
tracefree-curvature excursions. Neither a delta square nor the uncomputed
scheme-dependent `Box R` is assigned a numerical value.

With one normal pointing from cap to exterior, put

$$
t_\pm=\beta_{\pm,R},\quad b_\pm=\beta_\pm/R,\quad
A=\ddot R/R,\quad B=(1+\dot R^2)/R^2.
$$

The induced worldtube has `Ghat^i_j=diag(-B,-A,-A)`. Its Chern scalar and
integrated Euler coefficient are

$$
\mathcal B_3=4(J-2\widehat G_{ij}K^{ij})
=-8tb^2+8Bt+16Ab,\qquad J=-2tb^2,
$$

$$
\mathcal E_\Sigma=-[\mathcal B_3]
=8[tb^2]-8B[t]-16A[b].
$$

A direct smooth Lorentz curvature contraction checks the normal derivative
terms before taking any limit:

$$
(E_4)_{\rm normal}=-8t_n(B-b^2)-16b_n(A-tb)
=-\partial_n\mathcal B_3\big|_{\rm leading},
\qquad E_4=Riemann^2-4Ricci^2+R^2.
$$

Finite tangential/Codazzi terms do not contribute a leading normal
distribution. Pure-trace singular normal slopes give coefficients
`12*k^2`, `12*k^2`, `36*k^2` in the three squared curvatures, respectively.
They cancel in both `E4` and `Weyl^2`; the Euler normal term integrates
through its transgression instead. The four
[controls](../tests/test_nsc_atlas_dirac_euler_wall.py) contract the full
Lorentz-index curvature tensor, check the flat-ball boundary orientation
and static Schwarzschild Euler flux, and independently reduce the actual
signed wall equations. Reversing the common normal and side labels leaves
the surface coefficient unchanged.

For this quiet junction `[t]=[b]=-k_wall`, while
`A=t_minus*b_minus+H^2`, `B=b_minus^2+H^2` and
`2*b_minus+t_minus=3*(k_wall^2-H^2)/(2*k_wall)`. Consequently,

$$
\mathcal E_\Sigma=4k_{\rm wall}(k_{\rm wall}^2+3H^2)
=56/\ell^3.
$$

Exterior mass and instantaneous wall radius cancel. A smooth continuation
of the cap has no such surface coefficient. The corresponding Euler part
of the local Dirac trace is

$$
\frac{\langle T^\mu{}_{\mu}\rangle_{\Sigma,E}}{\hbar c}
=-\frac{11\mathcal E_\Sigma}{360(4\pi)^2}\delta_\Sigma
=-\frac{77}{720\pi^2\ell^3}\delta_\Sigma,
$$

where `delta_Sigma` is the invariant delta in proper normal distance.
This specifies one local interface effect of the attachment. The
scheme-dependent total derivative, finite interface terms, quantum state,
traceless stress and nonlocal response still have to be included to obtain
the full wall stress or `Psi/Phi`. The formal trace-only projection has
`ell*integral dn DeltaLambda_trace=77/(360*pi)*(lP^2/ell^2)`; it is an
integrated surface term, not a bulk cosmological constant or an assertion
of vacuum stress.

The [frozen wall evaluation](../results/summaries/nsc-atlas-dirac-euler-wall-v1.json)
at `Rwall=.45*ell` gives `B3_minus=-486.98504145775826/ell^3` and
`B3_plus=-542.9850414577583/ell^3`, hence `E_Sigma=56/ell^3`.
The trace coefficient is `-.010835737695083347/ell^3`. With the same
supplied `ell=10000 m`, the trace-only integrated projection is
`ell*integral dn DeltaLambda_trace=1.778517437787688e-79`.
The two independently formed junction/trace comparisons differ by zero
in reference arithmetic. Evaluation takes `.03125` process-CPU seconds,
with no field evolution, smoothing scan or new state. This is the Euler
surface contribution within the declared class; it is not a total
physical error estimate for the uncomputed quantum metric response.

## The curved metric-response calculation and its state dependence

The exact [wall-normal collar](nsc-atlas-conformal-action.md#a-shared-finite-field-action-on-a-quiet-wall-normal-collar)
provides a regular local time generator throughout the source strip. It
does not provide a global quantum state. In units `hbar=c=1`, the same
radial Dirac reduction, now with `a=1` and density spinor `xi=r*psi`, gives

$$
H_k=-i\alpha(N\partial_n+N_n/2)+\sigma_1 Nk/r.
$$

For a prescribed Gaussian-coordinate variation with the source off on the
initial surface, the operator insertion is

$$
V_k=\delta H_k=-i\alpha(\delta N\partial_n+\delta N_n/2)
+\sigma_1 k(\delta N/r-N\delta r/r^2).
$$

These coefficients use the actual cap/exterior `N,r`, not a flat
polarization multiplied by a new profile. A full covariant stress
variation also varies the radial metric and shift before fixing gauge.
The collar alone cannot impose a reflecting-box vacuum: the state and
causal propagator must include the regions influencing the measured
history, with their actual interface transmission and initial data.

At a finite canonical fermion regulator let
`C_ij=<c_j^dagger*c_i>`, `0<=C<=I`, and let `U` be the one-particle
propagator of this background. For an observable `c^dagger*O*c`, the CAR
identity `[c^dagger*A*c,c^dagger*B*c]=c^dagger*[A,B]*c` gives

$$
i\dot C=[H,C],\qquad
i\dot{\delta C}=[H,\delta C]+[V,C],\qquad \delta C(t_0)=0,
$$

$$
\delta\langle O(t)\rangle
=\operatorname{Tr}(C_0\,(\delta O)_I(t))
-i\int_{t_0}^t ds\,
\operatorname{Tr}\{C_0[O_I(t),V_I(s)]\},
\qquad O_I(t)=U(t,t_0)^\dagger O(t)U(t,t_0).
$$

Here `(delta O)_I=U^dagger*(delta O)*U` varies the observable at fixed
background propagator; the commutator term accounts for its dynamical
response. The first term is essential when the physical stress/observer operator
depends on the varied metric. Removing the regulator additionally needs
the covariant local and interface counterterms and their contact
variations, with the gravitational/junction Ward identities checked.
The formula is an exact finite-regulator identity and a specification of
the missing response calculation; no continuum polarization has yet been
evaluated. Its kernel depends on both the curved propagator and `C0`.

This dependence cannot be inferred from the local trace. For example, if
an admissible preparation has an unoccupied normalized smooth mode `u`,
then `C0+nu*|u><u|`, `0<=nu<=1`, is another allowed fermion occupation.
Where the original preparation has the required short-distance
singularity, a smooth transmitting finite-rank addition preserves it.
The corresponding finite state-dependent change is

$$
\Delta\delta\langle O(t)\rangle
=\nu\langle u,(\delta O)_I(t)u\rangle
-i\nu\int_{t_0}^t ds\,
\langle u,[O_I(t),V_I(s)]u\rangle.
$$

It can change the response while leaving the state-independent anomaly
above unchanged. The fixed local trace is therefore insufficient to
determine the anisotropic stress or a unique nested `Psi/Phi`.

The next physical calculation must specify the preparation `C0`, hold it
fixed under a declared later source, propagate its response on sufficient
actual chart coverage, and solve the same gravitational and wall
perturbation equations for a declared clock or light measurement. The
two flat static potentials cannot be assigned to a dynamical nested
geometry without that observer/gauge connection. These requirements are
missing model connections, not a claim that the already calculated
curvature response vanishes. No state, counterterm value or slip number
is silently selected here.
