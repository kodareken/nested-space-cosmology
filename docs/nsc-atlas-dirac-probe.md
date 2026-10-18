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

## Proper-wall radar coordinates give a finite causal propagation chart

The [radar component](../src/recursive_horizons/nsc_atlas_dirac_radar.py)
uses the same quiet trajectory and `.45<=Rwall/ell<=.55`. On each static
side let `u_i=t_i-rstar_i`, `v_i=t_i+rstar_i`, with
`rstar_minus=atanh(r)` and `rstar_plus=r+.8*log(1-r/.8)` on this
white-hole interior branch. Along the physical wall,

$$
\frac{du_i}{d\tau}=\frac1{\beta_i+\dot R},\qquad
\frac{dv_i}{d\tau}=\frac1{\beta_i-\dot R}.
$$

The cap maps both increase. The exterior `u` increases and `v` decreases;
keeping that sign is essential inside the horizon. Their inverses label
each null ray by the proper time of its wall intersection, `U,V`. Put

$$
T=(U+V)/2,\qquad x=(V-U)/2,\qquad
ds_2^2=C_i(T,x)(-dT^2+dx^2),
$$

$$
r_i^*(r)=\frac{v_i(V)-u_i(U)}2,\qquad
C_i=f_i(r)u_i'(U)v_i'(V)>0.
$$

The material wall is `U=V=tau`. Use the cap for `V<=U`, exterior for
`V>=U`. The null square `0<=U,V<=S`, `S=tau(.55)-tau(.45)`, is exactly
the causal diamond `0<=T-x,T+x<=S`. These are wall radar coordinates:
`T` is the mean emission/return wall clock and `|x|` its radar distance.
The construction has no spatial cutoff profile and no timelike box ends.
It preserves the failed wide-chart result as a different coordinate
construction; it does not tune that chart or change the sources.

Exterior optical monotonicity bounds its radius by `.45<=r/ell<=.55` on
the entire square. Conservative cap and lapse bounds follow before
sampling: `v>=sqrt(.8/.55-1)`, `v<max(beta_minus)` bound the wall duration
and each optical increment. Applying those increments to the two tortoise
inverses gives positive cap radius and positive `C`. The component retains
these analytic bounds separately from sampled minima.

At the wall `C=1`, `r_T=Rdot`, `r_x=beta_i`, `C_T=0` and
`C_x=2*beta_i,R`; both Israel jumps are retained. The canonical density
spinor and its actual time generator simplify to

$$
\xi=C^{1/4}r\psi,\qquad
H_\kappa=-i\sigma_2\partial_x+\sigma_1 W_\kappa,
\qquad W_\kappa=\kappa\sqrt C/r.
$$

`W` and its first derivatives are continuous at the pure-tension wall.
In this frame `xi` and `xi_x` are continuous. Undoing its density factor
still gives the original physical parallel-normal spinor jump
`[psi_n]=3*k_wall*psi/2`. The radial Dirac characteristics now have constant
speed. With `B=(I+i*sigma1)/sqrt(2)` and `eta=B^dagger*xi`,

$$
2\partial_V\eta_+=-iW_\kappa\eta_-,\qquad
2\partial_U\eta_-=-iW_\kappa\eta_+.
$$

Thus incoming data `g_plus(U)=eta_plus(U,0)` and
`g_minus(V)=eta_minus(0,V)` define the characteristic propagation.
Conservation of the Dirac current gives their norm
`integral_0^S |g_plus|^2 dU + integral_0^S |g_minus|^2 dV`
on the middle Cauchy segment `T=S/2,|x|<=S/2`. An inherited incoming
covariance must include both branches and their cross-correlations.
The existing [Gaussian boundary owner](../src/recursive_horizons/nsc_boundary_state.py)
can retain these state terms once supplied with the actual transfer; it
does not choose that covariance or its angular ultraviolet completion.

Six [controls](../tests/test_nsc_atlas_dirac_radar.py) compare the optical
maps with existing independent quadrature, pull back both physical metrics,
check tortoise inverses, derivatives, curvature and mass, and recover the
wall jumps and constant-speed generator. The finite clock protocol emits
at wall `U=tau(.475)` and receives at `(U,V)=(tau(.475),tau(.525))` on a
fixed-radar-position observer. Its radial optical frequency ratio is
`omega_received/omega_wall=1/sqrt(C_receiver)`, from the null phase `-U`
and the receiver's proper clock. This is a defined finite background
measurement; it is not a raw Fermi-origin metric ratio.

The diamond is a local physical chart rather than a whole global Cauchy
surface. Its middle segment differs from the old complete source-anchor
slice. No old state is copied into it, and no state, Dirac field evolution,
renormalized polarization or nested `Psi/Phi` is supplied by the geometry
recorder itself.

The [first frozen reference](../results/summaries/nsc-atlas-dirac-radar-v1.json)
covers `S/ell=.1256873632901346`. Its conservative analytic lapse bound is
`.3851891297`; the sampled minimum on the declared 81 points is
`.7481617071`. The sampled cap radius reaches `.4240479593*ell`, with
the exterior bounded by the original wall radii. The two wall maps differ
from independent quadrature by `2.11e-14` or less. The largest matching
check is `2.03e-12` for the Dirac-potential derivative. Evaluation and the
[causal figure](../paper/figures/dirac-wall-radar.png) use 1.21875 process-CPU
seconds and 104 wall-ODE calls; no field evolution runs.

At the predeclared receiver the nested geometry gives
`omega_rec/omega_wall=1.1536257442`, versus `1.0806856390` for the smooth
cap continuation. `W_1*ell` is `1.7100955965` versus `1.7247878777`.
These use the same wall-radar clock protocol. Their difference is a
finite conditional background prediction, not a calculated quantum slip.

This makes inheritance a concrete input/output problem: a prepared
two-branch covariance `C_in` is transported by the characteristic map
`Vprop` to `C_mid=Vprop*C_in*Vprop^dagger` on the middle segment. Positivity
and the CAR bounds follow from the current-normalized isometry. Neither
the wall equation nor `T^*D=D` supplies the root preparation or permits
initial cross-correlations to be discarded.

One possible additional axiom would compress the line Hardy kernels
`integral_(omega<0) exp(-i*omega*(s-s')) d omega/(2*pi)` onto each incoming
leg for the occupied block, and use its complement for the empty block.
That compression is a positive contraction, rather than an endpoint-box
projector. Setting the cross-leg block to zero is a further preparation
assumption. No such axiom is selected by the current recorder. Its
four-dimensional short-distance behavior requires a uniform angular-mode
and cross-leg analysis after actual propagation; per-mode CAR positivity
alone does not establish a renormalized stress.

Null-coordinate changes also carry the canonical half-density: under
`z=phi(s)`, the line kernel transforms with
`sqrt(phi'(s)*phi'(s'))/(phi(s)-phi(s')-i0)`. Its leading coincident
`1/(s-s'-i0)` coefficient is unchanged. A coefficient computed without
these weights cannot decide the state question. The remaining test is
the full four-dimensional symbol and interface behavior.

[Gerard and Stoskopf](https://arxiv.org/abs/2108.11630) construct Dirac
Hadamard states using Cauchy-surface pseudodifferential calculus and
smooth bounded-geometry methods. [Fewster and Lang](https://arxiv.org/abs/1408.1645)
show that an unsoftened fermionic-projector prescription can fail on
ultrastatic slabs, and establish a softened construction there for positive
mass and compact spatial sections. That mass assumption is material for
this massless probe; the preparation check below treats it explicitly. Those
results supply methods to inspect; neither is asserted here as a proof
for this finite two-branch, piecewise-smooth wall preparation.

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
For this free CAR bilinear commutator the dependence on `C0` is affine:
quadratic covariance terms from separate Wick contractions cancel. The
exact two-mode CAR contraction and an explicit first-order Fock unitary
check both the commutator and the `-i[O,V]` insertion sign.

A local Hadamard parametrix is a subtraction kernel, not by itself a
positive state covariance. Setting its smooth remainder to zero does not
select a physical vacuum. Renormalization also needs conditions on the
finite gravitational couplings, including curvature-squared terms and
their interface completions; fixing the bulk trace alone does not fix
those transverse metric responses.

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

In particular, raw Fermi-coordinate metric coefficients at the observer
origin do not define that ratio: proper-clock and orthonormal-frame
normalization can set their perturbations to zero there. A finite clock,
light-travel or tidal response must be specified before connecting it to
potentials in a declared gauge. A short source-to-observer separation
also does not alone prove that the full causal past or quantum initial
covariance lies inside the `.02*ell` collar. That coverage must be checked
using the actual null characteristics and initial segment.

## Dirac stress and the actual gravitational metric vertices

The [stress interface](../src/recursive_horizons/nsc_atlas_dirac_stress.py)
now connects the same radial Dirac action to gravitational sources. Write
`xi=sqrt(a)*r*psi` on
`ds2=-N^2dt^2+a^2(dy+beta*dt)^2`, with angular eigenvalue `k`.
For a whole magnetic multiplet the degeneracy is `d=2*abs(k)`. For any
two radial modes `u,v`, define the matrix insertions

$$
q_A(u,v)=-\frac i2(u^\dagger A v_y-u_y^\dagger A v),\qquad
s_1(u,v)=u^\dagger\sigma_1v,\qquad f=d/(4\pi).
$$

The real canonical Hamiltonian density, per coordinate `dy` and including
the whole multiplet, is

$$
h=d\{N(q_\alpha/a+ks_1/r)-\beta q_I\}.
$$

Integrating the unsymmetrized first-order generator differs by the endpoint
`-i*d*[(N/a)*xi^dagger*alpha*xi-beta*xi^dagger*xi]/2`.
An open physical patch retains that port; it is not a reflecting vacuum
boundary. Varying the metric at fixed canonical half-density gives the
following on-shell orthonormal stress matrices:

$$
\rho=\frac f{ar^2}(q_\alpha/a+ks_1/r),\quad
p_r=\frac{f q_\alpha}{a^2r^2},\quad
p_t=\frac{fks_1}{2ar^3},\quad
j=T^{\hat0\hat1}=\frac{fq_I}{a^2r^2}.
$$

Here `j` is outward energy flux; the charge current instead uses
`f*xi^dagger*alpha*xi/(a*r^2)`. Their bilinears are different.
The finite-mode stress has `-rho+p_r+2*p_t=0`; the anomalous trace calculated
above is a separate renormalized quantum contribution.

The four metric derivatives are fixed by this same action:

| Variation | Canonical derivative | Physical stress relation |
| --- | --- | --- |
| `N` | `d*(q_alpha/a+k*s1/r)` | `4*pi*a*r^2*rho` |
| `a` | `-d*N*q_alpha/a^2` | `-4*pi*N*r^2*p_r` |
| `r` | `-d*N*k*s1/r^2` | `-8*pi*N*a*r*p_t` |
| `beta` | `-d*q_I` | `-4*pi*a^2*r^2*j` |

The nonzero bare second derivatives are
`h_Na=-d*q_alpha/a^2`, `h_Nr=-d*k*s1/r^2`,
`h_aa=2*d*N*q_alpha/a^3`, `h_rr=2*d*N*k*s1/r^3`,
with symmetric mixed entries. They supply local canonical metric contacts;
physical operator/measure variation and the continuum/interface counterterms
must still be included in the final polarization.
The first-order time-end variation is
`i*d/2*integral dy (xi^dagger*delta_xi-delta_xi^dagger*xi)` at the two
time faces. It also remains part of an integrated response/variation;
holding the initial preparation fixed does not erase the final face.

On the actual normal chart `a=1,beta=0`, the time equation and its spatial
derivative give the local conservation identities

$$
\partial_\tau(r^2\rho)+N^{-1}\partial_n(N^2r^2j)
+2rr_\tau p_t=0,
$$

$$
\partial_\tau(r^2j)+\partial_n(Nr^2p_r)
+r^2\rho N_n-2Nrr_n p_t=0.
$$

These hold for cross-mode matrices as well as diagonal densities. Contract
with a fixed radial covariance to obtain an expectation, or with its
tangent/difference for a response. Complete angular normalization, global
radial mode normalization and the admissibility of that covariance remain
explicit inputs.

At the pure-tension wall, `xi` is continuous and the Dirac time equation
requires `[xi_n]=k_wall*xi/2`. That real scalar shift cancels in every
antisymmetrized `q_A`. All four finite-mode stress components above are
therefore continuous in the common wall frame. This supplies compatible
mixed/normal matter flux to the unchanged material wall law; it is not a
claim that the anomalous vacuum surface stress vanishes.

Five [controls](../tests/test_nsc_atlas_dirac_stress.py) compare normal-time
Hilbert energy/flux on the actual exterior collar, build the full
four-dimensional covariant divergence from exact on-shell local jets,
differentiate the canonical action and its contacts independently, and
check transmission, orientation, physical scaling and matrix/sample domains.
The matrices are computed directly from the actual curved field; no flat
polarization or new source is substituted. No modes, state or coupled
gravitational response has been evolved by this component.

## Causal transfer, inherited correlations and nonlocal metric response

The [characteristic component](../src/recursive_horizons/nsc_atlas_dirac_scattering.py)
supplies the finite propagator on the actual radar diamond.
Stagger `p=eta_plus` on `(U_mid,V_node)` and `m=eta_minus` on
`(U_node,V_mid)`. For `P=sqrt(du)*p`, `M=sqrt(dv)*m` and
`gamma=W*sqrt(du*dv)/4`, the characteristic midpoint equations give

$$
\begin{pmatrix}P\\M\end{pmatrix}_{out}
=\frac{(1-\gamma^2)I-2i\gamma\sigma_1}{1+\gamma^2}
\begin{pmatrix}P\\M\end{pmatrix}_{in}.
$$

Real W makes every gate unitary. Their causal product U maps both complete
incoming legs to the future legs, without reflecting ends or duplicated
corner weights. For a supplied correlated parent, `C_out=U*C_in*U^dagger`.
Future `V=S` is exterior and future `U=S` is cap. The existing
[Gaussian owner](../src/recursive_horizons/nsc_boundary_state.py) independently
reconstructs the four terms in the inherited cap block. Conditional
inheritance is determined; a root quantum state is not selected.

For `delta ln N=delta ln a=lambda*f, delta r=0`, kinetic variations cancel
and `W -> W*exp(lambda*f)`. Per angular copy, the action insertion is
`-integral W*f*xi^dagger*sigma_1*xi dT dx`. Differentiating the same gates gives

$$
J_f=\sum_{\rm cells}\frac{2\gamma f}{1+\gamma^2}
U_{\rm before}^\dagger\sigma_{1,cell}U_{\rm before}
=iU^\dagger\partial_\lambda U.
$$

For disjoint B strictly before A in both null labels,
`R_AB=-i*[J_A,J_B]` is the derivative of the later insertion. The cross
contact vanishes for these supports. This is an integrated matter
metric-vertex response, not the gravitational/junction or observer solve.

A conditional occupation difference avoids choosing a baseline vacuum.
For a normalized smooth incoming mode f, set `P=I-|f><f|`,
`C_empty=P*C_parent*P` and `C_nu=C_empty+nu*|f><f|`, `0<=nu<=1`.
The complement's correlations remain, with the mode's cross terms removed
by this declared preparation. CAR positivity is exact. The response change
is `2*abs(k)*nu*f^dagger*R_AB*f` in units hbar for a whole magnetic multiplet.
Degeneracy is applied once after the single-copy commutator; k and -k
are separate multiplets. Given an admissible parent mapping f smoothly,
the difference is smooth and finite rank. Local state-independent
counterterms cancel; the root vacuum, full angular sum and general quantum
thin-wall completion remain uncomputed.

Seven [controls](../tests/test_nsc_atlas_dirac_scattering.py) pass: exact
plane waves and a C1 wall-crossing solution, flux, units/angular sign,
direct action and retarded derivatives, causality and spacelike
commutation, finite CAR differences and correlated inheritance.
No P001 evolution is invoked.

The first application completed both Dirac resolutions (17.515625 CPU
seconds at its last checkpoint), then stopped at the 24000-byte summary
gate. No result or full matrices were retained. Aggregate matrix hashes
replace repeated sector hashes for one completion execution of the
unchanged mathematics. Charging 30 CPU seconds for that first attempt
keeps both executions within the original 180-second allowance.

### Evaluated nested-versus-cap response

The frozen `7567c2a` [reference](../results/summaries/nsc-atlas-dirac-scattering-v1.json)
evaluates the unchanged quiet geometry at 64 and 128 characteristic cells
per direction. The metric source B is exterior at
`(U,V)/S=(.35,.55)`; the later cap insertion A is at `(.80,.70)`.
Each profile is a product of peak-one smooth bumps of half-width `.04*S`.
They are off the wall and null edges, and every B point precedes every A
point. Both geometries use the same clocks, profiles and canonical input.
The smooth cap is a geometric control, not another matched material-wall
solution.

The state-independent operator difference is nonzero:
`||R_nested-R_cap||=8.388462619e-10` in k=1 and `4.168642516e-7` in
k=8 at N=128. Their respective 64/128 norm spreads are
`5.23797e-12` and `1.91885e-9`. These are norms on one radial angular
copy, without multiplying a commutator by angular degeneracy.

For the separate conditional occupation difference, choose the exterior
incoming mode `g_minus(V)=bump(V,.6*S,.3*S)*exp(-12i*V/S)`,
`g_plus=0`, normalized in canonical null flux. It is smooth and compact
on that leg, without a frequency cutoff. For each positive-k multiplet,
let chi denote the retarded response change per unit occupation in units
hbar. Then:

| k | chi nested, N=128 | chi smooth cap, N=128 | Nesting change | Absolute 64/128 spread of change |
| --- | --- | --- | --- | --- |
| 1 | -2.387787593e-9 | -2.396832494e-9 | 9.044900946e-12 | 3.066849359e-13 |
| 4 | -2.373505609e-6 | -2.382039025e-6 | 8.533415365e-9 | 3.023552295e-10 |
| 8 | -6.911853302e-5 | -6.932408164e-5 | 2.055486263e-7 | 8.551697259e-9 |

The response magnitude decreases by about .3774% for k=1 and .2965% for
k=8. This is one specified prepared response component. The spread is a
correlated two-resolution numerical diagnostic, not a continuum bound or
observational error. The positive-k cases are different preparations;
they are not an angular truncation of a vacuum calculation.

The maximum flux-unitarity defect over both resolutions, both geometries
and all eight k values is `7.78e-15`; the assembled metric vertices are
Hermitian. Execution uses eight Numba workers and one BLAS thread,
10.234375 process-CPU seconds, and a 22783-byte summary. The first attempt's
conservative charge is 30 seconds. No matrices, new P001 grid, inherited
baseline vacuum response or Psi/Phi value is retained or claimed.

## Einstein, moving-wall and proper-clock response

The [gravitational interface](../src/recursive_horizons/nsc_atlas_dirac_gravity.py)
connects the same stress to the actual orbit metric
`ds2=-C*dU*dV+r^2*dOmega2`. Put `L=log C`, `kappa=8*pi*G` and keep the
regional constants `Lambda_minus=3/ell^2`, `Lambda_plus=0`.
The nonvacuum stress t excludes those constants. Direct four-dimensional
curvature contraction gives

$$
r_{UU}-L_Ur_U=-\kappa r t_{UU}/2,\qquad
r_{VV}-L_Vr_V=-\kappa r t_{VV}/2,
$$

$$
r_{UV}=-\frac{r_Ur_V}{r}-\frac C{4r}
+\frac{\Lambda Cr}{4}+\frac{\kappa r t_{UV}}2,
$$

$$
L_{UV}=\frac{2r_Ur_V}{r^2}+\frac C{2r^2}
-\kappa t_{UV}-\frac{\kappa C p_t}{2}.
$$

The existing orthonormal Dirac stress transforms as
`t_UU=C*(rho+p_r+2*j)/4`, `t_VV=C*(rho+p_r-2*j)/4` and
`t_UV=C*(rho-p_r)/4`. Angular multiplicity is already included there.
The component differentiates these equations in `s=delta r,q=delta L`;
it does not replace their constraints by a separate projected field trace.
For example, if

$$
E_U=r_{UU}-L_Ur_U+\kappa rt_{UU}/2,
$$

$$
Q_U=\partial_Vt_{UU}+\partial_Ut_{UV}
+2r_Vt_{UU}/r+(2r_U/r-L_U)t_{UV}+Cr_Up_t/r,
$$

then the mixed equations imply
`partial_V E_U=-(r_V/r)*E_U+kappa*r*Q_U/2`.
The V identity exchanges U,V. Conserved matter has Q=0; compatible
initial null constraints propagate. A source that fails this identity
cannot be repaired by assigning arbitrary independent metric data.

The wall may move to `x=w(T)` in the conformal-coordinate gauge.
At the unperturbed wall, `C=1,r=R,r_T=v,r_x=beta_i`,
`L_T=0,L_x=2*beta_i_R`. The two common surface perturbations are
`Z=s_i+beta_i*w` and `Q=q_i+Lx_i*w`. Define time and angular extrinsic
eigenvalues with the normal directed from cap to exterior. Their variations are

$$
\delta K^T{}_T=\ddot w+\frac{q_x}{2}
+\frac{wL_{xx}}2-\frac{QL_x}{4},
$$

$$
\delta K^\theta{}_\theta=
\frac{s_x+wr_{xx}+\dot wv-Q\beta/2}{R}-\frac{\beta Z}{R^2}.
$$

Continuity of Z,Q and zero jumps of both variations are the fixed-tension
junction laws. The signed Gauss identity is
`kappa*[t_nn]-[Lambda]=-2*k_wall*average(K^T_T+2*K^theta_theta)`,
where t excludes the regional vacuum terms. Including those terms in a
total T gives the equivalent `kappa*[T_nn]=-2*k_wall*average(K)`.
Its perturbation is a compatibility condition from the bulk and Israel
equations, not another independent wall datum. The finite Dirac stress
transmits normally; its zero normal-pressure jump leaves zero variation
of that mean curvature. The physical proper-clock rate changes by Q/2.

To read the same radar protocol used above, transform to perturbed
proper-wall labels. If `delta_tau(T)` is the integrated wall-clock shift,
`Uhat=U+A(U),Vhat=V+B(V)` with `A=delta_tau+w`, `B=delta_tau-w`.
At fixed physical radar labels,

$$
\delta L_{\rm radar}=q-AL_U-BL_V-A'-B',\qquad
\delta r_{\rm radar}=s-Ar_U-Br_V.
$$

A radial phase emitted by the unit wall clock gives
`delta ln(omega_rec/omega_wall)=-delta L_radar/2` for the fixed
radar-position receiver. Both quantities vanish under a pure residual null
reparametrization. This avoids interpreting an arbitrary coordinate lapse
as an observable. A numerical clock correction still requires the actual
constrained gravitational solution.

For the retarded quantum calculation, conservation needs one more connection.
Varying a nonzero mean stress under a metric source obeys
`nabla(delta T)+delta(nabla)*T_mean=0`, not `nabla(delta T)=0`.
Let h_q solve the mean gravitational backreaction and h_B the classical
metric source. At the mixed order the Einstein source includes the
second variation `-E_second[h_B,h_q]` and the physical external source's
variation, together with the complete quantum stress response. Feeding
only the metric commutator into the vacuum linear equations would omit
those Ward/tadpole terms. The quantum state, source model and renormalized
local terms therefore remain explicit. A finite occupation difference
cancels the unchosen baseline, but still needs its own mean backreaction.

Six [controls](../tests/test_nsc_atlas_dirac_gravity.py) pass: full
four-dimensional curvature with independent metric jets, exact nonvacuum
constraint transport, linearization on both actual regions, direct
moving-normal acceleration, the signed wall Gauss law and residual gauge
removal in the physical clock. They use 4.390625 CPU seconds across the
initial suite and two focused follow-ups, including one scalar-versus-array
fixture correction. No gravitational or P001 evolution is run.

## Mean stress and a direct spherical metric reconstruction

The [finite mean-source interface](../src/recursive_horizons/nsc_atlas_dirac_mean.py)
propagates the same declared exterior mode with the shared Cayley gate.
Its compact envelope and phase are unchanged, with continuum rather than
mesh-dependent normalization. For P=eta_plus, M=eta_minus, analytic
characteristic differentiation gives

$$
\partial_V(P_U)=-iW_UM/2-W^2P/4,\qquad
\partial_U(M_V)=-iW_VP/2-W^2M/4.
$$

The bounded tensor-Hermite representation supplies the field and first
derivatives to the existing stress owner. No stored waveform is differentiated
and no extrapolation beyond the represented cell-center box is permitted.
The analytic quiet region V<=.3S is zero. Five field controls pass, including
wall-crossing derivatives, flux weights, continuum normalization, causal zero
and one-sided stress transmission. A uniform subcell derivative-error envelope
decreases `.00292686 -> .000741822 -> .000186756` at 8/16/32 cells.
The earlier ratio at one changing fractional cell position was not a uniform
order estimator; the envelope supplies the relevant check.

For the supplied conserved mean stress, the
[metric reconstruction](../src/recursive_horizons/nsc_atlas_dirac_backreaction.py)
uses spherical Einstein constraints directly. Write
`m_total=M_i+Lambda_i*r^3/6+m_i`, subtracting the cap vacuum mass at the
actual areal radius. To first order in the supplied stress,

$$
(m_U,m_V)=\frac{\kappa r^2}{C}
(t_{UV}r_U-t_{UU}r_V,\ t_{UV}r_V-t_{VV}r_U).
$$

This one-form is closed for a conserved source on the vacuum background.
Agreement of two integration paths is a necessary numerical check; it
does not by itself establish all stress-conservation equations.
At the wall,
`mu_i_dot=-kappa*R^2*(beta_i*j+v*p_r)/2`.
The existing **geometric** `wall_response_coefficients` supplies
rho_dot and delta beta from those transported exterior/cap masses.
No scalar mass functional is reused for the Dirac source. The scalar
certificates retain their original source and finite-domain claims;
this Dirac application has the independent stress, mass and metric checks here.

Use outgoing EF coordinates with
`ds2=-a^2*F*dU^2-2*a*dU*dr+r^2*dOmega2`. The background
`a0=C/(2*r_V)=1/(beta_i+v)` stays finite through beta_plus=0.
For psi=delta ln a at fixed U,r,

$$
\psi_V=\frac{\kappa r t_{VV}}{2r_V},\qquad
\psi_{\rm wall}=-\frac{\delta\beta_i+\dot\rho}{\beta_i+v}.
$$

Mass and psi follow by quadrature from the wall along constant U, on the
appropriate side. A nonzero homogeneous corner mass must be retained.
Both corner masses vanish for the prepared Dirac difference, whose input
vanishes there; the charged-family control correctly has nonzero corner mass.

The physical radar radius s=delta r obeys the linearized ingoing null ray
equation, with its wall radius rho as boundary data:

$$
s_U=bs+g,\quad b=-a_0f'/2,\quad
g=a_0m/r-a_0f\psi/2.
$$

Differentiating this equation analytically supplies z=s_V and hence
`q=delta ln C=psi+z/r_V`. The wall value q=0 follows from proper-time
normalization. The specified radial light protocol gives
`delta ln(omega_rec/omega_wall)=-q/2`. Thus no independent bulk metric
PDE grid is introduced for this first-order supplied source.

Five [metric controls](../tests/test_nsc_atlas_dirac_backreaction.py) pass.
They recover the nonvacuum mass one-form, the full nonlinear equal-charge
Einstein-Maxwell wall/ray derivative on both sides, regularity at beta_plus=0,
the zero-source limit, and a smooth-cap worldtube comparison. A separate
exact static-fluid family has nonzero T_rr and verifies the radial lapse
and resulting clock derivative, which the Maxwell family alone cannot test.
These manufactured sources are controls and are never inserted into the
physical Dirac calculation. The combined development controls use
23.03125 process-CPU seconds, including the bounded interpolation diagnosis.

The smooth cap contains no membrane. For the coming clock comparison its
observer worldtube is explicitly driven along the computed nested wall's
areal-radius history. It isolates geometry at a common observer motion;
it is not another independently evolving wall model or an observational fit.
The root vacuum, baseline polarization and full mixed retarded
source/tadpole response remain outside this mean calculation.

Physical conversion restores `8*pi*hbar*G/(c^3*ell^2)`, using the
independently supplied nominal ell=10 km and the pinned runtime's
[CODATA 2022 constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt).
G's measured uncertainty is separate from numerical resolution spread.
The incoming clock-generator energy is `2*hbar*c*12/(S*ell)` for this
complete k=1 multiplet; the smooth compact mode is not frequency limited.
No sky observation sets this preparation, length, or numerical result.

### Evaluated prepared mean-clock correction

The frozen `be1a11e` [reference](../results/summaries/nsc-atlas-dirac-mean-clock-v1.json)
connects the field, analytic derivatives, mean stress, both transported
masses, material-wall displacement, bulk lapse and physical radar clock.
It uses the previous proper-clock labels U=tau(.475), V=tau(.525);
the radii .475 and .525 label the quiet reference, not a constraint fixing
the perturbed wall radii. Put `eta_q=8*pi*hbar*G/(c^3*ell^2)` and define
`C_clock=delta ln(omega_rec/omega_wall)/(eta_q*nu)` for the declared
empty-to-occupied difference. At quadrature order 12:

| Geometry | C_clock at 128 cells | Absolute 64/128 spread | Normalized mass-path gap, 64 to 128 |
| --- | --- | --- | --- |
| Nested cap-wall-exterior | -.4177997648 | 1.117879e-5 | 1.150808e-6 to 2.924683e-7 |
| Smooth cap, shared observer motion | -.3992877906 | 9.947292e-6 | 3.252355e-6 to 5.845798e-8 |

The coefficient difference is `-.01851197417`, with a correlated 64/128
spread `1.231497e-6`. Nesting changes the magnitude of this logarithmic
mean-clock response by about 4.63625% against the specified cap control.
The nested mass-path discrepancy contracts by about a factor of four.
Quadrature orders 8/12 differ by at most `5.33e-15` in the clock coefficient.
The finest incoming flux is `1.00000000000964` with propagation defect
below `7e-16`. No validation condition fails.

For the independently supplied ell=10 km and nu=1, the fractional
frequency shift is `-2.743012753e-77`, versus `-2.621474673e-77` in the
smooth control. The nominal nesting difference is `-1.215380800e-78`.
The nested resolution spread converts to `7.33930e-82`; G's relative
standard uncertainty is `2.24743e-5` separately. These remain numerical
resolution diagnostics and a constant uncertainty, not a rigorous
continuum bound or observational confidence interval.

The reference uses an eight-process pool with one BLAS/Numba thread per
worker: 11.3722 wall seconds, 40.265625 observed process-CPU seconds and
48.265625 charged CPU seconds including overhead. Its summary is 13053
bytes; no field or metric arrays are retained. SciPy emitted a scalar-divide
warning in its local error estimator; all accepted solutions and the saved
results are finite and all prescribed mass/clock checks pass.

This is a calculated mean response for an explicit finite prepared-state
difference, conditional on a compatible parent. It supplies the h_q
prerequisite for the later mixed retarded source/tadpole response. It does
not compute the baseline vacuum, full Psi/Phi or an effective cosmological
Lambda, and it does not close the unconverged finite-amplitude scalar P001.

## Conserved compact-source mass tangent

The next physical source is a constant central excess mass b, present
before the finite window, whose body or inner horizon stays below the
diamond's innermost cap radius. Its exterior cap is
`f_minus=1-2*b/r-r^2`; the parent vacuum is `f_plus=1-2*M_parent/r`.
The body interior and a regular formation history are not supplied by this
finite annular construction. H=1, tension k=2, initial R0=.45 and initial
physical proper-wall speed are fixed. The signed Israel equations require

```text
beta_minus0(b) = sqrt(beta_minus0(0)^2-2*b/R0)
beta_plus0(b) = beta_minus0(b)-k*R0
M_parent(b) = R0*(1+v0^2-beta_plus0(b)^2)/2
dM_parent/db|0 = beta_plus0(0)/beta_minus0(0) = 2755/8587.
```

Holding both regional masses fixed would not preserve that initial speed.
On the annulus the nonvacuum bulk stress vanishes; the two mass parameters
are constants. The existing mean GR engine therefore supplies this
linearization with start=0, corner masses `(2755/8587,1)` and zero
initial displacement. Its initial displacement rate also vanishes.
There is no source switching event inside the window, so this is a
parametric susceptibility, not a retarded switched-mass response.

The [compact-source interface](../src/recursive_horizons/nsc_atlas_dirac_compact_source.py)
adds anchored optical sensitivities to that same wall solution:

```text
u_B' = -(delta_beta+rho_dot)/(beta+v)^2
v_B' = -(delta_beta-rho_dot)/(beta-v)^2
u_B(0)=v_B(0)=0
s_B = f*((v_B(V)-u_B(U))/2-delta_M*K(r))
q_B = (f'*s_B-2*delta_M/r)/f
      -(delta_beta(U)+rho_dot(U))/(beta(U)+v(U))
      -(delta_beta(V)-rho_dot(V))/(beta(V)-v(V))
delta_W = W*(q_B/2-s_B/r).
```

Here `delta_beta=beta_R*rho+(mu-nu)/(k*R^2)` includes displacement.
`K=partial_M(r_star)` is anchored to K(R0)=0. Its primitives are
`2*log(r)-log(1-r^2)+1/(1-r^2)` in the cap and
`2*log(.8-r)+1.6/(.8-r)` outside. Both have derivative `2/(r*f^2)`.
The interface supplies analytic first derivatives of s_B, q_B and delta_W,
including the optical rates, for a later differentiated Dirac calculation.
It never divides by beta_plus, which crosses zero regularly.

| Compute admission | Fixed scope |
| --- | --- |
| Borrowed and missing connection | Existing signed Israel, Misner-Sharp/EF and Dirac operators; connect a conserved physical mass source to the actual propagator, rather than an arbitrary compact metric strain. |
| Fixed assumptions and domain | The same quiet finite proper-wall diamond; infinitesimal b at zero, fixed H/k/R0/v0, compatible inner source below the domain, parent mass fixed by the junction. No baseline quantum state is selected. |
| Budget | At most 120 process-CPU seconds shared by development controls and diagnosis. Wall/optical ODEs and independent exact-family ray controls only; no P001 evolution, new physical field run, source scan or stored arrays. |
| Independent check | Full nonlinear SdS--wall--Schwarzschild family differentiated by complex step, independently integrated radial rays; wall clock/radius/velocity, cap and exterior, beta_plus=0, first potential derivatives and mass constraints. |
| Stop or redirect | A failed identity or exhausted budget stops this bridge for diagnosis. A validated tangent unlocks a separately admitted mixed prepared-state response; it does not itself calculate that response or full Psi/Phi. |

A radar clock/radius pair `(q_B,s_B)` is not the Newtonian pair `(Phi,Psi)`.
At the wall q_B=0 while s_B generally differs from zero. A full readout
requires a stated scalar-potential decomposition and its radial zero-mode
boundary data, or an explicitly named operational curvature observable.
Differentiating the later mean response must vary the background Einstein
operator, stress, source preparation, wall and readout together. Changing
only W on the fixed quiet GR solver would omit the mixed Einstein term.

Five [independent controls](../tests/test_nsc_atlas_dirac_compact_source.py)
pass in 3.84375 process-CPU seconds. They compare the optical metric,
radius gradients, lapse gradients and first Dirac-potential derivatives
with a full nonlinear SdS--wall--Schwarzschild family differentiated by
complex step. They also recover the existing EF metric reconstruction,
the initial physical speed, one-sided wall transmission and the regular
beta_plus=0 crossing. The anchored mass primitives agree with independent
quadrature. This validates the source tangent; no physical mixed quantum
response has been evaluated in these controls, and no P001 grids were run.

The next extension uses the exact real mass family on `0<=b<=2e-4`,
holding the original proper-clock duration fixed. Its computational
annulus is `.3<r<.65`; the inner body must stay below that annulus.
`CompactMassChart` integrates the actual wall and optical maps and
inverts its anchored tortoise by fixed 12-point Gaussian quadrature.
The mean GR engine reads this family's wall geometry, cap mass term
and signed beta derivatives. It no longer freezes those coefficients
when b changes. The zero-mass path remains the existing quiet engine.
The same 120-CPU-second development admission also covers the quiet
limit, exact finite-b rays/junctions, and an independent charged-family
check of mean metric reconstruction on a nonzero-b background. No
physical mixed-response execution is admitted by these controls alone.

The expanded source/family and existing GR controls pass: 13 tests in
15.3125 CPU seconds, bringing this source bridge's total development
charge to 19.15625 CPU seconds. This includes the finite-b and nonzero-b
charged-family checks above. No additional field resolution was run.

### Local curvature readout and its relation to scalar potentials

An independent observer readout uses actual Riemann components in the
fixed proper-wall radar observer's orthonormal frame. Define
`A=R_0r0r-R_0theta0theta` and
`B=R_rtheta_rtheta-R_thetaphi_thetaphi`. In the convention with
vacuum `R_0r0r=-2*M/r^3-Lambda/3`, direct Ricci/Weyl decomposition gives

```text
w = m_excess/r^3-kappa_g*(rho-p_r+p_t)/6
A = -3*w-kappa_g*(p_r-p_t)/2
B = -3*w+kappa_g*(p_r-p_t)/2.
```

The fixed vacuum term cancels. w and rho-p_r are invariant under radial
boosts; A/B and p_r-p_t depend on the specified observer. Flux is retained
in the stress transformation even though it does not appear explicitly
in this frame's two sectional contrasts. At fixed physical radar labels,
the quantum mean correction retains
`w_q=m_q/r^3-3*b*s_q/r^4-kappa_g*(rho_q-p_rq+p_tq)/6`.
For b>0 the classical cap has A0=B0=-3*b/r^3, so the leading occupation
correction to B/A is `kappa_g*(p_rq-p_tq)/A0`.
At b=0 there is no classical tidal denominator and that ratio is undefined.

For comparison, in a specified flat-FLRW Newtonian metric the first-order
contrasts in its *comoving* frame are `A=(Phi_chichi-Phi_chi/chi)/a^2`
and `B=(Psi_chichi-Psi_chi/chi)/a^2`. Their ratio compares these tidal
potential combinations; it is not the potentials' full ratio Psi/Phi.
The operator has radial zero modes `F(eta)+G(eta)*chi^2`, which require
boundary/background data beyond a local anisotropic-stress equation.
In particular, vanishing anisotropic stress alone does not remove those
zero modes on a finite annulus. The exterior is below its Schwarzschild
horizon and is not assigned this cap FLRW decomposition.

The [curvature interface](../src/recursive_horizons/nsc_atlas_dirac_tidal.py)
uses the same physical mass/stress. Its development admission is 60 shared
process-CPU seconds: direct four-dimensional Christoffel/Riemann contraction,
a changed radial observer with nonzero flux, the de Sitter Newtonian-metric
limit and an exact charged-family mean/ratio check. These are controls,
not additional physical field evolutions or a baseline vacuum calculation.

Three [readout controls](../tests/test_nsc_atlas_dirac_tidal.py) pass in
1.8125 CPU seconds. They independently contract the full four-dimensional
Riemann tensor with arbitrary radial flux and a changed observer, recover
the de Sitter Newtonian-metric tidal combinations, and differentiate an
exact charged source including its radius term. No potential zero mode
is fixed or inferred by these local curvature controls.

### Prepared compact-source response admission

The [bounded producer](../scripts/measure_nsc_atlas_dirac_mass_response.py)
will now evaluate the finite prepared mean *on* the exact mass family.
This differentiates the background Einstein reconstruction as well as
the Dirac propagation, stress, wall and observer. It retains the mixed
Einstein contribution that changing W on a frozen GR background would omit.

| Admission | Fixed decision |
| --- | --- |
| State and source | Same continuum unit-flux exterior mode, carrier and k=1 multiplet as the mean-clock reference. The occupation difference and initial quantum corner masses stay fixed. The compact classical mass is present before the window. A compatible parent is assumed; no baseline vacuum or angular sum is chosen. |
| Domain and observer | Proper-clock diamond with the original S, fixed H/k/R0/v0, .3<r<.65 and the body's interior below it. The cap receiver is U=.8S,V=.7S, fixed radar position; an inward phase -V is emitted at the earlier wall time V. The smooth SdS continuation shares classical and quantum areal-worldtube motion. |
| Parameters and numerics | Exactly b=0,5e-5,1e-4, with field cells 64/128 and piecewise quadrature orders 8/12. These give a fixed one-sided source derivative with h/h2 comparison and first-order extrapolation. They are not a source search or a P001 grid campaign. |
| Outputs | Mean clock/radius/mass, the clock's parametric mass derivative, and local tidal contrasts A/B. The finite-b state difference of B/A is recorded only for b>0. It is not full Psi/Phi. |
| Compute and storage | One physical reference creation, 480 process-CPU seconds including worker/import costs and eight seconds overhead; setup and each clock have 30-second CPU guards. Up to eight workers, one BLAS/Numba thread each. Summary at most 32 KiB, no saved arrays. |
| Checks and stop | Input norm within 1e-5, charge-flux defect below 2e-12, normalized mass-path gap below 1e-3 and contracting at 64/128 when above 1e-7, wall-clock defect below 2e-9; 8/12 clock spread below 2e-6 times max(1,clock magnitude). Source-derivative step and field-resolution spreads must each be below .02 times max(1,derivative magnitude). Failure stops for a diagnosis, without another mass choice or field grid. |
| Physical scale and uncertainty | Independently supplied ell=10 km; b=1e-4 means a geometric mass length of 1 m. Restore epsilon_q=8*pi*hbar*G/(c^3*ell^2); curvature contrasts have units ell^-2 and physical acceleration gradients add c^2. Grid, quadrature and one-sided derivative spreads are separate diagnostics, not rigorous remainders or observational confidence intervals. |

The likely failures are omission of b-dependence in the GR operator,
preparation or worldtube; unresolved differencing of the source response;
and mistaking a tidal ratio or state difference for a full vacuum potential
ratio. The exact nonzero-mass controls, fixed shared preparation, recorded
derivative spreads and explicitly defined local Riemann readout address
those failures. No parameter is selected using an observed target value.

### First mass-response result and precision diagnosis

The frozen `76a02c9`
[reference](../results/summaries/nsc-atlas-dirac-mass-response-v1.json)
finishes all 24 readouts in 62.171875 charged CPU seconds and a 26247-byte
summary. Its status is **unresolved**: the clock mass derivative fails
three prescribed step/resolution checks. The extrapolated nested
derivative moves from .01264 to .11540 between field resolutions; it
is not an accepted physical susceptibility. No tolerance or source
has been changed to turn this record into a pass.

The local tidal ratio uses the propagated stress directly and does not
depend on differencing the adaptive mean clock. At b=1e-4 its leading
occupation coefficients are -540762.8710915 (nested) and -540765.5712499
(smooth cap), with individual 64/128 spreads about 14.59. The nesting
difference is +2.700158436 with a correlated spread .001185480, about
five parts per million of the prepared correction. At the nominal
ell=10 km these are ratio changes -3.55031184e-71 and -3.55032957e-71.
These are local radar-frame state-difference diagnostics with numerical
spreads, not full Psi/Phi or observational uncertainties. The containing
mass-response record remains unresolved because of the clock derivative.

Before any further physical reference, a bounded
[same-field diagnosis](../scripts/diagnose_nsc_atlas_dirac_source_precision.py)
compares only N=128 at b=0 and b=5e-5. It regenerates each existing field
once, then reuses it for the original solver tolerances, a tighter wall
solve, and a tighter wall plus ray solve (2e-11 relative, 2e-13 absolute).
Independent wall-mass Gauss integration checks the ODE masses. The
diagnosis admits 120 shared process-CPU seconds, with two 50-second
worker guards, a four-second overhead charge and a summary no larger
than 12 KiB. No new field resolution, source choice or P001 grid is
admitted. A clock shift comparable to the small mass-difference signal
identifies an integration precision floor; if it does not, the source
derivative remains unexplained and this branch stops for a new diagnosis.

The frozen `eaaceb9`
[diagnosis](../results/summaries/nsc-atlas-dirac-source-precision-v1.json)
settles this precision question in 41.4375 charged CPU seconds. At b=0,
direct wall-mass gaps fall from (-1.176e-6,+9.460e-7) to
(-1.622e-8,+1.246e-8); at b=5e-5 they fall from
(+1.604e-6,-7.204e-7) to (-1.029e-8,+1.013e-8).
Tightening the ray alone adds no change after tightening the wall.
The same-field one-sided clock derivative moves from .07300170 to
-.002188022. The original subtraction was dominated by its wall ODE
precision floor, rather than by the ray or a different mass source.

One v2 reference is therefore admitted with the *same* b values,
64/128 fields, 8/12 quadrature, physical data, receiver and acceptance
criteria, changing only wall tolerances to relative 2e-11 and absolute
2e-13. The same 480-process-CPU-second and 32-KiB limits apply; the
per-setup and clock limits remain 30 CPU seconds. The original ray
tolerance is retained because the diagnosis found no ray shift.
Another failed check stops v2 for diagnosis; no additional field grid,
mass choice, P001 evolution or tolerance ladder is admitted.

For clarity, a cap scalar-potential map can be stated without identifying
it with the local tidal ratio. On the b=0 de Sitter background, use
`eta=-(exp(-u)+exp(-v))/2`, `chi=(exp(-u)-exp(-v))/2` and `a=-1/eta`.
Pulling the radar perturbation into this fixed background chart gives

```text
phi=q/2, shift=0
E_chichi-E_chi/chi=q/2-s/r
psi=E_chi/chi-s/r
Phi=phi-Hcal*E_eta-E_etaeta
Psi=psi+Hcal*E_eta.
```

The scalar decomposition permits `E -> E+F(eta)*chi^2/2+G(eta)`.
If E and E_chi are anchored to zero on a moving wall, then E_eta=0
there but `E_etaeta=chi_wall_prime^2*E_chichi`. Proper wall clocks give
q=0, hence `Phi_wall=chi_wall_prime^2*s/r` and `Psi_wall=-s/r`.
They generally differ even for a vacuum mass perturbation with A=B.
This explicit boundary convention is distinct from discarding the
radial zero modes; a full potential ratio must state which one it uses.
No numerical full-Psi/Phi result is inferred from this unevaluated map.

### Wall integrals after the v2 setup limit

The `bc4d39f`
[v2 execution](../results/summaries/nsc-atlas-dirac-mass-response-v2.json)
hits a setup's 30-CPU-second guard. Five preparation pairs finish, but
the N=128,b=5e-5 pair does not; no complete observation is collected.
The 133.90625 completed-worker CPU seconds exclude interrupted/uncollected
worker work and are not the total charge. This is a numerical execution
exit, not an additional physical exclusion. Its record is preserved.
The allowance is not increased, and no additional field grid is run.

The wall equations have a direct triangular solution. Integrate the
two specified mass-flux laws for mu and nu, then use
`(log v)_dot=acc/v` to obtain

```text
rho(t)/v(t) = rho(t0)/v(t0)
  + integral_t0^t [beta_plus*(mu-nu)/(k*R^2)+mu/R]/v^2 dt.
```

The [wall primitive](../src/recursive_horizons/nsc_atlas_dirac_wall_integral.py)
uses piecewise Gauss quadrature on existing field-cell breakpoints.
Masses at an outer quadrature node use direct partial-cell quadrature;
all inner nodes are evaluated together. It stores cumulative masses and
rho/v at edges, then directly integrates only the partial cell for an
arbitrary query. It differentiates neither a stored waveform nor a
mass interpolant and adds no physical boundary condition. The smooth
control only integrates its cap mass; its radius is still supplied by
the same nested observer worldtube.

This reformulation admits 120 development process-CPU seconds for
exact nonconstant-speed/corner integrating-factor controls and the
existing independent charged/fluid GR controls on both wall methods,
including a nonzero compact mass. The physical field data and initial
conditions remain fixed. A failed identity stops this formulation for
diagnosis before a new physical reference is admitted.

The wall-integral and GR controls pass: 21 tests in 70.84375 CPU seconds.
Seven selected callback/cache and process-pickling controls then pass in
31.65625 CPU seconds, for a shared development charge of 102.5 CPU seconds.
The exact clock inversions are shared through the existing field-point
cache; the calculation is unchanged. Bound method callbacks preserve
the mean metric across the spawned worker handoff.

One v3 physical reference is admitted under the same 480-CPU-second,
30-second setup/clock and 32-KiB limits. It replaces the wall ODE with
the triangular primitive, using wall order 8 together with bulk order 8
and wall order 12 together with bulk order 12. Both orders reuse the
same propagated field. The three masses, 64/128 field resolutions,
state, observer, physical scale and acceptance conditions stay fixed.
An eight-process pool still carries the independent cases; no new
P001 or Dirac field resolution is admitted. Another failed requirement
stops this reference for a recorded diagnosis rather than another grid.

### Evaluated compact-source response

The frozen `358bf25`
[v3 reference](../results/summaries/nsc-atlas-dirac-mass-response-v3.json)
passes every prescribed validation. Its 26964-byte summary records
6.97495 calculation wall seconds, 42.71875 observed process-CPU seconds
and 50.71875 charged CPU seconds, using eight independent workers.
No field arrays are retained, and no P001 grid or time step is run.

Let `C_clock(b)=delta ln(omega_rec/omega_wall)/(epsilon_q*n_occ)` and
differentiate in the dimensionless central mass b/ell. The fixed
one-sided h/h2 extrapolation at zero gives:

| Geometry and observer | Clock mass derivative, 128 cells | Absolute 64/128 spread |
| --- | --- | --- |
| Nested cap-wall-exterior | -.002329215549 | 2.002398e-7 |
| Smooth SdS cap, shared observer motion | -.002329519102 | 2.028155e-7 |

The nesting difference is +3.035527385e-7, with a correlated grid spread
2.575717417e-9. The correlated h/h2 step indicator is 5.950795412e-10;
the individual derivative-step indicators are about 1.655e-6. These are
numerical diagnostics, not rigorous error bars. The full mass-dependent
Einstein reconstruction, field stress, wall and observer vary together;
this is a prepared-state parametric mixed response to the conserved
source family, rather than a matter-only metric insertion or a switched
retarded source.

At b=1e-4 the independently propagated local tidal ratio keeps the
coefficients -540762.8710915 (nested) and -540765.5712499 (smooth cap).
The nesting difference +2.700158436 has correlated spread .001185480,
or about 4.9932 parts per million of the prepared correction. At the
supplied ell=10 km and unit occupation, the occupied-minus-empty ratio
changes are -3.550311840e-71 and -3.550329568e-71. Their difference is
+1.772755671e-76, with a grid-spread conversion 7.783122016e-80.
The empty state's baseline vacuum tide remains uncomputed; these are
not predictions of its absolute B/A or of the full potentials Psi/Phi.

The maximum normalized nested mass-path gap contracts from 8.075638e-7
to 1.930186e-8; the smooth control contracts from 8.077223e-7 to
1.934896e-8. Joint wall/bulk quadrature orders 8/12 differ by at most
4.45e-15 in the clock coefficient. The failed v1 derivative and v2
execution stay preserved; v3 changes only the numerical solution of
the same admitted wall equations. This finite prepared-source branch
is complete at the stated diagnostic precision. The next full-vacuum
calculation requires a specified joined state, ultraviolet angular
completion, renormalization and a chosen potential/readout convention.

### Vacuum tidal reduction and the remaining preparation rule

For the *local tidal channel*, a useful first-order cancellation reduces
the next vacuum calculation. Compare states on the same cap SdS
background and the same zeroth-order radar observer, with b>0. From
the exact sectional identities above,

```text
B-A = kappa_g*pi,  pi=p_r-p_t
(B/A)-1 = kappa_g*pi/A0 + O(hbar^2),  A0=-3*b/r^3.
```

The mean mass/radius contribution w_q enters both contrasts equally.
It cancels from their difference and enters the denominator only at
second quantum order. A first-order observer or evaluation-point shift
acting on the zeroth-order vacuum nonvacuum stress also contributes
zero to pi. This does not discard those terms from the separate
contrasts, clock response, or a full scalar-potential reconstruction.

Use a finite detector with a common zeroth-order test density f whose
support lies strictly inside the cap, separated from the wall. Define
its readout as the ratio of averaged contrasts,
`eta_f=(integral f*B dV)/(integral f*A dV)`, with a nonzero classical
denominator. Its first-order nesting difference is

```text
Delta eta_f = kappa_g * integral f*(pi_nested-pi_smooth) dV
             / integral f*A0 dV.
```

Under a common local renormalization prescription, bulk counterterm
tensors agree because the two cap metric jets agree on the detector.
Gravitational interface counterterms have support outside that detector;
their first-order metric effect cancels from this particular ratio as
above. This statement assumes perturbative Einstein response and fixed
renormalized classical source/couplings. It does not cover a detector
touching the wall, different cap backgrounds, direct spinor/interface
interactions, or nonperturbative higher-derivative gravity. It does not
identify two states or remove their vacuum stress.

For two states locally Hadamard on that common smooth cap neighborhood,
the difference of their two-point kernels is smooth there. With the
point-split stress operator D_pi from the same Dirac action,
`Delta pi(X)=lim_{Y->X} D_pi*(G_nested-G_smooth)(X,Y)` then needs no
additional local subtraction. Proving that the *propagated* state
families have this property remains an actual requirement; per-mode
CAR positivity or the conformal wall first-jet cancellation does not
prove the needed combined angular/frequency limit.

At a finite field/angular cutoff, the existing interfaces would evaluate
this through `Tr(C_nested*P_nested[f]-C_smooth*P_smooth[f])`, where
P[f] is the detector-integrated radial-minus-angular stress matrix
pulled back by the actual characteristic transfer. Full magnetic
degeneracy is included once in that matrix. A rotationally invariant
state permits the sum over both signed nonzero angular sectors; no
parity symmetry or cross-leg independence is supplied automatically.
The controlled cutoff limit, not one finite matrix, is the vacuum result.

The fixed-b quantity includes environment anisotropy already present
without the compact source. To call it *source-induced polarization*,
use a declared b=0 subtraction of `Delta pi_b-Delta pi_0`, or its
mass-family derivative, with the same preparation rule varied along
both families. Dividing a pre-existing wall anisotropy by b would not
be that source response.

The remaining physical datum is therefore an explicit common preparation
rule assigning admissible parent covariances to both geometries and to
their mass families, including their cross-leg correlations and local
short-distance structure. The parent metrics differ outside the cap;
copying one literal covariance matrix between them need not preserve
their distinct ultraviolet symbols. An instantaneous intrinsic negative
projector and a factorized Hardy prescription are unproved choices here,
not excluded by a blanket statement about every nonzero extrinsic
curvature. No such root rule is selected by the current calculation.
The finite occupation calculation itself shows why the inherited
operator relation alone does not fix it: changing an allowed covariance
at the same operator changes the computed response. This local tidal
reduction narrows the next calculation; full vacuum Psi/Phi remains open.

## A candidate parent preparation from clock energy

The next missing object is a state rule, not another prepared-packet
refinement. Three routes have different costs and assumptions:

| Route | What it supplies | Unresolved connection here |
| --- | --- | --- |
| Smooth ultrastatic deformation | A reference ground covariance transported by an explicitly supplied preparation geometry. | A common parent-time extension and control of the material-wall limit are not yet supplied. |
| Softened fermionic signature projector | A state selected by an averaged spinor bilinear. | The cited positive-mass theorem cannot be substituted for a massless calculation; its massless two-frequency block is tested below. |
| Minimum averaged clock energy | A covariance selected from the actual field's energy quadratic form in the wall clock. | This is a new state axiom. Its continuum/local short-distance admissibility remains to be established. |

Use the third route as a finite candidate to investigate, keeping the
deformation construction available for a parent-state completion. It uses
the physical clock and existing Dirac action, fixes both incoming branches
together, and assigns a geometry-dependent covariance under one common
rule. It does not assert that the operator identity determines that rule.

Let b(s) be the existing peak-one smooth bump centered at S/2 with
half-width S/2. For the actual propagated incoming modes, define

```text
A_k = (1/S) integral_0^S integral_0^S b(U)*b(V)*h_k(U,V) dU dV
h_k = 4*pi*a*r^2*(N*rho_k-beta*a*j_k)
    = d*[N*(q_alpha/a+k*sigma1/r)-beta*q_I], d=2*abs(k).
```

In the radar chart N=a=sqrt(C) and beta=0. This is the complete canonical
clock-energy form, including its kinetic term. The existing W-only metric
insertion is not this energy vertex. The positive weighting and full
angular normalization are fixed before evaluation. The rule fills the
negative spectral subspace of A_k. In a finite orthonormal flux basis this
minimizes Tr(C*A_k) over 0<=C<=I; unresolved zero/gap cases stop the
numerical selection rather than silently assigning occupations. The
existing ground-covariance implementation supplies the projector after
positive scale normalization.

Use smooth incoming trial functions on each null leg:
`y=log(s/(S-s))`, `J=dy/ds=S/(s*(S-s))`, and
`f_n(s)=sqrt(J)*h_n(y)`, where h_n are normalized Hermite functions on
the real line. Their continuum flux Gram is the identity. Every finite
mode has a smooth zero extension at the two ends, since its Gaussian in
log(s) dominates every inverse power. Analytic derivatives, not sampled
waveform differences, feed the current stress interface. This is a
Galerkin representation and supplies no angular ultraviolet completion.

Before a physical reference, a 60-process-CPU-second admission covers
only exact/small-matrix and one-dimensional analytic controls: minimization,
CAR and basis/unit covariance; the independent lapse/shift action vertex;
flux orthogonality, analytic derivatives and endpoint behavior; the
two-frequency signature counterexample. No Dirac evolution, P001 grid,
parameter search, reference JSON or Hadamard claim is admitted by these
controls. A failed identity stops the candidate for diagnosis. The likely
failures are confusing signature with energy, selecting a numerically
unresolved eigenspace, and mistaking finite CAR positivity for a full
four-dimensional state. All remain explicit in the next decision.

The five controls pass in 1.5 process-CPU seconds. In the softened
signature construction of [Fewster and Lang](https://arxiv.org/html/1408.1645v2),
the positive-mass assumption matters. Setting m=0 in its two-frequency
block gives B=[[0,q],[conj(q),0]]. For nonzero q the negative projector
has both diagonal entries 1/2, and its difference from diag(0,1) squares
to I/2, regardless of how small q becomes. A zero block leaves the
occupation undetermined. This exact calculation blocks importing that
positive-mass theorem into this massless candidate; it does not reject
every possible massless state construction.

The product b(U)b(V) defines an explicit spacetime window, including radial
weighting in T,x. It is not a time-only average over a whole Cauchy slice,
nor a causal switched preparation. Its physical duration S_prep is a
declared state parameter, held fixed if the computational domain changes.
Changing this physical window changes the state; it is not a numerical
cutoff to tune against observations. A smooth ultrastatic-deformation
completion remains a separate route with its own geometry assumptions.

For W=0 each chiral energy form reduces, up to positive normalization,
to i*(w*d/ds+w'/2), with w proportional to the positive bump. The
half-density change z=integral ds/w maps its complete flow to the real
line and its generator to i*d/dz. Its negative spectral subspace supplies
the free Hardy reference. No reflecting or periodic endpoint is imposed.
This one-dimensional check is not a four-dimensional state theorem.

Extend the same 60-CPU-second control admission to small manufactured
field grids of at most 16 cells: exact free chiral propagation, the
joint outgoing flux Gram, superposition, and signed-k energy/pressure
parity on a static cylinder. Recheck the five existing mean-field tests
after permitting its unchanged Numba kernel to release the GIL. Only
independent incoming modes run in parallel, with at most eight workers.
These are interface controls; no actual nested-state reference or P001
evolution is admitted until they pass.

All twelve preparation and existing mean-field controls pass in 4.03125
additional process-CPU seconds (5.53125 total charged to these controls).
The manufactured cylinder has full cross-leg covariance and exact signed-k
parity; the zero-potential control checks analytic chiral derivatives,
superposition and the full flux Gram. They validate the finite interface,
not the candidate's ultraviolet limit.

### Actual-geometry clock-energy candidate admission

The next bounded calculation supplies the missing finite parent covariance
on the actual cap-wall geometry, with the smooth cap continuation as its
control. It borrows the already checked Dirac action, energy-minimizing
fermionic covariance and characteristic transfer. The new connection is
applying one stated energy-window rule to both geometries and both members
of the conserved compact-mass family, rather than copying one flat state.

Freeze the method and controls before the reference. Use only N=64/128,
b=0 and 1e-4, signed k=-1,+1, and 2/4 smooth radial trial modes per null
leg. All sixteen cases propagate the four-mode basis once and reuse its
two-mode subspace. The two signs are calculated separately, including each
full magnetic multiplet once. The physical preparation duration is the
existing proper-wall duration .1256873632901346*ell, independent of grid.
The detector event is the existing cap event U=.8*S,V=.7*S. The state rule
is window-nonlocal and parametric; no causal source switch is asserted.

Record spectra/gaps, CAR/purity, raw flux-Gram conditioning, transfer flux,
cross-leg correlations, and grid/mode-count changes. Compute the finite
source-induced pressure difference
`[pi_nested(b)-pi_cap(b)]-[pi_nested(0)-pi_cap(0)]`, summing both signed
sectors. Divide by the common classical A0=-3*b/r^3 only after this
subtraction. Its tidal-ratio coefficient carries 8*pi*(l_P/ell)^2;
the independent illustrative length remains ell=10 km. No sky value sets
the state, length, masses or mode count.

The allowance is 240 charged process-CPU seconds including startup and an
8-CPU-second exit reserve, at most 30 CPU seconds per case, eight process
workers with one BLAS thread each, at most 2 GiB working arrays and 48000
summary bytes. No arrays are retained. Independent modes can use Numba's
GIL-free interface, but the reference instead parallelizes independent
cases to avoid nested worker pools. A failed gap, CAR/flux control, domain
condition or budget stops execution for diagnosis. A nonshrinking field
spread or a large change under radial trial count prevents a continuum
claim; it does not authorize more grids automatically.

The immediate decision is whether this common state rule yields a stable
finite nested-versus-cap response and what controls its cutoff dependence.
The next analytical test concerns the actual joint angular/frequency
short-distance symbol, including cross-leg covariance. These low sectors
alone cannot supply a Hadamard state, renormalized vacuum stress, an
observational error interval or the full scalar-potential ratio Psi/Phi.

### Evaluated finite clock-energy state

The [frozen reference](../results/summaries/nsc-atlas-dirac-preparation-v1.json)
uses method commit `0892dbc` and completes all sixteen cases. The physical
window, masses, detector event and state rule remain unchanged. The raw
joint flux-Gram defect is at most 1.18e-15; signed-sector pressure parity
differs by at most 8.53e-14. The largest 64/128 generator difference divided
by its fine spectral gap is 1.37e-5, and the covariance grid spread is at
most 3.08e-6. Thus these finite eigenspaces are numerically separated.

After summing the two signed multiplets and subtracting the b=0 nesting
anisotropy, the point tidal diagnostic is:

| Radial modes per leg | 128-cell coefficient of (l_P/ell)^2 | Absolute 64/128 spread |
| --- | --- | --- |
| 2 | -.001735011967 | 4.412917e-6 |
| 4 | -.022765039840 | 2.727998e-5 |

The mode-count difference .021030027873 is much larger than either grid
spread. At ell=10 km the finite coefficients convert to -4.53234e-81 and
-5.94687e-80. Those are cutoff-dependent state diagnostics, not vacuum
predictions or observational error ranges. The four-mode source pressure
difference is 1.920573595e-6 in the declared ell units, at r/ell=.5210822567;
its background denominator is A0*ell^2=-.002120324828. The b=0 contribution
has been explicitly subtracted, so this is not a pre-existing wall
anisotropy divided by the compact mass.

The eight-worker run uses 19.828125 observed CPU seconds, 27.828125 charged
including the reserve, and 2.81714 calculation wall seconds. It retains
47918 summary bytes and no full arrays. Field discretization is not the
dominant unresolved issue. The two/four-mode state truncation is; no
additional grid is scheduled. An exact free-energy spectral-projector
comparison is the next useful analytical control before changing that
cutoff. The renormalized vacuum and full Psi/Phi remain uncomputed.

The remaining part of the existing 60-CPU-second control allowance admits
one exact free-projector comparison, with matrices of size at most 64 and
one-dimensional quadratures of order at most 48. No physical field/grid,
angular extension or new reference is admitted. In a manufactured free
clock flow w=s*(S-s)/S the logit coordinate is exact, so Hermite/Fourier
identities give the continuum projector independently of the finite
energy minimizer. This directly tests the truncation exposed above.

### Exact free projector diagnoses the finite-state cutoff

That control passes in 1.53125 process-CPU seconds (7.0625 total for this
preparation's controls). For w=s*(S-s)/S, z=log(s/(S-s)) and
phi(z)=sqrt(w)*f(s), the free generator is i*d/dz. With normalized Hermite
functions h_n, its matrix and exact negative-projector compression are

```text
A_mn = i*[sqrt(n/2)*delta_(m,n-1)-sqrt((n+1)/2)*delta_(m,n+1)]
C_mn = i^m*(-i)^n * integral_0^infinity h_m(p)*h_n(p) dp.
```

The Fourier convention is exp(-i*p*z): i*d/dz has eigenvalue -p, so p>0
is occupied. In particular, C_01=-i/sqrt(2*pi), while the two-mode energy
minimizer has C_01=-i/2. Their entry difference is .1010577196. Independent
half-line quadratures at orders 16/48 agree; even matrix sizes 2 through
64 approach the exact entry. The compressed continuum covariance has
eigenvalues strictly between zero and one and is not itself a projector.
This is compatible with a pure state on the complete algebra: restricting
it to a proper trial subspace generally produces a mixed covariance.

Thus spectral projection and subspace compression do not commute:
`P_M*1_(-infinity,0)(A)*P_M != 1_(-infinity,0)(P_M*A*P_M)`.
The finite minimizer is valid on its own truncated CAR algebra, but purity
there is not a test that it equals the restriction of the continuum state.
This exact free result diagnoses an intrinsic low-mode truncation effect.
It neither proves that the physical bump/curved state diverges nor repairs
the actual joint angular/frequency limit. More field cells cannot fix it.

### Leading joint angular-frequency sign and the remaining UV equation

The first useful ultraviolet check keeps radial momentum p and angular
integer k large together. In the canonical radar frame, the local
principal symbol is

```text
H1 = alpha*p + sigma1*k*N/r
omega = sqrt(p^2+(k*N/r)^2)
P_sigma = (I+sigma*H1/omega)/2, sigma=+1,-1.
```

Since N,r are positive in this finite domain, omega is positive away
from the zero covector. In T,x coordinates the energy-window weight is
chi=2*b(T-x)*b(T+x)/S. The corresponding symmetric local energy insertion
is (chi*H+H*chi)/2; its derivative-of-weight term is essential. Transport
it to a common interior Cauchy slice using the actual field evolution.
On a geometrical-optics branch its leading expectation is
`a_sigma=sigma*d*integral chi(T,X_sigma(T))*omega_sigma(T) dT`.
For rays meeting a compact detector inside the positive window, that
integral is strictly positive. The two branch magnitudes need not be
equal: the rays generally cross different weights. Nevertheless the
negative spectral projector of this leading symbol has the usual
negative-frequency polarization P_minus. No particle label is needed
for this sign calculation.

This establishes a leading-symbol compatibility condition, not a
four-dimensional state theorem. Every lower symbol and the propagation
across the moving material interface remain relevant. A completed local
state requires the transported covariance to satisfy its projector and
Dirac evolution equations to all ultraviolet orders, modulo smooth
kernels. Matching P_minus only at leading order is insufficient.
The transparent conformal potential is C1 at the wall; a theorem that
assumes globally smooth coefficients cannot silently replace its
piecewise-smooth transmission analysis.

The next physical connection is a local two-point construction whose
common cap singular part is removed before the joint frequency/angular
limit. For the energy-state route this entails a self-adjoint continuum
energy operator, control of its spectral projector near zero and at the
null endpoints, and all-order transmitted asymptotics on the detector.
Only then can a resolvent representation of sign(A), or a transported
Hadamard two-point function, replace the present low-mode minimizers.
The finite reference does not determine which construction meets those
conditions. No new grid or vacuum number follows from this leading sign.

### Reusing the fourth-order symbol on the actual radar geometry

The existing `nsc_spatial_reference_symbol.reference_projector` already
solves the formal Weyl projector and time-transport identities through
order four. Its Taylor algebra is reusable; its old KS metric and state
are not inputs to this calculation. Supply instead the actual radar
lapse, radial length, radius and zero shift, with analytic derivatives
through order four. A fixed unitary Pauli cycle changes its massless
symbol p*sigma3+(k*N/r)*sigma2 to the current p*sigma2+(k*N/r)*sigma1.

Construct those jets by Taylor solving the existing proper-wall and
optical ODEs, then composing the exact one-sided optical inversion.
No waveform, stored field or finite-difference derivative is used.
At the wall, take separate cap/exterior jets; a single smooth jet there
would erase the actual curvature jump. Local subtraction at a detector
strictly inside the common cap uses that cap's jets.

This new connection has 120 process-CPU seconds for analytic-interface
controls: the known metric/curvature derivatives, the independent wall
second-jet identity, static and angular/unitary limits, the independent
first time-transport correction and joint momentum/angular homogeneity
through order four. Include the three existing Weyl-core tests. No
physical state evolution, P001 grid, new cutoff reference or full arrays
are admitted. Stop on a derivative, Pauli-frame or transport mismatch.

The decision is whether the existing finite algebra supplies the actual
local subtraction terms. It cannot select the missing global covariance
or prove its remainder. Likely failures are losing the proper-clock
derivatives in the optical composition, importing the KS Pauli frame,
and treating the formal subtraction symbol as a positive state.

The [radar adapter](../src/recursive_horizons/nsc_atlas_dirac_symbol.py)
now supplies those jets and all five local symbol orders, reusing the
existing Weyl core. Seven independent new/core controls pass across
the development runs, using 7.875 process-CPU seconds. The initial
fourth wall-derivative comparison missed an absolute threshold by
2.25e-10 at a derivative magnitude 12426.94 (1.81e-14 relative). A
scale-aware floating-point tolerance corrects that test; no geometry or
symbol equation changed. Curvature and null-constraint equations checked
through two derivatives independently control the order-four metric jets.

For H=p*sigma2+W*sigma1, the first transported correction is explicitly

```text
c0 = (I-H/omega)/2
c1 = -p*W_T*sigma3/(4*omega^3), omega=sqrt(p^2+W^2).
```

The first spatial Weyl terms cancel in this massless planar case. Solving
`i*d_T c0=[H,c1]` fixes the displayed sign. The actual geometry's nonzero
W_T therefore matters even at the first subleading order; a bare
instantaneous projector omits it. The controls compare this independently
derived expression, signed-sector unitary parity and the joint scaling
`c_n(lambda*p,lambda*k)=lambda^(-n)*c_n(p,k)` through n=4.

### The wall's curvature enters the second symbol order

Let d=M_parent-b and use the actual wall radius R(T), speed v and
acceleration Rddot. On each side the radar-wall values are

```text
N=1, N_T=N_TT=0, r=R, r_x=beta, N_x=beta'
r_xx=Rddot+F', N_xx=(beta')^2+F''/2.
```

Here primes are derivatives with respect to areal radius; x is the radar
coordinate, not a proper-normal coordinate away from the wall. Israel
matching gives common beta'-beta/R. Since W=k*N/r, its first derivatives
and all purely tangential jets agree. Its second radial radar jet obeys

```text
J=[W_xx] = k*(Delta F''/(2*R)-Delta F'/R^2+Delta F/R^3)
          = -6*k*(M_parent-b)/R^4.
[W_TT]=[W_Tx]=0.
```

The conformal radial curvature jump is 12*d/R, consistently giving
J=-k*Delta R2_conformal/(2*R^3). Thus making a single C-infinity wall
jet would delete a known part of the geometry.

With all lower and tangential second jets fixed, the first local
subtraction-symbol discontinuity is

```text
[c2] = -J/(16*omega^7) *
       (5*W^2*p^2*sigma1 + W*p*(2*p^2-3*W^2)*sigma2).
```

An independent Pauli/Sylvester calculation checks this against the two
one-sided fourth-order calculations. In detail, writing P=c0,Q=I-P,
c_W=d_W c0 and c_pp=d_p^2 c0, the two separate defects are
`Delta G2=-J*{c_W,c_pp}/8` and `Delta F2=J*[sigma1,c_pp]/8`.
Then `Delta c2=-P*Delta G2*P+Q*Delta G2*Q+[H,Delta F2]/(4*omega^2)`
reduces to the displayed expression. This is a jump of the formal local
subtraction data, not a new boundary condition on the physical state.

### Sufficient remainder for the current mean-stress comparison

The all-order local Hadamard target is stronger than the immediate
mean-stress calculation needs. If an actual covariance and its local
orders c0 through c4 differ by a symbol of order -5, uniformly on the
compact cap detector, its one-derivative stress remainder is integrable.
The joint radial/angular measure in three spatial dimensions scales as
lambda^2*d_lambda; stress contributes one power, so the remainder tail
is bounded by a constant times `integral_Lambda^infinity lambda^(-2)
d_lambda = 1/Lambda`. Order -4 would instead allow a logarithmic tail.
This conclusion assumes the stated symbol bounds and local operator
description; a five-term formal series supplies neither on its own.

At fixed b the nested and smooth models have identical cap jets, so their
local subtraction terms agree before the coincidence limit. A proved
order-minus-five bound for the actual state difference would suffice
for this mean-stress comparison even before an all-Wick-powers theorem.
The new adapter makes those local terms computable. The missing object
remains the joined covariance with a controlled remainder after actual
transmission; no bound is inferred by fitting the two/four-mode data.

Boundary-state existence results need their exact hypotheses. The
[published MIT-Dirac paper](https://ems.press/journals/dm/articles/13131883)
explicitly assumes an appropriate singularity-propagation result; it does
not prove transparent transmission for this wall. The
[2026 timelike-boundary construction](https://arxiv.org/html/2609.13358v1)
treats a smooth scalar Dirichlet problem. Its positive-energy splitting
argument is a useful route to inspect, rather than a theorem already
applicable to this piecewise-smooth Dirac geometry.

### A concrete static-product parent and deformation route

Develop the separately listed ultrastatic route as an explicit alternative
to the energy-window state. The root is the massless Dirac ground state on
R^(1,1) x S^2 with sphere radius R0=.45*ell. Its signed angular sector has
the exact radial Hamiltonian p*sigma3+mu*sigma1 in the chiral frame,
mu=k/R0, and omega=sqrt(p^2+mu^2). Filling its negative band specifies the
entire radial spectrum, with no fit to a sky observation or finite-mode
energy minimization. This is an additional state prescription, not a
consequence of the inheritance relation or the previous candidate's result.

For negative-energy plane modes set a=(omega-p)/2,b=(omega+p)/2. The
incoming components on V=0 and U=0 have phases exp(i*a*U),exp(i*b*V).
Their exact covariances are equal line Hardy kernels on each leg, plus
the nonzero cross-leg kernel
`C_plus_minus(U,V)=-mu*K0(abs(mu)*sqrt(U*V))/(2*pi)` for positive U,V.
Here K0 is the modified Bessel function. The full covariance is a
restriction of the global cylinder projector and is generally mixed on
the finite diamond. Taking a projector after finite compression would
change this declared state.

This route has 60 process-CPU seconds for bounded root-state controls:
one-dimensional Fourier/Bessel comparison, finite smooth-mode compressions
of sizes at most 4 per leg with Gauss orders at most 128, CAR and nested
subspace consistency, analytic incoming plane-wave derivatives on at most
16 manufactured cells, and folded-interface flux algebra. No new physical
field reference, UV cutoff campaign or P001 grid is admitted.

The preparation geometry is auxiliary state data, not a claimed physical
formation history or a sourced Einstein solution. Use a fixed smooth
switch chi(s), zero outside (.1*S,S), equal to one on [.2*S,.9*S]. Set
q_aux=1/R0+chi(U)*chi(V)*(N_physical/r_physical-1/R0). Blending this q,
instead of N and r separately, preserves the actual C1 Dirac potential
at the wall. The auxiliary product metric has r_aux=R0 and
N_aux=a_aux=R0*q_aux. Where both switches are one it is conformal to the
physical metric, with unchanged massless canonical spinor. The stress
and clocks must still be evaluated in the physical metric; no quantum
conformal invariance or anomaly cancellation within one model is assumed.
The previous cap detector (.8*S,.7*S) has a future Cauchy neighborhood
wholly inside this unchanged-operator region.

This specifies the missing parent data and a concrete preparation operator.
The bounded-potential argument below establishes transparent unitary
evolution on an admitted regular patch. Sufficient transmitted regularity
remains to be established before its full stress can be claimed finite.
The fixed support and root radius are state parameters. No result from
the energy-window state is transferred into this alternative.

### Exact correlated root data and its finite compressions

The [implemented root](../src/recursive_horizons/nsc_atlas_dirac_preparation.py)
uses the Fourier measure dp/(2*pi). For the normalized negative spinor,
`|u_plus|^2=a/omega=-da/dp`, `|u_minus|^2=b/omega=db/dp` and
`u_plus*conj(u_minus)=-mu/(2*omega)`. Changing variables gives the two
Hardy kernels
`C_diag(s,t)=delta(s-t)/2+i*PV(1/(s-t))/(2*pi)`.
The cross-leg integral is

```text
C_plus_minus(U,V) = -mu/(4*pi) * integral_0^infinity
                    exp(i*a*U-i*mu^2*V/(4*a)) da/a
                 = -mu*K0(abs(mu)*sqrt(U*V))/(2*pi).
```

At U=V=s it independently matches the equal-time cylinder Fourier kernel
`integral dp*(-mu/(2*omega))*exp(-i*p*s)/(2*pi)`. The two positive-energy
complements supply the other Hardy half and the opposite cross kernel,
so the full incoming CAR is the identity. Equivalently, restrict the
complete cylinder projector to a diamond Cauchy segment and use its
conserved flux map to the two incoming legs. Positivity 0<=C<=I follows
from that restriction. Its finite compressions need not be projectors.

The real smooth-mode PV pairing is evaluated after antisymmetrization;
its diagonal is removable and has limit
`(f_n*f'_m-f_m*f'_n)/2`. This contribution is retained at the quadrature
diagonal. Independent exact line-Hermite data check the normalization
C_01=-i/sqrt(2*pi), not -i/2. The K0 cross block is not discarded.

Six root/interface controls pass in 3.546875 process-CPU seconds across
development. They compare the Bessel kernel with an independent Fourier
integral, exact Hardy compression, 96/128 quadrature, CAR and mixedness,
consistent two/four-mode subspaces, units and signed sectors. Analytic
negative-frequency plane data propagate as the same known cylinder
solution on both legs in the manufactured field control. The auxiliary
switch preserves both first wall jets and agrees with the actual
canonical operator on the future detector neighborhood. No physical
Dirac cutoff reference or P001 evolution has run for this alternative.
An underresolved root quadrature returns no covariance; spectra are never
clipped or refilled to pass CAR. The exact analytic parent is unchanged.

### Folding the transparent wall and constructing unitary evolution

Fold exterior x=y and cap x=-y onto y>=0, with side order exterior,cap.
Let tau_i act on that side index, alpha=sigma2,beta=sigma1. The normal
and angular matrices are A=tau3 tensor alpha and B=I tensor beta.
Gamma=tau1 tensor I swaps the two sides. Gamma*Psi=Psi is exactly the
original equality of canonical spinors, rather than a reflecting wall.
The [time-generator owner](../src/recursive_horizons/nsc_atlas_dirac_evolution.py)
now exposes these algebraic data.

Gamma is a Hermitian involution, anticommutes with A and commutes with B.
Its allowed rank-two projector P obeys P*A*P=0. The folded interface
therefore has maximal zero flux. For the frozen imaginary-energy problem,
lambda=sqrt(W^2+zeta^2)>0, a decaying mode with allowed trace J*v,
J=(I,I)^T/sqrt(2), obeys

```text
M = (i*lambda*A+W*B-i*zeta*I)*J
M^dagger*M = 2*lambda^2*I.
```

It has no nonzero allowed decaying kernel. This verifies a uniform frozen
elliptic compatibility condition; it is not a full wavefront propagation
theorem. In particular, same-energy reflection p->-p is not a
positive/negative-energy mixing counterexample. The explicit local
W_xx discontinuity enters c2, not c1, as checked above.

For existence of the actual auxiliary evolution, unfold to the full
radial line. In each signed angular sector its canonical Hamiltonian is
`H_k(T)=-i*alpha*d_x+beta*k*q_aux(T,x)` on H^1(R,C^2). On a regular
admitted patch q_aux is real, bounded, C1 in x and continuous in T, and
equals 1/R0 outside a compact set. Its difference from the cylinder
Hamiltonian is a bounded self-adjoint multiplication operator.

In the interaction picture, the time-ordered Volterra series converges
in operator norm: its n-th term is bounded by A_k^n/n!, where
`A_k=abs(k)*integral ||q_aux(T)-1/R0||_infinity dT` is finite. Uniqueness
and differentiating the norm give unitary evolution, with the inverse
constructed backward. Continuity across the unfolded line implements
the same transparent domain. Taking the direct sum over signed k and
their 2*abs(k) magnetic states gives a strongly continuous unitary:
every sector has norm one, so finite-sector approximation controls the
tail. The exponential Volterra bound need not be uniform in k for this
existence argument; it would be insufficient for a stress-tail bound.

Consequently `C(T)=U(T,T0)*C_cylinder*U(T,T0)^dagger` defines an actual
quasifree CAR state for this declared preparation, including all angular
and radial modes. Restriction to the future matched region transfers it
to the physical massless field by the regular classical conformal map.
This establishes state data and evolution, not its renormalized stress.

There is also a useful exact gap estimate. If q_aux>=q_min>0, then
`{H_k,beta}=2*k*q_aux`. For any domain vector,
`abs(k)*q_min*||psi||^2 <= abs(Re <H_k psi,beta psi>)`
`<= ||H_k psi||*||psi||`. Thus the instantaneous operator has a gap
of at least abs(k)*q_min. Positivity and compact support of the blend
give q_min>0 on the admitted regular patch. Time derivatives of the
potential are bounded multiplication operators. These are concrete
inputs for a finite-order operator adiabatic estimate, including the
large-angular sector, without assuming a globally smooth glued metric.

The next unresolved estimate must also control large radial momentum and
local stress derivatives uniformly, rather than use the exponential
Volterra bound or the angular gap alone. The fourth-order cap symbol is
available, and the parent state is now explicit in this alternative.
What is still missing is an actual order-minus-five bound, or another
controlled sufficient coincidence estimate, for its propagated interior
two-point function. The gap estimate does not supply that bound by itself.

### Real-frequency boundary ellipticity for the next regularity estimate

The previous imaginary-energy check is not yet the real-frequency
boundary condition used in a microlocal argument. Check it directly,
using at most 30 process-CPU seconds for exact/frozen matrix identities.
No physical field, numerical grid or cutoff reference is admitted.

For real temporal frequency tau and signed angular potential W, a
decaying normal root exists in the spacelike region abs(tau)<abs(W),
with lambda=sqrt(W^2-tau^2). Its allowed boundary symbol is
`M=(i*lambda*A+W*B-tau*I)*J`. The opposite normal signs cancel its
cross terms, giving

```text
M^dagger*M = 2*(W^2*I-tau*W*beta)
min singular_value(M)^2 = 2*abs(W)*(abs(W)-abs(tau)).
```

For abs(tau)<=(1-epsilon)*abs(W), the lower bound is
`sqrt(2*epsilon)*abs(W)`. Both signs of k and either angular polarization
obey the same bound. At glancing the bound tends to zero; the calculation
does not discard that characteristic set or assert ellipticity there.
This supplies one actual input to the boundary regularity argument.

Extend the same 30-CPU-second allowance to bounded-band inequalities and
the exactly solvable constant-mass line resolvent. These controls concern
the following global estimates, not a new approximation to the state.

### Balanced energy estimates and an explicit limitation

Let A,B be positive spectral-band operators with lower bounds m_A,m_B.
For bounded Y the off-band Sylvester inverse is
`L^(-1)(Y)=integral_0^infinity exp(-s*A)*Y*exp(-s*B) ds`.
For a positive integer n, its n-fold iterate has the same integrand
weighted by s^(n-1)/Gamma(n). Operator Cauchy--Schwarz, applied before
taking the time integral, gives the exact balanced endpoint

```text
|| A^(n/2) * L^(-n)(Y) * B^(n/2) || <= 2^(-n) * ||Y||.
```

Indeed the two square integrals are Gamma(n)*I/2^n. Equal scalar bands
saturate the constant. Thus the divergent integral of separate supremum
norms at this endpoint is a failed estimate, not proof that the endpoint
bound fails. The same argument for the double inverse gives

```text
||A*L^(-2)(Y)|| <= ||Y||/(4*m_B)
||L^(-2)(Y)*B|| <= ||Y||/(4*m_A).
```

For an actual instantaneous projector P_minus, its first transported
off-band correction is -i*L^(-2)(H_T)_plus_minus. These inequalities
therefore control its full one-sided energy graph norm, uniformly in
the angular sector when ||H_T||/abs(k) is bounded. This goes beyond an
angular-only norm estimate and uses no spatial commutator across the wall.
It still does not express the full evolved covariance remainder as an
n-fold inverse of one bounded source. Products, time derivatives and
detector localization in that remainder need their own estimate.

For example, a CAR state with excitation fraction a_k=exp(-k^2), constant
in radial momentum, obeys every angular operator-norm decay bound. Its
energy excess per magnetic state up to radial cutoff L is
`a_k*(L*sqrt(L^2+mu^2)+mu^2*asinh(L/abs(mu)))/pi`, which diverges as
L^2 for every fixed k. This is a counterexample to the angular-only
inference, not a claim that the prescribed cylinder preparation has
those occupations.

### Frozen resolvent locality from the actual Dirac gap

For the instantaneous H_k at imaginary energy i*zeta, put
rho=sqrt(m^2+zeta^2), m=abs(k)*q_min. For a bounded real Lipschitz
weight F, exponential conjugation gives
`exp(F)*H_k*exp(-F)=H_k+i*alpha*F'`.
Choose abs(F')<=theta*rho, 0<theta<1, and F differing by theta*rho*d
between separated source and detector supports. The bounded perturbation
inverse then gives

```text
||chi_D*(H_k-i*zeta)^(-1)*chi_E||
    <= exp(-theta*rho*d)/((1-theta)*rho).
```

This uses the actual constant radial principal matrix and the proven
gap. It remains valid for a bounded potential; no globally smooth
spatial metric is assumed. Bounded weights and a limiting smooth
approximation avoid an unbounded-conjugation domain assumption.

If the nested and smooth q agree around D and differ only on E at
radial distance at least d, the resolvent identity supplies two such
factors. Let q_min be a common positive lower bound and delta_q a
supremum bound on their difference. Integrating their resolvent
difference in the sign-projector formula gives the coarse estimate

```text
||chi_D*(P_minus,N-P_minus,S)*chi_D||
 <= delta_q/(2*(1-theta)^2*q_min)
    * exp(-2*theta*d*q_min*abs(k)).
```

The intermediate integrand is bounded by abs(k)*delta_q times
`exp(-2*theta*d*sqrt(m^2+zeta^2))/((1-theta)^2*(m^2+zeta^2))`.
The remaining integral is at most pi*exp(-2*theta*d*m)/m; the
projector formula contributes 1/(2*pi). The radial spectral variable
has therefore been integrated, rather than replaced by an angular norm.

These are frozen operator bounds implemented in the
[bound calculator](../src/recursive_horizons/nsc_atlas_dirac_bounds.py).
Its q_min and delta_q inputs must be genuine uniform bounds, not sampled
minima or fitted values. No numerical bound on the actual full stress
is reported. The prepared state is time-evolved, so replacing it with
the instantaneous projector would change the prescribed calculation.

Six new real-boundary/operator controls pass in 1.953125 process-CPU seconds
across development. They include the full angular principal matrix before
diagonalizing k, independent Sylvester solves and sharp equal-band limits,
and the exact constant-mass whole-line resolvent between separated half
lines. That last norm is exp(-rho*d)/(2*rho), independently below the
general conjugation bound. Very small positive bounds are evaluated in
logarithms and do not silently become exact zero by underflow. No field
or P001 evolution ran.

### Normal regularity of the prepared one-particle field

One remaining microlocal hypothesis can be established directly from
the constructed unitary field. For an initial L2 vector v, let
psi_v=U*C_root*v on either folded side. Unitarity bounds its spacetime
L2 norm on any compact time interval by a constant times ||v||.
For a smooth compact tangential test phi(T,Omega), define
`h_v(y)=integral phi^dagger*psi_v dT dOmega`.
It is L2 in y. The first-order equation, with invertible constant A,
expresses each normal derivative as tangential T and angular derivatives
plus smooth one-sided multiplication. Integrating those tangential
derivatives onto phi and repeating gives
`||d_y^n h_v||_L2(I) <= C_(phi,n,I)*||v||` for every finite n on a
compact one-sided interval, including its boundary endpoint.

For rigor this is an operator estimate from the initial Hilbert space
to H^n(I), not an assumption that a point field has a finite L2 density.
The adjoint applied to delta_y and its derivatives, continuous in
H^(-n) for n sufficiently large, gives a Hilbert-valued function smooth
in y with values in tangential distributions. This proves normal
regularity. The transparent trace relation extends from a dense smooth
initial domain using the same bounds. The argument uses one-sided
smooth coefficients and transfers derivatives; it never differentiates
the glued coefficient through its W_xx jump.

There is also smooth well-posedness for smooth past-compact folded source
and boundary data. Lift the forbidden boundary trace to a collar, use
unitary Duhamel evolution for the homogeneous trace, and differentiate
in T and the commuting angular operator K. Energy estimates bound all
such derivatives using finitely many smooth source derivatives. Recover
normal derivatives from invertible A. Green's identity and the vanishing
allowed boundary flux give distributional uniqueness in the normally
regular class by duality against smooth advanced solutions.

<a id="the-precise-remaining-boundary-calculus-application"></a>

### Local Hadamard property for the declared preparation

**Proposition.** On an admitted regular patch with the positive auxiliary
q defined above, the cylinder-prepared state is locally Hadamard in every
relatively compact smooth interior neighborhood inside the matched
operator region. The nested and smooth-cap preparations consequently
have a smooth physical two-point difference on a common cap neighborhood.
This statement concerns the declared probe and preparation, with no
direct spinor surface interaction.

The standard local calculus input is the finite-bundle first-order
microlocal elliptic boundary-parametrix theorem with a non-characteristic
normal matrix and a complementing local boundary projection. Its
hypotheses can be checked here without assuming propagation across a
nonsmooth glued metric. The one-particle field is normally regular by
the energy argument above. On each folded side coefficients are smooth
up to the boundary. The real spacelike decaying map is bijective, with
matching rank and the uniform bounds just calculated.

More explicitly, in a signed angular eigenspace let
`N0=i*alpha*(W*beta-tau*I)/lambda` and `S=(I+N0)/2`.
Then N0^2=I. S and I-S are complementary stable projections, generally
oblique. The [implemented Poisson symbol](../src/recursive_horizons/nsc_atlas_dirac_evolution.py)
for the unnormalized boundary jump g is
`exp(-lambda*y)*(S*g,-(I-S)*g)`.
It solves the frozen operator and reproduces exactly
psi_exterior-psi_cap=g. On closed spacelike conic subsets lambda is
comparable to the tangential covector norm. The displayed inverse has
the smooth symbol estimates required for the standard boundary-parametrix
recursion; successive errors use the same inverse and smooth one-sided
jets. The calculation before angular diagonalization verifies the same
complementing condition in every tangential direction. Thus local boundary
elliptic regularity confines the b-wavefront set to the compressed causal
cone. Glancing belongs to that cone and has not been declared elliptic.

To propagate orientation, use a small-b microlocal partition into the
two disjoint oriented cones and the elliptic complement. The latter is
smooth. Applying D to the two remaining components gives opposite-cone
distributions whose sum is smooth, so each source is smooth. The same
argument for the forbidden boundary trace gives smooth boundary data;
the partition need not itself preserve the boundary condition. These
microlocal and trace properties extend componentwise to the finite bundle
in the normally regular class. A time cutoff of the initially absent
orientation has smooth source and boundary data and vanishes sufficiently
far in the past. The smooth well-posedness and distributional uniqueness
proved above make this component smooth everywhere. The two oriented
cones stay disjoint at glancing because their nonzero covectors have
opposite signs of tau and q_min>0.

The initial occupied and unoccupied cylinder fields supply the two
orientations. On a small smooth cap neighborhood separated from the
wall, the local CAR propagator is the ordinary smooth Dirac propagator.
Its splitting by the two opposite orientations is the local Hadamard
splitting. The physical conformal map is smooth and nonvanishing there,
so it preserves this property. Both states have the same cap geometry
at fixed b; subtracting their common local singular part proves the
stated smoothness of their difference. This completes the local argument.

The [source study](arxiv-studies.md#actual-cap--wall-dirac-evolution-18-october-2026)
records the normally-regular microlocal framework inspected. The scalar
Dirichlet existence theorem is not applied to this spinor system; the
coupled projection and first-order boundary symbols are calculated here.
A separate native Grok review checked this adaptation and its
normal-regularity estimate.

For a single covariance, subtracting c0 through c4 leaves a qualitative
order-minus-five local symbol remainder, not a remainder of every order.
For the two prepared states on the same cap, all local singular terms
cancel, so their difference is smooth. Neither conclusion supplies a
computable numerical tail constant. No full vacuum stress value or
Psi/Phi has been evaluated. Surface-coincident stress, full wall
renormalization, physical formation of the auxiliary preparation and
global horizon/asymptotic states are outside this local proposition.

### Parent modes and a common equal-time kernel

The existing 60-CPU-second root-control allowance also covers composition
of its plane data, auxiliary propagation and a common kernel. Only the
same manufactured 8/16-cell cylinder limit is used; no actual new field
reference is admitted. `PreparedCylinderModes` now propagates coherent
negative-energy parent columns and returns canonical xi and analytic xi_x.
Its equal-time point split is x+eta/2,x-eta/2, with incoming dp/(2*pi)
integrated before comparison with the independently integrated local Weyl
momentum reference. The existing kernel data type and spectral declaration
are reused without importing their old KS geometry or vertex frame.

The interface requires an explicit spectrum declaration and positive
quadrature weights. A partial panel remains partial. It checks that the
physical endpoints lie in the matched-operator region; auxiliary metric
factors are not substituted into physical stress. The kernel derivative
is `(xi_x,plus*xi_minus^dagger-xi_plus*xi_x,minus^dagger)/2`.
Incoming root momentum and local Weyl momentum are different integration
variables. Comparing their cutoffs mode by mode would be incorrect.

The independent conformal-cylinder limit compares this propagated kernel
and its split derivative with the exact negative-band Fourier projector,
using the same manufactured 8/16-cell control. The fine error is below
2e-4 and the ratio of coarse to fine errors exceeds three; coincident
Hermiticity and positivity also pass. This seventh root/interface control
used 1.28125 process-CPU seconds in the recovered check. It introduces no
physical reference result, spectral-tail estimate or new P001 evolution.

### Pairing the complete root covariance in null coordinates

A direct Volterra decomposition offers a radial-complete alternative to
cutting off incoming p. Set A=-i W/2 and write L_X=E_X+G_X, where
E_X f=(f_plus(U),f_minus(V)). The kernels obey the original integral field
equations, not a second propagation law. In particular,

\[
G_{-+}(u,v;s)=A(s,v)+\int_s^u A(t,v)G_{++}(t,v;s)dt,
\quad
G_{++}(u,v;s)=\int_0^v A(u,t)G_{-+}(u,t;s)dt.
\]

Exchange u,v and the two spin labels for the other column. Their endpoint
values include

\[
G_{--}(u,v;v)=-\frac14\int_0^u W(t,v)^2dt,
\qquad G_{+-}(u,v;v)=-\frac i2 W(u,v).
\]

Thus extension by zero beyond a leg's causal endpoint need not be H1.
Strong convergence of a compressed root covariance alone does not prove
convergence of its point stress. An H1 bound for an unmodified
logit-Hermite truncation fails at precisely this endpoint.

Instead retain E and the one-sided G explicitly and evaluate L_X C0 L_Y^dagger
against the analytic Hardy/Bessel root kernel. The direct/direct block
is common and cancels between the two geometries. Direct/G and G/G blocks
still require their Hilbert-transform and endpoint contributions before
the relative coincidence limit; they are not separately discarded.

**Bounded method control.** Reuse the existing root covariance and
characteristic law. The missing connection is its evaluation on supplied
one-sided polynomial Volterra maps, including their endpoint evaluation.
Use the exact constant-cylinder Riemann kernels and compare the full
off-diagonal covariance with its independent Bessel Fourier transform.
Only polynomial orders 12/24 on this manufactured limit are admitted,
with 60 process-CPU seconds, no physical field or P001 grid, no arrays or
reference JSON retained. Stop on a failed identity, an error outside its
numerical floor that does not contract, or budget. The same allowance
covers an analytic spatial response map, retaining both its moving upper
endpoint and its direct f-prime evaluation. Compare its split derivative
with the independent differentiated Fourier/Bessel formula. This is a
full radial root pairing on the supplied maps; their approximation error,
the relative coincidence limit and the angular sum remain separate
calculations.

The [pairing component](../src/recursive_horizons/nsc_atlas_dirac_root_pairing.py)
now evaluates those supplied polynomial maps on their own causal leg
intervals. A spatial-response map retains its direct f-prime term and the
upper-endpoint term `(-G_plus(U),+G_minus(V))`; only its remaining integral
kernel is differentiated. All same-leg delta and PV terms are kept. For
Legendre coefficients g_n on [0,L], the exact real Hilbert transform is

\[
\operatorname{PV}\int_0^L\frac{G(t)}{s-t}dt
=2\sum_n g_n Q_n(2s/L-1).
\]

Its derivative supplies the direct-derivative finite part. The direct
endpoints are required to be distinct; the component does not silently
assign a value at coincidence. The cross-leg Bessel kernel is evaluated
through the convergent identity, with z=abs(mu)*sqrt(L_left*L_right),

\[
K_0(z\sqrt{xy})=\sum_{j=0}^{\infty}
\frac{(z^2xy/4)^j}{(j!)^2}
\left[H_j-\gamma-\log(z/2)-\frac{\log x+\log y}{2}\right].
\]

Polynomial and log moments are integrated analytically. This includes
endpoint evaluation and its derivative; it never truncates the root
Fourier spectrum. The positive omitted-series majorant uses
`b_(j+1)/b_j <= z^2/(4*(J+1)^2)` beyond J and
`H_(J+m) <= H_J+m/(J+1)`. It is reported separately from the adaptive PV
quadrature estimate. The current stable series route is restricted to
z<=4; larger arguments are refused. Neither quantity bounds the supplied
maps' interpolation error, floating arithmetic, an angular tail or a
relative coincidence extrapolation.

Four [independent controls](../tests/test_nsc_atlas_dirac_root_pairing.py)
pass. The nonzero constant-cylinder Riemann kernels reproduce the full
off-diagonal vacuum Fourier/Bessel kernel within 2e-9 at both polynomial
orders; the analytic spatial derivative agrees within 3e-8. Endpoint
logarithms and the cross-leg series majorant are also checked directly.
Development and derivative regression together used 13.09375 process-CPU
seconds. No physical field reference or P001 grid ran. The next physical
input is the actual nested and continued-cap retarded maps, with their
analytic derivatives and controlled representation error. The full
renormalized stress and Psi/Phi remain unevaluated.

### Extracting retarded maps from the existing field kernel

The next missing connection is numerical response data for the complete
root pairing. Reuse `MeanModeGrid`, including its existing analytic
characteristic jets, on a rectangle whose last cell center is the detector.
Independent incoming sample columns give its field and spatial-response
rows. The endpoint value and f-prime functionals are retained explicitly;
only the remaining ordinary kernel is represented on causal source cells
and projected onto Legendre polynomials. No new propagation equation or
root state is introduced.

Before any physical application, compare with the exact nonzero cylinder
Riemann kernels and their spatial derivatives. Admit only N=32/64,
polynomial orders 12/24 from the same field rows, the existing nonzero split .008 and the same T=.11,
x=.006 cylinder control; 60 process-CPU seconds with at most eight workers
and 256 MiB of live field arrays. Independent columns are streamed and
discarded. Both field and derivative errors must contract; a failure stops
the physical application and identifies which row or endpoint term is
responsible. No P001 evolution, arrays or reference JSON are retained.

The first field/derivative row control contracts by about four. With
unconstrained order-24 cell projection, the full-kernel error contracts
only 2.1755 (1.26198e-5 to 5.80091e-6), failing the factor-three control;
the physical run is held. This identifies point evaluation of a projected
kernel as a separate requirement from its L2 projection. Enforce its exact
Volterra endpoint trace while preserving its integral. For p Legendre
coefficients, the smallest L2 correction of zero mean has
`delta g_n=(2n+1)*delta trace/(p^2-1)`, n=1,...,p-1; g0 is unchanged.
For the spatial kernel these traces follow by differentiating its integral
equations at the endpoint. With I_V=int W^2 dV, I_U=int W^2 dU,

\[
G_{++,x}(U,V;U)=-W^2/4+\tfrac14\int_0^V W W_U dt-I_V^2/32,
\quad G_{-+,x}=-i W_V/2-i W I_V/8,
\]

\[
G_{+-,x}=i W_U/2+i W I_U/8,
\quad G_{--,x}=W^2/4-\tfrac14\int_0^U W W_V dt+I_U^2/32.
\]

The same 32/64 control is repeated after imposing these mathematical
traces; no extra grid or covariance target selects the correction.

Endpoint constraints alone do not cure the step representation: its
order-24 derivative error grows 6.29197e-4 to 8.51322e-4. The step-kernel
L2 residual halves with h, rather than falling quadratically. Replace
the step reconstruction by continuous linear interpolation of its source
densities, with the exact causal endpoint and a linear extension to zero
source coordinate. Integrate each segment against Legendre polynomials
exactly before imposing the same endpoint constraint. This changes the
ordinary-kernel representation; the field solver, state, domain, two
grids and acceptance test stay fixed.

The corrected control passes at both polynomial orders. At order 24,
the complete cylinder kernel error falls 5.64820e-8 to 1.36737e-8 and
its spatial derivative error falls 2.10088e-6 to 5.17789e-7. The order-12
errors also contract by more than four. All development, including the
two failed representations and the variable-coefficient control below,
used 46.84375 process-CPU seconds of the
60-second allowance. Native Grok separately checks the direct-value and
direct-derivative split; the continuum Bessel formulas check the actual
numerical extraction. This control is closed; no finer cylinder grid is
needed unless the physical application exposes a failure.

The mixed endpoint terms also have a variable-coefficient control on the
same 32/64 rectangles: W=f(U)g(V), f=1+.5U, g=1-.3V. The exact change of
null variables Ubar=int f^2, Vbar=int g^2 and spinor factors (f,g) reduces
its propagation to unit W. A transported plane solution checks both maps
and their nonzero W_U/W_V endpoint terms. Its field and derivative errors
fall from 2.86251e-8/2.59639e-8 to 7.10273e-9/6.44708e-9, using an
additional 1.390625 process-CPU seconds within the same allowance.

### First actual full-radial finite-split comparison

Use the [bounded producer](../scripts/measure_nsc_atlas_dirac_root_response.py)
to evaluate the same declared parent state on the actual nested geometry
and the continued-cap control. This supplies the actual retarded maps to
the completed root pairing. Fix b=0, kappa=+1 with its magnetic degeneracy
2 counted once, root radius .45 ell, preparation duration S equal to the
quiet wall-clock strip, and detector (U,V)=(.8S,.7S). The switches and
parent state are identical in both geometries.

Admit only field grids 64/128, polynomial orders 12/24 extracted from the
same response rows, and nonzero spatial splits eta/S=.02/.01. Both split
endpoints stay in the common matched cap. Evaluate their complete radial
covariance difference before contracting the canonical symmetric
pressure-shaped bilocal, using the physical midpoint cap metric. The
reference ell=10 km converts its pressure units and actual proper split
distance; it is supplied independently of the result. This is a finite
bilocal observable for the first signed angular sector, not a coincident
stress, full angular sum, compact-source response or Psi/Phi.

Use eight streamed GIL-free field columns and single-threaded BLAS;
allow 240 process-CPU seconds and 256 MiB of live field arrays. Freeze
the method first. Save one summary below 24000 bytes, retaining numbers,
spreads and response/source hashes, never field arrays. Endpoint integrals
come from the potential and its analytic derivatives. Cancel common direct
terms algebraically before numerical covariance pairing.

Compare field, polynomial and separation variations separately. A spread
dominating the signal stops extrapolation and identifies the unresolved
representation. Persistent separation dependence keeps the observable
bilocal. Budget exhaustion saves a diagnosed partial result before any
new calculation is considered. These checks do not supply the uncomputed
angular tail or turn numerical spreads into a continuum error certificate.

### Evaluated full-radial finite-split reference

The [v1 summary](../results/summaries/nsc-atlas-dirac-root-response-v1.json),
frozen at method `dbbfba7`, completes this admitted comparison. It uses the
same complete correlated root covariance in both geometries, without a
root-momentum cutoff or finite occupation minimizer. The physical detector
has r=.5210778462 ell and radar lapse a=1.0124117332. For kappa=+1 only,
with degeneracy 2, the reported symmetrized canonical bilocal is

\[
\Delta\mathcal P_1(\eta)=
\frac{2}{4\pi a^2r^2}\operatorname{Re}[-i\operatorname{tr}(\sigma_3\partial_\eta\Delta K)]
-\frac{2}{8\pi a r^3}\operatorname{Re}[\operatorname{tr}(\sigma_1\Delta K)].
\]

Here a and r are the common physical midpoint metric, not the auxiliary
preparation metric. The real part gives the symmetric eta/-eta insertion
in the fixed canonical radar frame. Its pressure units are hbar*c/ell^4.
It becomes the local pressure difference only if its relative coincidence
limit is controlled; that limit is not taken in this reference.

| Field cells | Polynomial order | eta/S=.02, proper split 25.449495 m | eta/S=.01, proper split 12.724739 m |
|---|---|---|---|
| 64 | 12 | .002008494063 | .001911668141 |
| 64 | 24 | .001861644423 | .001785988527 |
| 128 | 12 | .002007087765 | .001909738536 |
| 128 | 24 | .001842679064 | .001765398991 |

Distances use the supplied ell=10 km. The last value corresponds to
5.58135618e-45 Pa. At this setting the grid spread is 2.05895358e-5
(1.17%), polynomial spread 1.44339544e-4 (8.18%) and separation spread
7.72800732e-5 (4.38%). These are separate numerical indicators, not a
combined error bar or an observational confidence interval. The PV
quadrature estimate is below 3.1e-13 in pressure units; representation and
separation effects dominate, so further quadrature refinement is unnecessary.

The run uses 47 process-CPU seconds and 42.47898 wall seconds, retaining
22899 bytes of observables, diagnostics and hashes. Eight GIL-free field
columns are used, but serial covariance quadrature dominates overall
throughput. A later multi-sector calculation should parallelize independent
cases with one field worker per process, rather than nest two eight-way
pools. No large field arrays are retained.

All admitted settings give a positive finite-split change for the first
signed angular sector under the declared preparation. This does not establish the
angular sum, the local vacuum stress, the compact-source response or full
Psi/Phi. The next error to resolve is polynomial representation, using the
same fixed field grid before considering any finer grid; the angular
series and source derivative remain necessary for the requested slip.
P001 remains unconverged and no new P001 evolution ran.

### Fixed-grid polynomial comparison and useful process parallelism

The v1 polynomial variation is the next unresolved numerical contribution.
Keep the actual field grid fixed at N=128, kappa=+1, the same parent,
detector and two splits. Compare polynomial orders 24/48 from the same
retarded rows. First control order 48 against the exact cylinder kernel
and derivative on the existing N=64 grid, within 30 process-CPU seconds;
no new control grid is needed.

Extend the existing producer with an explicit polynomial stage. Prepare
the small response maps once, then evaluate independent covariance
contractions in eight processes, single-threaded within each process.
Do not overlap two eight-way pools. The actual stage has 240 process-CPU
seconds, including reported child CPU, 256 MiB live field arrays and a
2 GiB process-runtime allowance. Retain only a compact summary and hashes.
The unchanged order-24 cases must reproduce v1 within 1e-9 in the
dimensionless bilocal before the new representation is interpreted.

A baseline mismatch stops interpretation and identifies a refactoring
error. A noncontracting polynomial change stops continuum extrapolation
and identifies representation error. Persistent split dependence keeps
the result bilocal; it is not handled by scanning smaller separations.
Budget exhaustion yields a diagnosed partial result. The angular sum
and compact-source response remain prerequisites for full Psi/Phi.

The fixed-grid order-48 cylinder control passes: kernel error
1.37221e-8 and spatial-derivative error 5.15858e-7, using 8.125
process-CPU seconds. The field grid and propagation equations are unchanged.
The explicit `--stage polynomial` route runs the new representation test;
the default initial route retains the v1 case selection.

### Polynomial refinement exposes an unresolved amplitude

The [v2 summary](../results/summaries/nsc-atlas-dirac-root-response-v2.json),
method `e6bd1bb`, completes the fixed-N128 comparison. The parallel stage
reproduces v1's order-24 values within 2.024e-14, well inside the 1e-9
baseline tolerance. Process refactoring does not explain the new spread.

| Polynomial order | eta/S=.02 | eta/S=.01 |
|---|---|---|
| 24 | .001842679064 | .001765398991 |
| 48 | .002173544722 | .002245374289 |

The new order difference at eta/S=.01 is 4.79975298e-4, 3.325 times the
earlier 12/24 difference, or a 27.19% change from order 24. Polynomial
refinement does not contract. The sign of the separation difference also
reverses. No continuum or coincidence amplitude is inferred from either
pilot. The earlier positive values remain numerical finite-representation
results; they are not converged predictions of the declared continuum state.

All eight covariance processes are observed. The run uses 54.625 CPU
seconds including 44.734375 child CPU seconds, and 12.77376 wall seconds
(4.28 effective cores averaged over both phases). It retains 13567 bytes.
No new field grid ran. The next diagnostic is the same detector's
endpoint source jets and their possible relative logarithmic residue;
more field refinement or finer PV quadrature would not settle that identity.

### Source-endpoint jets and a relative logarithmic diagnostic

For A=-iW/2 and B_plus=int A^2 dV, differentiation of the Volterra lower
limit gives `G_-+,s=A_U-A*B_plus` and
`G_++,s=int A*A_U dV-B_plus^2/2`. Exchange the two legs for the other
column. Thus the diagonal source jets are
`G_++,s=-J_V/4-I_V^2/32`, `G_--,s=-J_U/4-I_U^2/32`, distinct from the
observation-x derivatives already retained. The endpoint-jets component
now reports both derivatives from the same potential and integrals.

Let g_a, g_as and g_ax be a leg's endpoint value, source derivative and
observation derivative, with e_a the direct unit column. In the chiral
frame D0x_a=lambda_a*g_a, D1_a=lambda_a*e_a, lambda=(-1,+1). The
coefficient of log(1/abs(eta)) in q_alpha is

\[
L_q=\frac1{2\pi}\operatorname{Re}\operatorname{tr}\left[
\sigma_3\sum_a(D0x_a g_a^\dagger+D1_a g_{as}^\dagger-g_{ax}e_a^\dagger)\right].
\]

The direct/G Hilbert transform gives g_a*log(1/abs(eta)); its derivative
gives g_as*log(1/abs(eta)). Ordinary G_x*C0*G pairings have finite L2
limits and remain essential to the finite constant. Exact endpoint
identities reduce each leg's bracketed trace to W^2/2, giving
`L_q=W^2/(2*pi)`. The analogous angular insertion has `L_S=-W/pi`.
They are common to the two geometries. Merely fixing g_a and g_ax while
leaving g_as unconstrained need not cancel the relative stress logarithm.
This coefficient calculation does not yet supply a coincident finite part.

Control this identity against the independent constant-cylinder Fourier
asymptotics and the separable Riemann kernel's source derivative. A
zero-mean, zero-endpoint polynomial perturbation with a changed source
slope is a counterexample to endpoint-value sufficiency. Admit 30 CPU
seconds for these algebraic controls and one N128 central-detector
diagnosis in the two existing geometries, orders24/48. No new grid or
smaller separation is admitted. Print the small coefficient table and
retain its interpretation here; a separate large data record is unnecessary.

Three controls pass in 1.109375 process-CPU seconds: the independent
Fourier logarithms for both signed cylinder sectors, the separable
Riemann kernel's source derivative, and a zero-mean/zero-value polynomial
counterexample with a changed stress logarithm. The diagnosis route
`--diagnose-source-logs` evaluates only the same N128 central detector.
It reports the represented coefficients and those obtained by substituting
the exact source jets. The substitution changes no covariance or finite
stress; it isolates the missing identity before a reconstruction is changed.

The native Grok algebra audit confirms the lower-limit signs, the complete
leading-log coefficient and its reduction to the two universal local
values. Its scope is this coefficient, not the uncomputed finite part.
The first physical coefficient printout stopped on a scalar-shape bug:
`NullResponse.value` returned shape (1,2) rather than (2,) at a scalar
source coordinate. The pairing's outer products flattened it, so the
previous bilocal values are unaffected; indexing the diagnostic exposed it.
Correct the scalar contract and its control before repeating the diagnosis.
No coefficient table was obtained from that attempt. Charge its full
unreported 30-CPU allowance conservatively; the corrected fixed-domain
diagnosis has a fresh 30-CPU allowance after this implementation diagnosis.

The corrected diagnostic, method `8bad871`, completes at the same central
detector in 2.75 calculation CPU seconds / 1.57345 wall seconds; including
imports its wrapper records 3.59375 CPU seconds. The source fingerprint is
`1c4a40fbe52e2e14b3794d2d3726e716c2c5c9e53315da3e73b86d02e90a4357`.

| Order | Represented q-log, nested | Represented q-log, cap | Relative pressure-log coefficient |
|---|---|---|---|
| 24 | .445347440091 | .448849764909 | -.002002889525 |
| 48 | .484628661337 | .486418401135 | -.001023506179 |

The exact local q-log is .600799070856 in both geometries; the common
angular-log coefficient is -.618450133581. Substituting the exact source
jets cancels the relative coefficient to working arithmetic precision.
This confirms a missing cancellation in the finite representation. It
does not calculate the finite stress or prove every approximation error
has been identified.

The leading-log estimates for halving the separation are -.00138830
and -.00070944, whereas the observed v2 changes are -7.72801e-5 and
+7.18296e-5. The two existing splits therefore cannot be treated as a
pure leading-log asymptotic regime. No finite amplitude is extracted from
this diagnostic and no smaller split is scanned.

The next reconstruction must retain c0 and enforce, for each ordinary
field column on [0,L], both identities

\[
\sum_n c_n=G(L),\qquad
\frac1L\sum_n n(n+1)c_n=G_s(L).
\]

Their targets are the Volterra values and source derivatives above. A
minimum-L2 polynomial correction with these two constraints is fixed by
the equations before covariance evaluation. Recompute only after this
source compatibility is restored; finer field grids and further PV
quadrature do not repair the missing identity.

### Enforcing the Volterra source derivative before recomputation

For n=1,...,p-1 set w_n=2n+1 and d_n=n(n+1)/(p(p-1)). Holding c0 fixed,
the minimum-L2 correction has `delta c_n=w_n*(lambda0+d_n*lambda1)`.
The two multipliers solve the weighted 2x2 Gram system for the endpoint
value residual and the scaled source-derivative residual. This scaling
avoids conditioning the Gram by powers of p. Apply it to both spinor
components of the field kernel; the observation-spatial kernel keeps its
existing value/mean constraint. The old reconstruction remains available
explicitly for baseline comparisons.

Control mean preservation, both endpoint constraints and weighted
orthogonality against independent feasible variations. Check the full
cylinder kernel, analytic spatial insertion and universal logarithms on
the existing N64 grid at orders24/48. Admit 60 CPU seconds for these
controls, no new field grid and no retained arrays. A failed constraint
or cylinder control stops the physical calculation.

Then compare old and source-compatible reconstructions from the SAME
N128 response rows, state, clock, orders24/48 and splits .02S/.01S.
The old cases must reproduce v2 within 1e-9. Use the existing separated
eight-worker field phase and eight-process covariance phase, 240 total
process-CPU seconds including children, and a summary below 24000 bytes.
The corrected logarithmic identity must hold before a new amplitude is
interpreted. Noncontracting polynomial or split effects keep the result
unresolved. This correction does not by itself establish the finite
coincidence part, an angular sum, conservation at finite projection or
full Psi/Phi. No P001 grid or smaller separation is admitted.

Both controls pass in 11.75 process-CPU seconds. Weighted optimality is
checked against an independently constructed SVD nullspace of the three
constraints. On the unchanged N64 cylinder, the corrected field at orders
24/48 gives kernel errors 1.43270e-8/1.37257e-8 and spatial-insertion
errors 5.18031e-7/5.15725e-7, with the universal logarithmic coefficients.
The producer's explicit `--stage source-jets` compares both reconstructions;
legacy modes retain their original defaults and serve as numerical baselines.

### Source-compatible comparison: cancellation fixed, amplitude still open

The [v3 summary](../results/summaries/nsc-atlas-dirac-root-response-v3.json),
method `a097d42`, completes the admitted reconstruction comparison. All four
legacy values reproduce v2 exactly. The source-compatible maps satisfy the
universal endpoint logarithms at both polynomial orders and both splits;
their largest absolute coefficient discrepancy is below 4.5e-16. This
verifies the imposed algebraic identity, not a continuum error estimate.

| Order | Legacy, eta/S=.02 | Source-compatible, eta/S=.02 | Legacy, eta/S=.01 | Source-compatible, eta/S=.01 |
|---|---|---|---|---|
| 24 | .001842679064 | .002277828221 | .001765398991 | .001601037456 |
| 48 | .002173544722 | .002196411488 | .002245374289 | .002187989727 |

At order 48 the two source-compatible splits differ by 8.42176e-6,
about .385% of the smaller-split value. The 24/48 difference at the
smaller split remains 5.86952e-4, about 26.8%; the representation is
therefore not converged. The larger-split polynomial difference is
8.14167e-5. No finite local vacuum amplitude, angular sum or Psi/Phi is
inferred from these values.

The run uses all eight covariance processes: 92.296875 total CPU seconds,
including 83 child CPU seconds, and 17.42970 wall seconds (5.30 effective
cores over both phases). It retains 20761 bytes of numbers and hashes,
with no field arrays. No additional field resolution or smaller split ran.

### Observation derivative compatibility after constrained projection

Native Grok confirms the weighted projection and the source-chain rule.
For a moving source endpoint L_x=lambda and gamma(x)=G(x,L(x)),

\[
G_x(L)=\partial_x\gamma-\lambda G_s(L).
\]

The imposed field value/source jet and the separate spatial endpoint
therefore agree locally. This does not make their interiors identical:
differentiating the constrained field projection need not equal projecting
the physical observation derivative separately. A difference vanishing
at the endpoint can still enter finite Hilbert and G*C0*G insertions.
Log cancellation alone does not control that finite contribution.

One necessary check uses the maps at each existing split, without a
smaller separation. For I=L*c0 of the field kernel, J=L*h0 of its
spatial kernel, and side0/side1 at x+eta/2,x-eta/2, the mean residual is

\[
\mu=\frac{I_0-I_1}{\eta}
-\lambda\frac{\gamma_0+\gamma_1}{2}
-\frac{J_0+J_1}{2}.
\]

For an exact differentiated kernel this has the ordinary O(eta^2)
central-difference error. This mean check is necessary, not sufficient.
For example, the mean/value/source-slope preserving quadratic projection
of G(s)=s^3 on [0,L] is
`h=L^3/4-3*L^2*s/2+9*L*s^2/4`. Although G has no observation dependence,
`partial_L h=3*L^2/4-3*L*s+9*s^2/4` is nonzero. Its mean and endpoint
value vanish, while `int (partial_L h)/(s-L) ds=3*L^2/8` is nonzero.
Thus checking only the mean cannot bound the finite insertion. Nonconstant
source moments or the actual derivative commutator are the next required
compatibility evidence before a coincident amplitude is evaluated.
