# Scalar response to the computed metric

The [exterior and wall response](nsc-atlas-wall-evolution.md) and
[smooth cap](nsc-atlas-cap-response.md) supply the formal quadratic metric,
wall displacement and clock embeddings. This note derives the next scalar
coefficient on that geometry. It extends the same
[action and junction laws](nsc-atlas-cap-wall-response.md), with the
[analytic admission](prior-art-and-open-claim.md#scalar-response-to-the-computed-metric-15-october-2026)
recorded before the calculation. The equations below define a coupled linear
inhomogeneous transmission problem. They have not yet been numerically evolved.

The continuation chosen here keeps the prescribed reduced-field input
$\chi(U_0,V)=\varepsilon g(V)$ exactly in the original anchored source gauge.
This is a declared continuation of the preparation, not a derivation of a
physical source mechanism. Holding $\phi$ instead would define a different
next-order input despite giving the same leading field.

## Order, frames and the field being propagated

All coordinates and radii below are dimensionless; curvature denotes
$\mathcal R=\ell^2R_4$. Keep
$\alpha=\kappa_4\varepsilon^2/\ell^2$, $h=H\ell=1$ in the cap,
$h_+=0$ outside, $b=k\ell=2$ and $m=M/\ell=0.4$.
The Jordan scalar and invariant reduced field are expanded as

$$
\frac{\ell\phi_J}{\varepsilon}=p+\alpha p_{1,J}+O(\alpha^2),
\qquad c=\frac{\chi}{\varepsilon}=c_0+\alpha c_1+O(\alpha^2),
\qquad c_0=rp.
$$

The coefficient $c_1$ depends on the coordinate identification, which is
specified for each equation below. Introduce

$$
\mathcal F=1-\frac{\kappa_4\phi_J^2}{6}
=1-\frac{\alpha p^2}{6}+O(\alpha^2),\qquad
g_E=\mathcal Fg_J,\qquad r_E=\sqrt{\mathcal F}\,r_J.
$$

Conformal covariance uses the field
$\phi_{\mathrm{conf},E}=\mathcal F^{-1/2}\phi_J$. Consequently

$$
\chi=\ell r_J\phi_J=\ell r_E\phi_{\mathrm{conf},E}
$$

is invariant on the $\mathcal F>0$ branch, and $\phi_{\mathrm{conf},E}$ obeys
$(\Box_E-R_E/6)\phi_{\mathrm{conf},E}=0$. It is different from the
canonical Einstein scalar $\psi$, for which

$$
\psi=\sqrt{\frac6{\kappa_4}}
\operatorname{atanh}\!\left(\sqrt{\frac{\kappa_4}6}\phi_J\right)
=\phi_J+\frac{\kappa_4\phi_J^3}{18}+O(\phi_J^5),
$$
$$
\phi_{\mathrm{conf},E}
=\phi_J+\frac{\kappa_4\phi_J^3}{12}+O(\phi_J^5).
$$

These fields agree at leading order, which is sufficient for the existing
quadratic Einstein constraints. Their cubic difference matters here.

In each smooth Jordan region the on-shell improved trace gives
$\mathcal R_J=12h^2$ exactly. Thus a direct Jordan calculation has

$$
(\Box_0-2h^2)p_{1,J}=-\delta\Box_J\,p,
$$

with no independently adjustable curvature perturbation. In the Einstein
frame the same statement takes the form

$$
\mathcal R_E=12h^2+\alpha\mathcal R_{E,1}+O(\alpha^2),\qquad
\mathcal R_{E,1}=(\nabla_0p)^2+4h^2p^2.
$$

The $4h^2p^2$ term follows from the transformed fixed cap vacuum potential.
The distributional wall contribution is handled by the junction conditions,
separately from these smooth bulk equations. The correction $c_1$ changes
the gravitational source at order $\alpha^2$, so the already computed
order-$\alpha$ metric is sufficient to determine this scalar coefficient.

## Exterior equation and its source pullback

Retain the original exterior null label $V$ and use Einstein areal radius
$r$. The metric and its coefficients are

$$
\frac{ds_E^2}{\ell^2}=-Z^2A\,dV^2-2Z\,dr\,dV+r^2d\Omega^2,
\qquad A=A_0-\frac{2\alpha\mu_E}{r}+O(\alpha^2),
$$
$$
A_0=1-\frac{2m}{r},\qquad Z_0=-\frac{2a}{V}>0,\qquad a=2m,
\qquad Z=Z_0(1+\alpha\zeta_E)+O(\alpha^2).
$$

The negative $dr\,dV$ cross term follows from
$Z=C/(2r_U)$ in the retained sector $r_{0,U}>0$. Set

$$
P=ZA=P_0+\alpha P_1+O(\alpha^2),\qquad
P_0=Z_0A_0,\qquad
P_1=Z_0\left(A_0\zeta_E-\frac{2\mu_E}{r}\right).
$$

Direct spherical reduction gives the exact conformal equation

$$
-2c_{Vr}+\partial_r(Pc_r)
-\left(\frac{P_r}{r}+\frac{Z\mathcal R_E}{6}\right)c=0.
$$

Define $L_+c=-2c_{Vr}+\partial_r(P_0c_r)-P_{0,r}c/r$. Its correction is

$$
\boxed{L_+c_1=-\partial_r(P_1c_{0,r})
+\left(\frac{P_{1,r}}r+\frac{Z_0\mathcal R_{E,1}}6\right)c_0},
\qquad
\mathcal R_{E,1}=A_0p_r^2-\frac{2p_Vp_r}{Z_0}.
$$

Here $p_V$ is at fixed $r$. The existing null-coordinate jets convert by

$$
p_r=\frac{p_U}{r_{0,U}},\qquad
p_V\big|_r=p_V\big|_U-\frac{r_{0,V}}{r_{0,U}}p_U.
$$

The [radial metric constraints](nsc-atlas-wall-evolution.md) imply

$$
\mu_{E,r}=\frac{r^2A_0}{4}p_r^2,\qquad
\zeta_{E,r}=\frac r2p_r^2,\qquad
P_{1,r}=Z_0\left(A_{0,r}\zeta_E+\frac{2\mu_E}{r^2}\right).
$$

The last cancellation removes differentiated mass/lapse data from the
forcing. The divergence term can be treated as a known flux. A strong-form
implementation instead needs $c_{0,rr}$ from a continuum derivative.

The source curve in these coordinates moves:
$r_{E,\mathrm{src}}(V)=r_{0,\mathrm{src}}(V)+\alpha H_{E,\mathrm{src}}(V)$.
Exact input in the original source gauge therefore gives

$$
\boxed{c_1(r_{0,\mathrm{src}}(V),V)
=-H_{E,\mathrm{src}}(V)c_{0,r}(r_{0,\mathrm{src}}(V),V)},
\qquad H_{E,\mathrm{src}}=h_J-\frac{r_{0,\mathrm{src}}p^2}{12}.
$$

The existing [source traces](nsc-atlas-bulk-response.md#source-data-and-the-unchanged-affine-normalization)
fix this coefficient in the original source coordinate, including its
quiet integration constants. The correction is zero on the quiet bottom
$V=V_a$ and compatible at the source corner. Imposing zero correction on
the background source curve in areal coordinates would change the declared
input. No new incident waveform is introduced.

### Conservative exterior forcing on the existing wedge

Pull back the areal coefficient to the reference exterior coordinates:
$y(U,V)=c_1(r_0(U,V),V)$. Here $r_{0,V}=-P_0/2$,
$r_{0,U}>0$ and $C_0=2Z_0r_{0,U}$. The equation becomes

$$
y_{UV}+K_0y=\partial_U F-G,\qquad
K_0=\frac{C_0m}{2r_0^3},\qquad
F=\frac{Z_0}{2}\left(A_0\zeta_E-\frac{2\mu_E}{r_0}\right)
\frac{c_{0,U}}{r_{0,U}},
$$
$$
G=\left[\frac{C_0(m\zeta_E+\mu_E)}{2r_0^3}
-\frac{p_Up_V}{6}\right]c_0.
$$

The derivatives in the last expression are in the original null chart:
$\mathcal R_{E,1}=-4p_Up_V/C_0$. This formulation uses first scalar
jets and a known flux. At each fixed $V$, the radial mass and lapse can
be obtained on the active grid by one cumulative quadrature, then reused
in $F$ and $G$. Integrating separately from the source to every grid point
would repeat the same work.

Let $g_1(V)=-H_{E,\mathrm{src}}c_{0,r}|_{\mathrm{src}}$ and let
$\mathfrak b(V)=r_{0,U}c_{1,r}|_W$ be the boundary derivative derived below.
For $L(U)=U_W^{-1}(U)$ on newly exposed wall coordinates, define

$$
B_{\mathrm{lift}}(U)=\begin{cases}
0,&U\le U_a,\\
\displaystyle\int_{U_a}^{U}\mathfrak b(L(q))\,dq,&U>U_a.
\end{cases}
$$

Then $z=y-g_1(V)-B_{\mathrm{lift}}(U)$ has zero source and bottom values and
$z_U|_W=0$. Its evolution is

$$
z_V=F(U,V)-F(U_0,V)
-\int_{U_0}^{U}\left[K_0(z+g_1+B_{\mathrm{lift}})+G\right](q,V)\,dq.
$$

This has the existing clipped Volterra geometry. Grid values beyond the
moving wall are bookkeeping continuations of the wall value of $z$;
their reconstructed $B_{\mathrm{lift}}$ is evaluated at the clipped coordinate as well.
The lift removes a moving-wall derivative term without differentiating
$\mathfrak b$. It introduces no adjustable boundary datum. At the wall $F_W=0$
because $c_{0,U}|_W=0$, whereas $p_U|_W$ is generally nonzero. The
wall readout and its continuum derivative are

$$
y_W=z_W+g_1+B_{\mathrm{lift}}(U_W),\qquad
\frac{dy_W}{dV}=g_1'+U_W'\mathfrak b-F(U_0,V)
-\int_{U_0}^{U_W}(K_0y+G)\,dU,
$$
$$
g_1'=-\left(H_{E,\mathrm{src}}'
+\frac{H_{E,\mathrm{src}}P_{0,r}}2\right)c_{0,r}
+\frac{H_{E,\mathrm{src}}K_0c_0}{r_{0,U}}
\quad\text{on the source.}
$$

The last identity follows by differentiating the source pullback and using
$c_{0,UV}=-K_0c_0$. It avoids differences of stored waveforms. With
$K_0=F=G=0$, the manufactured solution $y=g_1+B_{\mathrm{lift}}$ checks the wall
transport term exactly. A numerical boundary interpolant and an exterior
evolution have not yet been computed.

## Cap equation in the central proper-clock gauge

Use the computed polar-areal cap metric

$$
\frac{ds_E^2}{\ell^2}=-e^{2\alpha\delta_E}A\,dT^2
+A^{-1}dr^2+r^2d\Omega^2,\qquad
A=A_0-\frac{2\alpha\nu_E}{r}+O(\alpha^2),\qquad A_0=1-h^2r^2.
$$

The time $T$ is Jordan central proper time through this order:
$\delta_E(T,0)=-p_0(T)^2/12$, where $p_0=p(T,0)$.
Put $\xi=r_*(r)=\operatorname{atanh}(hr)/h$ and

$$
B=e^{\alpha\delta_E}A=A_0+\alpha B_1+O(\alpha^2),\qquad
B_1=A_0\delta_E-\frac{2\nu_E}r,\qquad
D=\frac{B_1}{A_0}=\delta_E-\frac{2\nu_E}{rA_0}.
$$

The exact reduced equation is

$$
-\partial_T(B^{-1}c_T)+\partial_r(Bc_r)
-\left(\frac{B_r}r+\frac{e^{\alpha\delta_E}\mathcal R_E}{6}\right)c=0.
$$

The reference potential cancels because $A_{0,r}/r+2h^2=0$. Using the
cap mass and lapse constraints, its first-order potential reduces to

$$
W_1=\frac{B_{1,r}}r+2h^2\delta_E+\frac{\mathcal R_{E,1}}6
=\frac{2\nu_E}{r^3}-\frac{h^2p^2}{3}
-\frac{p_T^2}{6A_0}+\frac{A_0p_r^2}{6}.
$$

The cancellation uses both the cap scalar potential and its kinetic
contribution to the mass. With
$L_-c=-c_{TT}/A_0+\partial_r(A_0c_r)$, the forced equations are

$$
L_-c_1=-\partial_T\!\left(\frac{B_1}{A_0^2}c_{0,T}\right)
-\partial_r(B_1c_{0,r})+W_1c_0,
$$
$$
\boxed{(-\partial_T^2+\partial_\xi^2)c_1
=-\partial_T(Dc_{0,T})-\partial_\xi(Dc_{0,\xi})+A_0W_1c_0}.
$$

This conservative form uses the existing radial metric and leading scalar
jets. It permits integration of the forcing fluxes without differentiating
stored metric or waveform samples.

## Central-clock null coordinates

A further coordinate transformation removes the differentiated forcing.
Write $q_-=T-\xi$, $q_+=T+\xi$, so
$c_0=f(q_+)-f(q_-)$, retaining the reflected term. Define

$$
\lambda_-(T,\xi)=\int_0^\xi D(q_-+y,r(y))\,dy,\qquad
\lambda_+(T,\xi)=-\int_0^\xi D(q_+-y,r(y))\,dy.
$$

Then $\partial_+\lambda_-=D/2$ and $\partial_-\lambda_+=D/2$.
The coordinates $\mathsf q_\pm=q_\pm+\alpha\lambda_\pm$ are null
through this order and equal $T$ at the centre. The incoming correction
$\lambda_+$ is the existing fixed-event null travel integral. The outgoing
correction $\lambda_-$ is a second specified geometric integral. Portions
before the quiet initial slice use its zero field/metric extension, not
new initial data.

Expressing the reduced field in these coordinates gives

$$
c_{1,\mathrm{static}}=\widetilde c_1
+\lambda_-\partial_-c_0+\lambda_+\partial_+c_0,
\qquad
\boxed{\partial_-\partial_+\widetilde c_1=-K_1c_0},
\qquad K_1=\frac{A_0W_1}{4}.
$$

All coefficients on the right are evaluated in the reference coordinates.
The vanishing background cap potential is why no $K_0\widetilde c_1$ or
coordinate variation of $K_0$ remains. This is a forced characteristic
equation with regular reflection at the centre.

Regularity gives

$$
\nu_E=r^3\left(\frac{p_{T,0}^2}{12}
+\frac{h^2p_0^2}{6}\right)+O(r^5),\qquad K_1=O(r^2).
$$

In particular, the terms displayed separately in $W_1$ cancel at the
centre. Impose $c_1(T,0)=0$, equivalently
$\widetilde c_1(T,T)=0$, with a regular odd reduced field. The initial cap
slice $T=T_a$ is quiet, with $c_1=c_{1,T}=0$.

## Transmission at the displaced wall

Use the original signed cap-to-exterior normal and $[X]=X_+-X_-$.
The exact conformal junction has

$$
[\phi_J]=0,\qquad [n_J\phi_J]=k\phi_J,\qquad
[n_Jr_J]=-kr_J.
$$

Here the normal derivative is in physical units. The dimensionless normal
components below use $\ell n$ and $s=\tau_W/\ell$; multiplying these
junctions by $\ell$ replaces $k$ by $b=k\ell$.

It follows that $[\chi]=[n_J\chi]=0$. The conformal factor is continuous,
so $[n_E\chi]=0$ as well. These exact reduced-field laws specify the
moving transmission problem.

Keep $s$ as Jordan proper wall time and $v=x_s$. At fixed $s$ the
Einstein wall radius, exterior label and cap time are

$$
r_{E,W}=x+\alpha e,\quad e=\rho-\frac{xp_W^2}{12},\qquad
V_W=V_0+\alpha V_1,\qquad T_W=T_0+\alpha T_1.
$$

Their displacement vectors in the two bulk charts are
$\Xi_+=(V_1,e)$ and $\Xi_-=(T_1,e)$. Scalar continuity becomes

$$
c_{1,+}+V_1c_{0,+,V}+e c_{0,+,r}
=c_{1,-}+T_1c_{0,-,T}+e c_{0,-,r},
$$

with all terms evaluated at the reference wall. In central-clock cap null
coordinates the wall label shifts are

$$
\Delta\mathsf q_-=T_1-e/A_0+\lambda_{-,W},\qquad
\Delta\mathsf q_+=T_1+e/A_0+\lambda_{+,W}.
$$

The latter is the previously computed central-arrival correction at fixed
$s$. It supplies a boundary pullback, not the complete scalar correction.

Let $\delta n_{E,\mathrm{mat}}^a$ denote the variation of the normal's
components **at fixed material time $s$**. On either side, the correction
to its scalar normal trace is

$$
N_1=n_0^a\partial_a c_1
+\delta n_{E,\mathrm{mat}}^a\partial_a c_0
+n_0^a\Xi^b\partial_b\partial_a c_0,\qquad [N_1]=0.
$$

The following material variations already include the changed evaluation
point of the normal. Adding another $\Xi\cdot\partial n_0$ would count it
twice. The exterior components are

$$
n_{0,+}^V=V_{0,s},\quad n_{0,+}^r=\beta_+,\qquad
\delta n_{E,+}^V=V_{1,s}+\frac{p_W^2}{12}V_{0,s},\quad
\delta n_{E,+}^r=\delta\beta_+-\frac{xp_WB_+}{6}.
$$

The sign $n_{0,+}^V=+V_{0,s}$ agrees with the negative exterior cross term.
For the cap,

$$
n_{0,-}^T=\frac v{A_0},\quad n_{0,-}^r=\beta_-,\qquad
\delta n_{E,-}^T=
\frac{e_s+vp_W^2/12-v\delta_{E,W}}{A_0}
-\frac{vA_1}{A_0^2},\quad
\delta n_{E,-}^r=\delta\beta_--\frac{xp_WB_-}{6},
$$
$$
A_1=-2h^2xe-\frac{2\nu_{E,W}}x,\qquad
e_s=\rho_s-\frac{vp_W^2+2xp_WA_W}{12}.
$$

Here $A_0=A_0(x)$, $A_W=dp_W/ds$, and $B_\pm=n_{0,\pm}p$ are
the existing one-sided leading derivatives. The quantities $\rho$,
$\delta\beta_\pm$, $V_1$ and $T_1$ retain their original quiet anchors.
These formulas divide by $A_0$ in the regular static cap, not by $\beta_+$;
the exterior zero of $\beta_+$ creates no singular transmission condition.

The leading no-return solution has $c_{0,U}=0$ at the wall. The changed cap
metric generally has $K_1\ne0$ and can scatter a field back toward the wall
before any centre-reflected return. The correction must therefore use both
transmission equations; setting $c_{1,U}=0$ would remove a response the
inhomogeneous cap equation can generate. No independent outgoing waveform
is supplied in its place.

For an explicit analytic example, the static regular profile
$p=2\operatorname{atanh}(r)/r$ at $h=1$ solves the leading cap equation.
Its regular mass integral gives

$$
K_1=\frac{4r^2}{135}+\frac{4r^4}{189}+O(r^6),\qquad
-K_1c_0=-\frac{8r^3}{135}+O(r^5).
$$

This manufactured profile demonstrates nonzero forcing; it is not a new
physical preparation or a numerical result for the prescribed pulse.

## A triangular solution on the retained no-return interval

The transmission problem need not be solved by iteration in the current
finite domain. In null charts, differentiating scalar continuity along the
wall and combining it with normal continuity gives exactly

$$
U_s c_U=(\mathsf q_-)_s c_{\mathsf q_-},\qquad
V_s c_V=(\mathsf q_+)_s c_{\mathsf q_+}.
$$

The common wall conformal factor cancels. The established no-return
certificate has $q_{-,W}<q_a$, so the outgoing homogeneous correction at
the corresponding earlier centre event is quiet. Its wall value is therefore
the known integral

$$
J_-(s)=\partial_-\widetilde c_1(q_{-,W},q_{+,W})
=-\int_{q_{-,W}}^{q_{+,W}}
K_1(q_{-,W},q)c_0(q_{-,W},q)\,dq.
$$

The part with $q<q_a$ vanishes. This integral uses only the computed
quadratic metric and leading field. It is independent of the unknown
incoming correction, whose derivative with respect to $q_-$ vanishes.

Let $d=v-\beta_+>0$ and
$q_{-,0,s}=(\beta_--v)/A_0(x)$. In exterior areal coordinates the first
flux identity becomes

$$
\frac d2\left(c_{1,r}+V_1c_{0,Vr}+e c_{0,rr}\right)_W
=q_{-,0,s}J_-.
$$

The leading $c_{0,r}=0$ removes variations of this flux prefactor.
Differentiating that zero along the reference wall and using the leading
exterior equation gives, without a Hessian query,

$$
c_{0,Vr}=-\frac v{V_{0,s}}c_{0,rr},\qquad
c_{0,rr}=\frac{2m w}{x^3d^2}.
$$

Thus the explicit exterior correction boundary condition is

$$
\boxed{c_{1,r}\big|_W=
\frac{2q_{-,0,s}}d J_-
-\left(e-\frac{vV_1}{V_{0,s}}\right)\frac{2m w}{x^3d^2}}.
$$

After the exterior forced equation is solved with this boundary condition,
its wall trace fixes the cap incoming correction. Write $q=q_{+,W}$ and
use the existing arrival correction at fixed exterior $V$:

$$
\boxed{\widetilde c_{1,W}=c_{1,+,W}
-\Delta C_{\mathrm{fixed}\,V}\,f'(q)},\qquad
\Delta C_{\mathrm{fixed}\,V}
=\Delta\mathsf q_+-\frac{q_{+,0,s}}{V_{0,s}}V_1.
$$

Here $c_{1,+,W}$ is the exterior areal-gauge coefficient at the reference
wall; its leading radial derivative vanishes. This ordering retains cap
backscatter while determining it before the exterior correction is solved.
Once a centre-reflected correction can return to the wall, the quiet
outgoing argument no longer applies and the full transmission problem
must be revisited on its continued domain.

### Evaluating the outgoing boundary input

For the retained wall cuts set $u=q_W-2\mathrm{atanh}(x)$. Along the
active part of the outgoing ray,

$$
J_-=-\int_{q_a}^{q_W}K_1(T,r)f(q)\,dq,\qquad
T=\frac{u+q}{2},\qquad r=\tanh\frac{q-u}{2}.
$$

These points have $T\le T_W\le T_b<q_a$ and $u<q_a$. Thus their
reflected profile is quiet, and the metric queries remain within the
existing cap. Substitution of $p=f(q)/r$, $p_T=f'(q)/r$ and
$p_r=f'(q)/(rA_0)-f(q)/r^2$ into the general potential cancels the
terms quadratic in $f'$:

$$
K_1=\frac{A_0\nu_E}{2r^3}
-\frac{A_0ff'}{12r^3}
+\frac{A_0(1-3r^2)f^2}{24r^4}.
$$

This gives an independent expression to check against the general
first-jet potential. A second quadrature uses $q=u+2\mathrm{atanh}(r)$
and $dq=2\,dr/A_0$. Both integrals split at the prescribed profile edges.
The [boundary admission](prior-art-and-open-claim.md#scalar-backscatter-boundary-admission-16-october-2026)
fixes the comparison domain, settings and shared compute budget. Its
numerical outputs belong in the result record; sampled values alone
do not supply the boundary interpolant required for an exterior solve.

The implementation is `ScalarBoundary` in the
[boundary adapter](../src/recursive_horizons/nsc_atlas_scalar_boundary.py),
with a [comparison producer](../scripts/record_nsc_atlas_scalar_boundary.py).
The [boundary tests](../tests/test_nsc_atlas_scalar_boundary.py) check
manufactured integration and transmission; the
[flux tests](../tests/test_nsc_atlas_scalar_flux.py) independently check
the conservative exterior and clock-derivative identities. Compact
numerical outcomes and their reproduction commands are maintained in
the [current result](current-result.md) and [reproduction guide](reproducing.md).

## Central waveform and the remaining numerical connection

In the static cap coordinates, reconstruction at fixed Jordan central
proper time gives

$$
\frac{\ell\phi_C(T)}{\varepsilon}
=p_0(T)+\alpha\left[\partial_r c_{1,\mathrm{static}}(T,0)
-\frac{p_0(T)^3}{12}\right]+O(\alpha^2),\qquad p_0=2f'(T).
$$

At the centre $D_0=-p_0^2/12$, and the null-coordinate pullback gives
$\partial_r c_{1,\mathrm{static}}
=(\partial_+-\partial_-)\widetilde c_1+p_0^3/12$. Thus the same
Jordan waveform in central-clock null coordinates is

$$
\boxed{\frac{\ell\phi_C(T)}{\varepsilon}
=p_0(T)+\alpha(\partial_+-\partial_-)\widetilde c_1(T,T)
+O(\alpha^2)}.
$$

The frame factor and central lapse normalization cancel in this expression.
A solution evaluated at fixed $T$ already includes its propagation-time
change. Adding $-\Delta C\,p_0'(T)$ afterwards would count that change
twice. A waveform first parametrized by emission label instead needs its
coordinate pullback exactly once. This observable is the local classical
scalar field; a detector transition probability also needs a specified
detector coupling and state.

The known exterior field, radial mass/lapse and cap queries supply the
forcing. The existing
[wedge evaluator](../src/recursive_horizons/nsc_atlas_wedge_field.py),
[wall solution](../src/recursive_horizons/nsc_atlas_wall_evolution.py),
[cap field](../src/recursive_horizons/nsc_atlas_cap_field.py) and
[cap metric](../src/recursive_horizons/nsc_atlas_cap_response.py) remain
their computational owners. A new characteristic or conservative solver
must propagate $c_1$ with the derived exterior boundary forcing, then the
cap response, and check both moving-wall laws. Source/bottom data and the
regular centre stay fixed. The finite no-return certificate makes this
ordering possible; it does not remove the cap forcing.

If a strong-form exterior forcing or the explicit wall Hessian is used,
differentiate the continuum equation. With the existing lower limit
$L(U)=V_a$ for $U\le U_a$ and $L(U)=U_W^{-1}(U)$ otherwise,

$$
c_{0,UU}=L'(U)K(U,L)c_0(U,L)
-\int_{L(U)}^V\left(K_Uc_0+Kc_{0,U}\right)dq.
$$

The moving-lower-limit term is essential. At the reference wall,

$$
c_{0,UU}=K_Wc_{0,W}/U'_W,\qquad
c_{0,UV}=-K_Wc_{0,W},\qquad
c_{0,VV}=w_{VV}+U'_WK_Wc_{0,W}.
$$

The cap field interface already obtains $w_{VV}$ and $f''$ from the
continuum integral equation. No stored waveform is differentiated. The
cap central-null equation avoids new differentiated metric data; its
outgoing-label integral is an additional geometric query, not a new source.

For the later derivative of the incoming clock pullback, the optical
constraint and mass-flux law also give a first-jet formula:

$$
D_T(T,r)=-\frac{p_0p_{0,T}}6
+\int_0^r\left[
\frac{(1-3h^2s^2)p_Tp_r}{1-h^2s^2}
-\frac{2h^2spp_T}{1-h^2s^2}\right]ds.
$$

It follows from
$D_r=2\nu_E(1-3h^2r^2)/(r^2A_0^2)-h^2rp^2/A_0$ and
$\nu_{E,T}=r^2A_0p_Tp_r/2$, with the central clock anchor.
It specifies a route to differentiating the ray integral without
differentiating saved clock samples; that numerical derivative is not
part of the boundary calculation above.

The finite calculation uses the retained wedge and cap with
$q_+\le q_b$, including its established quiet earlier data. Initial
fixed-central-clock comparisons should lie on a compact subinterval of
$(q_a,q_b)$. The physical final arrival moves from $q_b$ by an
order-$\alpha$ amount. An endpoint coefficient can be stated as a one-sided
formal limit or at fixed terminal emission label; a finite-amplitude
fixed-clock endpoint assertion requires the corresponding continued domain.
The prepared tail is not set to zero at $q_b$.

The next numerical evidence must check the source and wall pullbacks,
both transmission residuals, regular-centre cancellation and agreement of
the two central waveform reconstructions under independent refinement.
The [independent analytic tests](../tests/test_nsc_atlas_scalar_metric_response.py)
derive the spherical operators from four-dimensional divergence forms,
check their perturbations and coordinate transformations, and evaluate
the manufactured cap example. They do not evolve the prepared scalar.
The analytic equations do not provide a computed $c_1$, a nonzero amplitude
range, a remainder bound, or formation and stability results. Those remain
subsequent calculations in the [atlas programme](nsc-cosmological-atlas.md).
