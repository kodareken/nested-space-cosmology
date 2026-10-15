# Smooth exterior metric response to the prepared field

The [cap and wall reduction](nsc-atlas-cap-wall-response.md) fixes conditional
mass and displacement functionals. This calculation propagates the smooth
exterior radius and lapse coefficients from their existing characteristic
data. Its first domain is the rectangle

$$
\mathcal D_0=[U_0,U_a]\times[V_a,V_b]
$$

inside the retained exterior wedge. Both past null sides of this rectangle
are already specified. The later strip $U>U_a$ instead meets the timelike
wall in its past; extending the quiet bottom through the excised cap would
supply different data. The continuation described below is an analytic
route, not a numerical result for that strip.

The [admission](prior-art-and-open-claim.md#smooth-exterior-metric-propagation-admission-15-october-2026)
fixes the original source, corner, finite domain and compute limits. The
borrowed equations and frame map are those already checked in the
[two-sided owner](nsc-atlas-cap-wall-response.md). The leading scalar remains
the [prepared response](nsc-atlas-prepared-response.md), with its
[characteristic initialization](nsc-atlas-characteristic-initialization.md)
and [source mass](nsc-atlas-source-mass.md). No new scalar preparation or
finite physical amplitude is selected.

The [bulk implementation](../src/recursive_horizons/nsc_atlas_bulk_response.py)
provides `source_traces`, `solve_bulk_rectangle`, and bounded
`BulkSolution.sample`/`summary` queries. Its
[focused tests](../tests/test_nsc_atlas_bulk_response.py) check the source,
equation and domain contracts. The
[prepared-field evaluator](../src/recursive_horizons/nsc_atlas_prepared_response.py)
supplies the scalar coefficients under the shared finite compute budget.

## Metric coefficients in the original null coordinates

Use the dimensionless background areal radius $r_0=R_0/\ell$,
$a=2m=0.8$, and the original Kruskal labels:

$$
UV=(1-r_0/a)e^{r_0/a},\qquad
C_0=\frac{4a^3}{r_0}e^{-r_0/a},\qquad
ds_0^2/\ell^2=-C_0\,dU\,dV+r_0^2d\Omega^2.
$$

In this retained white-hole sector $C_0,r_{0,U},r_{0,V}$ are positive.
The known leading scalar is $\phi=(\varepsilon/\ell)p+O(\varepsilon^3)$,
where $p=\chi/(\varepsilon r_0)$. Its first derivatives are obtained from
the continuum characteristic integral equations, not differences of saved
waveform samples.

Set $\alpha=\kappa_4\varepsilon^2/\ell^2$ and write the smooth
Einstein-frame response at fixed original $(U,V)$ as

$$
r_E=r_0+\alpha H+O(\alpha^2),\qquad
\log C_E=\log C_0+\alpha c+O(\alpha^2).
$$

The exterior vacuum potential is zero. Expanding the spherical Einstein
equations gives

$$
H_{UV}=-\frac{r_{0,U}H_V+r_{0,V}H_U}{r_0}
+\left(\frac{r_{0,U}r_{0,V}}{r_0^2}
+\frac{C_0}{4r_0^2}\right)H-\frac{C_0}{4r_0}c,
$$
$$
c_{UV}=\frac{2(r_{0,U}H_V+r_{0,V}H_U)}{r_0^2}
-\left(\frac{4r_{0,U}r_{0,V}}{r_0^3}
+\frac{C_0}{r_0^3}\right)H+\frac{C_0}{2r_0^2}c-p_Up_V.
$$

For evolution use

$$
y=r_0H,\qquad z=c+H/r_0.
$$

The background mixed radius equation cancels the first derivatives:

$$
\boxed{y_{UV}=\frac{C_0}{4r_0^2}y-\frac{C_0}{4}z,}
\qquad
\boxed{z_{UV}=-\frac{3C_0}{4r_0^4}y
+\frac{C_0}{4r_0^2}z-p_Up_V.}
$$

This is a forced linear characteristic system. If its right sides are
$F_y,F_z$, its rectangle formulation is

$$
y_V(U,V)=y_{0,V}(V)+\int_{U_0}^{U}F_y(q,V)\,dq,
\qquad
z_V(U,V)=z_{0,V}(V)+\int_{U_0}^{U}F_z(q,V)\,dq.
$$

The zero bottom initializes $y=z=0$. Auxiliary derivatives obey
$(y_U)_V=F_y$ and $(z_U)_V=F_z$ with zero bottom values. They supply
metric gradients without spatial differentiation of the output arrays.

## Source data and the unchanged affine normalization

On $U=U_0$, let $x=r_0(U_0,V)$ and use the existing Jordan source radius
coefficient $h_J(x)$, $I_0(x)=\int_{0.4}^{x}p^2dt$, and analytic source
jet. The fixed original Jordan lapse means $c_J=0$ there. Consequently

$$
H_0=h_J-\frac{xp^2}{12},\qquad c_0=-\frac{p^2}{6},
$$
$$
\boxed{y_0=xh_J-\frac{x^2p^2}{12},\qquad
z_0=\frac{h_J}{x}-\frac{p^2}{4}.}
$$

All source derivatives follow by the analytic $x(V)$ chain rule and the
existing $p,p_x,p_{xx},h_J,h_{J,x},h_{J,xx}$ formulas. The quiet source
corner makes these data agree with the zero bottom. No new mass, expansion
or lapse integration constant is introduced.

The initialization API's `h_einstein` is expressed in the transformed
Einstein affine parameter. It is **not** $H_0$ at the original fixed
coordinate. Here

$$
H_0=-\frac12\int_{0.4}^{x}(x-t)t p_x(t)^2dt-\frac16I_0(x).
$$

Dropping the last term would silently change the source-coordinate gauge.
The evolving variables are dimensionless coefficients; the physical Jordan
radius correction is $\alpha\ell H_J$, not $H_J$ itself.

## Independent radial mass and lapse reconstruction

At fixed original labels, the geometric Einstein-frame mass coefficient is

$$
\mu_E^{\rm geom}=\frac{mH}{r_0}
+\frac{2r_0}{C_0}
\left(r_{0,V}H_U+r_{0,U}H_V-c\,r_{0,U}r_{0,V}\right).
$$

It can be compared with the independent minimal-scalar mass transport:

$$
\mu_E^{\rm flux}(U,V)=\mu_{E,0}(V)
-\int_{U_0}^{U}\frac{r_0^2r_{0,V}}{C_0}p_U^2\,dq,
$$
$$
\mu_{E,0}(V)=-\frac14\int_{0.4}^{x(V)}t(a-t)p_x(t)^2dt.
$$

On the rectangle there is also the independently initialized second path,
$\mu_E(U,V)=-\int_{V_a}^{V}r_0^2r_{0,U}p_V^2/C_0\,dq$.
Its constant is zero on the quiet bottom.

There is a separate lapse diagnostic. In areal-null coordinates define

$$
Z=\frac{C_E}{2r_{E,U}},\qquad
ds_E^2/\ell^2=-Z^2\left(1-\frac{2\mathfrak m_E/\ell}{r_E}\right)dV^2
-2Z\,dV\,dr_E.
$$

The exterior background has $Z_0=-2a/V>0$, independent of radius.
For $\log Z=\log Z_0+\alpha\zeta_E+O(\alpha^2)$, the $UU$
constraint gives

$$
\zeta_E^{\rm geom}=c-\frac{H_U}{r_{0,U}},\qquad
\zeta_E^{\rm radial}(U,V)=\zeta_{E,0}(V)
+\int_{U_0}^{U}\frac{r_0p_U^2}{2r_{0,U}}\,dq.
$$

The source constant follows from the existing transverse radius data.
In the source's $(U,x)$ coordinates put
$C_x=-4a/U_0$, $d=(a/x-1)/4$, $p_U=C_x s_\phi$ and
$H_{J,U}=C_x k$. Then

$$
\zeta_{E,0}=-\frac{k}{d}-\frac{p^2}{12}
+\frac{xps_\phi}{6d}
=\boxed{\frac{h_J}{x}-\frac{p^2}{4}
-\frac{I_0}{6(a-x)}}.
$$

The constants vanish at the same quiet corner. The denominator $a-x$
stays positive on the original source leg. These diagnostics share the
known field but integrate different gravitational equations from the mixed
solver; they do not fit new source data.

Recover the original metric by

$$
H_J=H+\frac{r_0p^2}{12},\qquad c_J=c+\frac{p^2}{6}.
$$

Because both the areal radius and its derivative transform, the areal-null
lapse has a different correction:

$$
\zeta_J=\zeta_E+\frac{p^2}{12}
-\frac{r_0p p_U}{6r_{0,U}}.
$$

This follows by expanding $Z_J=C_J/(2r_{J,U})$. Moving between the two
areal-radius coordinates adds no background evaluation term here because
$\partial_{r_0}Z_0=0$. At the source it gives $\zeta_{J,0}=-k/d$, as
required by the fixed Jordan lapse and corrected transverse radius.
The full Jordan mass conversion, for a separate frame check, is

$$
\mu_J=\mu_E+\frac{mp^2}{12}
+\frac{r_0^2p}{3C_0}(r_{0,V}p_U+r_{0,U}p_V).
$$

## Constraint checks and measured calculation

The two linear null constraints are

$$
H_{UU}-(\log C_0)_U H_U-r_{0,U}c_U+\frac{r_0}{2}p_U^2=0,
$$
$$
H_{VV}-(\log C_0)_V H_V-r_{0,V}c_V+\frac{r_0}{2}p_V^2=0.
$$

Their derivatives can be evaluated from the evolution equations. For
example, if $A_0=C_0/(4r_0^2)$,
$(y_{UU})_V=A_{0,U}y+A_0y_U-C_{0,U}z/4-C_0z_U/4$.
The analogous $y_{VV}$ follows by differentiating the $y_V$ integral,
using the source's analytic second derivative. No scalar Hessian in the
bulk or derivative of sampled gravitational curves is required for this
formulation.

The implementation currently checks the equivalent integrated constraints.
Define $k_E=(H_U-r_{0,U}c)/C_0$ and
$l_E=(H_V-r_{0,V}c)/C_0$. Then

$$
k_E(U,V)=k_{E,0}(V)-\int_{U_0}^{U}\frac{r_0p_U^2}{2C_0}\,dq,
\qquad
l_E(U,V)=-\int_{V_a}^{V}\frac{r_0p_V^2}{2C_0}\,dq,
$$

where $k_{E,0}=-(r_{0,U}/C_0)\zeta_{E,0}$. The $V$ integral is
evolved separately from the metric variables. The areal-null lapse check
is the same $U$ constraint expressed in different variables; it is not an
additional independent field equation. Spatial prefix quadrature and
time evolution supply the mixed solution, with analytic first derivative
states and optional distinct endpoint quadrature for the mass comparison.

The [comparison producer](../scripts/record_nsc_atlas_bulk_response.py)
records separate changes of metric grid, leading-field resolution,
field-derivative quadrature, source quadrature, background geometry and
time-integration tolerance. Each comparison starts with a fresh cumulative
120 process-CPU-second, 30-million-scalar-kernel budget, includes two
metric solves, and saves nine transverse-leg cuts and maximum defects;
dense trial grids remain transient. It also compares the metric coefficients
on shared interior nodes. Zero source reproduces the exact reference. Source/bottom
conditions, the two null constraints, the mixed equations, both independent
mass expressions and both lapse expressions check different parts of the
calculation. A refinement difference is a numerical diagnostic rather than
a finite-amplitude remainder bound. The [current results](current-result.md)
record the measured values and refinement status. To create a fresh
comparison or authenticate it without running the evolution:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --comparison metric --output .build/bulk-metric-review
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --check --output .build/bulk-metric-review
```

Other comparison names are `field`, `geometry`, `derivative`, `source` and
`time`. The earlier v8 prepared-response record retains its original source
hashes. Its historical replay remains `scripts/replay_history.py atlas-v8`;
the added field interface does not rewrite that evidence.

## Continuation to the wall: an unevaluated analytic route

The radial mass and lapse integrals above follow constant-$V$ segments
from the source to the wall entirely inside the retained wedge. Extending
the bounded field evaluator to those segments would therefore supply the
smooth exterior metric in areal-null coordinates without assuming vacuum
data on an artificial extension of the bottom. This route is proposed for
subsequent calculation; the rectangle result alone does not establish it.

For a wall with known formal displacement $\rho(s)$ and signed
$\delta\beta_+(s)$ from the two-sided mass functionals, write
$V_W(s)=V_0(s)+\alpha V_1(s)$. The exact Jordan normalization and
signed normal give

$$
V_{W,s}=\frac{1}{Z_J(r_{W,s}-\beta_+)},
$$
$$
V_{1,s}=-V_{0,s}\left[
\zeta_J+\frac{Z_{0,V}}{Z_0}V_1
+\frac{\rho_s-\delta\beta_+}{v-\beta_+}\right],
\qquad V_1(s_a)=0.
$$

The denominator $v-\beta_+$ is positive in this white-hole strip,
including the zero of $\beta_+$. The residual $U$ reparametrization can
keep the wall coordinate curve equal to its reference $U_W(V)$ while
remaining the identity on the old bottom. At fixed $V$, define the wall
proper-time shift $S=-V_1/V_{0,s}$. Its Einstein radius boundary value is

$$
H_W=\rho+vS-\frac{xp^2}{12}.
$$

One may then recover the original double-null radius throughout the wedge
from the first-order null transport equation

$$
H_V=-\frac{Z_0}{2}\left[
\frac{2m}{r_0^2}H+\left(1-\frac{2m}{r_0}\right)\zeta_E
-\frac{2\mu_E}{r_0}\right],
\qquad c=\zeta_E+H_U/r_{0,U}.
$$

For $U\le U_a$ its initial value is the quiet bottom. For $U>U_a$ it is
the specified wall value at the inverse wall label $V=L(U)$. This names
the missing timelike input rather than setting it to zero. Differentiating
this transport equation supplies $H_U$ with the corresponding moving
lower-endpoint term.

The cap still needs its own smooth metric reconstruction with its quiet
initial geometry, regular centre and fixed central clock normalization.
Its leading regular field is already supplied by the transmitted
generating function and its reflected counterpart. Its mass is
vacuum-subtracted at the actual radius, and its Einstein-frame scalar
potential includes $h^2p^2$ at this order; the exterior zero-potential
equations cannot be copied into it. Subsequent observer arrival and field
corrections require the reconstructed cap and wall embeddings. Neither
this rectangle calculation nor these analytic continuation equations
establish finite-amplitude validity, formation, stability or a cosmological
prediction.
