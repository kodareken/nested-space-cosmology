# The smooth cap and its central clock

The [moving-wall calculation](nsc-atlas-wall-evolution.md) supplies both
one-sided masses, the wall displacement and its exterior null embedding.
The remaining interior connection is the metric between that wall and the
regular centre. It determines the central proper time at which an inward
null ray arrives. This note uses the same action, prepared scalar and
finite no-return interval; its [admission](prior-art-and-open-claim.md#smooth-cap-and-central-clock-admission-15-october-2026)
fixes the calculation and budget.

All response quantities are dimensionless coefficients of
$\alpha=\kappa_4\varepsilon^2/\ell^2$. The reference has $h=H\ell=1$,
$b=k\ell=2$, $m=M/\ell=0.4$ and wall radius $x\in[0.45,0.55]$.
No finite physical amplitude is selected here.

## The same scalar inside the cap

Use static cap time $T$ and reference areal radius $r$, with

$$
A_0(r)=1-h^2r^2,\qquad r_*(r)=\frac{\operatorname{atanh}(hr)}h.
$$

At the centre, $T$ is the proper time of the reference observer. Let
$q(V)$ be the existing arrival map and $w(V)=\chi_W(V)/\varepsilon$.
Define $f(q(V))=w(V)$, and use the established quiet preparation for
earlier labels. The regular leading cap field is

$$
p(T,r)=\frac{f(T+r_*)-f(T-r_*)}{r},
$$
$$
p_T=\frac{f'(T+r_*)-f'(T-r_*)}{r},\qquad
p_r=\frac{f'(T+r_*)+f'(T-r_*)}{rA_0}-\frac p r.
$$

The reflected term remains part of the interior solution even when the
no-return certificate sets it to zero at the wall. In particular,
$p(T,0)=2f'(T)$ and $p_r(T,0)=0$. The source ending does not set the
transmitted tail to zero. Labels later than the supplied $q_b$ require
new data and are outside this calculation.

The cap vacuum term transforms to
$V_E=3H^2/\kappa_4+H^2\phi^2+O(\phi^4)$.
Consequently the effective scalar mass squared is $2H^2$. The reduced
field $rp$ obeys the free equation in $(T,r_*)$ because
$A_0(A_{0,r}/r+2h^2)=0$. This uses the transformed fixed vacuum term;
the exterior zero-potential shortcut would give a different equation.

The numerical field is queried from the existing dense Volterra solution.
Its first derivative is $f'=w_V/q_V$, with a bounded monotone inverse for
$q$. Where a regular-centre evaluation needs a second derivative, it comes
from the differentiated continuum equation, not saved sample differences:

$$
w_{VV}=g''(V)-U'_W K_Ww
-\int_{U_0}^{U_W}\left(K_V\chi+K\chi_V\right)dU,
\qquad
f''=\frac{w_{VV}q_V-w_Vq_{VV}}{q_V^3}.
$$

Here $K_V=-(1/a+4/r)K r_V$ and $a=2m$. The compact polynomial source
supplies $g''$ analytically. Geometry interpolation and field discretization
remain separate numerical errors in this evaluation.

For small $r$, a stable form avoids subtracting almost equal endpoint
values. Set $d=r_*(r)$ and integrate over $z\in[-1,1]$:

$$
I_1=\int_{-1}^1f'(T+dz)\,dz,\quad
I_T=\int_{-1}^1f''(T+dz)\,dz,\quad
I_2=\int_{-1}^1z f''(T+dz)\,dz.
$$
$$
p=\frac d r I_1,\qquad p_T=\frac d r I_T,\qquad
p_r=\frac{r/A_0-d}{r^2}I_1+\frac d{rA_0}I_2.
$$

The first radial prefactor has the regular expansion
$2h^2r/3+4h^4r^3/5+6h^6r^5/7+\cdots$.
This calculation uses derivatives of the continuum field equation, with
its numerical representation errors retained in the refinement account.

## Radial metric and an independent mass path

Write the Einstein-frame polar-areal metric as

$$
\frac{ds_E^2}{\ell^2}=-e^{2\alpha\delta_E}A\,dT^2
+A^{-1}dr^2+r^2d\Omega^2,\qquad
A=A_0-\frac{2\alpha\nu_E}{r}+O(\alpha^2).
$$

The quadratic constraints and mass flux are

$$
\partial_r\nu_E=\frac{r^2}{4}
\left(\frac{p_T^2}{A_0}+A_0p_r^2\right)
+\frac{h^2r^2p^2}{2},
\qquad
\partial_T\nu_E=\frac{r^2A_0}{2}p_Tp_r,
$$
$$
\partial_r\delta_E=\frac r2
\left(p_r^2+\frac{p_T^2}{A_0^2}\right).
$$

Regularity fixes $\nu_E(T,0)=0$. Choose $T$ to remain **Jordan central
proper time** at this order. Since $g_J=(1+\alpha p^2/6)g_E+O(\alpha^2)$,
this fixes the time-dependent integration constant
$\delta_E(T,0)=-p(T,0)^2/12$. It is a clock normalization, not an
additional physical source.

The scalar equation makes
$\partial_T(\partial_r\nu_E)=\partial_r(\partial_T\nu_E)$.
This checks conservation independently of the radial quadrature. The
Jordan vacuum-subtracted mass at a reference point is

$$
\nu_J=\nu_E-\frac{r^2A_0p p_r}{6}
-\frac{h^2r^3p^2}{12}.
$$

At the wall, use the already defined tangential and normal derivatives
$A_W,B_-$:

$$
A_W=\frac{\beta_-}{A_0}p_T+vp_r,\qquad
B_-=\frac v{A_0}p_T+\beta_-p_r,
\qquad vA_W-\beta_-B_-=-A_0p_r.
$$

The last two terms in $\nu_J$ are exactly the prior $Q_-$.
The independently integrated radial $\nu_E$ must therefore agree with
the prior timelike $\int L_-ds$. This comparison includes the fixed cap
potential and checks the actual interior field against the boundary law.

## Wall time and central arrival

Keep $s$ as Jordan proper wall time. Write
$T_W=T_0+\alpha T_1$, $T_{0,s}=\beta_-/A_0(x)$ and
$r_{E,W}=x+\alpha e$, where

$$
e=\rho-\frac{xp_W^2}{12},\qquad
A_1=-2h^2xe-\frac{2\nu_{E,W}}x.
$$

The Jordan-to-Einstein change of the wall normal and proper clock gives

$$
T_{1,s}=\frac{\delta\beta_- -xp_WB_-/6-\beta_-p_W^2/12}{A_0(x)}
-\frac{\beta_-\delta_{E,W}}{A_0(x)}
-\frac{\beta_-A_1}{A_0(x)^2},\qquad T_1(s_a)=0.
$$

The quiet event anchors the relative time origins. The cap mass in $A_1$
is the Einstein excess $\nu_E$, not the Jordan mass or the total vacuum
mass. This distinction is measurable in the radial/timelike comparison.

There is also a direct induced-metric check. With
$e_s=\rho_s-(vp_W^2+2xp_WA_W)/12$, expand the wall's Jordan proper
normalization without using its cap time equation:

$$
N=2\beta_-T_{1,s}+\frac{2\delta_{E,W}\beta_-^2}{A_0}
+\frac{A_1(\beta_-^2+v^2)}{A_0^2}
-\frac{2ve_s}{A_0}+\frac{p_W^2}{6}.
$$

The prior angular constraint and frame term reduce this residual to
$N=2(\nu_{E,W}^{\rm radial}-\int L_-ds)/(xA_0)$.
It exposes the same mass-path discrepancy in the induced metric;
it is not an independent total error estimate.

At fixed proper emission time $s$, the reference central arrival is
$q_0=T_0+r_*(x)$. Expanding the inward null equation gives

$$
T_C=q_0+\alpha\Delta_C+O(\alpha^2),
$$
$$
\Delta_C=T_1+\frac e{A_0(x)}
+\int_0^x\left[
\frac{2\nu_E(q_0-r_*(r),r)}{rA_0(r)^2}
-\frac{\delta_E(q_0-r_*(r),r)}{A_0(r)}
\right]dr.
$$

The integrand is regular at the centre. For emission at a **fixed exterior
null label** $V$, the same observable needs the already computed wall
embedding correction $V_1$:

$$
\left.\Delta_C\right|_V
=\Delta_C-q_{0,s}\frac{V_1}{V_{0,s}}.
$$

The two coefficients answer different clock questions. Their units are
$\ell$ after multiplying the dimensionless time by $\ell$ (and by $1/c$
if physical seconds are used). A correction to a null arrival map does
not yet compute the order-$\varepsilon\alpha$ scalar amplitude; that
requires its perturbed wave equation, moving boundary and stated source
continuation.

An independent combination of the radial constraints explains the optical
integrand. For $h=1$, define $D=\delta_E-2\nu_E/(rA_0)$. Then

$$
D_r=\frac{2\nu_E(1-3r^2)}{r^2A_0^2}-\frac{rp^2}{A_0},
\qquad D(T,0)=-\frac{p(T,0)^2}{12}.
$$

The kinetic terms cancel in this identity, and the ray integrand is
$-D/A_0$. Positive mass alone therefore does not select its sign. Signed
quadrature contributions and sampled extrema locate the sign in the
computed cap; they do not establish a separate continuum error bound.

## Finite domain and verification

Put $T_a=q_a-r_*(x_a)$, $T_b=q_b-r_*(x_b)$ and
$q_{\min}=T_a-r_*(x_a)$. Before $T_b$, the retained cap is
$T_a\le T\le T_b$, $0\le r\le x_W(T)$. After the last wall event it
closes on the inward ray: $T_b\le T\le q_b$ and
$0\le r\le\tanh[h(q_b-T)]/h$. All required profile labels lie between
$q_{\min}$ and $q_b$. Constant-$T$ radial metric integrals on each
arrival ray remain inside that domain and require no later wall history.

The existing preparation certificate supplies the quiet outgoing wall
labels. Verification compares the radial and timelike cap masses, the
wall's scalar derivatives, the central Jordan normalization, zero-source
recovery and analytic manufactured fields. Field, geometry, centre
evaluation, radial, wall-time and ray quadratures are refined separately.
Their observed differences are numerical sensitivity estimates; they are
not a proven finite-amplitude remainder.

## Implementation and bounded comparisons

The [cap field evaluator](../src/recursive_horizons/nsc_atlas_cap_field.py)
reuses the original dense field. Its bounded arrival inverse uses at most
56 floating bisections; the reported bracket width is conditional on the
reference spline and is not an interval enclosure of geometry error. The
small-radius integral formulas apply at $r\le10^{-4}$ by default, with
$10^{-3}$ retained as a separate numerical control. Their Gauss order is
8, compared with 16. The differentiated exterior integral uses 2 points per
active field cell, compared with 4. Every kernel query shares the original
field's cumulative compute budget.

The [cap metric and clock integrator](../src/recursive_horizons/nsc_atlas_cap_response.py)
splits constant-time radial quadratures at both incoming and reflected
profile edges. Ray quadratures also split where the central-clock
normalization crosses a profile edge. The cap wall clock uses 4-point
Gauss panels, with partial-panel values from their integrated cubic.
Tests of [field regularity](../tests/test_nsc_atlas_cap_field.py) and
[metric/clock identities](../tests/test_nsc_atlas_cap_response.py) include
independent polynomial profiles, closed radial and ray integrals, exact
conservation and optical identities, quiet recovery and clock conversion.

The [comparison producer](../scripts/record_nsc_atlas_cap_response.py)
uses a primary setting of 257 field nodes, 257 geometry knots, radial and
ray orders 16, and 16 cap clock panels. It reuses the wall integrator with
radial order 2, field-derivative order 32, 8 exterior clock panels, source
order 64 and ODE tolerance $2\times10^{-9}$; the separate wall evidence
records the sensitivity of those inputs.

Nine comparisons vary one numerical setting within each pair: field
129/257, radial 8/16, ray 8/16, cap time panels 8/16, geometry 129/257,
cap derivative 2/4, centre Gauss order 8/16, centre threshold
$10^{-4}/10^{-3}$, and a further ray 16/32 comparison. The geometry,
derivative, centre and further-ray comparisons fix 129 field nodes to fit
the complete pair within the same 120 process-CPU-second and 30-million
scalar-kernel budget. These fixed overrides are recorded explicitly;
each pair compares the same physical preparation and otherwise identical
settings. The main field comparison checks 129 against 257 nodes at the
primary cap quadrature settings.

Each record retains nine common wall cuts, mass/clock defects, signed ray
contributions, sampled integrand extrema and cumulative execution counts.
Dense pilot arrays stay under the ignored build directory. Authentication
checks the declared settings and bound source bytes without another
physical solve. The [current result](current-result.md) owns measured
outcomes and links the completed evidence. The higher-order scalar and
finite-amplitude remainder remain the next dynamical connection.
