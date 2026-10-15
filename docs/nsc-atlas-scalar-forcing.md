# Scalar boundary lift and cumulative exterior forcing

The [scalar transmission equations](nsc-atlas-scalar-metric-response.md)
specify the next reduced-field coefficient. This successor implements its
boundary and bulk inputs without changing the source-bound calculations
already recorded there. The [admission](prior-art-and-open-claim.md#scalar-lift-and-forcing-admission-16-october-2026)
fixes the model, finite domain and compute budget.

All quantities are dimensionless coefficients on the original reference
wedge, with $h=1$, $b=2$, $m=0.4$, $a=2m$ and the prescribed original
reduced-field input held exactly fixed. The expansion parameter remains
$\alpha=\kappa_4\varepsilon^2/\ell^2$. No finite value of $\alpha$ is
selected by this calculation.

## Boundary derivative from first scalar jets

Write $u=q_{-,W}$, $q_W=q_{+,W}$ and $A=1-r^2$. Before the first reflected
return, $c_0(u,q)=f(q)$ on the active outgoing ray, and

$$
J_-=-\int_{q_a}^{q_W}K_1(u,q)f(q)\,dq,\qquad
T=(u+q)/2,\quad r=\tanh((q-u)/2).
$$

At fixed $q$, $r_u=-A/2$. The radial mass constraint and the mass flux
$\nu_{E,T}=r^2A p_Tp_r/2$ give

$$
\nu_{E,u}=\frac{\nu_{E,T}-A\nu_{E,r}}2
=-\frac{A(1+r^2)f^2}{8r^2}.
$$

Differentiating the quiet-profile potential with $f(q)$ fixed therefore gives

$$
K_{1,u}=
\frac{A(3-r^2)\nu_E}{4r^4}
-\frac{A(3-r^2)ff'}{24r^4}
+\frac{A(1-8r^2+3r^4)f^2}{48r^5}.
$$

Thus, for reference proper wall time $s$,

$$
J_{-,s}=-q_{+,s}K_{1,W}w
-q_{-,s}\int_{q_a}^{q_W}K_{1,u}f\,dq,
\qquad q_{\pm,s}=\frac{\beta_-\pm v}{1-x^2}.
$$

The extra integral uses the existing radial metric and first scalar jets;
it does not require $f''$ or a finite difference of a stored waveform.
The moving upper limit is retained. At the wall $w=f(q_W)$.

Set $d=v-\beta_+$, $E=e-vV_1/V_s$, $e=\rho-xp_W^2/12$, and
$C=2mw/(x^3d^2)$. The established boundary formula and its derivative are

$$
c_{1,r}|_W=\frac{2q_{-,s}}dJ_--EC,
$$
$$
(c_{1,r}|_W)_s=
\left(\frac{2q_{-,s}}d\right)_sJ_-
+\frac{2q_{-,s}}dJ_{-,s}-E_sC-EC_s.
$$

Here $v_s$ and $\beta_{\pm,x}$ are analytic reference coefficients, and

$$
d_s=v_s-\beta_{+,x}v,\qquad
e_s=\rho_s-\frac{vp_W^2+2xp_W A_W}{12},
$$
$$
C_s=\frac{2m}{x^3d^2}
\left[j-w\left(\frac{3v}{x}+\frac{2d_s}{d}\right)\right],
\qquad j=w_s,\quad A_W=(p_W)_s,
$$
$$
E_s=e_s-\frac{v_sV_1+vV_{1,s}}{V_s}
+\frac{vV_1V_{ss}}{V_s^2},\qquad
\frac{V_{ss}}{V_s}=\frac{V_s}{V}-\frac{d_s}{d}.
$$

The alternative identity $q_{-,s}=1/(\beta_-+v)$ supplies its derivative
without subtraction of nearly equal quantities. The exterior pullback
$\mathfrak b(V)=r_Uc_{1,r}|_W$ has

$$
\frac{(r_U)_s}{r_U}=\frac{V_s}{V}-v\left(\frac1a+\frac1x\right),\qquad
\mathfrak b_V=
\frac{(r_U)_s c_{1,r}|_W+r_U(c_{1,r}|_W)_s}{V_s}.
$$

### A pre-arrival profile moment

There is also a simpler exact representation on this quiet domain. Define

$$
\mathcal C(q)=\frac12\int_{q_a}^q[f'(z)]^2\,dz.
$$

For $T<q_a$ and $q=T+\mathrm{atanh}(r)$, regularity and the radial
constraint give

$$
\nu_E(T,r)=\mathcal C(q)-\frac{(1-r^2)f(q)^2}{4r}.
$$

To check it, differentiate at fixed $T$, using $q_r=1/(1-r^2)$:

$$
\nu_{E,r}=\frac{f'^2}{2(1-r^2)}-\frac{ff'}{2r}
+\frac{(1+r^2)f^2}{4r^2}.
$$

Substituting $p=f/r$, $p_T=f'/r$ and
$p_r=f'/[r(1-r^2)]-f/r^2$ into the full kinetic-plus-potential mass
constraint gives the same expression. Its centre anchor is zero because
$T<q_a$. Consequently the potential is

$$
K_1=\frac{(1-r^2)\mathcal C(q)}{2r^3}
-\frac{(1-r^2)ff'}{12r^3}
-\frac{(1-r^2)f^2}{12r^4}.
$$

This identity supplies an independent manufactured mass for the derivative
tests and a possible later replacement for repeated radial integrations.
The current physical boundary evaluator retains the existing radial metric
queries. Once $T\ge q_a$, the centre-reflected term is active and this
one-profile formula is insufficient; the full reflected field remains the
cap model there.

## Hermite boundary and an exact polynomial lift

At wall knots $V_i$, store $U_i=U_W(V_i)$, $b_i=\mathfrak b(V_i)$ and
$b_{U,i}=\mathfrak b_V(V_i)/U'_W(V_i)$. Piecewise cubic Hermite interpolation
in $U$ defines $b_H(U)$. Its exact polynomial antiderivative is

$$
B_H(U)=\begin{cases}
0,&U\le U_a,\\
\displaystyle\int_{U_a}^U b_H(\eta)\,d\eta,&U>U_a.
\end{cases}
$$

The finite domain is $U_0\le U\le U_b$. This choice satisfies
$dB_H(U_W(V))/dV=U'_W(V)b_H(U_W(V))$ algebraically, with no inverse-wall
query or extra quadrature. Initial quiet data remain exactly zero.

Knots at the source edges are retained. A 9/17-knot comparison uses the
16 additional midpoints of the 17-knot mesh only for withheld checks.
Report errors in $\mathfrak b$, local-width-scaled slopes and the lift.
Continuum slopes and numerically interpolated values have different input
errors; the withheld comparison checks their practical compatibility.
These sampled discrepancies are numerical sensitivity estimates, not a
uniform enclosure between all queries. With a uniform boundary error
$d_b$, the conditional lift bound is
$\|\delta B\|_\infty\le(U_b-U_a)d_b$.

## Cumulative exterior constraints

At a fixed original null label $V$, the Einstein constraints are

$$
\mu_{E,U}=-\frac{r^2r_Vp_U^2}{C_0},\qquad
\zeta_{E,U}=\frac{rp_U^2}{2r_U}.
$$

Integrate these once across the sorted union of the original field-grid
cells, requested active coordinates and moving wall. Gauss quadrature
within each cell and prefix sums supply all requested masses and lapses.
One source-trace query fixes both integration constants. The field evaluator
supplies first jets in batches within its existing finite query limit.
Inactive grid coordinates are clipped to the wall.

The scalar evolution uses the known flux and potential input

$$
F=\frac{Z_0}{2}\left(A_0\zeta_E-\frac{2\mu_E}{r}\right)
\frac{c_{0,U}}{r_U},\qquad Z_0=-\frac{2a}{V},\quad A_0=1-\frac{2m}{r},
$$
$$
G=\left[\frac{C_0(m\zeta_E+\mu_E)}{2r^3}
-\frac{p_Up_V}{6}\right]c_0,\qquad K_0=\frac{C_0m}{2r^3}.
$$

Both $p_U$ and $p_V$ here belong to the original null chart. The wall
condition $c_{0,U}|_W=0$ makes $F_W=0$, but does not remove the curvature
term in $G$. Source compatibility checks compare the field value and its
source-tangent derivative with the independently integrated source data.
The source pullback and its continuum derivative are

$$
g_1=-H_{E,\mathrm{src}}c_{0,r},\qquad
g_1'=-\left(H_{E,\mathrm{src}}'
+\frac{H_{E,\mathrm{src}}P_{0,r}}2\right)c_{0,r}
+\frac{H_{E,\mathrm{src}}K_0c_0}{r_U}.
$$

The transverse derivative can persist after the prescribed source trace
ends, so $g_1$ is not truncated at the end of that pulse. The independent
wall mass transport checks the cumulative metric through
$\mu_{E,W}+Q_{+,W}=\mu_{J,W}$. Repeating an endpoint radial quadrature is
only an implementation comparison, since it uses the same constraint.

## Subsequent evolution and conditional error propagation

The lifted field $z=y-g_1-B_H$ obeys

$$
z_V=F(U,V)-F(U_0,V)
-\int_{U_0}^{U}\{K_0[z+g_1+B_H]+G\}\,dU.
$$

Its zero wall derivative permits the existing clipped continuation.
Forcing should be cached on a checked time mesh before evolution;
recomputing the full radial integration at every ODE stage is unnecessary.
The present input comparison does not evolve this equation.

A later time cache can tabulate the known $F,G$ on
$\eta=(U-U_0)/(U_W(V)-U_0)$, then evaluate that interpolant at each solver
coordinate. This removes clipping transitions from the cache coordinates.
It changes neither the evolved coordinates nor the unknown, so it introduces
no transport or Jacobian term in the displayed equation. The kernel can
still be evaluated directly on the solver grid. Source edges and numerical
field-cell transitions still need checks; mapped targets can nearly double
the cumulative integration's panel count. That cost and interpolation error
must be measured before selecting the cache mesh.

On the same fixed wedge, kernel, quiet anchors and wall map, the existing
analytic operator bound is $Q=0.13882099234357653<1$. If the input errors
have genuine uniform bounds $d_g,d_b,d_F,d_G$, then

$$
\|\delta y\|_\infty\le
\frac{d_g+(U_b-U_a)d_b+2\Delta Vd_F+\Delta U\Delta Vd_G}{1-Q}.
$$

This follows by writing the clipped integral as either order of integration
on its moving domain. The zero wall derivative gives no additional transport
term. An integral-equation residual adds its sup norm to the numerator.
Changing the kernel, geometry or quiet datum needs the corresponding
additional error term. Sampled refinements alone do not provide the uniform
bounds in this statement, nor does this reduced-field bound control detector
derivatives or the central amplitude.

## Implementation

The [boundary lift](../src/recursive_horizons/nsc_atlas_scalar_lift.py)
provides joint outgoing values and derivatives, followed by a Hermite
antiderivative. The [forcing adapter](../src/recursive_horizons/nsc_atlas_scalar_forcing.py)
supplies one cumulative exterior slice. Their controls live in the
[lift tests](../tests/test_nsc_atlas_scalar_lift.py) and
[forcing tests](../tests/test_nsc_atlas_scalar_forcing.py). The
[record producer](../scripts/record_nsc_atlas_scalar_inputs.py) binds each
comparison to its inputs and producer bytes; authentication does not repeat
the physical calculation. Numerical outcomes are maintained in the
[current result](current-result.md) and commands in the
[reproduction guide](reproducing.md).
