# The finite-amplitude conformal cap and material wall

The [owned action](nsc-atlas-cap-wall-response.md) fixes an exact nonlinear
model, while the existing metric, wall and scalar results are perturbative
coefficients of that model. This note identifies its finite-amplitude source
data, smooth equations and interface laws. It distinguishes an elementary
source-leg continuation bound from the still unestablished existence and
remainder estimate for the full transmitted spacetime.

The physical inputs are unchanged: a canonical Jordan scalar with coupling
`xi=1/6`, fixed vacuum densities `V0,-=3H^2/kappa_4` and `V0,+=0`, constant
Jordan material tension `sigma`, and the independently prescribed reduced
field on the anchored exterior source. The reference has `h=H ell=1`,
`b=k ell=2`, `k=kappa_4 sigma/2`, and `M/ell=0.4`. The source profile below
is called `g_src`; `b` always denotes the wall parameter.

## Exact frame and potential

Start from the same split Jordan action, including its `2fK` boundary terms,

$$
S_J=\sum_\pm\int_{\mathcal M_\pm}\sqrt{-g_J}
 \left[\left(\frac1{2\kappa_4}-\frac{\phi_J^2}{12}\right)R_J
 -\frac12(\nabla\phi_J)^2-V_{0,\pm}\right]d^4x
 -\sigma\int_\Sigma\sqrt{-\gamma_J}\,d^3y.
$$

On the branch $\mathcal F=1-\kappa_4\phi_J^2/6>0$, define

$$
z=\sqrt{\kappa_4}\phi_J,\qquad
\Psi=\sqrt{\kappa_4}\psi=\sqrt6\operatorname{atanh}(z/\sqrt6),
\qquad g_E=\mathcal F g_J.
$$

Thus $z=\sqrt6\tanh(\Psi/\sqrt6)$,
$\mathcal F=\operatorname{sech}^2(\Psi/\sqrt6)$, and
$d\Psi/dz=\mathcal F^{-1}$. The smooth Einstein-frame action has
canonical scalar kinetic term. In dimensionless coordinates and curvature
units `ell`, its potential is

$$
\widehat V_\pm(\Psi)=\kappa_4\ell^2 V_{E,\pm}
=3h_\pm^2\mathcal F^{-2}
=3h_\pm^2\cosh^4(\Psi/\sqrt6),\qquad h_+=0,\quad h_-=1.
$$

The fixed Jordan vacuum density has transformed; it has not been multiplied
by the Jordan gravitational coupling. In particular,
$\widehat V_-=3+\Psi^2+O(\Psi^4)$ reproduces the cap potential used in
the quadratic metric calculation. The canonical scalar `Psi` is distinct
from the conformally transformed scalar used for invariant reduced-field
propagation in the [scalar response](nsc-atlas-scalar-metric-response.md).
For example, the exact invariant reduced field is
$\chi/\varepsilon=(\sqrt6/\sqrt\alpha)r_E\sinh(\Psi/\sqrt6)$,
with its continuous normalized limit at `alpha=0`; it is not
$r_E\Psi/\sqrt\alpha$ at finite amplitude.

## Exact smooth double-null equations

Use dimensionless Einstein areal radius `r` and metric

$$
ds_E^2/\ell^2=-C\,dU\,dV+r^2d\Omega^2,\qquad C>0.
$$

The two null constraints, radius equation, mixed lapse equation and scalar
equation are

$$
r_{aa}-(\log C)_a r_a=-\frac r2\Psi_a^2,
\qquad a=U,V,
$$
$$
r_{UV}=-\frac{r_Ur_V}{r}-\frac C{4r}
       +\frac{Cr}{4}\widehat V_\pm(\Psi),
$$
$$
(\log C)_{UV}=\frac{2r_Ur_V}{r^2}+
\frac C{2r^2}-\Psi_U\Psi_V,
$$
$$
\Psi_{UV}=-\frac{r_U\Psi_V+r_V\Psi_U}{r}
           -\frac C4\widehat V_{\pm,\Psi}.
$$

The potential cancels from the displayed mixed lapse equation, but remains
in the radius and scalar equations. Expanding
$\Psi=\sqrt\alpha\,p+O(\alpha^{3/2})$ and
$\alpha=\kappa_4\varepsilon^2/\ell^2$ recovers the existing retained-order
equations. The exact equations require the evolved canonical field, not
the fixed leading `p` in a finite-amplitude source term.

These coordinates require a regular formulation at `r=0`. The centre has
zero reduced field, even regular scalar and metric profiles, and enclosed
mass of order `r^3`; the apparent inverse-radius terms must be combined or
treated with weighted regular variables. In polar-areal cap coordinates,
keeping `T` equal to Jordan central proper time fixes exactly

$$
e^{2\delta_E(T,0)}=\mathcal F(T,0),\qquad
\delta_E(T,0)=\tfrac12\log\mathcal F(T,0).
$$

This reduces to the previously used $-\alpha p_0^2/12$ anchor. It is not
a free lapse function or additional source.

## Exact material interface

Keep the same signed normal from cap to exterior. Write
$\widehat n=\ell n$ and $\widehat K^a{}_b=\ell K^a{}_b$.
The conformal Jordan junction laws on the positive-coupling branch are

$$
[\gamma_J]=0,\quad [\phi_J]=0,\quad
[\widehat K_J{}^a{}_b]=-b\delta^a{}_b,\quad
[\widehat n_J\phi_J]=b\phi_J.
$$

Because the conformal factor is continuous, so are the Einstein induced
metric and `Psi`. The same wall action becomes

$$
S_\Sigma=-\int_\Sigma\sqrt{-\gamma_E}\,
\sigma_E(\Psi)\,d^3y,\qquad
\sigma_E=\sigma\mathcal F^{-3/2}
=\sigma\cosh^3(\Psi/\sqrt6).
$$

This scalar dependence is forced by the frame change. It is not a new
surface interaction. To see the extrinsic jump explicitly, use

$$
\widehat K_E{}^a{}_b=\mathcal F^{-1/2}
 \left(\widehat K_J{}^a{}_b+
 \delta^a{}_b\widehat n_J\log\sqrt{\mathcal F}\right),
\qquad
[\widehat n_J\log\sqrt{\mathcal F}]
=-\frac{bz^2}{6\mathcal F}.
$$

The exact Einstein-frame conditions are therefore

$$
[\widehat K_E{}^a{}_b]
=-b\mathcal F^{-3/2}\delta^a{}_b,
\qquad
[\widehat n_E\Psi]=b z\mathcal F^{-3/2}
=b\sqrt6\sinh(\Psi/\sqrt6)\cosh^2(\Psi/\sqrt6).
$$

The second equation also equals
$\kappa_4\ell\,d\sigma_E/d\Psi$ from the transformed action.
Keeping constant Einstein tension or imposing a continuous normal
canonical-scalar derivative would change the original model. Induced-metric
matching, proper wall normalization and the signed angular/temporal
conditions must be used together. The reference zero of `beta_plus` is not
a singularity; a reduction dividing by that quantity is unsuitable.

## Exact invariant-field reduction and regular readout

There is an exact reduced equation for the same invariant field propagated
in the perturbative calculation. Set

$$
c=\chi/\varepsilon,\qquad u=c/r,\qquad
w=\sqrt6\sinh(\Psi/\sqrt6)=\sqrt\alpha\,u.
$$

Here `u` is a scalar ratio, not a null coordinate. The continuous normalized
equations at `alpha=0` follow from

$$
\Psi=\sqrt6\operatorname{asinh}(\sqrt\alpha\,u/\sqrt6),\qquad
\mathcal F=(1+\alpha u^2/6)^{-1},\qquad
\Psi_a=\frac{\sqrt\alpha\,u_a}{\sqrt{1+\alpha u^2/6}},\qquad
\widehat V=3h^2(1+\alpha u^2/6)^2.
$$

Indeed, $w_{,\Psi\Psi}=w/6$ and this potential satisfies
$w_{,\Psi}\widehat V_{,\Psi}=2\widehat Vw/3$.
The canonical equation and Einstein trace
$R_E=(\nabla\Psi)^2+4\widehat V$ therefore give
$(\Box_E-R_E/6)w=0$. On either smooth side of the wall,

$$
\boxed{c_{UV}+\mathcal Kc=0,\qquad
\mathcal K=-\frac{r_{UV}}r+\frac{CR_E}{24}
=\frac{Cm}{2r^3}-\frac{C\widehat V}{12}
-\frac{\Psi_U\Psi_V}{6},}
\qquad
m=\frac r2\left(1+\frac{4r_Ur_V}{C}\right).
$$

The mass `m` includes the cap vacuum contribution. Substitution of the
relations above makes all scalar backreaction terms explicitly proportional
to `alpha`. The vacuum exterior gives its original Schwarzschild kernel;
the vacuum cap has `m=h^2 r^3/2` and `K=0`. Differentiating this mass using
the null constraints and mixed radius equation also gives the exact laws

$$
m_U=\frac{r^2\widehat V}{2}r_U-
\frac{r^2r_V}{C}\Psi_U^2,\qquad
m_V=\frac{r^2\widehat V}{2}r_V-
\frac{r^2r_U}{C}\Psi_V^2.
$$

These are first-derivative identities, not replacements for proving that
the constraints propagate in a chosen independent evolution system.

The reduced field has homogeneous transmission conditions exactly.
The angular Jordan jump gives
$[\widehat n_J r_J]=-br_J$, which cancels
$[\widehat n_J\phi_J]=b\phi_J$ in their product. Consequently

$$
[c]=[\widehat n_Jc]=[\widehat n_Ec]=0.
$$

In exterior null labels `(U,V)` and cap null labels `(Q_-,Q_+)`, continuity
of tangential and normal derivatives is equivalently

$$
\dot U_+c_{U,+}=\dot Q_-c_{Q_-,-},\qquad
\dot V_+c_{V,+}=\dot Q_+c_{Q_+,-},
$$

where a dot uses one common proper-time parameter on the wall. Its timelike
label rates, rather than `beta_plus`, are the denominators of this
transmission map. This statement includes the original wall source:
the singular conformal-curvature term and the jump of the conformal
scalar's normal derivative cancel. The displayed smooth formula for `K`
must not be applied as a distributional curvature formula at the wall.
For a continuous timelike multiplier tangent to the wall, the scalar energy
flux is proportional to $(n c)(Xc)$ and is continuous. This cancels the
basic interface flux; it does not estimate differentiated metric or moving
interface terms.

Choose actual cap null labels with `Q_-=Q_+=T` on the axis, where `T` is
Jordan central proper time. At a regular centre,

$$
c=0,\quad c_{Q_-}+c_{Q_+}=0,\qquad
r_{Q_-}=-r_{Q_+},\qquad
C=4r_{Q_+}^2=\mathcal F,\qquad r_{Q_+}=\sqrt{\mathcal F}/2.
$$

The normalized Jordan central scalar is therefore exactly

$$
\boxed{p_{J,C}=\lim_{r\to0}\sqrt{\mathcal F}\frac cr
=2c_{Q_+}\big|_{Q_-=Q_+=T}.}
$$

In this gauge a uniform first-jet remainder controls the physical central
readout without an additional frame or clock correction. At the axis
$\mathcal F=1-2\alpha c_{Q_+}^2/3$; positivity would follow from a proved
uniform derivative bound and a corresponding small amplitude. Sampled
derivatives do not establish that bound. Regularity gives

$$
\frac m{r^3}\longrightarrow\frac{\widehat V_C}{6}
+\frac{\Psi_T^2}{12\mathcal F_C},\qquad
\Psi_{Q_-}=\Psi_{Q_+}=\Psi_T/2.
$$

Thus the constant terms in `K` cancel at the centre. With sufficient even
radial regularity, `K=O(r^2)`. A numerical or analytic estimate must retain
this cancellation or use the regular mass-density variable; independent
bounds on singular-looking terms do not prove the regular limit.

### Temporal junction on the expanding Einstein branch

One junction redundancy can be resolved exactly. In dimensionless Einstein
proper wall time `s_E`, let

$$
\omega=\frac{dr_W}{ds_E},\qquad
\beta_\pm=\widehat n_{E,\pm}r,\qquad
B=b\mathcal F^{-3/2}.
$$

The angular condition is $[\beta]=-Br_W$.
The mixed smooth Einstein equation gives
$\operatorname{Hess}r(\dot\Sigma,n)=-r\Psi_{s_E}\Psi_n/2$,
and hence

$$
\frac{d\beta_\pm}{ds_E}
=\widehat K_{E,\pm}{}^{s_E}{}_{s_E}\,\omega
-\frac r2\Psi_{s_E}\widehat n_{E,\pm}\Psi.
$$

Since $B_{,\Psi}=Bz/2$ and
$[\widehat n_E\Psi]=Bz=2B_{,\Psi}$, differentiating the angular
condition cancels the derivative of `B` and yields

$$
\omega\left([\widehat K_E{}^{s_E}{}_{s_E}]+B\right)=0.
$$

On a tube preserving `omega>0`, the temporal Israel condition follows
from the angular condition, scalar jump and smooth Einstein equation.
This requires Einstein areal expansion; Jordan expansion alone is not
the finite-amplitude condition. No division by `beta_plus` occurs.
Once the full extrinsic jump holds, the Gauss constraint gives
$[\widehat{\mathcal T}_{nn}]=-2B\overline{\widehat K_E}$, with
dimensionless Einstein stress defined by
$\widehat G_E=\widehat{\mathcal T}$ and `K` here denoting the extrinsic
trace. This identifies the normal-force condition as a consistency
consequence when the bulk constraints hold. It does not prove their
propagation or the existence of that expanding tube.

## Exact source family in the original affine gauge

The [initialization owner](nsc-atlas-characteristic-initialization.md)
keeps the Jordan null normalization fixed on `U=U0`. Its original
$x=R_0(U_0,V)/\ell$ remains an affine parameter on that source leg. Preserve
the specified continuation

$$
\chi(U_0,V)=\varepsilon g(V),\qquad
g_{\rm src}(x)=g(V(x)),\qquad r_J(x)=R_J(U_0,V(x))/\ell.
$$

The normalized compact profile is the existing `C3` function, including
its quiet initial interval. Set

$$
z(x)=\sqrt\alpha\frac{g_{\rm src}(x)}{r_J(x)},\qquad
r_J(x_a)=x_a=0.4,\qquad r_J'(x_a)=1.
$$

In this affine coordinate the exact Jordan null constraint is

$$
-\frac{2\mathcal F r_J''}{r_J}
=\frac23 z'^2-\frac13 zz''.
$$

Multiplication by $-3r_J$ gives
$6\mathcal F r_J''=r_Jzz''-2r_Jz'^2$. Substituting the prescribed
reduced field produces a term $-z^2r_J''$ on this right side. Since
$6\mathcal F+z^2=6$, the second derivative is isolated without a residual
conformal denominator:

$$
\boxed{r_J''=\frac{\alpha}{6r_J}
 \left(g_{\rm src}g_{\rm src}''-2(g_{\rm src}')^2\right)
 +\frac{\alpha}{3r_J^2}g_{\rm src}g_{\rm src}'r_J'.}
$$

Writing $r_J=x+\alpha h+O(\alpha^2)$ gives

$$
h''=\frac{g_{\rm src}g_{\rm src}''-2(g_{\rm src}')^2}{6x}
+\frac{g_{\rm src}g_{\rm src}'}{3x^2}
=-\frac{x}{3}p'^2+\frac{x}{6}pp'',\qquad p=g_{\rm src}/x,
$$

which is exactly the owned quadratic source constraint. The finite-amplitude
Einstein initial traces follow by transforming the same source solution:
$r_E=\sqrt{\mathcal F}r_J$, $C_E=\mathcal F C_{J,0}$ and the canonical
field above; explicitly $r_E=\sqrt{r_J^2-\alpha g_{\rm src}^2/6}$.
After the compact source ends, the exact source ODE has $r_J''=0$,
while its radius and slope retain their evolved values.
The bottom and cap initial data remain quiet reference data;
the fixed corner expansion and lapse are unchanged. A separate Einstein
affine reparametrization must not be mistaken for the original source gauge.
There is no additional freely selectable source mass function. On the
exterior source, its transformed `r,C` traces and the original corner
value determine the transverse expansion through
$(r_U)_V=-r_Ur_V/r-C/(4r)$; the mass then follows from its geometric
definition. The specified source continuation, quiet data, corner and
wall branch identify the exact family whose construction is still needed.

### What a source-leg bound establishes

There is an elementary finite-interval continuation argument for this ODE.
Let `Lx` be the finite affine interval length, choose `0<r_min<x_a`,
`1<v_max<2` and `0<F_min<1`, and define exact profile bounds

$$
M_0=\|g_{\rm src}\|_\infty,\quad
M_1=\|g_{\rm src}g_{\rm src}'\|_\infty,\quad
M_2=\|g_{\rm src}g_{\rm src}''-2(g_{\rm src}')^2\|_\infty,
\qquad L=\frac{M_2}{6r_{\min}}+
\frac{M_1v_{\max}}{3r_{\min}^2}.
$$

While $r_J\ge r_{\min}$ and $|r_J'|\le v_{\max}$, the ODE gives
$|r_J''|\le\alpha L$. Integration from the fixed corner gives

$$
|r_J'-1|\le\alpha L L_x,\qquad
|r_J-x|\le\tfrac12\alpha L L_x^2.
$$

Consequently the strict conditions

$$
\alpha L L_x<v_{\max}-1,\qquad
\tfrac12\alpha L L_x^2<x_a-r_{\min},\qquad
\frac{\alpha M_0^2}{6r_{\min}^2}<1-\mathcal F_{\min}
$$

close those bounds and keep the conformal factor above `F_min` throughout
the source leg. They also give $r_J'>2-v_{\max}>0$.
The vector ODE is locally Lipschitz for positive radius;
these bounds keep its solution in a bounded nonsingular region, so the
usual finite-interval ODE continuation applies. Finite profile derivatives
make these conditions nonempty for sufficiently small positive `alpha`.
No numerical values for the constants or interval are asserted here.

The same ODE permits a source-only remainder estimate. On a specified
nonsingular tube, the exact vector field and its derivatives have finite
Lipschitz bounds. Substitution of the exact quadratic coefficient
`(x+alpha h, 1+alpha h')` gives a residual of order `alpha^2`; a Gronwall
estimate then bounds the source radius/slope error with an explicit tube
constant. Numerical approximations to `h` would add their own enclosed
error. This argument concerns the initial leg only: it proves neither a
two-sided evolved solution nor a central-waveform remainder.

## Missing full-domain estimate

A justified finite-amplitude spacetime interval needs a quantitatively
well-posed coupled problem for the smooth equations, the timelike interface
and the regular centre. The equations above identify its data; they are not
such a theorem. In particular, the following remain to be established.

- Choose an independent characteristic/interface system with anchored gauge
  and a fixed-domain identification. Show that the remaining null Einstein
  constraints propagate from compatible initial and interface data. The
  temporal and normal consequences derived above apply only when their
  smooth equations, constraints and positive Einstein expansion hold.
  Distributional stress conservation must use the transformed
  scalar-dependent tension if the Einstein frame is used.
- Prove a regular-centre estimate and a boundary energy or characteristic
  estimate for the moving transmission law, without derivative loss or
  division by `beta_plus`. Standard smooth-bulk scalar evolution does not
  alone establish this moving-interface estimate.
- Control a tube with positive conformal factor and lapse, positive exterior
  radius, timelike wall, monotone label maps and the no-return margin. The
  central comparison must remain on a compact interior arrival interval,
  or the terminal domain must be continued consistently.
- Obtain uniform bounds for a continuous approximate solution, its bulk
  and boundary residuals and the inverse linearized evolution. Stored cuts,
  interpolation checks and refinement differences are not those bounds.
  The `C3` compact preparation and its source edges must fit the chosen
  function space; unlimited source differentiability is unavailable.

The linear preparation's known Volterra norm bound `q_*<1` supplies an
inverse estimate for the scalar on the **fixed reference** wedge. It does
not bound the full metric/interface/centre inverse. Its reduced-field
supremum bound alone does not control the central limit of `chi/r`;
regular first-derivative bounds are also needed. The quadratic metric
and cubic scalar calculations evaluate particular forced directions, not
that inverse's operator norm. Older bounds for other actions or domains
cannot fill this gap.

A concrete next estimate is the forced reference transmission resolvent in
a norm that controls one scalar derivative at the centre, together with
the metric-constraint and wall maps in the same regularity class. On a
regular cap radial coordinate `x` with `r_J/x>=a_min>0`, the centre condition
`c(T,0)=0` gives
$\|c/r_J\|_\infty\le a_{\min}^{-1}\|c_x\|_\infty$.
A first-commuted scalar energy controlling `c` in spatial `H2` and its
time derivative in `H1` can therefore control this readout; uncommuted
energy alone cannot. In a wall-adapted Gaussian-normal collar, coordinate
tangential differentiation commutes with the normal derivative, so both
homogeneous scalar transmission conditions commute. The bulk commutator,
moving coordinate coefficients and matching to a regular centre still
need bounds. The compact preparation is `C3` and piecewise `W4,infinity`;
unlimited differentiability must not be used to close that estimate.

The local commuted operator can be written explicitly. In that collar take
the quotient metric to be
$h=-a(s,n)^2ds^2+dn^2$, with the wall at `n=0`, common wall proper time
`s`, and `a(s,0)=1`. Define the scalar quotient potential
$W=(\Box_h r)/r+R_E/6$, which equals `4K/C` in a double-null chart.
Then

$$
Lc=-a^{-2}c_{ss}+\frac{a_s}{a^3}c_s+c_{nn}
 +\frac{a_n}{a}c_n-Wc.
$$

For `Lc=0` and `v=c_s`, direct differentiation gives

$$
Lv=-\frac{2a_s}{a^3}v_s
-\left(\frac{a_{ss}}{a^3}-\frac{3a_s^2}{a^4}\right)v
-\left(\frac{a_{ns}}a-\frac{a_na_s}{a^2}\right)c_n+W_sc.
$$

A prescribed forcing adds its `s` derivative to this right side. There
is no `c_nn` commutator term. Since the coordinate fields commute,
`[c]=[c_n]=0` implies `[v]=[v_n]=0`; the energy fluxes of both equations
cancel at the wall. This is a local derivative-preserving scalar estimate
on a prescribed collar with `a` bounded above and away from zero and
bounded one-sided `a_s,a_n,a_ss,a_ns,W,W_s`. An energy with a positive
added zeroth-order scalar term handles a potential of either sign.
These hypotheses are to be bounded, not inferred from the transmission
identities. The collar need not reach the centre: a compatible regular
centre patch, patch transition estimates and bounds for coefficients
generated by the coupled Einstein and wall equations remain necessary.

This route need not infer constants from the saved central cuts. For
example, after constructing a bounded reference resolvent and a
constraint-compatible map on a specified regular tube, a preconditioned
fixed-point map with defect at most `alpha A` and Lipschitz constant at most
`alpha L` would preserve a radius-`R` ball and contract if
$\alpha A\le R/2$ and $\alpha L\le1/2$, together with the geometric
margin inequalities. Finite proved `A,L` would then supply a nonzero
interval by explicit inequalities. The first-jet kernel and homogeneous
transmission laws above identify the operators to estimate; they do not
establish these constants, a derivative-preserving coupled map, or its
constraint compatibility. A bounded scalar interface estimate is a useful
next component, rather than a complete nonlinear certificate.

One concrete sufficient route is a residual estimate on this finite domain.
Use regular weighted variables and the normalized field `chi/epsilon`,
whose equations have a smooth `alpha -> 0` limit after the canonical frame
change. The approximate metric and wall retain their order-`alpha` terms;
the normalized scalar retains `c0+alpha c1`. If the coupled solution map and
its inverse estimates are controlled in that function space, their exact
residual starts at order `alpha^2`. A bounded inverse and a locally bounded
nonlinear derivative can then support a contraction or continuation
estimate on a stated tube. These inverse and nonlinear constants have not
been proved or enclosed for this model. If a whole-domain contraction is
too large, finitely many compatible characteristic patches are a possible
construction, not an already justified certificate.

The intended form of a remainder statement would be

$$
\|g-g_0-\alpha g_1\|+\|r_W-x-\alpha\rho\|
\le C_g\alpha^2,
\qquad
\|\chi/\varepsilon-c_0-\alpha c_1\|_{\rm reg,1}
\le C_c\alpha^2,
$$

on a named nonzero interval and domain. All quantities require an explicitly
fixed pullback to the reference domain, with the source
labels, quiet wall event and central proper clock anchored. In this display
`g` and `r_W` are Jordan-frame quantities, so `rho` is the existing Jordan
wall displacement. The scalar's regular first-derivative norm
must control the central scalar readout. The corresponding dimensional
scalar remainder is of order `epsilon alpha^2 / ell`. Numerical errors in
the leading and correction fields enter separately; they cannot be hidden
inside an asserted perturbative remainder.

Small sampled values of `alpha g1`, `alpha rho` or `alpha p1,C`, and a
positive truncated conformal factor, are useful coefficient-size
diagnostics. They do not imply existence, constrain the omitted coefficient
or bound the error between the truncated and exact solution. No nonzero
finite-amplitude validity interval for the complete cap-wall-exterior
spacetime is claimed by this note. Independent review of the exact
nonlinear reduction precedes its implementation.
