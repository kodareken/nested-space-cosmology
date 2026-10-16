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

These are first-derivative identities. The smooth constraint transport
and source compatibility are derived below; their interface and centre
closure must still belong to a consistent independent coupled formulation.

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
The used Einstein equation in the wall tangent-normal frame includes null
constraints. This shortcut cannot be used circularly to establish those
constraints from the angular jump alone.
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

## Exact constraint transport and boundary input

The smooth mixed equations determine the propagation of the null
constraint residuals without assuming that those residuals vanish. Define

$$
\mathcal E_U=r_{UU}-(\log C)_U r_U+\frac r2\Psi_U^2,
\qquad
\mathcal E_V=r_{VV}-(\log C)_V r_V+\frac r2\Psi_V^2.
$$

Differentiating and substituting the mixed radius, lapse and scalar
equations cancels the potential and its derivative, giving

$$
\boxed{(r\mathcal E_U)_V=0,\qquad
(r\mathcal E_V)_U=0.}
$$

This proves smooth propagation along each null segment with an appropriate
constraint-satisfying endpoint. It does not supply an endpoint on a newly
entering wall characteristic. At a regular centre, bounded residuals imply
`r E -> 0`; the constant transported along a segment ending there is then
zero. Boundedness of the residuals is part of the regular-centre
construction, not a consequence of the displayed inverse-radius formula.

Source compatibility also holds exactly. In the original affine source
coordinate `x`, write primes for `x` derivatives. Its fixed Jordan
normalization has constant `C_J` along the leg, while

$$
r_E=\sqrt{\mathcal F}\,r_J,\qquad C_E=\mathcal F C_J,
\qquad \Psi'=z'/\mathcal F.
$$

The transformed constraint is

$$
\mathcal E_x^E
=\sqrt{\mathcal F}\,r_J''+
\frac{r_J(2z'^2-zz'')}{6\sqrt{\mathcal F}}=0,
$$

by the exact Jordan source ODE. Hence the exterior source supplies
`E_V=0`. The quiet reference bottom supplies `E_U=0` only for
$U_0\le U\le U_a$. Exterior characteristics with `U>Ua` first enter
through the wall and need a boundary constraint argument.

That argument must distinguish Codazzi from the normal force. Let
$\mathcal R_{ab}=G_{E,ab}-\mathcal T_{ab}$ denote the dimensionless
Einstein residual, and put

$$
S=b\mathcal F^{-3/2},\qquad N=\widehat n_E\Psi.
$$

If the full Israel and scalar laws are imposed independently,
$[\widehat K^a{}_b]=-S\delta^a{}_b$ and
$[N]=2S_{,\Psi}$, Codazzi gives
$[\mathcal R_{ns}]=0$. Gauss instead gives

$$
[\mathcal R_{nn}]
=-2S\overline{\widehat K}-2S_{,\Psi}\overline N
+[\widehat V].
$$

This expression does not vanish merely from the tangential constraint.
Its vanishing is the normal-force condition to impose independently in
this constraint-propagation argument, or to obtain from another
independently justified incoming-constraint prescription. It is a law of
the fixed model, not a new freely selectable incident datum.

Indeed, with the mixed evolution equations satisfied, on either side

$$
\mathcal R_{ns}=\frac2r
\left(\dot U^2\mathcal E_U-\dot V^2\mathcal E_V\right),\qquad
\mathcal R_{nn}=-\frac2r
\left(\dot U^2\mathcal E_U+\dot V^2\mathcal E_V\right),
$$

using that side's null labels and a common Einstein proper wall time.
Suppose the exterior incoming source has already supplied
$\mathcal E_{V,+}=0$ and the quiet or regular cap side supplies
$\mathcal E_{Q_-,-}=0$. Set

$$
\mathcal A_{\rm in}=\dot U_+^2\mathcal E_{U,+},\qquad
\mathcal B_{\rm in}=\dot Q_+^2\mathcal E_{Q_+,-}.
$$

Codazzi only gives
$\mathcal A_{\rm in}=-\mathcal B_{\rm in}$: one common incoming
constraint mode remains. An independently imposed zero normal-force
residual also gives
$\mathcal A_{\rm in}=\mathcal B_{\rm in}$, so both modes vanish.
The wall is timelike, so its null label rates are nonzero. This identifies
the boundary condition needed by the transport identities. It does not
prove that a proposed reduced wall evolution imposes the conditions
independently or has a bounded solution map. In particular the temporal
shortcut above already uses Einstein constraints and cannot supply this
missing condition circularly.

### One noncircular constraint-boundary selection

The incoming mode can be closed without assuming angular Israel during
the propagation argument. Assume a shared induced metric, common
Einstein proper tangent and radius `r(s)`, continuous scalar trace, and
the scalar jump $[N]=2S_{,\Psi}$. Write
$v=r_s$, $\beta_\pm=\widehat n_{E,\pm}r$ and
$K_{\tau,\pm}=\widehat K_{E,\pm}{}^s{}_s$.
Define the three mismatches

$$
\mathcal A_\Sigma=[\beta]+Sr,\qquad
\mathcal T_\Sigma=[K_\tau]+S,\qquad
\mathcal F_\Sigma=2S\left(\overline K_\tau+
\frac{2\overline\beta}{r}\right)
+2S_{,\Psi}\overline N-[\widehat V].
$$

Retaining the Einstein residual explicitly, the spherical tangent-normal
identity reads

$$
\beta_{\pm,s}=K_{\tau,\pm}v-
\frac r2\left(\Psi_sN_\pm+\mathcal R_{ns,\pm}\right).
$$

The common scalar tangent and $S_s=S_{,\Psi}\Psi_s$ then give

$$
\mathcal A_{\Sigma,s}=v\mathcal T_\Sigma-
\frac r2[\mathcal R_{ns}].
$$

Gauss, without setting either extrinsic mismatch to zero, gives

$$
[\mathcal R_{nn}]=-\mathcal F_\Sigma+
\frac2r\left(\overline K_\tau+
\frac{\overline\beta}{r}\right)\mathcal A_\Sigma
+\frac{2\overline\beta}{r}\mathcal T_\Sigma.
$$

Under the already anchored outgoing constraints
$\mathcal E_{V,+}=\mathcal E_{Q_-,-}=0$, put
$J_+=\mathcal A_{\rm in}$ and $J_-=\mathcal B_{\rm in}$ as defined
above. The two residual identities become

$$
\mathcal A_{\Sigma,s}=v\mathcal T_\Sigma-(J_++J_-),
$$
$$
J_--J_+=-\frac r2\mathcal F_\Sigma+
\left(\overline K_\tau+\frac{\overline\beta}{r}\right)
\mathcal A_\Sigma+\overline\beta\mathcal T_\Sigma.
$$

An independent boundary choice is therefore
$\mathcal E_{U,+}=0$, $\mathcal T_\Sigma=0$ and
$\mathcal F_\Sigma=0$, together with the initial angular match
$\mathcal A_\Sigma(s_0)=0$. Any sufficiently regular solution of this
boundary system satisfies

$$
\mathcal A_{\Sigma,s}=-
\left(\overline K_\tau+\frac{\overline\beta}{r}\right)
\mathcal A_\Sigma,
\qquad\text{hence}\qquad
\mathcal A_\Sigma=J_-=0.
$$

The conclusion uses a locally integrable coefficient and positive wall
radius; it divides by neither `v` nor `beta_plus`. Thus angular Israel
and the cap incoming constraint propagate from the initial match. The
argument is noncircular because the residuals were retained until the
independent boundary choice was imposed. It proves this implication for
solutions of that system, not that its boundary equations are solvable,
complete or possess a bounded coupled solution map.

## Missing full-domain estimate

A justified finite-amplitude spacetime interval needs a quantitatively
well-posed coupled problem for the smooth equations, the timelike interface
and the regular centre. The equations above identify its data; they are not
such a theorem. In particular, the following remain to be established.

- Choose an independent characteristic/interface system with anchored gauge
  and a fixed-domain identification. The smooth constraint transport and
  exact source compatibility are established above; supply the new wall
  characteristics through an independently closed interface system and
  preserve the bounded regular-centre residuals. Distributional stress
  conservation must use the transformed scalar-dependent tension if the
  Einstein frame is used.
- Extend the bounded reference scalar collar and centre estimates below
  to compatible coupled coefficients and domain maps, without derivative
  loss or division by `beta_plus`. The reference scalar bounds alone do
  not establish the nonlinear moving-interface estimate.
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

The reference estimates below control scalar transmission on a wall
collar and one scalar derivative at the centre. A compatible coupled
inverse additionally needs metric-constraint and wall maps in the same
regularity class. On a
regular cap radial coordinate `x` with `r_J/x>=a_min>0`, the centre condition
`c(T,0)=0` gives
$\|c/r_J\|_\infty\le a_{\min}^{-1}\|c_x\|_\infty$.
A first-commuted scalar energy controlling `c` in spatial `H2` and its
time derivative in `H1` can therefore control this readout; uncommuted
energy alone cannot. In a wall-adapted Gaussian-normal collar, coordinate
tangential differentiation commutes with the normal derivative, so both
homogeneous scalar transmission conditions commute. Bounds for the
commutator and centre patch are given below on the reference geometry;
their coupled coefficient and patch-matching bounds remain open.
The compact preparation is `C3` and piecewise `W4,infinity`; unlimited
differentiability must not be used to close that estimate.

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
These hypotheses must be bounded rather than inferred from transmission.
The following reference bounds do so locally. The collar need not reach
the centre: matching to the separate centre patch and bounding coefficients
generated by the coupled Einstein and wall equations remain necessary.

### A bounded vacuum reference collar

The coefficient hypotheses above can be made explicit for the vacuum
reference `alpha=0`. This result concerns a Gaussian wall collar, not the
sourced metric or the centre. Put
$\delta=1/20$, $9/20\le x(s)\le11/20$, and orient `n` from cap to
exterior. The two collar pieces are $-\delta\le n\le0$ and
$0\le n\le\delta$. For the reference static functions

$$
f_\pm(r)=1-\frac{2m_{0,\pm}}r-h_\pm^2r^2,
\qquad (m_{0,-},h_-)=(0,1),\quad (m_{0,+},h_+)=(2/5,0),
$$

the normal geodesic and Jacobi equations are

$$
r_{nn}=\frac12f_\pm'(r),\qquad
a_{nn}=\frac12f_\pm''(r)a,\qquad
r(s,0)=x(s),\quad r_n(s,0)=\beta_\pm(x),\quad
a(s,0)=1,\quad a_n(s,0)=\kappa_\pm(x).
$$

Here `kappa` denotes the reference normal acceleration, not the gravitational
coupling. The signed reference formulas are

$$
\beta_- =\frac1{5x^2}+\frac{3x}4,\qquad
\beta_+ =\frac1{5x^2}-\frac{5x}4,\qquad
\kappa_\pm=\frac{d\beta_\pm}{dx},\qquad
v=x_s>0,\quad v^2=\beta_-^2-1+x^2.
$$

They give the deliberately coarse rational bounds

$$
0<v<\frac43,\quad |v_s|<5,\quad
0<\beta_-<\frac43,\quad |\beta_+|<\frac12,\quad
-4<\kappa_-<0,\quad -6<\kappa_+<0,\quad
|\kappa_{\pm,s}|<40,\quad |\kappa_{\pm,ss}|<700.
$$

For example,
$v_s=-2/(25x^5)-3/(20x^2)+25x/16$,
$\kappa_s=6v/(5x^4)$ and
$\kappa_{ss}=-24v^2/(5x^5)+6v_s/(5x^4)$, so these bounds follow on
the stated compact radius interval without sampling a computed solution.
The cap collar is explicit:

$$
r=x\cos n+\beta_-\sin n,\qquad
a=\cos n+\kappa_-\sin n.
$$

Consequently
$r\ge(9/20)(1-\delta^2/2)-(4/3)\delta>3/8$,
$1-\delta^2/2\le a\le6/5$,
$|a_s|\le2$, $|a_n|\le81/20$,
$|a_{sn}|\le40$ and $|a_{ss}|\le35$.

For the exterior, bootstrap $2/5\le r\le3/5$. Since
$0<r_{nn}\le5/2$, integration gives the strictly interior bounds

$$
r\ge\frac9{20}-\frac\delta2=\frac{17}{40}>\frac25,
\qquad
r\le\frac{11}{20}+\frac\delta2+\frac54\delta^2
=\frac{37}{64}<\frac35.
$$

Set $Q=4/(5r^3)\le25/2$. As long as `a` is positive,
$a_{nn}=-Qa$, $\kappa_+<0$ imply $a\le1$, while its Volterra
equation gives

$$
a\ge1-6\delta-\frac{25}{4}\delta^2
=\frac{219}{320}>\frac23,\qquad
|a_n|\le6+\frac{25}{2}\delta=\frac{53}{8}.
$$

These strict inequalities close the radius and lapse bootstraps on each
normal segment. The Jacobi identity $r_s=va$ gives
$|Q_s|\le125$. With $A_1=\|a_s\|_\infty$ on the exterior collar,
differentiating the Volterra equation yields

$$
A_1\le40\delta+
\frac{\delta^2}{2}\left(125+\frac{25}{2}A_1\right),
\qquad
A_1\le\frac{46}{21}<\frac94,
$$
$$
|a_{sn}|\le40+\delta\left(125+\frac{25}{2}\frac94\right)
=\frac{1525}{32}<48.
$$

Next $|r_{ss}|\le|v_s|a+|v a_s|<8$ gives
$|Q_{ss}|\le12Q r_s^2/r^2+3Q|r_{ss}|/r\le7250/3$.
Thus $A_2=\|a_{ss}\|_\infty$ obeys

$$
A_2\le700\delta+\frac{\delta^2}{2}
\left(\frac{7250}{3}+2\cdot125\cdot\frac94+
\frac{25}{2}A_2\right),\qquad
A_2\le\frac{7435}{189}<40.
$$

The smooth reference quotient potential is `W=0` in the cap and
$W=4/(5r^3)$ outside. Both sides therefore satisfy the usable bounds

$$
\frac23<a\le\frac65,\qquad
|a_s|<\frac94,\quad |a_n|\le\frac{53}{8},\quad
|a_{sn}|<48,\quad |a_{ss}|<40,\quad
0\le W\le\frac{25}{2},\quad |W_s|\le125.
$$

Positive `a` excludes a focal degeneracy of this local normal chart. It
does not prove that the entire finite characteristic domain fits one
injective Gaussian chart. Apply the estimates only on compatible collar
patches lying inside the available field domain; no field extension
through an artificial null endpoint is supplied by these geometry bounds.

For a field or residual equation `Lc=j` on fixed normal-coordinate
intervals in such a patch, with `[c]=[c_n]=0` at `n=0`, define the sum
of the two one-sided energies

$$
E[c]=\frac12\sum_\pm\int_{I_\pm}
\left(\frac{c_s^2}{a}+a c_n^2+a c^2\right)dn.
$$

Multiplication by $-a c_s$ and integration give

$$
\frac{dE}{ds}=\mathcal B_{\rm out}
+\sum_\pm\int_{I_\pm}
\left[\frac{a_s}{2a^2}c_s^2+
\frac{a_s}{2}(c_n^2+c^2)+a(1-W)cc_s-a j c_s\right]dn.
$$

Here $\mathcal B_{\rm out}$ is the sum of the outward endpoint terms
$[a c_s c_n]_{\partial I_\pm}$ after cancelling the internal wall terms.
It must be retained or bounded from supplied boundary traces. The regional
bounds are $|a_s/a|\le27/8$ and $a|1-W|\le23/2$; the latter uses
`a<=1` on the exterior and `W=0` on the cap. Hence

$$
\boxed{\frac{dE}{ds}\le\mathcal B_{\rm out}+\frac{119}{8}E
+\sqrt{2E}\,\|a^{3/2}j\|_{L^2(I_-\cup I_+)}.}
$$

This is the energy inequality for one scalar equation. Its use for `c_s`
also includes the displayed commutator in that equation's forcing.
Moving or clipped outer endpoints require their Reynolds flux terms;
their omission is not justified by cancellation at the material wall.
The reference collar estimate bounds neither the central patch nor the
sourced metric coefficients, full coupled inverse, or nonlinear remainder.

### Closed first-commuted reference estimate

The commutator and coefficient bounds give a closed estimate for
$\mathcal H=E[c]+E[v]$, $v=c_s$. This uses the same homogeneous
transmission conditions and compatible fixed collar patches. Let

$$
J^2=\|a^{3/2}j\|_2^2+\|a^{3/2}j_s\|_2^2,
\qquad R_a=\frac{27}{8},\qquad
X=(v/\sqrt a,\sqrt a\,c_n,\sqrt a\,c),\quad
Y=v_s/\sqrt a.
$$

Then $\|X\|_2^2=2E[c]$ and $\|Y\|_2^2\le2E[v]$.
The contribution of the commutator to the differentiated energy is

$$
\int\left[2\frac{a_s}{a}Y^2+
\left(\frac{a_{ss}}a-3\frac{a_s^2}{a^2}\right)X_1Y+
\left(a_{sn}-\frac{a_na_s}a\right)X_2Y-aW_sX_3Y\right]dn.
$$

The three cross coefficients are bounded respectively by
$6027/64$, $4503/64$ and `125`. The last bound uses `W_s=0` in the
cap and `a<=1` outside. Define

$$
M^2=\left(\frac{6027}{64}\right)^2+
\left(\frac{4503}{64}\right)^2+125^2.
$$

Cauchy--Schwarz bounds the displayed contribution by
$4R_aE[v]+2M\sqrt{E[c]E[v]}$. In particular, the self term is
$4R_aE[v]=(27/2)E[v]$ before combining the energies. The associated
matrix on $(\sqrt{E[c]},\sqrt{E[v]})$ is

$$
\begin{pmatrix}0&M\\M&4R_a\end{pmatrix},\qquad
\lambda_{\max}=2R_a+\sqrt{4R_a^2+M^2},\qquad
4R_a^2+M^2=\frac{120788362}{4096}<172^2.
$$

Adding the two original energy inequalities and combining the prescribed
forcing terms therefore gives the closed bound

$$
\boxed{\mathcal H'\le\mathcal B+
\frac{1549}{8}\mathcal H+\sqrt{2\mathcal H}\,J,\qquad
\mathcal B=\sum_\pm[a(c_sc_n+v_sv_n)]_{\partial I_\pm}.}
$$

The internal wall terms cancel. With enough piecewise regularity,
`[c]=[c_n]=0` implies `[v]=[v_n]=0` by tangential differentiation;
the forcing needs its weak time derivative `j_s` on each spatial side,
not a continuous cross-wall value. A temporal jump cannot be omitted
from that derivative. No distributional interface source is included in
this homogeneous statement.

For clarity about the remaining boundary input, write the fixed outer
endpoints as `n_L,n_R`. Their incoming characteristic traces give

$$
\mathcal B\le\mathcal I(s)=\frac14\sum_{z=c,v}
\left[(z_s-a z_n)^2\big|_{n_L}+
(z_s+a z_n)^2\big|_{n_R}\right].
$$

These traces must be supplied or bounded by compatible source or adjacent
patch data; they are not set to zero. The traces for `v` must agree with
`v=c_s` and are not a second freely chosen incident wave. Young's
inequality yields, with $\gamma=1557/8$,

$$
\mathcal H(s)\le e^{\gamma(s-s_0)}\mathcal H(s_0)
+\int_{s_0}^s e^{\gamma(s-\sigma)}
\left[\mathcal I(\sigma)+\tfrac12J(\sigma)^2\right]d\sigma.
$$

This is a finite reference-collar amplification and forcing estimate.
It requires finite compatible initial energy, square-integrable forcing
and commuted forcing, and the displayed boundary-trace control. Moving
endpoints or nonzero interface residuals add their respective fluxes.
It supplies an a priori estimate and uniqueness for zero inputs in this
class, not a construction of the nonlinear coupled solution.

### Reference centre patch and central trace

The reference cap admits a separate regular patch. Put
$y=\operatorname{atanh}r$, $A=1-r^2=\operatorname{sech}^2y$ and
$h=A(-dT^2+dy^2)$. Its reduced potential is zero. For a normalized
reduced-field error or forced field `e`, use the explicit convention

$$
L_0e=j,\qquad e_{TT}-e_{yy}=F,\qquad F=-Aj,\qquad e(T,0)=0.
$$

The regular central scalar readout is
$p_e(T)=\lim_{r\to0}e/r=e_y(T,0)$. Regular spherical fields and
forcings have odd extensions at `y=0`; in particular `F(T,0)=0` where
classical traces are used.

The finite labels are those of the
[cap domain](nsc-atlas-cap-response.md#finite-domain-and-verification):
$q_a=q(V_a)$, $q_b=q(V_b)$,
$T_a=q_a-\operatorname{atanh}x_a$,
$T_b=q_b-\operatorname{atanh}x_b$, and
$q_{\min}=q_a-2\operatorname{atanh}x_a$, with
$x_a=9/20$, $x_b=11/20$. Choose

$$
0<Y<\operatorname{atanh}x_a,\qquad
q_a<T_1<q_b-Y,\qquad T_0=q_a-Y.
$$

The strip $T_0\le T\le T_1$, $0\le y\le Y$ lies inside the supplied
cap: `T0>Ta`, its radius is below `x_a` before `T_b`, and
`T+y<q_b` puts it inside the terminal inward cone afterward.
Its outgoing labels satisfy $T-y\ge q_a-2Y>q_{\min}$.
This restriction keeps the readout on a compact interior interval and
does not require a later wall history or a zero tail after `q_b`.

Let $R=e_T+e_y$. Since $(\partial_T-\partial_y)R=F$, the classical
incoming characteristic identity is

$$
p_e(T)=R(T_*,T-T_*)+\int_{T_*}^T F(s,T-s)\,ds,\qquad
T_*=\max(T_0,T-Y).
$$

The first term is either the initial trace
$R(T_0,y)=e_1(y)+e_0'(y)$, with $e_0=e(T_0,\cdot)$ and
$e_1=e_T(T_0,\cdot)$, or the incoming outer trace
$g(T-Y)=R(T-Y,Y)$. This `g` is a matched cap trace, not an additional
physical source and not an automatically vanishing boundary condition.

Define the two centre-patch energies

$$
E_0=\frac12\int_0^Y(e_T^2+e_y^2)dy,\qquad
E_1=\frac12\int_0^Y(e_{TT}^2+e_{Ty}^2)dy,\qquad E_C=E_0+E_1.
$$

The axis fluxes vanish, and direct integration gives

$$
E_C'=\mathcal B_Y+
\int_0^Y(F e_T+F_Te_{TT})dy,\qquad
\mathcal B_Y=(e_Te_y+e_{TT}e_{Ty})\big|_Y,
$$
$$
E_C'\le\mathcal B_Y+\sqrt{2E_C}
\left(\|F\|_2^2+\|F_T\|_2^2\right)^{1/2},\qquad
\mathcal B_Y\le\frac{g^2+g_T^2}{4}.
$$

Thus initial energy, `F,F_T` and compatible incoming traces `g,g_T`
control this patch. In particular, applying Young's inequality gives an
explicit Gronwall bound with energy coefficient `1` and forcing
$(g^2+g_T^2)/4+(\|F\|_2^2+\|F_T\|_2^2)/2$.
The initial commuted energy is

$$
E_1(T_0)=\frac12\int_0^Y
\left[(e_0''+F(T_0,y))^2+(e_1')^2\right]dy.
$$

The prepared leading field is quiet on the initial strip because its two
profile labels are at most `q_a`. For an error or residual field, however,
quiet `e_0,e_1` alone does not imply zero `E_1(T0)` if the residual forcing
there is nonzero. Full quiet initial energy requires compatible quiet
forcing as well.

Finally, for $w=e_y$ the exact one-dimensional identity

$$
w(0)=\frac1Y\int_0^Yw(z)dz-
\frac1Y\int_0^Y(Y-z)w'(z)dz
$$

and $e_{yy}=e_{TT}-F$ give the pointwise central bound

$$
\boxed{|p_e(T)|\le\sqrt{\frac{2E_0(T)}Y}
+\sqrt{\frac Y3}\left(\sqrt{2E_1(T)}+\|F(T,\cdot)\|_2\right).}
$$

This provides the missing central trace estimate for the reference forced
scalar patch. Matching its incoming traces to collar or intervening patches,
and bounding the metric-generated nonlinear forcing, remain separate
requirements. Neither this estimate nor the collar estimate substitutes
sampled coefficient sizes for those bounds.

### What remains for nonlinear validity

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
constraint compatibility. The bounded reference collar and central
readout estimates are established components; compatible coupled inverse,
coefficient and residual bounds remain necessary for a nonlinear
certificate.

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
