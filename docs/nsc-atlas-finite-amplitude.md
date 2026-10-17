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

### Normalized positive-radius bulk operator

The invariant field $c=\chi/\varepsilon$ retains the prepared scalar
information at zero amplitude. For $r,C>0$, put

$$
u=c/r,\quad u_a=(c_a-u r_a)/r,\quad
q=\alpha u^2/6,\quad D=1+q,\quad
P=3h^2D^2,\quad \sigma_{ab}=\alpha u_a u_b/D.
$$

Here $\Psi=\sqrt6\operatorname{asinh}(\sqrt\alpha u/\sqrt6)$,
so $\Psi_a\Psi_b=\sigma_{ab}$. The exact mixed system becomes

$$
r_{UV}=-r_Ur_V/r-C/(4r)+CrP/4,\qquad
(\log C)_{UV}=2r_Ur_V/r^2+C/(2r^2)-\sigma_{UV},
$$
$$
c_{UV}=-Kc,\qquad
K=\frac{Cm}{2r^3}-\frac{CP}{12}-\frac{\sigma_{UV}}6,\qquad
m=\frac r2\left(1+\frac{4r_Ur_V}{C}\right).
$$

An independent canonical-field chain rule gives

$$
u_{UV}=-\frac{r_Uu_V+r_Vu_U}{r}
-\frac{CPu}{6}+\frac{\alpha u\,u_Uu_V}{6D};
$$

differentiating $c=ru$ recovers the displayed reduced-field equation.
No expression divides by $\sqrt\alpha$. For numerical evaluation, form
$\Delta P=3h^2q(2+q)$ directly and use

$$
K=\frac{r_Ur_V+C(1-h^2r^2)/4}{r^2}
-\frac{C\Delta P}{12}-\frac{\sigma_{UV}}6.
$$

This exposes the vacuum cap cancellation and avoids obtaining the potential
increment by subtracting nearly equal values. It does not regularize arbitrary
first jets at the axis; the positive-radius operator excludes $r=0$.

First jets determine the right sides of the null constraints, but cannot
measure their residuals $E_U,E_V$, which contain pure second derivatives.
The often used mass rates

$$
M_U=\frac{r^2P}{2}r_U-\frac{r^2r_V}{C}\sigma_{UU},\qquad
M_V=\frac{r^2P}{2}r_V-\frac{r^2r_U}{C}\sigma_{VV}
$$

are derivatives of the geometric mass only on the constraint surface.
Direct differentiation using the mixed equation gives instead

$$
m_U=M_U+\frac{2rr_V}{C}E_U,\qquad
m_V=M_V+\frac{2rr_U}{C}E_V.
$$

The [bulk implementation](../src/recursive_horizons/nsc_atlas_nonlinear_bulk.py)
keeps these qualifications explicit in its outputs. Its finite first-jet
evaluation supplies the regional equations, not a spacetime trajectory or
a freely selectable mass history.

### Bounded characteristic cell

A local cell evolves $X=(r,L,c)$, $L=\log C$, from three known corners.
Let `SE` advance the first null coordinate by $\Delta U>0$, and `NW`
advance the second by $\Delta V>0$, relative to `SW`. With unknown
$z=X_{NE}$, the implicit-midpoint equation is

$$
\mathcal R(z)=z-X_{NW}-X_{SE}+X_{SW}
-\Delta U\Delta V\,\mathcal B(X_m,X_{U,m},X_{V,m})=0,
$$
$$
X_m=\frac{z+X_{NW}+X_{SE}+X_{SW}}4,\qquad
X_{U,m}=\frac{z+X_{SE}-X_{NW}-X_{SW}}{2\Delta U},\qquad
X_{V,m}=\frac{z+X_{NW}-X_{SE}-X_{SW}}{2\Delta V}.
$$

Here $\mathcal B$ is the three-component mixed bulk operator, with
$C_m=e^{L_m}$. The cell coordinates may be the source's affine null
coordinates below; they need not be the original $U,V$ labels.

The [cell implementation](../src/recursive_horizons/nsc_atlas_characteristic_cell.py)
uses a bounded damped Newton iteration and the analytic directional
derivative of these formulas. Each variation in the unknown gives
$dX_m=dz/4$, $dX_{U,m}=dz/(2\Delta U)$ and
$dX_{V,m}=dz/(2\Delta V)$, followed by the same normalized field and
potential chain rules. Its error scales are fixed from the supplied corners
and reference scales; a growing iterate cannot enlarge its own tolerance.
A declared compact state tube, positive radius and representable positive
lapse are checked at corners, midpoint and candidate endpoint.

The default limits are twelve Newton updates and eight trial step lengths
per update, with separately counted operator and Jacobian evaluations and
a shared finite CPU budget. Failure reports the iterate and residual
diagnosis; it does not silently change cell size, amplitude or corner data.
The final residual and correction control the local algebraic solve. They
do not themselves certify Einstein constraints, regional convergence or a
continuum error bound.

Minkowski geometry with constant canonical scalar gives an exact affine
control. Curved analytic Schwarzschild and zero-amplitude cap corners give
an independent local refinement test: balanced halving reduces the local
midpoint defect and one-cell endpoint error at fourth order before
roundoff. A complete smooth grid has second-order accuracy; the local
order is not a fourth-order regional claim. A later prepared-source march
must also measure propagated constraints and cumulative error.

### Original exterior rectangle and independent diagnostics

The characteristic input below determines the finite rectangle

$$
0\le\bar U\le\bar U_b=C_{J0}(U_a-U_0),\qquad
x_a=0.4\le x\le x_b=x_{\rm last},\qquad
\bar U=C_{J0}(U-U_0),\quad C_{J0}=-3.2/U_0.
$$

Its lower sides are the exact source at $\bar U=0$ and the quiet
Schwarzschild bottom at $x=x_a$. A future timelike wall starting at
$(\bar U_b,x_a)$, parametrized by Einstein proper wall time $s$, obeys
$C\bar U_s x_s=1$ with both rates positive. While
that continuation and chart exist, the wall lies beyond $\bar U_b$;
the fixed rectangle has no incoming characteristic from that later wall.
This causal statement allows the regional calculation before the coupled
wall march. It does not itself establish the latter's existence interval.

The [regional midpoint march](../src/recursive_horizons/nsc_atlas_nonlinear_rectangle.py)
stores $X_{ij}=(r,L,c)$ with $\bar U$ as
the first array axis and $x$ as the second. For the cell ending at
$(i,j)$, its three inputs are $X_{i-1,j-1}$, $X_{i-1,j}$ and
$X_{i,j-1}$. Boundary first jets are retained separately as independent
comparators; they are not fitted from the stored scalar waveform.

To measure null constraints independently of the mixed evolution equation,
use centered three-point first and pure second differences on interior
nodes. At the low endpoint of a uniform coordinate with spacing $h$, use

$$
D f_0=\frac{-3f_0+4f_1-f_2}{2h},\qquad
D^2f_0=\frac{2f_0-5f_1+4f_2-f_3}{h^2},
$$

with reflected formulas at the high endpoint. For each $a=\bar U,x$,
the nodal diagnostics are

$$
u=\frac cr,\qquad D=1+\frac{\alpha u^2}{6},\qquad
u_a^{\rm FD}=\frac{D_a c-uD_a r}{r},
$$
$$
E_a^{\rm FD}=D_{aa}r-(D_aL)(D_ar)
+\frac{\alpha r}{2D}(u_a^{\rm FD})^2,\qquad
m^{\rm FD}=\frac r2\left(1+4(D_{\bar U}r)(D_xr)e^{-L}\right).
$$

Pure second derivatives in these expressions are never replaced by the
solved mixed right-hand side. Polynomial data with a prescribed nonzero
constraint provide an omission control. The zero-amplitude geometry gives
an independent full-grid comparator. With
$p=(1-\bar U/3.2)(1-x/0.8)e^{x/0.8}$,

$$
r_0=0.8\,[1+W_0(-p/e)],\qquad
L_0=\log(x/r_0)+(x-r_0)/0.8,\qquad m_0=0.4.
$$

This principal Lambert branch is the positive trapped-radius branch of the
original null chart. The normalized scalar still propagates at $\alpha=0$;
zero amplitude removes its gravitational forcing, not its limiting field
equation. In particular, the source value can vanish after the pulse while
its transverse derivative $c_{\bar U}$ remains nonzero.

The first admitted calculation compares 16 and 32 cells per direction at
zero amplitude and half the source-only threshold. It retains raw fields,
raw constraint residuals, and finite-minus-zero responses as distinct
quantities. Common interior nodes supply the refinement comparison. Stencils
touching either fixed source support label are flagged separately, alongside
the all-interior diagnostic; boundary defects are compared with the exact
source/bottom first jets. The inherited $C^3$, piecewise $W^{4,\infty}$
profile warrants measuring derivative refinement separately from field
refinement. Two resolutions give sampled differences, not an established
asymptotic rate or a rigorous error enclosure.

The [creation-only recorder](../scripts/record_nsc_atlas_nonlinear_rectangle.py)
retains the exact boundaries, nodal states, local iterations, independent
diagnostics and cumulative source/cell accounting. Its authentication path
recomputes the saved-array diagnostics without a source integration or cell
solve. The [admission](prior-art-and-open-claim.md#nonlinear-exterior-rectangle-16-october-2026)
fixes the domain, four comparisons, tolerances and stopping budget.

The [first record](../results/development/nsc-atlas-nonlinear-rectangle-v1/rectangle.json)
completes all four grids. The [measured comparison](current-result.md#first-nonlinear-exterior-rectangle-16-october-2026)
shows refining metric and constraint diagnostics, with the largest raw
$E_x$ error already present when differencing the known source. To extract
jets for the next characteristic leg, write
$X(\bar U,x)=X_0(x)+\Delta X(\bar U,x)$ and retain the exact source jets:

$$
r_x^{\rm lift}=r_{0,x}+D_x\Delta r,\qquad
r_{xx}^{\rm lift}=r_{0,xx}+D_{xx}\Delta r,\qquad
c_x^{\rm lift}=g_x+D_x\Delta c,\qquad
L_x^{\rm lift}=L_{0,x}+D_x\Delta L.
$$

This differentiates the propagated increment and adds known initial
derivatives; it does not substitute a mixed equation for a pure constraint.
For values between nodes, keep the compact $g(x)$ analytic and interpolate
its propagated correction, so grid placement does not flatten the source
peak. The increment still needs its own refinement control. The supplied
source jets also retain their numerical initialization error; the split is
an exact algebraic identity, not a numerical enclosure. The terminal
leg $\bar U=\bar U_b$ supplies wall-born exterior cells; the hybrid interface
must be evaluated at actual evolved wall events.

### Regular polar-areal cap

Inside an untrapped cap, use the Einstein metric

$$
ds_E^2/\ell^2=-Ae^{2\delta}dT^2+A^{-1}dr^2+r^2d\Omega^2,
\qquad A=1-2m/r>0.
$$

The time $T$ is Jordan central proper time, anchored at the same quiet
event. Define $f=\Psi/\sqrt\alpha$, with its continuous normalized limit,
$\Phi=f_r$, $\Pi=e^{-\delta}f_T/A$ and $B=Ae^\delta$. With
$q_f=\sqrt{\alpha/6}\,f$ and $\operatorname{sinhc}(q)=\sinh(q)/q$
continued to one at zero, the potential and normalized derivative are

$$
P=3\cosh^4q_f,\qquad
U_\alpha(f)=P_{,\Psi}/\sqrt\alpha
=2f\operatorname{sinhc}(q_f)\cosh^3q_f.
$$

In particular $U_0=2f$; the normalized scalar still evolves when its
gravitating amplitude vanishes. The scalar system and radial constraints are

$$
f_T=B\Pi,\qquad \Phi_T=(B\Pi)_r,\qquad
\Pi_T=\frac1{r^2}(r^2B\Phi)_r-e^\delta U_\alpha(f),
$$
$$
m_r=\frac{r^2P}{2}+\frac{\alpha r^2A}{4}(\Pi^2+\Phi^2),\qquad
\delta_r=\frac{\alpha r}{2}(\Pi^2+\Phi^2),
$$
$$
m(T,0)=0,\qquad \delta(T,0)=-\log\cosh q_f(T,0).
$$

These constraints supply the cap wall metric from its field and central
clock. For example, with $W=\Pi^2+\Phi^2$ and
$I(r)=\int_0^r\alpha sW(s)/2\,ds$, the mass has the radial representation

$$
m(r)=e^{-I(r)}\int_0^r e^{I(s)}s^2
\left(\frac{P(s)}2+\frac{\alpha W(s)}4\right)ds,
\qquad \delta(r)=\delta(0)+I(r).
$$

The omitted temporal Einstein equation is a consistency check, not another
free evolution law. Differentiating the radial mass constraint and using
the scalar equations gives

$$
\partial_r\mathcal E_m=-\delta_r\mathcal E_m,\qquad
\mathcal E_m=m_T-\frac{\alpha r^2A^2e^\delta}{2}\Pi\Phi.
$$

Regularity makes $\mathcal E_m(0)=0$, hence it vanishes throughout the cap.
The remaining angular equation follows from the Einstein identity and the
scalar equation in a sufficiently regular solution. Also
$(\Phi-f_r)_T=0$, so compatible initial scalar derivatives stay compatible.

To evaluate the axis without dividing by its radius, use
$\rho=r^2$, $Q=\Phi/r$ and $\mu=m/r^3$, all continued by regular even
profiles. Then

$$
f_T=B\Pi,\qquad Q_T=2(B\Pi)_\rho,\qquad
\Pi_T=3BQ+2\rho(BQ)_\rho-e^\delta U_\alpha(f),
$$
$$
A=1-2\rho\mu,\qquad
3\mu+2\rho\mu_\rho=\frac P2+
\frac{\alpha A}{4}(\Pi^2+\rho Q^2),\qquad
\delta_\rho=\frac\alpha4(\Pi^2+\rho Q^2),
$$
$$
\mu_0=\frac{P_0}{6}+\frac{\alpha\Pi_0^2}{12},\qquad
B_0=e^{\delta_0},\qquad
(\Pi_T)_0=3B_0Q_0-e^{\delta_0}U_\alpha(f_0).
$$

The radial integral supplies $\mu$ away from the axis and this limit
supplies its centre value. The $Q$ equation uses the regular derivative
$(B\Pi)_\rho$, not a separately prescribed centre waveform. The central
Jordan scalar is $p_{J,C}=f_0\operatorname{tanhc}(q_{f,0})$, where
$\operatorname{tanhc}(q)=\tanh(q)/q$ has the continuous value one.
For implementation, the same mass integral has an explicit regular form.
Writing $W_\rho=\Pi^2+\rho Q^2$ and
$\gamma(\rho)=\alpha\int_0^\rho W_\xi\,d\xi/4$ gives

$$
\mu(\rho)=e^{-\gamma(\rho)}\int_0^1 t^2e^{\gamma(t^2\rho)}
\left[\frac P2+\frac{\alpha W_\rho}{4}\right]_{t^2\rho}dt,
\qquad B_\rho=e^\delta(\mu-P/2).
$$

The last identity follows by cancellation in the radial constraints and
avoids differentiating sampled metric values. Likewise,

$$
\frac{P-3}{\alpha}=\frac{f^2}{2}\operatorname{sinhc}^2q_f
\left(\cosh^2q_f+1\right)
$$

has a continuous value at zero amplitude and avoids subtracting two vacuum
potentials. The finite domain must retain $A>0$ and regular profiles. The
owned pulse is $C^3$ with piecewise bounded fourth derivatives: use the
equations piecewise or weakly at its fronts, and retain the weighted
$\rho Q_\rho$ term instead of assuming $Q_\rho$ alone is bounded there.
The characteristic speeds in $\rho$ degenerate at the axis; these identities
do not establish uniform Cartesian hyperbolicity. A stable discretization
and its error estimate remain to be constructed.

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

### Completing the source first jets

The exact source ODE also supplies all first jets required to start the
mixed regional system. In the original affine coordinate $x$, let
$C_0=C_{J,0}=-3.2/U_0>0$, and temporarily use
$\bar U=C_0(U-U_0)$. The source is at $\bar U=0$ and its Einstein
lapse is $C=\mathcal F$.
With $g=g_{\rm src}$, the known longitudinal traces are

$$
r=\sqrt{r_J^2-\alpha g^2/6},\qquad
r_x=\frac{r_Jr_J'-\alpha gg'/6}{r},\qquad
c=g,\quad c_x=g',\quad
\mathcal F=\frac{r^2}{r_J^2},\quad
L_x=2r_x/r-2r_J'/r_J.
$$

Introduce $a=r_{\bar U}$, $b_c=c_{\bar U}$ and
$\lambda=L_{\bar U}$; $b_c$ is unrelated to the wall parameter $b=2$.
The mixed exterior equations give three ordinary differential equations:

$$
a'=-ar_x/r-\mathcal F/(4r),\qquad
b_c'=-K g,\qquad
\lambda'=2ar_x/r^2+\mathcal F/(2r^2)
-\frac{\alpha u_{\bar U}u_x}{D},
$$
$$
u_{\bar U}=(b_c-u a)/r,\qquad u_x=(g'-u r_x)/r,\qquad
K=\frac{ar_x+\mathcal F/4}{r^2}
-\frac{\alpha u_{\bar U}u_x}{6D}.
$$

All three corner values follow from the same quiet Schwarzschild corner:

$$
a(x_a)=\frac{0.8/x_a-1}{4},\qquad b_c(x_a)=0,\qquad
\lambda(x_a)=-\left(\frac1{0.8}+\frac1{x_a}\right)a(x_a).
$$

Multiplying these transverse derivatives by $C_0$ restores the original
$U$ derivatives. The change from $V$ to $x$ has a Jacobian depending only
on $x$, so it introduces no additional transverse lapse term. The geometric
mass follows from these data; neither a mass function nor another scalar
history is free.

For stable evaluation near zero amplitude, the source radius can instead
be integrated through $r_J=x+\alpha h$, $r_J'=1+\alpha w$ with
$h'=w$ and

$$
w'=\frac{gg''-2(g')^2}{6r_J}
+\frac{gg'(1+\alpha w)}{3r_J^2}.
$$

At positive amplitude these are exact scaled changes. At zero amplitude
they give the quadratic sensitivity without subtracting two nearby radii.
The existing analytic source interface returns $p=g/x,p',p''$; recover
$g=xp$, $g'=p+xp'$ and $g''=2p'+xp''$, rather than differentiating stored
samples. The [exact source implementation](../src/recursive_horizons/nsc_atlas_nonlinear_source.py)
integrates this system with bounded DOP853 steps, split at the two support
endpoints. It supplies bounded dense first-jet queries and the exact quiet
bottom in the same shifted affine coordinates. For small-amplitude lapse
arithmetic it uses $L=\log(1-\alpha p_J^2/6)$ through `log1p` and
$L_x=-\alpha p_J(p_J)_x/(3\mathcal F)$, where $p_J=g/r_J$.
The analytical null-constraint outputs check frame consistency; they are
not an independent enclosure of the numerical ODE defect.

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
The owned compact profile gives explicit constants below.

The same ODE permits a source-only remainder estimate. On a specified
nonsingular tube, the exact vector field and its derivatives have finite
Lipschitz bounds. Substitution of the exact quadratic coefficient
`(x+alpha h, 1+alpha h')` gives a residual of order `alpha^2`; a Gronwall
estimate then bounds the source radius/slope error with an explicit tube
constant. Numerical approximations to `h` would add their own enclosed
error. This argument concerns the initial leg only: it proves neither a
two-sided evolved solution nor a central-waveform remainder.

### Continuous bounds for the prescribed source

For the unchanged compact profile $g=B(z(x))$, where
$B(z)=(1-z^2)^4$ on $|z|<1$ and zero elsewhere, direct differentiation gives

$$
BB_z=-8z(1-z^2)^7,\qquad
BB_{zz}-2B_z^2=-8(1-z^2)^6(1+9z^2).
$$

Their absolute maxima on $[-1,1]$ are

$$
A_1=\frac8{\sqrt{15}}\left(\frac{14}{15}\right)^7,\qquad
A_2=\frac{80}{7}\left(\frac{20}{21}\right)^6,
$$

attained at $z^2=1/15$ and $z^2=1/21$, respectively. These are continuous
polynomial bounds, not maxima estimated from a sampling grid.
Put $a=0.8$, $\Delta=|U_0|(V_b-V_a)>0$, $x_b=x_{\rm last}$, and

$$
Z_1=\frac{4x_b e^{x_b/a}}{a^2\Delta},\qquad
Z_2=\frac{4e^{x_b/a}(1+x_b/a)}{a^2\Delta}.
$$

Both $z_x$ and $z_{xx}$ are positive and increasing on the source leg,
so these endpoint values bound them on the compact support. The chain
rule therefore supplies

$$
M_0=1,\qquad M_1\le A_1Z_1,\qquad
M_2\le A_2Z_1^2+A_1Z_2.
$$

The zero extension contributes no extra first or second derivative at the
support endpoints. With $r_{\min}=1/5$, $v_{\max}=3/2$ and
$\mathcal F_{\min}=1/2$, define

$$
\overline L=\frac56(A_2Z_1^2+A_1Z_2)+\frac{25}{2}A_1Z_1,
\qquad L_x=x_b-0.4.
$$

The original validated interval has $0<L_x<3/20$. Thus the radius
condition is weaker than the slope condition, and the strict source-only
threshold is

$$
0\le\alpha<\alpha_*=
\min\left\{\frac1{2\overline L L_x},\frac3{25}\right\}.
$$

The deterministic initialization control $\alpha_{\rm test}=\alpha_*/2$
has $|r_J'-1|\le1/4$, $|r_J-x|\le L_x/8$ and
$\mathcal F\ge3/4$ from the same bootstrap argument. These bounds choose
the test using only its prescribed source and interval. Evaluation with
numerical reference coordinates inherits their numerical accuracy; it is
not an interval-arithmetic certificate for the geometric quadrature.
Neither this threshold nor a successful source integration establishes a
finite-amplitude interval for the evolved cap and wall.

The initialization control evaluates $\alpha_*=2.5966254285102885\times10^{-4}$
and $\alpha_{\rm test}=1.2983127142551443\times10^{-4}$ from the existing
order-64 reference geometry. Its shared-budget comparison tightens relative
and absolute ODE tolerances from $(2\times10^{-10},10^{-12})$ to
$(2\times10^{-12},10^{-14})$ at 33 source cuts. The maximum changes are
$1.44\times10^{-14}$ in Einstein radius, $9.82\times10^{-13}$ in its
longitudinal derivative, $4.72\times10^{-15}$ in its transverse derivative,
$4.02\times10^{-14}$ in $c_{\bar U}$ and $1.71\times10^{-13}$ in geometric
mass. The scaled Jordan slope response changes by $7.56\times10^{-9}$.
The two solves use 828 and 1,200 RHS evaluations; geometry, validation and
queries bring the shared count to 3,063. Their analytical constraint
consistency residual is below $1.06\times10^{-15}$. These are source-only
numerical controls with inherited geometric accuracy, not a rigorous ODE
error enclosure or a propagated spacetime result.

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

### An explicit outgoing-to-incoming wall map

The [boundary and clock helpers](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
implement the following algebra; their
[independent controls](../tests/test_nsc_atlas_nonlinear_boundary.py)
exercise the vacuum limit, canonical wall action and coordinate changes.

There is also a pointwise map that imposes angular, temporal and normal
matching independently at every wall time. It is distinct from the
preceding choice that imposes an incoming constraint and propagates the
angular match. Keep `r` as the common Einstein radius, `c` as the invariant
normalized field and `s` as Einstein proper wall time. Define the weighted
outgoing traces

$$
R_+=v+\beta_+=2\dot V_+r_{V,+},\qquad
R_-=v-\beta_-=2\dot Q_-r_{Q_-,-},
$$
$$
C_+=c_s+n c=2\dot V_+c_{V,+},\qquad
C_-=c_s-n c=2\dot Q_-c_{Q_-,-}.
$$

The notation `C_+,C_-` in this paragraph denotes scalar traces, not metric
lapses. Put

$$
u=c/r,\qquad D=1+\alpha u^2/6,\qquad
\mathcal F=D^{-1},\qquad S=bD^{3/2}.
$$

Angular Israel and homogeneous reduced-field transmission give

$$
v=\frac{R_++R_-+Sr}{2},\qquad
\overline\beta=\frac{R_+-R_-}{2},\qquad
\beta_+=\overline\beta-Sr/2,\quad
\beta_-=\overline\beta+Sr/2,
$$
$$
c_s=\frac{C_++C_-}{2},\qquad n c=\frac{C_+-C_-}{2}.
$$

The exact canonical derivatives, evaluated by their continuous normalized
form even at zero amplitude, are

$$
\frac{\Psi_s}{\sqrt\alpha}=\frac{c_s-uv}{r\sqrt D},\qquad
\frac{N_\pm}{\sqrt\alpha}=\frac{n c-u\beta_\pm}{r\sqrt D}.
$$

The scalar jump follows directly from the angular jump:
$[N]=\sqrt\alpha\,uS/\sqrt D=2S_{,\Psi}$.
Moreover $S_{,\Psi}/S=z/2$, and normal force determines the mean
normal acceleration without dividing by either signed `beta`:

$$
\overline K_\tau=
\frac{3(h_+^2-h_-^2)\sqrt D}{2b}
-\frac{\alpha u(n c-u\overline\beta)}{2rD}
-\frac{2\overline\beta}{r},\qquad
K_{\tau,+}=\overline K_\tau-S/2,\quad
K_{\tau,-}=\overline K_\tau+S/2.
$$

Here the model remains `b=2,h_+=0,h_-=1`. The incoming radius and scalar
traces are

$$
2\dot U_+r_{U,+}=R_-+Sr,\qquad
2\dot Q_+r_{Q_+,-}=R_++Sr,
$$
$$
2\dot U_+c_{U,+}=C_-,\qquad
2\dot Q_+c_{Q_+,-}=C_+.
$$

For compatible positive wall lapses `C_W,+,C_W,-`, parameterize the
future null rates on each side by
$\dot U=C_W^{-1/2}e^{-\eta}$ and
$\dot V=C_W^{-1/2}e^{\eta}$, using `Q_-,Q_+` on the cap.
Let $L=\log C$, $\Gamma_\pm=(\log C_{W,\pm})_s$,
$\ell_+=2\dot V_+L_{V,+}$ and
$\ell_-=2\dot Q_-L_{Q_-,-}$. Then

$$
2\dot U_+L_{U,+}=2\Gamma_+-\ell_+,\qquad
2\dot Q_+L_{Q_+,-}=2\Gamma_--\ell_-,
$$
$$
nL_+=\ell_+-\Gamma_+,\qquad nL_-=\Gamma_--\ell_-,\qquad
\eta_{\pm,s}=K_{\tau,\pm}-\tfrac12 nL_\pm.
$$

These formulas supply six incoming first derivatives and the common
radius, field, embedding and rapidity rates. Their denominators are
positive radius, lapse and null rates. They are smooth at `alpha=0` and
at the reference zero of `beta_plus`. They are an algebraic boundary
operator, not a construction of the regional evolution or its gauge.
Compatible corner data and the same physical source remain required.

In particular, fixing both cap null labels to Jordan proper time at the
axis already fixes their gauge. In that gauge the cap wall lapse and
`Gamma` must come from compatible regional evolution. Prescribing wall
lapses instead requires an observer-clock conversion; the two gauge
choices cannot be imposed independently on the same null labels.

### Closing a polar cap against the null exterior

The regular polar cap removes that competing lapse prescription. At the
wall its field and radial constraints supply $A_-,\delta_-$ and the
outgoing canonical combination $w_-=\Pi-\Phi$. Keep only the exterior
null lapse trace $C_{W,+}(s)>0$ and its compatible derivative
$\Gamma_+=d_s\log C_{W,+}$ as gauge input. It agrees with the reference
gauge on the quiet initial interval; the source $V$ labels and existing
bottom $U$ labels stay anchored. The remaining relabelling of wall-born
$U$ coordinates realizes this exterior gauge.

The exterior outgoing derivatives give $R_+=2\dot V_+r_{V,+}$ and
$C_+=2\dot V_+c_{V,+}$ as above. Set

$$
H=R_++Sr,\qquad A_->0,\qquad H>0.
$$

The polar cap has $\beta_-^2-v^2=A_-$ and $\beta_->0$.
Angular matching therefore solves the wall velocity algebraically:

$$
\boxed{v=\frac{H^2-A_-}{2H},\qquad
\beta_-=\frac{H^2+A_-}{2H},\qquad
\beta_+=R_+-v.}
$$

In particular $v+\beta_-=H$ and $\beta_--v=A_-/H$; no division
by $v$ or the sign-changing $\beta_+$ occurs. The common field value
is $c=r f\operatorname{sinhc}(q_f)$. With $u=c/r$ and $D$ as before,
transmission gives

$$
C_- = \frac{A_-}{H}\left(r\sqrt D\,w_--u\right),\qquad
w_+=\Pi+\Phi=\frac{C_+/H-u}{r\sqrt D}.
$$

Only $w_+$ is incoming to the cap. From these two combinations recover
$\Pi=(w_++w_-)/2$, $\Phi=(w_+-w_-)/2$, then use
$c_s=(C_++C_-)/2$ and $n c=(C_+-C_-)/2$ in the same normal-force
formula for $\overline K_\tau,K_{\tau,\pm}$.
The wall evolution and exterior incoming derivatives are

$$
r_s=v,\qquad T_s=\frac{\beta_-}{A_-e^{\delta_-}},\qquad
\eta_{+,s}=K_{\tau,+}-\frac{2\dot V_+L_{V,+}-\Gamma_+}{2},
$$
$$
2\dot U_+r_{U,+}=v-\beta_+,\qquad
2\dot U_+c_{U,+}=C_-,\qquad
\dot U_+L_{U,+}=\Gamma_+-\dot V_+L_{V,+}.
$$

The existing positive exterior null rates evolve $U,V$. The cap wall
lapse is reconstructed from its regular interior and central clock,
not prescribed again. The common scalar value is a compatibility condition:

$$
f_s=\beta_-\Pi+v\Phi
=\frac{c_s-uv}{r\sqrt D}.
$$

Thus its transport preserves $c=r f\operatorname{sinhc}(q_f)$ from
the common quiet anchor. Adding an independent Dirichlet waveform for
$f$ would overdetermine the single incoming scalar channel. An expanding
wall admits new fixed-radius worldlines: their scalar value is initialized
from this transported common wall value. Combining its tangent derivative
with $f_T=B\Pi$ preserves $\Phi=f_r$ there. This derived value is needed
alongside the incoming characteristic, but is not a second free source.
At zero wall speed no new fixed-radius worldline enters the cap.

This mixed-chart prescription has a conditional constraint implication.
When the temporal cap mass law holds as well as the radial constraints,
the cap normal acceleration computed from its metric and wall curve is

$$
K_{\tau,-}^{\rm actual}
=\frac{v_s+A_{-,r}/2}{\beta_-}
+\beta_-\delta_{-,r}+\alpha r v\Pi\Phi.
$$

Let $d=K_{\tau,-}^{\rm actual}-K_{\tau,-}^{\rm predicted}$, where
the prediction is the normal-force expression used above. The exterior
rapidity equation imposes its predicted acceleration. Assume the cap
constraints hold, including temporal mass consistency beyond radial
reconstruction alone. With the exterior outgoing
constraint zero by exact source transport, the earlier residual identities
have $\mathcal A_\Sigma=0$, $\mathcal T_\Sigma=-d$,
$\mathcal F_\Sigma=Sd$ and $J_-=0$. They give separately

$$
J_+=-vd,\qquad J_+=\beta_-d.
$$

Consequently $Hd=0$, hence $d=J_+=0$ on this branch. Both remaining
wall laws and the incoming exterior constraint follow without assuming
that incoming constraint first. This is an implication for a sufficiently
regular solution; it does not establish solvability or stability of a
discrete coupled iteration.

The zero-outgoing-constraint assumption matters for a source interpolant.
Retain the weighted null residuals
`J_+=Ubar_s^2 E_Ubar,+`, `O_+=x_s^2 E_x,+`,
`J_-=Q+_s^2 E_Q+,-` and `O_-=Q-_s^2 E_Q-,-`. Under the same smooth
mixed equations, common induced metric and scalar transmission, the
general residual identities are

$$
\mathcal A_{\Sigma,s}
=v\mathcal T_\Sigma-(J_++J_--O_+-O_-),
$$
$$
J_--J_++O_--O_+
=-\frac r2\mathcal F_\Sigma
+\left(\overline K_\tau+\frac{\overline\beta}{r}\right)
 \mathcal A_\Sigma+\overline\beta\mathcal T_\Sigma.
$$

If actual angular matching holds, exterior acceleration equals its
prescription, and `d=K_tau,-^actual-K_tau,-^predicted`, then
`T_Sigma=-d` and `F_Sigma=S d`. The identities give

$$
J_+=O_-+\frac{\beta_--v}{2}d,\qquad
J_-=O_+-\frac H2d.
$$

Since `beta_-^2-v^2=A_-` and `H=beta_-+v`, equivalently

$$
d=\frac{2(O_+-J_-)}H,\qquad
J_+=O_-+\frac{A_-}{H^2}(O_+-J_-).
$$

With vanishing cap constraints this reduces to `d=2O_+/H` and
`J_+=(A_-/H^2)O_+`. The actual Hermite source has nonzero `O_+`.
A nonzero incoming wall constraint alone therefore does not identify an
additional junction failure. Its source contribution and any cap
constraint contribution require evaluation with compatible actual
geometry; these formulas do not presume a discrete solution satisfies
their smooth hypotheses.

For an independent acceleration measurement, the full polar metric gives

$$
K_{\tau,-}^{\rm geom}
=\frac{v_s+A_{-,r}/2}{\beta_-}+\beta_-\delta_{-,r}
 -\frac{v A_{-,T}}{A_-^2 e^{\delta_-}}.
$$

The earlier compact expression substitutes the temporal mass law. If
`E_m=m_T-alpha*r^2*A_-^2*exp(delta_-)*Pi*Phi/2`, the difference between
the geometric expression and that substituted expression is
`2v E_m/(r*A_-^2*exp(delta_-))`. Radial reconstruction alone cannot
measure this temporal residual. An independent junction test needs the
actual acceleration and fixed-radius mass derivative, with their errors,
rather than substituting the desired mass law.

The [hybrid boundary helper](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
implements this map using the constrained cap metric and its outgoing
scalar characteristic. The regional update must use the same stage fields
for radial reconstruction and wall motion. On a moving cap grid it must
also include the coordinate advection: for $y=r/r_W(T)$,
$\partial_T|_y=\partial_T|_r+y(r_W)_T\partial_r$.
After the final exterior source label $V=V_b$, the central readout follows the
terminal inward ray $dr/dT=-B$ within the determined domain; no extra
incoming waveform is supplied to complete that ray.

### A simultaneous exterior, wall and cap stage

The [coupled operator](../src/recursive_horizons/nsc_atlas_coupled_evolution.py)
uses the full active exterior rather than identifying the completed
rectangle's terminal leg with the wall. Put

$$
0\le\eta=\bar U/W(x)\le1,\qquad
X=(r,L,c),\quad L=\log C,\quad Y=X-X_0(x).
$$

Here `X0` is the original source, and `W(.4)=Ubar_b>0`. The source
endpoint has `Y=0`; the other endpoint is the common wall state. Only
interior increments are evolved independently. For a simultaneous radial
slice, let `P=X_Ubar` and `Q=X_x`, with the latter derivative at fixed
`Ubar`. The exterior potential is zero. The exact bulk equations give

$$
(r r_x)_{\bar U}=-\frac C4,\qquad
(r f_x)_{\bar U}=-r_x f_{\bar U},\qquad
L_{\bar U x}=\frac{2r_{\bar U}r_x}{r^2}
+\frac C{2r^2}-\alpha f_{\bar U}f_x,
$$
$$
f_{\bar U}=\frac{c_{\bar U}-u r_{\bar U}}{r\sqrt D},\qquad
u=c/r,\quad D=1+\alpha u^2/6.
$$

Consequently three ordered radial prefixes determine `r_x`, `f_x` and
`L_x` from the original source derivatives, then
`c_x=u*r_x+r*sqrt(D)*f_x`. The
[prefix implementation](../src/recursive_horizons/nsc_atlas_exterior_characteristic_slice.py)
accepts spatial first jets explicitly. Its trapezoidal and Simpson
controls measure quadrature error separately from any later numerical
spatial derivative. No pure Einstein constraint is substituted into these
diagnostics.

The moving-coordinate update is

$$
Y_x\big|_\eta=Q-X_0'(x)+\eta\frac{W_x}{W}Y_\eta.
$$

For positive `W_x` its transport speed is negative, so the prototype uses
forward spatial upwinding towards the wall. The terminal backward jet is
compared with the incoming jet required by the hybrid law. It is not
replaced by that required jet when reporting compatibility.

The [source snapshot](../src/recursive_horizons/nsc_atlas_characteristic_trace.py)
retains the analytic compact `g,g_x,g_xx` and interpolates only the
integrated Jordan correction and transverse jets. Its `h,h_x,h_xx` are
derivatives of the same cubic Hermite polynomial. The source-ODE target
curvature is reported separately; inserting that target in place of the
interpolant's actual second derivative would conceal an input defect.
All live-source samples are obtained within that source's original budget;
later snapshot evaluations use the caller's cumulative allowance.

#### Fixed wall label and evolving lapse

Choose `W(x)=W_ref(x)` from the unchanged vacuum reference. This uses
the remaining outgoing-label freedom after the original source and quiet
bottom are anchored. It fixes a coordinate trajectory, not the physical
wall radius. The actual exterior wall lapse remains an evolved unknown.
Write `b` for the exterior null boost here, to distinguish it from the
radial grid coordinate:

$$
b=-\tfrac12\log W_x,\qquad
x_s=(C_W W_x)^{-1/2},\qquad \bar U_s=W_x x_s.
$$

The hybrid law supplies `r_s,c_s,T_s,k_plus`. Its acceleration equation
`b_s=k_plus-x_s*L_x+Gamma/2`, with `Gamma=d_s L_W`, then gives

$$
\frac{d r_W}{dx}=\frac{r_s}{x_s},\quad
\frac{d c_W}{dx}=\frac{c_s}{x_s},\quad
\frac{dT}{dx}=\frac{T_s}{x_s},\quad
\frac{ds}{dx}=\frac1{x_s},\qquad
\boxed{\frac{d L_W}{dx}=2b_x-\frac{2k_+}{x_s}+2L_x\big|_W.}
$$

The implementation first evaluates the boost-independent hybrid
quantities, then evaluates the incoming lapse jet using the actual
`Gamma=x_s*d_x L_W`. Prescribing `L_W` in addition to `W` would use the
same coordinate freedom twice. The measured endpoint mismatch is

$$
\mathcal R_W^X=\frac{dX_W}{dx}-Q_W-W_xP_W^{\rm numerical}
=W_x(P_W^{\rm required}-P_W^{\rm numerical}).
$$

These are two expressions for one residual, not independent error bounds.

The [reference gauge evaluator](../src/recursive_horizons/nsc_atlas_reference_wall_gauge.py)
inverts the reference's local advanced-time integral between radii `.45`
and `.55`, then differentiates its defining equations analytically.
It returns `W,W_x,W_xx` and `b,b_x` without differencing coarse wall
samples. `T=s=0` at the quiet anchor. The same evaluator supplies analytic
vacuum radius, lapse and clock comparisons, including the event where
the signed exterior `beta_plus` passes through zero.

#### Conservative cap update

The [cap operator](../src/recursive_horizons/nsc_atlas_polar_cap.py) uses
`y=r/R(x)` and the three regular cell variables `f,Pi,Q=Phi/r`.
Their averages carry weights `dr`, `r^2 dr` and `r dr`, respectively.
On a cell with fixed faces `a,b` in `y`, define
`d1=(b^2-a^2)/2` and `d2=(b^3-a^3)/3`. Reynolds transport gives

$$
\begin{aligned}
\Pi_x&=\frac{[y^2(T_xB\Phi+R_xy\Pi)]_a^b}{R d_2}
-\frac{3R_x}{R}\Pi-T_x\langle e^\delta U_\alpha\rangle_{r^2},\\
Q_x&=\frac{[T_xB\Pi+R_xy\Phi]_a^b}{R^2d_1}
-\frac{2R_x}{R}Q,\\
f_x&=T_x\langle B\Pi\rangle_{dr}
+\frac{R_x[yf]_a^b}{R(b-a)}-\frac{R_x}{R}f.
\end{aligned}
$$

At the axis the weighted `Pi` flux vanishes, but the `B*Pi` flux in
the `Q` equation remains. `Phi(0)=0` and the two characteristics reflect.
The scalar characteristics have speeds `-T_x B-R_x y` and
`T_x B-R_x y`; the timelike wall admits one incoming characteristic.
That value is the hybrid `w_plus`; the common transported `f_W` enters
only as the value of newly admitted material when the mesh expands.
It is derived from `c_W/r_W`, not an additional source waveform.

`Pi` is reconstructed linearly in `rho=r^2`, with limited slopes and
the exact `r^2`-weighted mean radius squared in each cell. This retains
the even quadratic axis term. A literal constant reconstruction of a
quadratic `Pi` would leave a nonvanishing first-cell error in `Q_x`.
`f` and `Q` are piecewise constant, so the overall prototype remains
first order. Smooth-axis manufactured controls explicitly test the
quadratic `Pi` case.

Every cap stage integrates the radial constraints using the same stage
fields, open Gauss nodes and a vacuum-subtracted mass,

$$
m=\frac{r^3}{2}+\alpha n,\qquad
n_r+\frac{\alpha r}{2}(\Pi^2+\Phi^2)n
=r^2\left[\frac{P-3}{2\alpha}
+\frac{1-r^2}{4}(\Pi^2+\Phi^2)\right].
$$

The continuous expression for `(P-3)/alpha` is used at zero or tiny
amplitude. The centre anchors `delta=-log cosh(sqrt(alpha/6)*f_C)`.
The incoming endpoint trace has zero radial measure and is not used as
an extra constraint input. Independently differenced mass/lapse residuals
and `Phi-f_r` are reported separately.

The complete common stage advances with fixed-step SSPRK2. Both its
predictor stages check the declared domain and characteristic Courant
limit; a first failure stops the calculation with its context. A CFL
diagnostic alone does not prove stability of the coupled metric, reaction
and boundary terms. The [admitted comparison](prior-art-and-open-claim.md#coupled-characteristic-exterior-and-polar-cap-16-october-2026)
therefore retains vacuum, spatial, temporal, source and constraint
measurements as separate evidence.

#### Independent constraints on a moving saved grid

The [history diagnostics](../src/recursive_horizons/nsc_atlas_coupled_diagnostics.py)
differentiate saved increments with centred stencils, independently of
the upwind stage or its mixed-equation integrals. Let
`b=W_x/W` and `a=eta*b` in this paragraph. Then

$$
\partial_x\big|_{\bar U}=\partial_x\big|_\eta-a\partial_\eta,
\qquad
r_{xx}\big|_{\bar U}=r_{0,xx}+Y_{r,xx}-2aY_{r,x\eta}
+a^2Y_{r,\eta\eta}
+\eta\left(2b^2-\frac{W_{xx}}W\right)Y_{r,\eta}.
$$

First derivatives use the analogous chain rule, with the source's actual
first jets, and `r_Ubar_Ubar=Y_r,etaeta/W^2`. This retains the steep
analytic source instead of differentiating its sampled waveform again.
The pure constraints and geometric mass are

$$
E_a=r_{aa}-L_a r_a+\frac{\alpha r}{2D}
\left(\frac{c_a-u r_a}{r}\right)^2,
\quad a\in\{\bar U,x\},\qquad
m=\frac r2\left(1+4e^{-L}r_{\bar U}r_x\right).
$$

All boundary stencils are excluded, and source-front stencils are
identified separately. Source constraint defects remain a separate input
diagnostic; lifting the known source does not remove its interpolation
error. Coarsening an already saved history measures diagnostic sampling
sensitivity. The separately evolved time-step comparison measures time
integration sensitivity. Neither replaces spatial refinement or a
continuum residual bound.

The [first five-case calculation](../results/development/nsc-atlas-coupled-evolution-v1/coupled.json)
completes this numerical window with the same source and action. Its
[measured responses and unresolved constraints](current-result.md#first-coupled-exterior-wall-and-cap-window-16-october-2026)
separate spatial, temporal and source sensitivity. In particular, the
early quiet wall layer retains a nondecreasing maximum `E_Ubar`, and
the cap's discrete field compatibility and scalar wall response are
unresolved. The completed trajectory is evidence for the implemented
coupling and its numerical behaviour, not the missing uniform estimate.

#### Saved-data diagnosis and the compatible cap connection

The flat maximum `E_Ubar` is present before the source pulse and moves
into a thinner wall-adjacent strip with refinement. Exact-vacuum
preservation at the discrete interface is therefore a direct control
for the next exterior scheme. Subtracting the vacuum diagnostic after
the calculation would not repair that boundary evolution.

The cap has a different local compatibility problem. Its separately
evolved `f` and `Q` have a nearly correct total jump but different local
profiles. On the 32-cell, 256-step nonzero-amplitude final slice, integrating the
piecewise `Phi=rQ` gives a wall-to-centre jump differing from the common
field jump by about `5.39e-8`. This diagnostic infers the regular centre
value as `f_C=f[0]-Q[0]*r1^2/6` from the first cell average and its outer
radius `r1`; the retained cap metric uses its original `f[0]` approximation.
The largest cell-average value
discrepancy from its primitive is about `.120849`. The latter is about
9.9% of the final peak `f`. The derivative diagnostic's largest defect
lies inside the pulse rather than at the axis. These comparisons use
the declared cell weights; a fine-cell restriction does not establish
decreasing maximum error yet.

A precise continuum connection is available without a new source law:

$$
f(T,r)=f_W(T)-\int_r^{R(T)}\Phi(T,z)\,dz.
$$

The existing wall chain rule and cap conservation equation give

$$
\frac{d f_W}{dT}=B_W\Pi_W+R_T\Phi_W,\qquad
\Phi_T=(B\Pi)_r
\quad\Longrightarrow\quad f_T=B\Pi,\qquad f_r=\Phi.
$$

Thus conservative `Q` plus the transported common wall value can
determine a compatible `f` reconstruction and its potential, instead
of carrying an independent locally inconsistent value profile.
For piecewise constant `Q`, the face increment is exactly
`Q*(r_right^2-r_left^2)/2`; the corresponding cell average is obtained
by integrating that same quadratic primitive. This is the next
numerical connection to implement and test, not an established
improvement to the retained run.

There is also an exact zero-amplitude transmission control. In the
de Sitter cap, `A=1-r^2`, `delta=0`, `r_star=atanh(r)` and `chi=r*f`
obeys the free reduced wave equation. Before an outgoing signal can
return from the centre, the quiet initial state implies

$$
(\partial_T-\partial_{r_\star})\chi=0
\quad\Longleftrightarrow\quad
R\,w_- - \frac{c_W}{R}=0.
$$

Every saved zero-amplitude wall event has `T<.198`, earlier than even
the initial inward light-travel interval `atanh(.45)>.484`. This
control therefore applies throughout the present window. Reconstructing
the declared outgoing last-cell trace from saved data gives peak
absolute defects about `8.72` and `9.84` on 16 and 32 cells. Resolving
this cancellation and local value/gradient compatibility is necessary
before interpreting the small finite-minus-zero scalar response.

#### Reduced canonical field and an exact incoming-wave control

An equivalent continuum system exposes that cancellation directly. Define
the reduced **canonical** field and its two first derivatives by

$$
\chi=rf,\qquad p=r\Pi,\qquad q=\chi_r=f+r\Phi.
$$

Here $q$ is not $Q=\Phi/r$. At finite amplitude $\chi$ also differs
from the earlier $c=ru$: with $a=\sqrt{\alpha/6}$,
$\chi=r\operatorname{arsinh}(ac/r)/a$, continued to $c$ at $a=0$.
The scalar equations give, without a new field assumption,

$$
\chi_T=Bp,\qquad
p_T=(Bq)_r-B_r f-r e^\delta U_\alpha(f),\qquad
q_T=(Bp)_r.
$$

The two radial metric constraints cancel their kinetic terms in
$B_r=e^\delta(A_r+A\delta_r)$, leaving the useful identity

$$
\boxed{B_r=e^\delta\left(\frac{2m}{r^2}-rP\right).}
$$

For $m=r^3/2+\alpha n$ and $P=3+\alpha P_1$, the combined source in
the $p$ equation can therefore be evaluated as

$$
B_r f+r e^\delta U_\alpha
=e^\delta\left[r(U_\alpha-2f)
+\alpha\left(\frac{2n}{r^2}-rP_1\right)f\right].
$$

Both bracketed terms vanish at zero amplitude. In that limit
$B=1-r^2$, $U_0=2f$ and the outgoing variable $d=p-q$ obeys

$$
\boxed{d_T+(Bd)_r=0,\qquad d_W=R w_- - f_W.}
$$

Thus a compatible reconstruction of $f$ and preservation of zero outgoing
data are distinct controls. The latter requires the numerical evolution
and its boundary flux to preserve this homogeneous equation; enforcing a
primitive identity alone does not establish it.

Regular-axis data have $\chi$ and $p$ odd, with $\chi(0)=p(0)=0$,
and $q$ even, with $q(0)=f_C$. Consequently the reduced characteristics
reflect with opposite signs: $p-q=-(p+q)$ at the axis. For a moving wall,

$$
\frac{d\chi_W}{dT}=B_Wp_W+R_Tq_W,\qquad
p_W+q_W=R w_++f_W.
$$

Integrating $q_T=(Bp)_r$ shows that
$\chi_W-\int_0^R q\,dr$ is constant when the regular axis and this
same wall law hold. A discrete analogue must use consistent wall fluxes
and time updates. It must also test propagation against an independent
solution, rather than count an imposed reconstruction identity as
evolution accuracy.

Such a solution is available before central arrival. For a compact
profile $F(\zeta)$, with $\zeta=T+\operatorname{atanh}r$ and support
above the tested times, set $\chi=F(\zeta)$. Then

$$
f=\frac F r,\quad \Pi=\frac{F'}{(1-r^2)r},\quad
\Phi=\frac{F'}{(1-r^2)r}-\frac F{r^2},\quad
p=q=\frac{F'}{1-r^2},\quad
w_+=\frac{2F'}{(1-r^2)r}-\frac F{r^2}.
$$

The central neighbourhood is exactly quiet, so these expressions need
no division at the axis. For cell faces $a<b$, write
$F_a=F(T+\operatorname{atanh}a)$ and likewise for $F_b$. Exact ordinary
radial means of the reduced variables are

$$
\langle p\rangle_{dr}=\langle q\rangle_{dr}
=\frac{F_b-F_a}{b-a}.
$$

The original weighted variables have independent targets

$$
\langle Q\rangle_{r\,dr}
=\frac{f(b)-f(a)}{(b^2-a^2)/2},\qquad
\langle\Pi\rangle_{r^2dr}
=\frac{bF_b-aF_a-\int_a^b F(T+\operatorname{atanh}r)\,dr}
{(b^3-a^3)/3}.
$$

For example, $R(T)=.45+.5T$ on $0\le T\le.2$ is strictly timelike,
and the compact $C^5$ profile $F=.01(1-z^2)^6$ for $|z|<1$, zero
otherwise, with $z=(\zeta-.65)/.10$, is initially quiet and cannot
reach the centre during this test. Its known wall value, incoming
characteristic, cell means and time-dependent profile separately test
transmission, compatibility and propagation. The successor's manufactured
controls now include this independent propagation comparison. Its prepared
physical-source comparison remains pending; these identities do not alter
the retained five-case result or resolve its measured errors.

#### Conservative reduced variables and exact vacuum departures

The [reduced-cap operator](../src/recursive_horizons/nsc_atlas_reduced_cap.py)
stores the physical cell integrals `Ip=integral(p dr)` and
`Iq=integral(q dr)`. On the moving grid `r=R(x)y`, their balances are

$$
\partial_x I_p=[T_xBq+R_xy p]_{\rm faces}
-T_x\int\left(B_r f+r e^\delta U_\alpha\right)dr,
\qquad
\partial_x I_q=[T_xBp+R_xy q]_{\rm faces}.
$$

There are no volume-dilution terms for these actual integrals. The regular
axis has zero `Iq` flux. Consequently `chi_W=sum(Iq)` evolves by the same
wall flux used for transmission. The field is the primitive of the
reconstructed `q`, starting from `chi(0)=0`. Away from the first cell,
both reduced characteristics use the same limited linear reconstruction.
The first cell uses odd `p=r(P0+P2*r^2)` and even `q=Q0+Q2*r^2`, with
coefficients preserving its integrals and regular limits. The radial
constraints evaluate this field, rather than a separate cell value, at
their quadrature nodes. The exact within-cell lapse integral and bounded
Gauss mass integration retain the zero-amplitude metric analytically.

The [coupled successor](../src/recursive_horizons/nsc_atlas_compatible_evolution.py)
therefore has no independent scalar wall state. It converts
`chi_W` to `c_W=R*sinh(a*chi_W/R)/a`, with its continuous `a=0` limit,
before applying the unchanged hybrid map. The resulting incoming
characteristic determines the cap face flux. The converted conservative
`c_x` is checked independently against the hybrid `c_s/x_s`. Its
numerical one-sided exterior jet remains a separately reported residual.

For the exterior, let `V=(r_v,L_v,0)` be the exact Schwarzschild fields in
the same `Ubar,x` coordinates. With `a_h=.8`, they satisfy

$$
r_v=a_h\left[1+W_0\left(-\frac{(1-\bar U/(4a_h))
(1-x/a_h)e^{x/a_h}}{e}\right)\right],
\qquad
C_v=\frac{x}{r_v}e^{(x-r_v)/a_h},\qquad L_v=\log C_v.
$$

Here `W0` is the principal Lambert function; this is unrelated to the
moving wall label `W(x)`. Write

$$
X=V+X_0(x)-V(0,x)+Z,\qquad Z(0,x)=0.
$$

The [departure-prefix operator](../src/recursive_horizons/nsc_atlas_vacuum_departure.py)
differentiates `V` analytically and only differences `Z`. Its radial
prefix is the same nonlinear equation rearranged as

$$
r_x=\frac{r_0r_{0x}-[x-r_v r_{vx}]
-\int_0^{\bar U}(C-C_v)/4\,d\bar U'}{r}.
$$

Likewise, with `F=2*r_Ubar*r_x/r^2+C/(2*r^2)-alpha*f_Ubar*f_x`,
it evaluates `L_x=L_vx+L_0x+integral(F-F_v)`. The canonical scalar
prefix remains the full equation. Stable difference formulas avoid
cancellation between almost equal background terms. The source's values
and first derivatives remain the original supplied data. Evolve

$$
Z_x\big|_\eta=Q-V_x-[X_0'-V_x(0,x)]
+\eta\frac{W_x}{W}Z_\eta.
$$

Wall radius, lapse and clock variables also store departures from their
reference trajectories. Fixed-step SSPRK2 therefore preserves the moving
vacuum to roundoff rather than repeatedly integrating its discretization
error. The removed background quadrature errors are retained separately.
Independent pure constraints still require derivatives of saved fields;
report both ordinary finite differences and analytic vacuum jets plus
finite differences of the departures. Neither diagnostic substitutes the
mixed evolution equations for the derivatives it is testing.

#### Compatible window outcome and the remaining cap transport

The [five-case record](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
uses the same preparation, amplitude and interval as the first coupled
window. Its 16/32-cell, 256-step cases and separate finite 32-cell,
512-step case all finish. Creation and authentication together cost
64.40625 process CPU seconds and 5,371,737 charged points. The
[current result](current-result.md#compatible-coupled-evolution-16-october-2026)
contains the response comparisons.

Vacuum geometry and wall clocks now agree with their exact trajectories
to approximately `1e-16`. The cap's zero-amplitude outgoing reduced trace
`a_minus=p-q` is exactly zero at every saved wall event. These remove the
predecessor's quiet geometric layer and its spurious outgoing cap signal.
The constructed primitive and wall-flux identities also hold to roundoff;
they do not measure propagation error.

The actual finite response remains nonzero. At 32 cells its endpoint
radius change is `-8.56927e-4`, central-clock-label change `-1.09815e-3`
and conformal wall-scalar change `+3.29131e-5`. The paired scalar spatial
change is `2.49582e-6`, about 6.9% of its peak response. The separate raw
finite scalar time comparison is `1.24653e-4`, with only `2.00485e-8`
at the endpoint. Without a corresponding zero-amplitude 512-step run,
this does not establish time convergence of the small nonlinear history.
Finite departure constraints still peak during the pulse at the final
interior spatial node: `abs(E_Ubar)=.00959384` and `abs(E_x)=.00748492`.
Their spatial changes improve; their values are not error enclosures.

There is a direct propagation control before any central return. At zero
amplitude the same equations give

$$
\chi(T,r)=F(T+\operatorname{atanh}r),\qquad p=q.
$$

The saved wall history supplies `F` on
`zeta_W=T+atanh(R)`, since `chi_W=c_W` in this limit. Quiet initial data
supply zero before the wall's first label. Pull this history back along
incoming characteristics and compare it with the independently evolved
cell-face primitive `cumsum(Iq)`. This comparison tests transport of the
actual computed wall input; it does not assume that input is exact.
Linear interpolation gives maximum primitive differences `.56115` and
`.38095` on 16 and 32 cells, at the final `y=.75` face. A shape-preserving
cubic interpolation changes the differences by about `.001`, much less
than the discrepancy. The incoming pulse has amplitude about `.995`.
Thus exact cancellation of the outgoing characteristic coexists with
substantial spreading of the incoming pulse.

The finite outgoing maximum changes from `1.60392e-4` to `3.44027e-4`,
while its separate time change is `4.00569e-7`. The fine cap's energy
source is consequently not established by the integral identities or by
the small wall time-step change. Its final mass and lapse endpoint-trap
defects are `.000498315` and `.00244719`; these compare one radial
quadrature with another and do not bound the Gauss constraint error.
The next numerical connection is accurate incoming characteristic
transport and its gradient-squared gravitational source. It must retain
the single common wall scalar and the original action. Extending the
domain or increasing every grid without addressing this identified
leading-wave error would not resolve the present uncertainty.

#### Carrying a characteristic background and evolving its correction

The [background adapter](../src/recursive_horizons/nsc_atlas_characteristic_background.py)
defines one finite C2 quintic Hermite profile `F` from its values and first
two derivatives. All evaluated jets differentiate that same interpolant.
Quiet continuation below its first label requires matching zero value and
both derivatives; no upper extrapolation is supplied. With
`xi=atanh(r)` and `A0=1-r^2`, define

$$
\chi_0=F(T+\xi)-F(T-\xi),\qquad
p_0=\frac{F'(T+\xi)-F'(T-\xi)}{A_0},\qquad
q_0=\frac{F'(T+\xi)+F'(T-\xi)}{A_0}.
$$

Both cell integrals have exact endpoint primitives:
`Iq0=[F(T+xi)-F(T-xi)]` and `Ip0=[F(T+xi)+F(T-xi)]`.
Reflection retains `chi0(0)=p0(0)=0`, with `f0(0)=2F'(T)` and
`Pi0(0)=2F''(T)`. Piecewise C2 interpolation supplies these finite axis
limits, not arbitrary smoothness: a jump in the third derivative can
leave a radial derivative cusp in `Pi0` at a knot-centre event. The present
physical window has a quiet centre; later centre-regularity claims must
state their required profile smoothness.

The [departure cap](../src/recursive_horizons/nsc_atlas_background_cap.py)
reuses the conservative reconstruction for `dp=p-p0` and `dq=q-q0`.
It evolves their physical cell integrals. Subtracting the exact background
balances gives

$$
\begin{aligned}
\partial_x\Delta I_p&=
 [T_x(B\,dq+(B-A_0)q_0)+R_xy\,dp]_{\rm faces}
 -T_x\int S\,dr,\\
\partial_x\Delta I_q&=
 [T_x(B\,dp+(B-A_0)p_0)+R_xy\,dq]_{\rm faces},\\
S&=\alpha e^\delta r\left[
 (2n/r^3-P_{\rm diff})f+(U_\alpha-2f)/\alpha\right],\qquad
P_{\rm diff}=(P-3)/\alpha.
\end{aligned}
$$

The last quotients use their continuous limits. The total scalar is
`chi=chi0+integral(dq dr)` everywhere, including the wall. The hybrid map
continues to supply the total incoming characteristic. Its departure is
`R*w_plus+f_wall-(p0+q0)_wall`, which generally is nonzero. Consequently
choosing `F` does not impose a second boundary datum or suppress a finite
waveform correction. The total wall rate restores
`T_x*A0*p0+R_x*q0` before converting the canonical field to `c`.

The metric uses total `f,Pi,Phi`, including their background/departure
cross terms. Split the radial integration at every cell edge and mapped
profile knot. On each panel integrate

$$
\gamma'=a,\quad n'+an=b,\qquad
a=\alpha r(\Pi^2+\Phi^2)/2,\quad
b=r^2[P_{\rm diff}/2+(1-r^2)(\Pi^2+\Phi^2)/4].
$$

Let `G_poly` be the Gauss polynomial antiderivative approximation to
`integral(a)` from the panel's left edge. At each target node and endpoint,
use a separate positive-weight Gauss partial quadrature `Q` of the actual
total fields:

$$
G_t=Q_{[r_l,r_t]}(a),\qquad
n_t=e^{-G_t}n_l+Q_{[r_l,r_t]}(e^{G_{\rm poly}-G_t}b),\qquad
\gamma_t=\gamma_l+G_t.
$$

The shifted exponential avoids an unnecessary large intermediate factor.
The factors restart at each panel without changing its prefix values.
Reconstruct `m=r^3/2+alpha*n`,
`delta=-log(cosh(sqrt(alpha/6)*f_center))+gamma` and `A=1-2m/r`.

The first implementation integrated the mass source through an
antiderivative matrix. A narrow-pulse control exposed a negative nodal
mass contribution despite positive energy. The positive-weight partial
quadrature above removes that defect without clipping: `a,b,n_l>=0`
give `G_t,n_t>=0`. Separate partial quadratures do not by themselves
guarantee monotonicity between nodes. The difference between `G_poly`
and direct `G_t` remains a reported quadrature-consistency diagnostic.
Independent radial integration and order-4/8/16 narrow-pulse comparisons
test accuracy; constructed node and primitive identities cannot do so.

The [input adapter](../src/recursive_horizons/nsc_atlas_background_profile.py)
maps supplied zero-case wall values and same-stage rates to the profile:

$$
\zeta_W=T+\operatorname{atanh}R,\qquad
F(\zeta_W)=c_W,\qquad
F'(\zeta_W)=\frac{c_{W,x}}{T_x+R_x/(1-R^2)}.
$$

It uses a stated second-order difference of these first jets to select
`F''` as computational interpolation curvature. It performs no physical
source solve and does not advertise those second jets as an independent
physical calculation. The finite profile range remains enforced at every
runtime query. Withheld wall values and first jets can check the chosen
interpolant separately from the evolved incoming correction.

The [coupled adapter](../src/recursive_horizons/nsc_atlas_background_evolution.py)
retains the original exterior and wall formulas and fixed SSPRK2. Its
cap state stores departure integrals, and every saved physical view also
stores the total integrals at that event's actual `T,R`. Its wall chain
uses the full background-plus-departure rate. Seventy-two controls pass across
these four operators, including actual moving-vacuum steps, an active
finite-amplitude background with nonzero incoming correction and the
regular-axis and radial-integration checks above. Exact terminal-label
evaluation passes, while a one-ulp excursion outside the profile is rejected.
The [first physical comparison](current-result.md#characteristic-background-coupled-evolution-16-october-2026)
now completes all six cases in the [calculation record](prior-art-and-open-claim.md#coupled-characteristic-background-comparison-16-october-2026).
Creation and authentication use 228.9375 process CPU seconds and
43,551,872 charged points. The independent zero-cap pullback reduces the
32-cell primitive discrepancy from about `0.38` to below `0.0009`;
remaining linear/PCHIP interpolation sensitivity is `0.00133`. This
demonstrates removal of the dominant incoming-pulse smearing in this window.
Geometry/clock paired spatial changes are below `0.063%` of response peaks,
but the scalar paired time change is `15.83%` of its fine peak. The three
frozen radial order-4/8 controls change `n_wall` by at most `3.00472e-8`
and `delta_wall` by `2.28539e-11`. Pointwise exterior constraints and
wall-slope matching still need resolution. The saved comparisons do not
provide an error enclosure or a new physical preparation mechanism.

#### Isolating the scalar time contribution

For the canonical wall primitive, write `chi=B(T,R)+Z`, where `B` is
the analytic background value and `Z` is the sum of departure `Iq`.
The cap flux gives `Z_x=J_chi-B_x`, with
`B_x=T_x*(1-R^2)*p0+R_x*q0`. Consequently, an SSPRK2 update gives

$$
\chi_{n+1}-\chi_n=\frac h2(J_{\chi,n}+J_{\chi,*})+\mathcal R_B,
\qquad
\mathcal R_B=B_{n+1}-B_n-\frac h2(B_{x,n}+B_{x,*}).
$$

`B_{n+1}` is evaluated at the accepted geometry, not an average of the
two stage background values. The conformal junction map separately gives
`J_c=Q_c+C_minus/(2*x_s)`. Thus the known-profile chain contribution can
be separated from the exterior and returning-cap rate components.
Neither identity identifies the observed error by itself.

The [saved-data diagnostic](../scripts/diagnose_nsc_atlas_background_time.py)
and [35 controls](../tests/test_diagnose_nsc_atlas_background_time.py)
implement these identities, signed polynomial quadrature oracles and
strict profile-knot partitioning. Controls use 10,956 points and 0.265625
process CPU seconds, plus a three-point intentional exhaustion control;
an independent 18-point cubic check confirms the sign and factor-four
two-half-step reduction. The [bounded audit](prior-art-and-open-claim.md#saved-profile-time-integration-diagnosis-16-october-2026)
uses the existing reference cache and four frozen stages. No new physical
evolution is required for that diagnosis.

The completed audit uses 57,152 points and 1.0625 process CPU seconds.
The reference-chain time difference peaks at `1.24476e-4`, versus the
actual zero-case difference `2.87645e-6`. At cut 82 the paired rate change
splits into exterior `4.64186e-6` and cap return `-6.65783e-9`, while
the independent identities close within `2.84e-14`. This demonstrates
cancellation but does not attribute the full-history error. The
[result](current-result.md#scalar-time-contribution-audit-16-october-2026)
retains the distinction between rates, local path defects and accumulated
field differences.

For the next actual-step comparison, let angle brackets denote the
two-stage trapezoidal increment and define

$$
\chi_{\rm RK}=\chi_n+\langle J_\chi\rangle,\qquad
R_{\rm RK}=R_n+\langle R_x\rangle,\qquad
C_\alpha(R,\chi)=\frac{R}{k}\sinh(k\chi/R),\quad
k=\sqrt{\alpha/6},\quad C_0=\chi.
$$

The accepted primitive is `chi_1=chi_RK+R_B`. The physical accepted
radius `R_1` includes the exact reference-gauge value at the endpoint;
since the code evolves departures, it need not equal `R_RK`. The
following five terms give the exact invariant scalar increment:

$$
\begin{aligned}
c_1-c_n={}&\langle Q_c\rangle+
\left\langle\frac{C_-}{2x_s}\right\rangle+
\mathcal R_{\rm frame}+\mathcal R_{\rm geometry}+\mathcal R_{\rm background},\\
\mathcal R_{\rm frame}={}&C_\alpha(R_{\rm RK},\chi_{\rm RK})-c_n-\langle J_c\rangle,\\
\mathcal R_{\rm geometry}={}&C_\alpha(R_1,\chi_{\rm RK})-C_\alpha(R_{\rm RK},\chi_{\rm RK}),\\
\mathcal R_{\rm background}={}&C_\alpha(R_1,\chi_{\rm RK}+\mathcal R_B)-C_\alpha(R_1,\chi_{\rm RK}).
\end{aligned}
$$

Both frame and geometry terms vanish at zero amplitude. Clock and radius
reference-increment defects are retained separately as diagnostics; their
effect on `B_1` already belongs to the background term. Each half-step
has its own accepted geometry and ledger. Adding the two half-step
ledgers before subtraction makes the local comparison exact without
asserting a propagated error estimate. The [admission](prior-art-and-open-claim.md#actual-scalar-step-balance-16-october-2026)
fixes cut 82, both amplitudes, twelve stages and the existing input cache.

The [completed local ledger](current-result.md#actual-scalar-step-balance-16-october-2026)
uses 138,891 points and 1.546875 CPU seconds. The paired fine-minus-coarse
increment is `8.03607226e-8`: background chain `8.32362063e-8`, exterior
flux `-2.65646853e-9`, frame conversion `-2.21866157e-10`, cap return
`2.85128618e-12` and radius conversion at roundoff. Its closure is
`5.10e-17`. The individual amplitudes contain much larger cancelling
terms. This is a measured local truncation sensitivity, not a reason to
drop the exact background contribution from the update.

#### Higher-order time integration with the same regional equations

The next adapter changes the time update alone. With `Y` denoting the
raw exterior, wall and cap departures, use SSPRK(3,3):

$$
\begin{aligned}
Y_1&=Y_n+hF(x_n,Y_n),\\
Y_2&=\tfrac34Y_n+\tfrac14[Y_1+hF(x_n+h,Y_1)],\\
Y_{n+1}&=\tfrac13Y_n+\tfrac23[Y_2+hF(x_n+h/2,Y_2)].
\end{aligned}
$$

The [standard method](https://math.umd.edu/~tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf)
has weights `b=(1/6,1/6,2/3)` at `c=(0,1,1/2)`. Its order conditions
are `sum(b)=1`, `sum(b*c)=1/2`, `sum(b*c^2)=1/3` and `b*A*c=1/6`.
The time coordinate can be included as an additional equation `x'=1`,
which fixes the nonautonomous stage coordinates. Source and reference
geometry are read at these coordinates; physical predictor geometry
is not formed by mixing reference values from different times.

In the five-part scalar ledger above, replace every trapezoidal increment
by `h*(J_0/6+J_1/6+2*J_2/3)`. The accepted canonical primitive still
satisfies `chi_1=chi_RK+R_B`, and the same telescoping frame/geometry/
background decomposition remains exact. The cap wall scalar is neither
prescribed nor independently evolved. Explicit midpoint nodes retain
exact access to the original cache despite floating-point rounding.

The [bounded local comparison](prior-art-and-open-claim.md#ssprk3-with-fixed-physical-inputs-16-october-2026)
holds every spatial and physical input fixed. Computational profile knots
and limiter changes can reduce the measured temporal order. Any measured
improvement leaves the spatial constraint and cap-refinement questions
separate; the Runge--Kutta coefficients do not resolve those questions.

The [first local result](current-result.md#local-ssprk3-comparison-16-october-2026)
uses six stages and 70,118 points. Its paired increment differs from the
SSPRK2 asymptotic extrapolant by `-2.71979e-9`, compared with
`-2.67869e-8` for the two SSPRK2 half steps. The weighted canonical and
invariant balances close at roundoff. This motivates the fixed-input
128/256-step comparison; the extrapolant is not an exact solution.

The [full comparison](current-result.md#fixed-input-ssprk3-coupled-window-16-october-2026)
now completes at the same two amplitudes and fixed 32-cell spatial grid.
Its paired radius and clock temporal differences are below `1.34e-9`,
while the paired scalar difference is `7.57345e-6` (21.54% of the fine
peak). The new fine scalar history differs from the earlier pointwise
SSPRK2 extrapolant by at most `5.60381e-8`. Agreement between those
approximations does not certify their common limit or uniform error.
The finite constraint residuals persist. A small pre-support scalar is
also present in both amplitudes and decreases with time resolution.
The actual source/profile regularity and stage sampling determine which
time comparisons are useful; a formal method order is insufficient.

#### Computational support and profile-knot sampling

The retained computational profile supplies a concrete pre-support defect.
Its knot 17 lies at reference `x=0.4216495553`, before the physical source
front `x=0.4226917906`, and has `F=F'=0` but `F''=0.1045440553`.
Knot 16 has all three zero jets. If their label spacing is `d` and
`t=(zeta-zeta_16)/d`, the unique quintic on that interval is

$$
F(\zeta)=\frac{F''(\zeta_{17})d^2}{2}\,t^3(1-t)^2.
$$

Its nonzero interior follows directly from the six endpoint constraints.
The estimated curvature borrows information from the later nonzero knot,
although the physical source `c,c_x` are exactly zero before its front.
At saved cut 35 the wall background is `7.63534e-6`; the evolved departure
nearly cancels it, leaving `-3.17616e-7` at 128 steps and `-1.98830e-8`
at 256 steps. The fine cap's two adjacent incoming integrals are about
`-1.57101e-7` and `1.39906e-7`. Pre-support outgoing `Ip-Iq`, outgoing
wall characteristic and central field are zero in all four histories.
This is an auxiliary-profile cancellation error, not a physical precursor.
Respecting the known support in a future profile construction is a separate
correction; it must not be mixed into the present fixed-profile time test.

The late paired discrepancy has a different location. Define from the
saved physical wall values

$$
e_i=(c_{\alpha,256,i}-c_{0,256,i})
 -(c_{\alpha,128,i}-c_{0,128,i}),\qquad
q_i=T_{\alpha,128,i}+\operatorname{atanh}r_{\alpha,128,i}.
$$

Let `I_i=1` when a retained profile label lies strictly between
`q_(i-1)` and `q_i`. Exactly 45 intervals satisfy this endpoint-bracket
test, and

$$
\frac{\sum_i I_i(e_i-e_{i-1})^2}{\sum_i(e_i-e_{i-1})^2}
=0.999885990677127.
$$

The 5% and 95% cumulative quantiles of `sum(e_i^2)` are cuts 89 and 98.
These are reproducible saved-array associations, not recovered internal
stage histories. Limiter branch disagreements peak elsewhere.

There is an exact quadrature reason to retain the crossing phase. For a
local cubic kink `B(x)=A*(x-x_*)_+^3/6`, write `theta=(x_*-x_0)/h`.
For `0<=theta<=1`, direct integration gives the SSPRK3/Simpson chain defect

$$
\begin{aligned}
R_B&=B(x_0+h)-B(x_0)
 -\frac h6\{B'(x_0)+4B'(x_0+h/2)+B'(x_0+h)\}\\
&=\frac{Ah^3}{12}(1-2\theta)
   \min(\theta,1-\theta)^2.
\end{aligned}
$$

For `theta<=1/2`, expansion gives `theta^2*(1-2*theta)`; for
`theta>=1/2`, it gives `(1-theta)^2*(1-2*theta)`. A cubic kink aligned
with an endpoint or the midpoint has zero defect. An interior kink has
a phase-dependent cubic defect, whose sign can change on refinement.
Its largest magnitude is `abs(A)*h^3/324`, at phases `1/3` and `2/3`.
For a monotone increasing smooth wall label `zeta(x)`, a jump of `F'''` contributes
this local leading term with `A=[F''']*zeta_x^3`; the remaining terms
depend on the smooth profile and trajectory. The full coupled stages
need not follow that prescribed trajectory, so this identity alone does
not assign the entire scalar error.

On the zero reference the knots align with step boundaries. Evaluating
`B=F(T+atanh r)-F(T-atanh r)` and its analytic first jet on the original
513 gauge rows gives cumulative Simpson defects `2.06210e-6` and
`1.28889e-7` at 128/256 steps. Their difference follows the actual zero
wall time difference with correlation `0.997555` and maximum residual
`2.23903e-7`. The saved-data analysis uses 727,310 charged arithmetic
points and 0.921875 CPU seconds, with no coupled stage or source solve.
The finite trajectory shifts knots into the step interiors. The admitted
[512-step comparison](prior-art-and-open-claim.md#midpoint-complete-ssprk3-comparison-16-october-2026)
now measures a paired scalar time difference `5.84100e-7`, or 1.67616%
of the fine response peak, while preserving that profile. Its difference
from the old SSPRK2 extrapolant increases to `5.92862e-7`, so the earlier
close cross-method agreement did not certify the limit. Finite constraint
and quiet cap residuals persist; no further time refinement is automatic.

A support-preserving profile would give a conditional scalar-zero
invariance: if all reached cap labels satisfy `T+atanh(r)<=zeta_front`,
the profile and its first two jets vanish there, all incoming source scalar
jets vanish, and the cap/exterior scalar states are zero, then the bulk,
cap and hybrid-wall scalar rates vanish. SSPRK3 preserves that zero state.
This does not prove exact finite vacuum geometry. The retained source
snapshot has prefront Hermite leakage, with `max(abs(h))=1.26817e-14`,
`max(abs(h_x))=5.08617e-11` and finite `max(abs(r_x-1))=6.66134e-15`,
even though `c,c_x,c_xx` are exactly zero. The scalar check uses 26,178
saved values and 0.140625 CPU seconds, with no source or stage evaluation.

Existing prefront profile values and first jets can be retained, but
their nonzero computational curvature must change if the profile is to
vanish there. A zero value/first/second jet at the exact reference front
permits a C2 connection to the later retained knots. Matching a C3 onset
requires an additional condition. The source-front label and actual-stage
coverage need explicit checking before that successor is constructed.

#### Inherited constraints and saved-cut differentiation

The exact transport identity above also applies to the numerical source
representation's nonzero constraint. In the separate null coordinates
`(Ubar,x)`, set `u=c/r`, `D=1+alpha*u^2/6` and

$$
E_a=r_{aa}-L_a r_a+\frac{\alpha r}{2D}u_a^2,
\qquad a\in\{\bar U,x\}.
$$

On a piecewise smooth solution of the mixed equations, a segment from
the source at `Ubar=0` carries

$$
E_x^{\rm inherited}(\bar U,x)
=\frac{r_{\rm src}(x)E_{\rm src}(x)}{r(\bar U,x)},
\qquad
\mathcal T_x=rE_x^{\rm FD}-r_{\rm src}E_{\rm src}.
$$

The source term uses its actual interpolated `r_xx`; substituting the
exact source ODE would erase a supplied residual. The
[saved-data calculation](prior-art-and-open-claim.md#saved-constraint-transport-and-sampling-16-october-2026)
retains the raw diagnostic, inherited prediction and subtracted remainder
as separate arrays. Radius is reconstructed from the saved physical
increment and source value, without evaluating a coupled stage.

Every one of the 129 saved cuts is a knot of the source Hermite
interpolant. At the interior cuts its second derivative uses the right polynomial trace,
whereas centred differences of propagated departures sample both sides.
The 65-cut diagnostic has wider temporal stencils at a subset of those
same knots. Comparing the two diagnostics at identical points measures
differentiation sensitivity, including this trace convention. It does
not by itself isolate evolution error or establish its order. Front-near
comparisons use the union of both stencils' front masks.

The trace effect has a conditional local formula. Suppose the fields and
first jets are continuous across a source knot, their pieces are C2, and
the mixed and tangential second jets have common traces. If
`J=[r_src,xx]`, exact transport and the departure decomposition imply

$$
[r_{xx}]_{\bar U}=\frac{r_{\rm src}}r J,
\qquad [Z_{xx}]_{\bar U}=\left(\frac{r_{\rm src}}r-1\right)J.
$$

A centred second difference tends to the mean of the two curvature
traces. With other derivatives evaluated exactly, adding the right source
trace therefore leaves the limiting contribution
`(1-r_src/r)*J/2` in `E_x^FD-E_x^inherited`. Smooth `W` and the common
mixed/tangential traces give the same jump in the moving eta description.
C1 regularity alone does not establish those second-jet assumptions.
This identity identifies a possible sampling contribution; no correction
is subtracted from the saved diagnostic on its authority alone.

For the other constraint, `Ubar=eta W(x)` determines the origin of its
constant-null-coordinate segment. Points with `Ubar<=W(.4)` inherit the
quiet bottom value. Points with `Ubar>W(.4)` enter through the wall and
instead satisfy

$$
rE_{\bar U}(\bar U,x)
=r_W(x_b)E_{\bar U,W}(x_b),\qquad W(x_b)=\bar U.
$$

The incoming wall value remains to be determined independently. A saved
first-jet tangential defect does not determine `r_Ubar_Ubar`; setting
that wall value to zero would assume the needed result. Partition the
three-point eta stencils into wholly bottom-origin, wholly wall-origin
and separator-crossing sets. The 129/65 comparison uses the same eta
stencil, so equality of `E_Ubar` at matching cuts is an indexing control,
not spatial convergence. General wall constraint identities must also
retain the outgoing source residual rather than assume it vanishes.

#### Independent incoming wall curvature

The [saved audit](current-result.md#saved-constraint-transport-and-sampling-16-october-2026)
finds an incoming constraint maximum `5.54765e-4` on wholly wall-origin
stencils, compared with `1.15980e-5` on wholly bottom-origin stencils.
An independent wall datum can be reconstructed without the hybrid wall
rates. At the wall, recover each departure from
`Z_W=X_W-X_vac(W,x)-source_offset`. The source offset is constant in
`Ubar`; the analytic vacuum supplies its own first and second jets.

Let `N=32`, `h=s/32` for stride `s=1` or `2`, and order the wall samples
as `Z_N,Z_(N-s),...`. The baseline one-sided differences are

$$
D_1 Z=\frac{3Z_N-4Z_{N-s}+Z_{N-2s}}{2h},\qquad
D_2 Z=\frac{2Z_N-5Z_{N-s}+4Z_{N-2s}-Z_{N-3s}}{h^2}.
$$

The five-point alternatives have coefficient rows
`(25,-48,36,-16,3)/(12h)` for `D_1` and
`(35,-104,114,-56,11)/(12h^2)` for `D_2`. Add `D_1 Z_X/W` to the vacuum
first jets and `D_2 Z_r/W^2` to its pure transverse radius curvature.
The actual saved `r_W,c_W` then give

$$
E_{\bar U,W}
=r_{\bar U\bar U,W}-L_{\bar U,W}r_{\bar U,W}
 +\frac{\alpha r_W}{2D_W}
 \left(\frac{c_{\bar U,W}-(c_W/r_W)r_{\bar U,W}}{r_W}\right)^2.
$$

This uses measured spatial curvature. It neither inserts a prescribed
extrinsic curvature nor infers it from first-jet tangential compatibility.
The formal stencil orders require sufficient smoothness and are not
assumed error estimates for these fields.

The monotone saved labels `W(x_k)` directly tabulate
`B(W(x_k))=r_W(x_k) E_Ubar,W(x_k)`. Thus no numerical inverse-gauge
solve is needed to compare `r E_Ubar` with the boundary value at an
interior `Ubar`. Piecewise-linear interpolation on 129 wall samples,
compared with the 65 alternating samples, tests interpolation sampling
separately from the four spatial derivative stencils. Both bracket ends
must precede or equal the observation cut. Unsupported points remain
masked; no future wall value or extrapolation supplies them.

The widest boundary stencil reaches eta index 24. A disjoint interior
three-point stencil therefore has centre index at most 22. Report these
points separately from overlapping stencils: shared values can correlate
derivative errors and make apparent transport agreement less informative.
Even disjoint stencils share one numerical evolution. The
[fixed diagnostic](prior-art-and-open-claim.md#incoming-wall-curvature-and-transported-constraint-16-october-2026)
tests the missing boundary connection; it does not prove continuous
constraint propagation or the complete junction system.

The [completed calculation](current-result.md#incoming-wall-curvature-and-transported-constraint-16-october-2026)
supports all 916 wall-origin points at both interpolation samplings, but
does not resolve this connection. On the 162 disjoint-stencil points,
every wall subtraction increases the raw weighted-constraint RMS
(`7.16444e-5` becomes `8.14450e-5`--`1.38386e-4`). Changing the wall
derivative stencil changes its transported value by up to
`1.08923e-4`--`2.67968e-4`, far more than the same-stencil interpolation
change (`2.53441e-6`--`3.86554e-6`). These disjoint points occur only
late in the window and exclude the earlier near-wall maximum. Spatial
derivative and bulk-transport accuracy therefore precede interpretation
as an actual wall or junction defect. The general identities above
retain the nonzero outgoing-source contribution throughout.

#### Separating spatial evolution and diagnostic width

The [fixed-input spatial comparison](prior-art-and-open-claim.md#fixed-input-spatial-refinement-16-october-2026)
changes the exterior and cap grids from 32 to 64 cells at the same 512
SSPRK3 steps. The source, analytic reference gauge, computational profile
and physical model remain identical. Compare exterior values at common
`eta=j/32`, and sum adjacent total physical cap integrals when restricting
the fine grid. These are common normalized coordinates; the evolved cap
radii need not be identical.

For a constraint diagnostic `D_h` and saved exterior solution `X_h`,
retain the separate differences

$$
D_{1/32}(\mathcal R X_{1/64})-D_{1/32}(X_{1/32}),\qquad
\mathcal R D_{1/64}(X_{1/64})-D_{1/32}(\mathcal R X_{1/64}).
$$

Here restriction `R` selects common nodes, with source and wall endpoints
retained. The first difference changes the evolved data at fixed
diagnostic width; the second changes diagnostic width on one evolved
history. At the wall, the stride-two fine-grid estimate has the same
width as the stride-one coarse-grid estimate. Retain the native fine-grid
estimate separately. Subtract inherited outgoing-source constraints only
after assembling the physical radius and actual source jets. Neither
difference supplies a continuum error bound, and neither assumes a zero
incoming constraint at the wall.

The [completed 32/64-cell comparison](current-result.md#fixed-input-spatial-refinement-16-october-2026)
finds a paired wall-scalar difference of `2.69431e-7`, or `0.774698%`
of the fine response peak. Geometry and clock differences are below
`0.071471%` of their respective fine peaks. At equal diagnostic width,
the incoming constraint maximum decreases from `5.54765e-4` to
`3.31042e-4`; its native fine-grid maximum is `3.53518e-4`.
On common nodes the equal-width evolution difference reaches
`2.75515e-4`, while the same-history derivative-width difference reaches
`4.97718e-5`. Thus the observed improvement involves the evolved data.
The corresponding source-subtracted outgoing maxima are `7.78011e-3`
(old), `5.15143e-3` (fine at equal width) and `4.93833e-3` (fine native).

The weighted wall estimates remain sensitive: equal-width evolution
differences reach `1.84444e-4` and `2.22042e-4`, while narrowing the
stencil on the fine history changes them by `1.14253e-4` and
`1.30794e-4`. These measurements support further spatial-operator work
but do not supply a converged wall datum. The unchanged profile's quiet
leakage and the source's curvature-trace convention remain separate
identified issues; the new grid does not remove them.

Small wall-observable changes do not bound the whole cap state. The
adjacent-sum paired cap `Iq` difference reaches `3.09202e-3`, about
7% of the old 32-cell paired peak. The own-trace cap comparison also
remains limited by its wall-history interpolation sensitivity.

The next implementation question is the first-order exterior departure
gradient in `background_stage`, used both for transverse prefix jets and
the moving-grid evolution rate. A successor needs one consistent
higher-order reconstruction for both uses, with the same source, profile
and cap laws. Its wall closure must be derived and checked on manufactured
inflow/advection data before a new coupled run. Substituting the desired
wall derivative into a diagnostic would not test this correction. The
existing frozen operator and comparisons remain the baseline; no further
grid doubling or changed profile is selected by this result alone.

#### A common second-order departure derivative

Let `Z` be a field departure on `eta_i=i h`, `h=1/N`, including its
physical source and wall endpoints. The source departure is zero by
construction; the wall departure comes from the cap primitive and evolved
wall geometry. Quadratic interpolation gives the candidate derivative

$$
G_i=\frac{-3Z_i+4Z_{i+1}-Z_{i+2}}{2h}\quad(0\le i\le N-2),\qquad
G_{N-1}=\frac{Z_N-Z_{N-2}}{2h},\qquad
G_N=\frac{3Z_N-4Z_{N-1}+Z_{N-2}}{2h}.
$$

These rows differentiate every quadratic exactly. For a sufficiently
smooth field their leading errors are `-h^2 Z'''/3`, `+h^2 Z'''/6`
and `-h^2 Z'''/3`, respectively. The centered penultimate row uses the
supplied wall value, without a ghost field outside the declared region.
The wall row reconstructs the gradient independently of the wall equation.

The same `G` enters both the prefix jets `delta X_Ubar=G/W` and the
moving-grid contribution `eta W_x G/W`. Replacing only one of these would
change their common discrete field interpretation. In particular, a desired
wall-compatibility derivative is not supplied as the reconstructed value.
The exact vacuum derivative remains analytic and is added separately.

For the isolated homogeneous-inflow equation `Z_x=a(eta) G`, with positive
`a` and only interior values evolved, the matrix is upper block triangular.
Its first diagonal blocks are `-3a_i/(2h)`; its last block is

$$
\frac1h\begin{pmatrix}-3a/2&2a\\-b/2&0\end{pmatrix},
\qquad a=a_{N-2}>0,\quad b=a_{N-1}>0.
$$

Its trace is negative and determinant `ab/h^2` positive. At constant speed
the two eigenvalues are `a(-3 +/- i sqrt(7))/(4h)`. This excludes a
positive-real eigenmode in that isolated matrix. It does not bound
nonnormal transient amplification, establish the explicit time step's
stability, or prove stability of the coupled cap, wall and exterior.
The [operator controls](prior-art-and-open-claim.md#second-order-exterior-departure-operator-16-october-2026)
therefore test known inflow solutions, finite amplification and complete
stage equivalence before a separately specified physical comparison.

The [completed same-grid comparison](current-result.md#fixed-grid-exterior-operator-comparison-16-october-2026)
reduces the final radius, lapse and scalar kinematic defects by about
21--33 times. The paired wall-scalar response changes by `0.524076%`
of its new peak; paired wall geometry and clock changes are below
`0.017798%`. At identical diagnostic widths, incoming-constraint RMS
decreases from `6.20959e-5` to `2.24521e-5`, and the source-subtracted
outgoing maximum from `4.93833e-3` to `3.21118e-3`.

This is a mixed accuracy result. The incoming maximum rises slightly
from `3.53518e-4` to `3.64417e-4`, at the last interior eta node.
Native weighted wall-curvature peaks increase, and varying the one-sided
width changes them by up to `1.80909e-4` or `3.02654e-4`. Thus improved
first-jet compatibility does not yet supply a resolved second jet.
The near-wall closure and independent curvature must be examined before
another physical run or a continuum junction conclusion. Cap spatial
sensitivity, profile support and source curvature sampling remain separate
issues. In particular, the general J/O relations above retain the cap
constraints; a future cap check must measure its temporal mass residual
instead of substituting the desired flux law into the acceleration.

The boundary-row diagnostic has two exact polynomial controls. Write
`Delta^3 Z_j=Z_(j+3)-3Z_(j+2)+3Z_(j+1)-Z_j`. On identical arrays,

$$
G^{\mathrm{center}}_{N-1}-G^{\mathrm{back}}_{N-1}
=\frac{\Delta^3Z_{N-3}}{2h},\qquad
G^{\mathrm{forward}}_{N-2}-G^{\mathrm{center}}_{N-2}
=-\frac{\Delta^3Z_{N-3}}{2h}.
$$

Thus the penultimate row changes the leading first-derivative truncation
coefficient by `h^2 Z'''/2`. This is a property of the stencils; it is
not an estimate of the accumulated coupled evolution error.

The wall diagnostic labelled `p3` uses three nodes for the first derivative
and four for the second; `p5` uses five for each. A perturbation `epsilon`
only at node `N-1` contributes `-5 epsilon/(h W)^2` and
`-26 epsilon/[3(h W)^2]` to their native wall curvatures, respectively.
The stride-two stencils omit that node entirely. Their second-derivative
Taylor errors on smooth exact samples start at orders two and three,
respectively. On evolved data these errors also contain the differentiated
state error. The [saved-row admission](prior-art-and-open-claim.md#saved-near-wall-row-audit-16-october-2026)
therefore separates stencil changes on identical arrays from changes in the
evolved arrays, retaining scalar stress and both background cross terms.
A concentrated last-row contribution alone cannot identify its physical or
numerical origin.

The [completed three-cut audit](current-result.md#saved-near-wall-row-audit-16-october-2026)
finds that localization. At cut 63 the new finite incoming constraint
changes from `9.28122e-5` at node 62 to `3.64417e-4` at node 63,
while the same-cut interior control is below `2.895e-10`. The node-63
value combines departure curvature `-0.00217152184`, background-radius
and lapse-departure coupling `+0.00252877654`, and smaller terms.
The new native/wide `p5` wall values at cut 68 are `6.38441e-4`
and `3.27035e-5`. These are unweighted constraints.

The radius third difference decreases between the old and new methods
at each selected cut, despite the unresolved local constraint. Therefore
localization alone does not select a corrected evolution stencil. The
same-grid diagnostic ends here. An independent cap temporal-mass and
radial-lapse check is the next selected calculation, retaining actual
reconstructed derivatives and incoming-field trace differences before
interpreting the J/O relations or replacing the boundary treatment.

#### Independent cap mass and wall-curve consistency

The next selected diagnostic reconstructs the cap from saved states, rather
than substituting its desired temporal mass equation into the wall force.
Write `m=r^3/2+alpha*n` and `P=3+alpha*Pdiff`, using the stable potential
quotient already returned by the cap implementation. At the moving wall
`r=R(x)`, the chain rule gives

$$
n_T\big|_r=\frac{D_x n_W-R_x D_r n_W}{T_x}.
$$

Here all temporal derivatives come from neighboring saved wall states and
reconstructed metrics; radial derivatives come from the returned cap face
arrays. The scalar fields below are the left reconstructed cap traces.
The normalized mass residuals and the lapse residual are

$$
e_m=n_T-\frac{R^2A^2e^\delta}{2}\Pi\Phi,\qquad
e_r=D_r n-\frac{R^2P_{\mathrm{diff}}}{2}
             -\frac{R^2A}{4}(\Pi^2+\Phi^2),
$$

$$
E_m=\alpha e_m,\qquad E_r=\alpha e_r,\qquad
E_\delta=D_r\delta-\frac{\alpha R}{2}(\Pi^2+\Phi^2).
$$

The brackets `e_m,e_r` have finite normalized-probe limits at `alpha=0`.
Their values in that limit are not nonzero vacuum Einstein residuals:
the actual mass residuals retain the factor `alpha`. The analytic vacuum
derivative `3r^2/2` is restored when reporting `m_r`; it is not replaced
by a finite difference of `r^3/2`.

An acceleration check can also use the actual saved curve. Define

$$
\sigma^2=Ae^{2\delta}T_x^2-\frac{R_x^2}{A},\qquad
v=\frac{R_x}{\sigma},\qquad
\beta=\frac{Ae^\delta T_x}{\sigma},\qquad
v_s=\frac{R_{xx}\sigma-R_x\sigma_x}{\sigma^3}.
$$

Require measured `T_x>0` and `sigma^2>0`. Compute `sigma_x` by
differentiating the displayed expression using
measured `A_W,x`, `delta_W,x`, `R_xx` and `T_xx`. Compare `sigma`
with the derivative of the saved proper clock instead of imposing their
equality. The independently reconstructed metric derivatives are

$$
A_r=-2R-2\alpha\left(\frac{n_r}{R}-\frac{n}{R^2}\right),\qquad
A_T=-\frac{2\alpha n_T}{R}.
$$

They enter the full `K_geom` expression above. Neither the hybrid predicted
acceleration nor the desired mass-flux law supplies these derivatives.
The incoming cap characteristic reconstructed on the left, the one
implied by the saved wall curve, and the one imposed by the hybrid map
remain separate reported quantities. Their associated flux differences
are measured explicitly, without replacing the primary left-trace residual.

The [completed audit](current-result.md#selected-cap-mass-and-wall-curve-audit-16-october-2026)
uses matched integration orders on all neighboring slices, two time and
radial widths, actual boundary-map reconstruction and the
[fixed allowance](prior-art-and-open-claim.md#selected-cap-mass-and-wall-curve-audit-16-october-2026).
Pulse-time pointwise residuals change strongly with radial width, while
quadrature-order sensitivity is small. At the terminal cut the finite
radial mass/lapse residuals are small, but the actual temporal mass
residual is `-8.319e-5` / `-7.626e-5` at the two time widths. The normalized
left-to-hybrid flux difference is only `-3.379e-5`, compared with the
normalized temporal residual `-0.640753` / `-0.587390`. That trace change
cannot account for the mass discrepancy. The next connection is the
reconstructed cap mass functional's balance under the existing discrete
field transport, including its moving radial cells. These observations
identify a numerical consistency question; they do not establish a new
physical source or a failed junction law.

#### Cap mass balance under field transport

The terminal discrepancy motivates a connection between the radial mass
functional and the scalar transport equations. On a smooth cap piece write
`S=Pi^2+Phi^2`, `B=A exp(delta)`, `U=dPdiff/df`, and

$$
a=\frac{\alpha r S}{2},\qquad
q=\frac{r^2P_{\mathrm{diff}}}{2}
       +\frac{r^2(1-r^2)S}{4},\qquad
J=\frac{r^2 A^2e^\delta\Pi\Phi}{2}.
$$

Then `e_r=n_r+a n-q`, `E_delta=delta_r-a`, and `e_m=n_T-J`.
Keep the field-equation defects explicit:

$$
\mathcal R_f=f_T-B\Pi,\qquad
\mathcal R_\Pi=\Pi_T-r^{-2}\partial_r(r^2B\Phi)+e^\delta U,
\qquad
\mathcal R_\Phi=\Phi_T-\partial_r(B\Pi).
$$

Differentiating the radial mass equation at fixed areal radius, using
`A_T=-2 alpha n_T/r`, and collecting the product derivatives gives

$$
\boxed{\partial_r e_m+a e_m
=\partial_T e_r+J E_\delta
 +\frac{r^2}{2}\left[U\mathcal R_f
       +A(\Pi\mathcal R_\Pi+\Phi\mathcal R_\Phi)\right].}
$$

In particular, the potential terms cancel when `f_T=B Pi` and the
`-exp(delta) U` term in `Pi_T` are inserted. The remaining scalar products
give `J_r+J delta_r`; moving `a J` to the right produces `J E_delta`.
This identity was derived independently twice. It explains why a small
radial defect at the wall alone does not establish the temporal mass law:
the latter also depends on field transport throughout the cap.

For piecewise reconstructed fields, apply the identity on each smooth
piece. Set `mu(r)=exp(integral_0^r a d rho)` and denote its right-hand side
by `F`. If `n` is continuous at a moving interface `r_j(T)`, differentiation
along that interface gives

$$
[n_T]_j=-r_{j,T}[n_r]_j,\qquad
[e_m]_j=-r_{j,T}[n_r]_j-[J]_j.
$$

Consequently, with regular-centre `e_m(0)=0`,

$$
\mu(R)e_m(R^-)
=\sum_{\mathrm{pieces}}\int\mu F\,dr
 +\sum_{0<r_j<R}\mu(r_j)[e_m]_j.
$$

This is piecewise calculus with explicit traces. It does not assign products
such as `Pi R_Pi` a meaning when the residual is distributional. Moving
cell-interface terms cannot simply be omitted from a numerical balance.

The actual finite-dimensional mass functional offers a direct next test.
At fixed `y=r/R`, let `N(y)=n(Ry)` and

$$
\mathcal A=\frac{\alpha R^2yS}{2},\qquad
\mathcal B=R^3y^2\left[\frac{P_{\mathrm{diff}}}{2}
                         +\frac{(1-R^2y^2)S}{4}\right].
$$

The exact radial problem is $N_y+\mathcal A N=\mathcal B$, $N(0)=0$.
Along the existing rates `(R_x,T_x,Ip_x,Iq_x)`, its tangent $M=\partial_x N$
satisfies

$$
M_y+\mathcal A M=\mathcal B_x-\mathcal A_x N,\qquad M(0)=0.
$$

The wall balance compares `M(1)` with `R_x n_r+T_x J`, retaining left and
imposed-face fluxes separately. Field tangents include the changing
background profile, cap reconstruction and moving radius. For the
implemented quadrature, its actual nodes, weights, split points and
recurrence must also be differentiated or tested by controlled directional
variation; the exact integral formula is not automatically its derivative.
Limiter ties and branch changes require separate one-sided treatment.
The zero-backreaction terminal profile reaches its upper supplied label,
so a forward finite perturbation can exit the domain. An analytic tangent
or explicitly backward probes must preserve that boundary.

The completed selected-state calculation below separates instantaneous spatial/flux compatibility from
time differentiation of saved histories, without supplying a continuum
error bound or changing the physical wall law.

For a direct finite-dimensional test, freeze the selected hybrid rates
`R_x=r_s/Vdot_plus`, `T_x=T_s/Vdot_plus` and its incoming `w_plus`. One
`background_cap_rhs` evaluation supplies the departure-integral rates.
Perturb `R,T,Ip,Iq` together along that direction, leaving the profile
fixed; reconstruct the mass from scratch at each admissible nearby input.
These are variations of one functional, not accepted evolution steps.
For the base value `N0` and backward values `N(-h)`, `N(-2h)`, use

$$
D_hN=\frac{3N_0-4N(-h)+N(-2h)}{2h}.
$$

Three widths `h,h/2,h/4` require four distinct backward probes. Compare
the resulting normalized temporal defect with the saved-history result,
and also report the global balance

$$
\mathcal B_{\mathrm{mass}}
=D_hN-R_x\rho_n-T_xJ,\qquad
\rho_n=\frac{R^2P_{\mathrm{diff}}}{2}
       +\frac{R^2 A(\Pi^2+\Phi^2)}{4}.
$$

If `e_m^dir=(D_hN-R_x D_rn)/T_x-J`, then
`B_mass-T_x e_m^dir=R_x e_r`. The global form avoids using a radial
difference to define its expected moving-boundary contribution, while
the displayed relation preserves the independently measured radial defect.
Both left and imposed-face flux comparisons remain explicit. Stable
nonzero instantaneous defects implicate spatial transport/reconstruction
compatibility; small ones redirect to time integration or saved-history
differentiation. Branch and width sensitivity can leave that attribution
unresolved.

The instantaneous direction uses hybrid wall/clock rates, whereas the
saved-history residual used measured rates. Keep that change separate:

$$
e_{\mathrm{hist,hyb}}
=\frac{D_x^{\mathrm{hist}}n_W-R_x^{\mathrm{hyb}}D_rn}
       {T_x^{\mathrm{hyb}}}-J,\qquad
e_m^{\mathrm{dir}}-e_{\mathrm{hist,hyb}}
=\frac{D_hN-D_x^{\mathrm{hist}}n_W}{T_x^{\mathrm{hyb}}}.
$$

Report `e_hist,hyb-e_hist` separately. The already recorded `n_T,T_x,R_x`
and radial derivative recover `D_x^hist n_W`; this bridge needs no new
metric or evolution. Without it, a directional-versus-history comparison
would also contain changed wall and clock rates.

The [completed directional calculation](current-result.md#terminal-cap-mass-functional-derivative-16-october-2026)
finds finest finite `e_m^dir=-0.6612504470`, with order sensitivity
`1.09497e-5`, compared with `e_hist,hyb=-0.6407403055`.
The last two widths differ by `0.000515001`; the finest stencil has no
recorded limiter branch changes. Thus the mass discrepancy survives
direct use of the implemented rates. The earlier small terminal radial
mass and lapse residuals do not imply temporal mass compatibility.
An exact discrete balance must now separate geometry/background changes,
boundary flux, interior flux and reaction contributions. This result
does not assign the discrepancy to a particular coupling term or establish
a physical junction failure. Further grid scans are not the next step.

#### Discrete mass derivative and departure-flux balance

The [diagnostic companion](../src/recursive_horizons/nsc_atlas_cap_mass_tangent.py)
differentiates the implemented mass-edge quadrature. It leaves the retained
metric and transport producers unchanged. At this stage its manufactured
controls and the four retained terminal cases pass, including creation
validation and a standalone check. The [completed ledger](current-result.md#discrete-cap-mass-balance-16-october-2026)
records the numerical contributions and the corrected accumulation bound.

On panel j, with half-width h and Gauss weights omega, the edge recurrence is

$$
N_{j+1}=e^{-A_j}N_j+S_j,\qquad
S_j=h_j\sum_k\omega_k e^{G_{jk}-A_j}b_{jk},
$$
$$
a=\alpha r(\Pi^2+\Phi^2)/2,\qquad
b=r^2P_{\rm diff}/2+r^2(1-r^2)(\Pi^2+\Phi^2)/4.
$$

Here A is the producer's edge partial-Gauss integral; G is its interpolated
primitive, evaluated at the edge quadrature nodes. Both are differentiated,
including their field samples, moving positions and half-width. Interior
target masses do not feed the final edge mass. A forward tangent and a
reverse sweep with weights `lambda_j=exp(-A_j)*lambda_(j+1)`, terminal weight
one, give separate contractions of the same local differentials. The
companion first checks every reconstructed primal face mass against the
unchanged producer.

One common physical direction selects the branch for negative parameter
increments before any component weight is taken. Zero minmod output can
open onto a nonzero slope derivative when both zero secants have compatible
directional rates. At a coincident cut, retain its preimages and order them
by their backward velocities; a born panel contributes
`width_rate*(b-a*N)`. Its regular-axis contribution vanishes because
`b-a*N` tends to zero there. Unresolved coincident-cut ordering is reported,
rather than silently dropping the panel. These are branch component weights,
not an assertion of an unrestricted Jacobian at a limiter tie.

Let `w_i^p` and `w_i^q` be the resulting derivatives with respect to the
departure integrals. The unchanged cap RHS is
`Ip_x[i]=Fp[i+1]-Fp[i]-T_x*reaction[i]` and
`Iq_x[i]=Fq[i+1]-Fq[i]`. A finite telescoping sum gives

$$
\begin{aligned}
D_xN={}&N_RR_x+N_TT_x+w^p_{n-1}F^p_n+w^q_{n-1}F^q_n
-w^p_0F^p_0-w^q_0F^q_0\\
&+\sum_{j=1}^{n-1}\left[(w^p_{j-1}-w^p_j)F^p_j
+(w^q_{j-1}-w^q_j)F^q_j\right]
-T_x\sum_iw^p_i\,\mathrm{reaction}_i .
\end{aligned}
$$

These are departure fluxes. The analytic background motion is already
present in N_R and N_T; adding its integral rates again would double-count
it. Only the centre q-flux is assigned zero. The centre p-flux remains in
the identity even though the reduced p field vanishes on the axis.

For characteristic speeds `v_+=-T_x*B-R_x*y` and
`v_-=T_x*B-R_x*y`, define each field jump as right minus left and
`Delta w=w_left-w_right`. The interior upwind contribution is exactly

$$
\sum_{j,\,\pm}\frac{|v_\pm|}{4}
(\Delta w^p_j\pm\Delta w^q_j)[p\pm q]_j .
$$

The rest is the central face contribution, including the metric change
acting on the known background. This upwind work need not be negative for
the nonlinear mass functional. At fixed state, metric and clock/wall rates,
the incoming trace sensitivity is

$$
\frac{\partial(D_xN)}{\partial w_+}
=(w^p_{n-1}+w^q_{n-1})\,R(T_xB_W+R_x)/2.
$$

The ledger must telescope and agree with the tangent before its total is
compared with `R_x*rho_n+T_x*J`. Large cancelling channels alone do not
identify a faulty equation. Attribution needs a matching mathematical
counterpart for the proposed term and an accounted-for complementary residual.

The [integrand derivatives](../src/recursive_horizons/nsc_atlas_cap_mass_derivatives.py)
also retain the producer's small-radius formulas. For `d=atanh(r)`, its
split integrals I1, I2 and Iw integrate F', F'' and z F'' on `[-1,1]`.
Their derivatives use the three moments `J_k=integral z^k F''' dz`:

$$
(I_{1,T},I_{1,d})=(I_2,I_w),\quad
(I_{2,T},I_{2,d})=(J_0,J_1),\quad
(I_{w,T},I_{w,d})=(J_1,J_2).
$$

C2 matching cancels internal moving-knot terms, despite possible jumps in
F'''. The actual truncated radial coefficient in Phi is differentiated as
written; it is not substituted for the derivative of `atanh(r)/r`.
Similarly, for `z=sqrt(alpha/6)*f`, `s=1+z^2/6+z^4/120`, and `c=cosh(z)`,
the potential's small-z branch has derivative

$$
\partial_fP_{\rm diff}
=f s\left[(c^2+1)(1+z^2/2+z^4/24)+z s c\sinh z\right].
$$

Replacing this by the exact-hyperbolic shortcut changes its rational z^6
coefficient from `37/45` to `33/40` after division by f. The production
difference is below ordinary binary64 rounding, so this distinction has
an exact symbolic control as well as high-precision branch comparisons.
The zero-alpha limit remains `Pdiff=f^2`, `Pdiff_f=2*f`.

The tangent's admitted domain excludes a quadrature sample exactly on
`r=0.01` or `abs(sqrt(alpha/6)*f)=1e-4`, where the implemented formulas
switch. This conservative restriction avoids assuming a derivative across
different truncated expressions. At the quiet-profile boundary, the
required zero value and first two jets make the first background
derivatives continuous at positive radius; an axis sample remains excluded.

Forward and reverse contraction use identical local gradient vectors. Their
floating summation check retains the absolute work
`W = sum_p attenuation_after_p * sum_j abs(local_gradient_pj * direction_j)`.
With P panels and K components, set `k=4*(P+K)+32` and
`gamma=k*eps/(1-k*eps)`. The allowance is
`gamma*W/(1-gamma)`, where `eps` is binary64 spacing at one.
This conservatively covers the two accumulation paths and rounding of the
positive work sum. Scaling by the small final contraction instead failed
on the terminal zero case because of intermediate cancellation. This is
an arithmetic check of two contractions of the same stored tape, not an
uncertainty estimate for the physical mass equation. Ledger, face-split
and saved-state tolerances keep their separate scales.

For the finite order-eight state, the implemented upwind pairing is
`-2.01783226502` of total balance `-2.12697432799`; the remaining pairing is
`-0.10914206297`. The next connection is the piecewise mass-transport
identity above: map its reconstructed field-equation defects, moving
interface jumps and quadrature defects onto this discrete remainder.
The radius and clock jets already differentiate all panel/node motion
and downstream attenuation. Adding a separate width-times-density term
would not recover a missing contribution. Similarly, a one-sided face
trace is not automatically the counterpart of the nonlocal mass weights.
No modified transport law is supplied by the one scalar mass constraint;
a proposed change must preserve its field equations and recheck the full
direction's limiter/cut branch. No further numerical run is admitted here.

#### Separating radial quadrature from reconstructed field transport

The ledger can be connected to the field equations without identifying its
nonlocal weights with pointwise face traces. Freeze the same full backward
direction `(R_x,T_x,Ip_x,Iq_x)` and the branch it selected. On the resulting
piecewise total fields, define a **shadow radial solution**

$$
n^E_r+a n^E=b,\quad n^E(0)=0,\qquad
\delta^E_r=a,\quad\delta^E(0)=\delta_0.
$$

Here `a,b` are the actual field/potential expressions in the preceding
recurrence, and `delta_0` is the same centre lapse. This exact linear radial
problem is a comparison for the existing reconstruction, not a different
cap candidate. Set `A^E=1-r^2-2 alpha n^E/r`, `B^E=A^E exp(delta^E)`.
The fields themselves still follow the implemented departure rates.

Let `N^Q` be the producer's wall mass, `D^Q` its checked full tangent, and
`E=N^Q-n^E(R)`. With left wall traces,

$$
\rho^E=b-a n^E,\quad\rho^Q=b-a N^Q,\qquad
J^E=\frac{r^2(A^E)^2e^{\delta^E}\Pi\Phi}{2},
\qquad \rho^E(R)-\rho^Q(R)=a(R)E.
$$

`J^Q` denotes the same left-trace flux used by the ledger. The imposed
incoming-face flux remains a separately reported comparison. Dots below
differentiate the complete moving-wall map; subscripts `x` on interior
fields mean the full selected direction **at fixed physical radius**.
They include dependence on `R`, the cell moments and the background clock.
They are not partial derivatives of the profile with frozen cell state.

Define the reconstruction's field-equation defects in that direction:

$$
\begin{aligned}
\mathcal R_f^x&=f_x-T_xB^E\Pi,\\
\mathcal R_\Pi^x&=\Pi_x-T_xr^{-2}\partial_r(r^2B^E\Phi)
                         +T_xe^{\delta^E}U,\\
\mathcal R_\Phi^x&=\Phi_x-T_x\partial_r(B^E\Pi),
\qquad U=\partial_fP_{\rm diff}.
\end{aligned}
$$

Use the derivative of the implemented potential branch in `U`. Substituting
a different small-field expression would add a different residual. With
`M=n^E_x-T_x J^E` and `W(r)=exp(-integral_r^R a ds)`, the earlier product
identity gives

$$
M_r+aM=\frac{r^2}{2}
\left[U\mathcal R_f^x+A^E(\Pi\mathcal R_\Pi^x+\Phi\mathcal R_\Phi^x)\right].
$$

The shadow mass is continuous. At an interior interface with velocity
`v_j=dr_j/dx`, right-minus-left jumps therefore obey
`[M]_j=-v_j[rho^E]_j-T_x[J^E]_j`. Integrating from the regular centre,
where `M(0)=0`, and converting the moving-wall derivative yields

$$
\boxed{\begin{aligned}
\mathcal B_{\rm mass}^Q={}&\dot E+R_xa(R)E+T_x(J^E_W-J^Q_W)\\
&+\sum_{\rm pieces}\int W\frac{r^2}{2}
  \left[U\mathcal R_f^x+A^E(\Pi\mathcal R_\Pi^x+\Phi\mathcal R_\Phi^x)\right]dr\\
&-\sum_{0<r_j<R}W(r_j)\left(v_j[\rho^E]_j+T_x[J^E]_j\right).
\end{aligned}}
$$

The first line isolates radial quadrature and its wall-flux conversion.
The remaining lines measure compatibility of the reconstructed field rates
with the shadow geometry, including its moving interfaces. The wall is
excluded from the jump sum because its contribution is already explicit.
Cell faces have `v_j=R_x*y_j`; the two profile-knot families have
`v_j=+(1-r_j^2)*T_x` and `v_j=-(1-r_j^2)*T_x`. The exact integration
partition must also account for formula switches, including `r=0.01`
(velocity zero) and the potential threshold where its level set supplies
a nonzero jump. Tangencies, coincidences and branch openings need their
one-sided limits; a zero-width panel cannot be discarded before taking a
directional derivative.

There is an independent route to the first channel. On each original
producer panel, let `L_j^Q` be its actual affine edge map and `L_j^E` the
exact radial flow of the same reconstructed fields. With the numerical
entering mass `N_j^Q`, set

$$
\varepsilon_j=(L_j^Q-L_j^E)(N_j^Q),\qquad
\sigma_{j+1}^E=\exp\left(-\int_{r_{j+1}}^R a\,dr\right).
$$

Affine telescoping, using the **exact** downstream attenuation, gives

$$
E=\sum_j\sigma_{j+1}^E\varepsilon_j,\qquad
\dot E=\sum_j\left(\dot\sigma_{j+1}^E\varepsilon_j
                         +\sigma_{j+1}^E\dot\varepsilon_j\right).
$$

The derivative of each local error includes its moving endpoints, field
samples and entering mass `N_j^Q`. Both exact and implemented panel maps
participate. Additional formula-switch splits may be needed inside an
exact flow even when the producer did not split there. The stored
`gamma_interpolant_defect_nodes` is not this panel error. Defining `dot E`
as `D^Q-d(n^E(R))/dx` would make the final comparison circular; that
difference can only be an independent check of the differentiated panel sum.

The [exact controls](../tests/test_nsc_atlas_cap_mass_transport_identity.py)
check the full-parameter product identity without imposing the field
equations, the moving-wall conversion, differentiated affine telescoping,
and a manufactured moving-interface example with an independently
integrated radial mass. Omitting attenuation derivatives, entering-mass
derivatives, the radius term, or the interface jump fails its control.
Four controls pass with zero field samples and 0.640625 active CPU seconds.
An independent Grok derivation agrees with the corrected full-direction
formula; its initial split into profile time and other state changes is
absorbed into the full fixed-radius derivative above.

`CapTransportFields` in the [field adapter](../src/recursive_horizons/nsc_atlas_cap_transport_fields.py)
implements the field inputs to this identity. Its first tangent column
contracts `(R_x,T_x,Ip_x,Iq_x)` at fixed physical radius; its second varies
only the query radius. The common backward direction selects the same
limiter branch as the mass tangent. A query names its owner cell, including
at a shared face, so the material trace derivative is explicitly
`X_x+v_j*X_r`. The adapter returns `f,Pi,Phi,a,b` and both derivatives,
plus the implemented `Pdiff` and `Pdiff_f`. It freezes the profile and
reconstruction inputs, and charges the caller's existing budget.

The [23 manufactured controls](../tests/test_nsc_atlas_cap_transport_fields.py)
cover regular polynomial coefficients, nonconstant field derivatives,
upstream moments, first-cell regularity, limiter openings, both material
face traces, profile endpoints, input immutability and budget rejection.
Axis jets and exact queries at `r=0.01` or the potential formula threshold
are excluded. The interface inventory and one-sided switch limits required
by the shadow integration are separate work; this adapter does not silently
select those limits or integrate a new metric.

No value of these new channels has been calculated for the retained cap.
Their evaluator needs controlled exact-flow approximation, positive `A^E`
and finite lapse on its declared domain, full reconstructed field jets,
branch/interface accounting and an independent error budget.
Once the complete balance closes through these independent channels, the
already measured upwind work can be subtracted to resolve the remaining
`-0.10914206297`. A frozen-branch central-flux subtraction is still an
algebraic comparison, not a validated replacement evolution.

#### Interfaces of the reconstructed radial comparison

The literal-potential comparison below is mathematically valid, but its
threshold-root inventory is not required by the selected smooth-potential
reference described after these interface identities.

The departure scalar is continuous at every cell face. In an ordinary cell,
`q=qc+qs*r` has mean `Iq/(b-a)` and `qc=Iq/(b-a)-qs*(a+b)/2`, so its
primitive increases by exactly `Iq`. Both traces of `f=chi/r` therefore
agree. The regular first cell also satisfies
`Q0*b+Q2*b^3/3=Iq[0]`. Its quadratic term is necessary for that identity.
The independent limited traces of `Pi` and `Phi` can still jump. At a pure
profile-knot preimage, the background fields are continuous because the
profile is C2; both mass-density and mass-flux jumps vanish there. Keep the
preimage as an integration breakpoint for the piecewise higher derivatives.

The implementation switches require separate accounting. At `r0=0.01`,
the background scalar and `Pi` formulas agree in real arithmetic, while
the right-minus-left `Phi` jump is `DeltaC*I1`, where `I1` is the common
split integral of `F'` and

$$
\Delta C=\frac{r/(1-r^2)-\operatorname{atanh}r}{r^2}
 -\sum_{k=1}^{5}\frac{2k}{2k+1}r^{2k-1}
 =\frac{12}{13}r^{11}+\frac{14}{15}r^{13}+\cdots.
$$

For `0<r<1`, this positive remainder lies between `(12/13)*r^11`
and `r^11/(1-r^2)`. The switch has velocity zero. With continuous shadow
mass and lapse, its density jump is
`[rho]=r^2*A^E*(2*Phi_left*[Phi]+[Phi]^2)/4`; its flux jump must also be
retained. Roundoff from the two evaluated formulas is a separate error.

For `alpha>0`, write `q=sqrt(alpha/6)*f`, `q0=1e-4`, and
`f_star=q0*sqrt(6/alpha)`. Potential switches occur at `f=+/-f_star`.
With `S4=1+q^2/6+q^4/120`, the real-arithmetic difference between the
outer hyperbolic branch and inner polynomial branch is

$$
\Delta P=f^2D(q),\qquad
D(q)=\frac{\cosh^2q+1}{2}
 \left[\left(\frac{\sinh q}{q}\right)^2-S_4(q)^2\right]
 =\frac{q^6}{2520}+O(q^8).
$$
$$
\Delta U=f\,[2D(q)+qD'(q)]=\frac{f q^6}{315}+O(f q^8).
$$

These differences are nonzero. At a simple root the full fixed-radius jets
give `v=-f_x/f_r`. With no other coincident switch, `a` and `J` are
continuous, while `[rho]=r^2*[P]/2`. The right-minus-left sign of `[P]`
follows whether `|f|` enters or leaves the polynomial branch. The term
`-W*v*[rho]` belongs in the differentiated balance. Its omission cannot be
justified from the small value jump alone: a nearly tangent crossing can
have large velocity. At `alpha=0` this threshold is absent.

A complete numerical inventory must isolate all threshold roots on each
analytic piece, or bound the omitted contribution independently. Endpoint
sign samples alone cannot establish absence of roots. A non-simple root,
unresolved root enclosure, or unresolved coincident motion stops that
comparison. At coincident stations retain their preimages, order the
backward ray by `(r,-v)`, and keep any born zero-width panel until its
directional derivative is included. A root on a cell face uses each side's
continuation and its owner admissibility; continuity of `f` does not imply
equality of the two radial slopes. Axis and wall remain endpoint terms.

Three [exact controls](../tests/test_nsc_atlas_cap_interface_identity.py)
prove the cell primitive identities, both switch leading orders, and a
moving potential-jump contribution against an independently integrated
manufactured mass. They pass with zero field samples and 1.71875 active CPU
seconds. No threshold roots or shadow solution have been computed for the
retained cap. The current field adapter intentionally rejects exact switch
queries; branch-limit evaluation and root completeness are still required
by the literal-branch comparison. The selected reference below removes its
potential-root requirement explicitly.

#### Selected smooth-potential radial reference

Use the smooth potential already specified by the model,

$$
P_*(f)=\frac{f^2}{2}\left(\frac{\sinh q}{q}\right)^2
(\cosh^2q+1),\qquad U_*=2f\frac{\sinh q}{q}\cosh^3q,
\qquad q=\sqrt{\alpha/6}\,f,
$$

with the continuous `q=0` limits. The implemented `P_Q` is still used at
every native quadrature node and differentiated on its selected branch.
Define the reference by `n^E_r+a*n^E=b_*`, with the same reconstructed
fields, `a`, and centre lapse as before; only replace the potential inside
`b_*`. Write `E=N^Q-n^E(R)`. Now

$$
\rho^E-\rho^Q=aE+(b_*-b_Q),
$$
$$
\boxed{\begin{aligned}
\mathcal B^Q_{\rm mass}={}&\dot E
 +R_x\bigl(a(R)E+(b_*-b_Q)_W\bigr)+T_x(J^E_W-J^Q_W)\\
&+\sum\int W\frac{r^2}{2}
 [U_*\mathcal R_f^x+A^E(\Pi\mathcal R_\Pi^x+\Phi\mathcal R_\Phi^x)]\,dr\\
&-\sum_jW_j\bigl(v_j[\rho^E]_j+T_x[J^E]_j\bigr).
\end{aligned}}
$$

Use `U_*` in both the local residuals and this bulk expression. The jump
sum uses the reference density `[b_*]-n^E[a]`. Potential thresholds are
absent from the reference; cell/profile interfaces and the fixed small-radius
switch remain. The exact panel-error sum and its derivative are unchanged
after replacing `L_j^E` by the smooth reference map. Consequently `dot E`
contains **both quadrature and potential implementation error**, with no
additional interior integral of `b_*-b_Q`. Counting that difference again
would double-count it. The wall conversion above still requires its one
explicit value, `r^2*(P_*-P_Q)/2`.

This is a diagnostic reference, not a change to the retained equations,
source or evolution. It avoids locating discontinuities introduced only by
the potential's numerical approximation. Native nodes exactly on that
branch threshold still require their previously declared branch guard.
The tiny gap can be evaluated from the positive sinh-series remainder,
without subtracting almost equal potentials. Grok independently confirmed
the revised balance, the wall-term sign and the absence of another bulk term.

For a smooth panel with moving endpoints `l(x),h(x)`, let `L=h-l` and
`r=l+u*L`, `0<=u<=1`. Its affine map is `n_out=S*n_in+T`. With
`A=L*a`, `B=L*b_*`, differentiate the map in this fixed panel coordinate:

$$
\begin{aligned}
S_u&=-AS,& T_u&=B-AT,\\
\dot S_u&=-A\dot S-\dot A S,&
\dot T_u&=\dot B-\dot A T-A\dot T,\\
\dot A&=\dot L a+L(a_x+\dot r a_r),&
\dot B&=\dot L b_*+L((b_*)_x+\dot r(b_*)_r).
\end{aligned}
$$

The initial map is `(S,T,dot S,dot T)=(1,0,0,0)`. At a born panel
`L=0` on an admissible backward ray, the derivative is
`dot S=-dot L*a`, `dot T=dot L*b_*`; its value remains the identity map.
An independent check integrates the fixed-radius variation
`nu_r+a*nu=(b_*)_x-a_x*n`, starting at
`nu(l)=dot n_in-rho(l)*dot l`, and converts back at the right endpoint by
`dot n_out=nu(h)+rho(h)*dot h`. Agreement with that check and manufactured
exact maps tests endpoint motion separately from the panel-coordinate flow.

The [reference component](../src/recursive_horizons/nsc_atlas_cap_radial_reference.py)
implements these four map equations. The caller supplies one smooth piece's
six coefficients `(a,b_*,a_x,(b_*)_x,a_r,(b_*)_r)` and its endpoint limits;
the component does not discover interfaces or select their owners. It uses
the same cumulative budget for callback work and map evolution. Requested
ODE tolerances and solver completion are not a uniform error enclosure.

Its potential gap uses positive series terms through `q^10`. The next tail
is bounded in real arithmetic by `|q|^12/13!/(1-q^2/210)`, and its
`q*d/dq` derivative by `12*|q|^12/13!/(1-q^2/180)`. Propagating these
positive tails through `D` and `2D+qD'` gives the returned series bounds.
Floating-point evaluation and the reference integration errors remain
separate; these series bounds alone do not certify either.

The adapter's `reference_sample()` uses this smooth potential. At `r=0.01`
it requires `small_radius_side="left"` or `"right"`, independent of the
explicit cell owner. The left trace reuses the polynomial-moment branch;
the right evaluates the ordinary formulas at the same radius, without a
displaced sample. `reference_axis_coefficients()` returns the regular limits
`a=b=a_x=b_x=b_r=0`, `a_r=alpha*Pi(0)^2/2`. No radial derivative of `Pi`
at an axis profile knot is asserted. These coefficient limits suffice at
the fixed centre endpoint of the affine map. Both switch traces have
independent cubic-profile controls, and the reference can pass a native
potential threshold without creating an artificial interface.

The [cap composer](../src/recursive_horizons/nsc_atlas_cap_radial_comparison.py)
preserves the original quadrature partition for each native map. Reference
maps additionally split at `r=0.01`, ordered by `(r,-r_x)` along the backward
ray. This retains zero-width pieces that open under endpoint motion. Each
reference piece inherits its native owner; no cell face lies inside that
native panel. A native born panel exactly at a formula switch remains outside
the original tangent's domain and is rejected, not assigned a new branch.

For consecutive maps, composition is `S=S2*S1`, `T=S2*T1+T2`, with the
full product rule for both derivatives. The local differences use the native
entering mass and rate, while the backward weighting uses the reference maps
and their derivatives. Thus the reported error derivative is constructed from
local map differences; subtracting the two final derivatives is an independent
algebra check. Positive absolute recurrences retain cancellation work. With
`p` native panels, `k=32*p+64` exceeds the operation count along each expanded
product path; `gamma=k*eps/(1-k*eps)` gives the recurrence comparison allowance
`gamma*absolute_work/(1-gamma)`. This allowance treats the computed maps as
inputs. Coefficient evaluation and radial integration errors are separate.

Seven manufactured controls include an exact zero-coupling constant-field
mass, independent fixed-radius variation on a nonconstant cap, moving faces
at the fixed switch, and a two-sided born map. They test the local-error
channel only. The complete balance still needs the reference lapse, bulk
field residuals, interface flux/density terms and controlled integration error.
No retained-state comparison or corrected evolution follows from these tests.

#### Reference lapse and field-residual channels

The [smooth-piece integrator](../src/recursive_horizons/nsc_atlas_cap_transport_reference.py)
uses six states: `n,delta,Z_f,Z_Pi,Z_Phi,Z_jump`. It integrates
`n_r=b_*-a*n`, `delta_r=a`, each `Z_i,r=F_i-a*Z_i`, and
`Z_jump,r=-a*Z_jump`, where

$$
F_f=\frac{r^2}{2}U_*\mathcal R_f^x,\qquad
F_\Pi=\frac{r^2}{2}A^E\Pi\mathcal R_\Pi^x,\qquad
F_\Phi=\frac{r^2}{2}A^E\Phi\mathcal R_\Phi^x.
$$

The radial product derivatives use
`A_r=-2*r-2*alpha*(rho/r-n/r^2)` and
`B_r=exp(delta)*(A_r+A*a)`. Multiplying the Pi residual by `r^2`
before evaluation avoids unnecessary singular divisions. The regular centre
has zero mass, accumulated channels and RHS, with its supplied centre lapse;
axis field jets are not queried. A zero-width piece is the identity map for
these states. This component never reads the native mass tangent or the
local radial-map error.

The complete caller must supply the same centre lapse as the existing cap,
split at every reference interface, and add `-v*[rho]-T_x*[J]` to `Z_jump`
between adjacent pieces. That interface assembly is not implemented by the
smooth-piece integrator. Select the global centre and wall by their ordered
event identities. An interior interface may coincide with `R` at the
evaluation time and still contribute in the backward limit; its term must
not disappear through an instantaneous `r<R` filter. Coincident interfaces
retain separate velocities and intermediate traces. At the centre, the
regular `rho,J` limits make their contributions vanish. Seven exact transport
controls include these wall/centre and coincident-cut cases.

Five manufactured numerical controls check the smooth-piece integration
against exact polynomial integrals and an independent 60-digit constant-
attenuation mass integral. The reported minimum `A` is over solver samples,
not a proof of positivity between them. Requested solver tolerances likewise
do not enclose integration error. Interface assembly, a controlled error
comparison and application to the retained cap remain open.

The [full cap assembler](../src/recursive_horizons/nsc_atlas_cap_transport_balance.py)
now connects those pieces. `reference_transport_cap()` authenticates the
centre values against the frozen reconstruction and integrates every ordered
piece, adding each interface event to the jump channel. It independently
reconstructs the weighted jump sum from the lapse increments. It never reads
the native mass tangent. `cap_transport_balance()` separately obtains the
affine local-map error and combines it with the wall, bulk and jump terms in
the selected-reference identity. The wall fields must share a limiting trace;
a wall exactly at `r=0.01` remains excluded from this full comparison.

The affine chain and residual ODE are distinct numerical solutions of the
reference radial problem. Their mass and material-rate differences, and the
jump-sum discrepancy, are returned alongside the full balance defect. They
are checks, never replacements for a missing channel. A returned defect is
not automatically accepted: retained-state interpretation requires the
explicit accuracy checks in the [calculation record](prior-art-and-open-claim.md#next-fixed-terminal-comparison).
Ten manufactured controls cover exact constant caps, nonconstant coupled
rates, coincident wall/axis/switch events and rejected input inconsistencies.
The existing coupled-rate ledger is recovered with defect about `8.1e-16`.
The decimal constant-cap control also tests native branch agreement when
directional rates cancel. The adapter now selects its limiter in the exact
same component basis and arithmetic order as `mass_tangent`, then contracts
each reconstruction jet to the full-direction column. The independent radial
query column remains zero in reconstruction data and one in the query radius.
This fixes branch-representative drift without changing the retained native
tangent, imposing a residual or weakening the branch check. The first retained
attempt stopped at that branch comparison. The corrected two-state terminal
comparison now passes all declared accuracy gates; its
[recorded channel attribution](current-result.md#terminal-cap-transport-attribution-16-october-2026)
locates the imbalance in reconstructed field transport and interfaces, with
radial-map contributions below the stated diagnostic resolution. It supplies
no corrected evolution or uniform finite-amplitude bound.

#### Variational mass response and flux projection

The terminal attribution identifies field transport and reconstructed
interfaces as the unresolved connection. The following exact calculation
relates the mass functional to that transport. It uses the same equations,
with an exact scalar primitive and smooth potential; it does not change the
saved evolution. The [fifteen symbolic controls](../tests/test_nsc_atlas_cap_flux_projection_identity.py)
use no retained field samples or numerical radial solves.

Write `P=Pdiff`, `U=P_f`, and use total fields

$$
\chi=rf,\qquad q=\chi_r,\qquad p=r\Pi,\qquad
A=1-r^2-2\alpha n/r,\qquad B=Ae^\delta.
$$

Take a regular centre, `n(0)=chi(0)=0`, fixed `R,T`, and piecewise smooth,
integrable variations with fixed interfaces. The exact radial constraint is

$$
n_r=\rho=\frac{r^2P}{2}+\frac A4[p^2+(q-f)^2],\qquad
a=\delta_r=\frac{\alpha}{2r}[p^2+(q-f)^2],\qquad
W(r)=\exp\left(-\int_r^R a(s)ds\right).
$$

Denote a variation by `d`, to distinguish it from the lapse. Differentiating
this constraint, including the dependence of `A` on `n`, gives

$$
(dn)_r+a\,dn=\frac{r^2U}{2}df+
\frac A2\{p\,dp+(q-f)(dq-df)\}.
$$

Since `W_r=aW`, integration to the wall gives
`dn(R)=integral(h_p dp+h_q dq+h_f df) dr`, where

$$
h_p=\frac{WAp}{2},\quad h_q=\frac{WA(q-f)}2,\quad
h_f=\frac W2[r^2U-A(q-f)].
$$

The primitive relation fixes `df(r)=integral_0^r dq(s)ds/r`. Changing the
order of integration therefore yields the functional derivatives

$$
\boxed{dn(R)=\int_0^R(g_p\,dp+g_q\,dq)dr,\qquad
g_p=h_p,\quad g_q=h_q+\int_r^R\frac{h_f(s)}s ds.}
$$

The integral in `g_q` retains the field-value/potential contribution. For
a branch-fixed reconstruction with moments `I`, its mass weight is
`w_i=integral(g_p p_{I_i}+g_q q_{I_i})dr`. Both shape derivatives enter:
characteristic limiting can make `q` depend on a `p` moment and conversely.
A face trace of `g` is not this integrated moment weight.

The radial constraints also give
`A_r+aA=-2r+2alpha n/r^2-alpha rP`. Thus, for the original reaction

$$
S=e^\delta r[(U-2f)+\alpha(2n/r^3-P)f],\qquad
H(r)=\int_r^R S(s)ds,
$$

the same derivative has the useful equivalent form

$$
g_p=\frac{e^{-\delta_W}}2Bp,\qquad
\boxed{g_q=\frac{e^{-\delta_W}}2(Bq+H)-\frac{A_Wf_W}{2}.}
$$

Indeed, `(g_q-WAq/2)_r=-WS/(2e^delta)`, and `W/e^delta=e^{-delta_W}`;
the wall value supplies the last constant. With `u=Bp`, `v=Bq`,
`c=e^{-delta_W}` and `C=A_Wf_W/2`, the original field equations
`p_T=v_r-S`, `q_T=u_r` then satisfy

$$
g_pp_T+g_qq_T
=\partial_r\left[\frac c2u(v+H)-Cu\right].
$$

The reaction cancels because `H_r=-S`. At the regular centre the primitive
vanishes; at the wall it equals
`J_W=A_W^2 exp(delta_W) p_W(q_W-f_W)/2`. This proves `n_T(R)=J_W`
for the stated continuum equations. No independent removal of the reaction
or field-value channel is compatible with this calculation.

At an ordinary fixed interface, `f,n,A,W` and the tail integral are
continuous, although `p,q` can jump. Consequently

$$
[g_p\pm g_q]=\frac{WA}{2}[p\pm q].
$$

Use right-minus-left jumps and left-minus-right cell-weight differences
`Delta w`. The native characteristic speeds are
`s_+=-T_x B-v_j`, `s_-=T_x B-v_j`, with radial face velocity `v_j`.
For `j_+=[p+q]`, `j_-=[p-q]`, the exact interior upwind work is

$$
\begin{aligned}
U_j&=\frac14\sum_{\sigma=\pm}|s_\sigma|
 (\Delta w^p+\sigma\Delta w^q)j_\sigma\\
&=-\frac{WA}{8}\sum_{\sigma=\pm}|s_\sigma|j_\sigma^2
 +\frac14\sum_{\sigma=\pm}|s_\sigma|\eta_\sigma j_\sigma,\\
\eta_\sigma&=\Delta w^p+\sigma\Delta w^q+\frac{WA}{2}j_\sigma.
\end{aligned}
$$

For `A>0`, the first term is nonpositive. The projection remainder has no
fixed sign. This is not merely freedom in arbitrary weights: at `alpha=0`,
take constant `p=1,9/8` in the interior cells `[1/4,1/2]`, `[1/2,3/4]`,
with `q=f=0`, `T_x=1`, and a stationary face. Differentiating their actual
mass integrals with respect to their integral moments gives upwind work
`67/16384>0`: trace damping `-9/4096` plus projection `103/16384`.
This counterexample concerns that piecewise constant representation, not
the native limited-linear state or a nonzero gravitational response at
`alpha=0`.

The formula may use native weights `w^Q`, but the reference functional just
derived, `N*`, requires `q=(rf)_r` exactly. It is not automatically identical
to the earlier diagnostic reference, which retains native field formulas.
In particular, the native small-radius approximation for `Phi` must be
accounted for before replacing it by `(q-f)/r`. The weight difference
`w^Q-w*` can contain field-representation and potential errors as well as
radial quadrature error. No value for that difference is supplied here.
Ordinary-face variations keep the background clock, radius, physical
endpoint fluxes and interface locations fixed. Their other contributions
remain in the full ledger; this local identity does not discard them or
merge coincident moving events.

#### A mass-gradient route and its geometric reference limit

The preceding adjoint gives a stronger structural consequence:

$$
\boxed{p_T=2e^{\delta_W}\partial_r g_q,\qquad
q_T=2e^{\delta_W}\partial_r g_p.}
$$

These are the original field equations, including the reaction through
`g_q`, expressed through derivatives of the constrained mass. At fixed
geometry a differentiable discrete mass with moment gradients `w` therefore
suggests the central gradient fluxes
`F_p=e^{delta_W}(w_L^q+w_R^q)` and
`F_q=e^{delta_W}(w_L^p+w_R^p)`. Their interior mass work is exactly
`2e^{delta_W}(w_L^p w_L^q-w_R^p w_R^q)`. For cells `0,...,N-1`,
the sum leaves the explicit boundary expression

$$
\dot N=
w_{N-1}\cdot F_W-2e^{\delta_W}w_{N-1}^p w_{N-1}^q
-w_0\cdot F_C+2e^{\delta_W}w_0^p w_0^q.
$$

This is an exact algebraic route for constructing both field equations,
not an added scalar correction to a computed mass rate. It is not yet a
consistent cap scheme: physical incoming data must give the required wall
term, reconstruction gradients must approximate the continuum derivatives,
and moving geometry, analytic background transport, limiter branches and
time integration still require their own compatible identities. A reaction
already included through the gradient cannot be added again as a separate
source. Arbitrary boundary fluxes do not satisfy the required wall identity.

At `alpha=0`, in the existing centre-normalized gauge `delta=0`, the
structure can be made local. Here `P=f^2`, `A0=1-r^2`, and

$$
\rho_0=\frac{A_0}{4}(p^2+q^2)
-\partial_r\left(\frac{A_0rf^2}{4}\right).
$$

Use the already established geometric coordinate `z=atanh(r)` and define
`u=A0 p`, `v=A0 q`. The moments are unchanged,
`integral p dr=integral u dz`, `integral q dr=integral v dz`, while

$$
u_T=v_z,\qquad v_T=u_z,\qquad
n(R)=E-K_W,\quad E=\int_0^{z_W}\frac{u^2+v^2}{4}dz,\quad
K_W=\frac{A_0(R)\chi_W^2}{4R}.
$$

This adds an energy/mass relation to the existing characteristic background;
the coordinate and wave solution themselves are not new. For a moving wall,
`z_{W,x}=R_x/A0(R)` and `chi_{W,x}=T_x u_W+z_{W,x}v_W`. Differentiating both
`E` and `K_W`, including their moving endpoints, gives exactly
`n(R)_x=R_x rho_W+T_x J_W`, with
`J_W=u_W(v_W-A0(R)f_W)/2`. Omitting `K_W` loses this identity.

For piecewise constant `u,v` on moving `z` cells of length `h`, the discrete
energy is `(I_u^2+I_v^2)/(4h)`. Its geometric derivative `-h_x(u^2+v^2)/4`
cancels the interior mesh-velocity term in central transport. Upwinding then
adds exactly

$$
-\frac18\left(|-T_x-z_{j,x}|[u+v]^2
+|T_x-z_{j,x}|[u-v]^2\right).
$$

There is no projection remainder for this particular quadratic energy and
constant reconstruction. Its physical boundary completion is also explicit.
Let `G=(u+v)_incoming=A0(R)(R w_plus+f_W)`, retain the outgoing last-cell
`u-v`, and set `u_hat=(G+u-v)/2`, `v_hat=(G-u+v)/2`. Use these two traces
in both conservative wall fluxes. At the centre, retain the outgoing `u+v`
and impose `u_hat=0`, giving `v_hat=u+v` and zero q flux. Since
`chi_W=sum I_v`, the same wall flux supplies
`chi_{W,x}=T_x u_hat+z_{W,x}v_hat`.

For `T_x>0` and timelike faces `|z_{j,x}|<T_x`, the resulting semidiscrete
mass `N_h=sum (I_u^2+I_v^2)/(4h)-K_W` obeys

$$
\dot N_h=R_x\rho(\widehat u,\widehat v,f_W)
+T_xJ(\widehat u,\widehat v,f_W)-D_{\rm int}-D_C-D_W,
$$

where `D_int` is the sum of the positive jump squares above,
`D_C=T_x u_0^2/2`, and
`D_W=(T_x+z_{W,x})(u_{N-1}+v_{N-1}-G)^2/8`. A three-cell symbolic
calculation checks the complete moving mass identity directly from the
conservative moment rates. It includes the wall scalar and all cell-width
derivatives. The interface and boundary losses are numerical, not extra
physical matter. This is a piecewise constant weak-field control, not a
pointwise regular centre reconstruction or a repaired nonzero-alpha cap.

Higher-order regular shapes, finite-backreaction gradients and the actual
wall dynamics remain to be constructed together. These exact results select
a compatible-gradient construction and geometric energy control as the next
implementation route. They do not supply a new evolved history, a nonlinear
stability theorem, or finite-amplitude cap-wall closure.

#### Regular-centre fixed-cap mass-gradient control

The [geometric reference component](../src/recursive_horizons/nsc_atlas_cap_geometric_reference.py)
constructs one instantaneous weak RHS at `alpha=0`, fixed `0<R<1` and
`delta=0`. It uses continuous piecewise linear `chi(z),u(z)` on
`0=z_0<...<z_N=atanh(R)`, with `chi_0=u_0=0`. Here `u=chi_T` and
`v=chi_z`; the first element therefore has `chi=c z`, `u=d z`, giving
the exact centre limits `f(0)=c`, `Pi(0)=d`, `Phi(0)=0`. It supplies the
regular-centre reference missing from the constant-cell control above.

For an element of length `h`, use the exact matrices

$$
M_e=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad
K_e=\frac1h\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
$$

Assemble them and eliminate the fixed centre degree of freedom. With
`e` selecting the wall, `kappa=A0(R)/R`, set `S=K-kappa e e^T`. The mass
and weak equations with the existing incoming `w_plus=Pi+Phi` are

$$
N_h=\frac14(u^TMu+\chi^TS\chi),\qquad
\chi_T=u,\qquad
M u_T=-S\chi+e\{A_0(R)R w_+-u_W\}.
$$

The boundary conversion is
`v_hat=A0(R)(R w_plus+chi_W/R)-u_W`. In the usual weak wave equation
`M u_T=-K chi+e v_hat`, its scalar term is precisely the wall subtraction
already included in `S`. The same scalar enters the field and mass equations.
Multiply both rates by the supplied positive `T_x` to obtain x-clock rates.

The mass is positive definite. The trace inequality
`chi_W^2 <= z_W integral chi_z^2 dz` and

$$
\gamma=\frac1{z_W}-\frac{A_0(R)}R>0
$$

give coercivity factor `beta=z_W gamma`. Positivity follows from
`R-(1-R^2)atanh(R)>0`: it vanishes at zero and its derivative is
`2R atanh(R)>0`. In particular,

$$
N_h\ge\frac14u^TMu+\frac\beta4\int_0^{z_W}\chi_z^2dz,\qquad
\chi^TS\chi=\sum_e h_e\left(\chi_{z,e}-\frac{\chi_W}{z_W}\right)^2
+\gamma\chi_W^2.
$$

The implementation uses the latter expression and its gradient, avoiding
subtraction of nearly equal energies at small radius. It evaluates
`gamma=1/z-2/sinh(2z)` through the positive sinh remainder for `2z<=1`.
The first omitted series contribution is below `2.4e-22` relative to its
leading term in real arithmetic; floating-point evaluation is checked
separately against high precision and is not covered by that tail bound.

Differentiating the quadratic mass against the weak RHS gives

$$
\boxed{\dot N_h=\frac{T_x}{2}u_W(d-u_W),\qquad
d=A_0(R)R w_+.}
$$

This is the physical mass flux evaluated with `u_W,v_hat,f_W`. Completing
the square proves `dot N_h <= T_x d^2/8`; homogeneous incoming data give
`dot N_h=-T_x u_W^2/2`. These are exact semidiscrete statements for this
fixed-radius reference. The last-element derivative need not equal `v_hat`.
Both that discrepancy and its converted incoming-characteristic residual
are returned, separately from the mass/weak-flux balance residual.

The [26 controls](../tests/test_nsc_atlas_cap_geometric_reference.py) include
an exact time-linear wave on two prescribed meshes, independent dense
element matrices, direct integration of the original radial mass density,
symbolic centre limits, the completed-square bound, small-radius arithmetic,
quiet data, input immutability and rejected data/budget cases. Development
uses **10,600 charged coefficient/manufactured values / 0.53125 active CPU
seconds** in 97 test executions. No saved field, radial ODE or physical time
evolution is evaluated. The production cap is unchanged. Moving radius,
finite backreaction and the background/departure decomposition remain
outside this reference component.

Grok's source review confirms the element indexing, gap remainder, mass
gradient, regular limits and weak wall flux. Its small-radius boundary
observation is addressed by evaluating `v_trace-A0 chi_W/R` through the
stable wall stiffness action. Both the incoming residual and derivative
discrepancy now retain their small nonzero limits without subtracting two
large inverse-radius terms. The kinked-field control also distinguishes the
true centre slope from the wall secant.

#### Finite-backreaction fixed-cap variational reference

The [variational component](../src/recursive_horizons/nsc_atlas_cap_variational_reference.py)
extends this reference to nonnegative alpha with the same smooth potential
and radial constraints. All radii and fields use the existing dimensionless
normalization; `N=n(R)` is the normalized cap mass, whose contribution to
the physical mass function is `alpha N`. Radius and alpha are held fixed
through this calculation.

Keep continuous linear `chi,k` in `z=atanh(r)`, with centre values zero, and
define `k=A0 p`. At nonzero alpha `k` is not the scalar time derivative.
The exact primitive supplies

$$
r=\tanh z,\quad A_0=\operatorname{sech}^2z,\qquad
f=\chi/r,\quad p=k/A_0,\quad q=\chi_z/A_0,\quad h=q-f.
$$

For one affine element `chi=a_e z+b_e`, the stable identity
`h=a_e z gamma(z)/A0-b_e/r`, with `gamma=1/z-A0/r`, reuses the positive
gap above. On the first element `b_e=0`; the centre limits are
`f=chi_z`, `Pi=k_z`, `Phi=0`. The radial mass, lapse and mass-sensitivity
RHS values all vanish at the centre. Each element retains its own
one-sided derivative at a shared node.

The inherited centre clock is
`delta(0)=-log cosh(sqrt(alpha/6) f(0))`, not zero when alpha and the
central scalar are nonzero. The mass and its nodal sensitivities have zero
central data; the lapse starts from this field-dependent value. The mass
functional is independent of that integration constant, but the field rates
and physical flux use its value through `exp(delta_W)`.

For a free coefficient `theta` of either field, the exact sensitivity obeys

$$
\begin{aligned}
n_z&=A_0\rho,\qquad \delta_z=A_0a,\\
(G_\theta)_z&=A_0(-aG_\theta+s_\theta),\qquad G_\theta(0)=0,\\
s_\theta&=\frac{r^2U}{2}f_\theta
+\frac A2\{p p_\theta+h(q_\theta-f_\theta)\}.
\end{aligned}
$$

Here `a=alpha(p^2+h^2)/(2r)`, `rho=r^2P/2+A(p^2+h^2)/4`,
`A=A0-2alpha n/r`, and `U=P_f`. The `-aG` term differentiates the mass
inside the metric. A second route computes

$$
G_\theta(R)=\int_0^{z_W}A_0(z)
 e^{\delta(z)-\delta_W}s_\theta(z)dz.
$$

This attenuation-weighted quadrature reads the forward mass and lapse but
does not read the forward sensitivity components. It is an independent
gradient calculation on that reference metric, not an independent metric
solution. The component stops before constructing a weak RHS if their
absolute discrepancy exceeds the admitted `2e-8` threshold. A deliberate
corruption of the forward sensitivities is rejected by this check.

An all-radius metric bound is available before integration. For exact
declared nodal data let `L_chi,L_k` bound the absolute element slopes,
`Z=z_W`, and set

$$
F=L_\chi Z/R,\quad p_{\max}=L_k Z/A_0(R),\quad
h_{\max}=L_\chi/A_0(R)+F,\quad
b_{\max}=\frac{R^2P(F)}2+\frac{p_{\max}^2+h_{\max}^2}{4}.
$$

The potential is even and increasing with `|f|`, and `z/r` increases with
radius, so these bound the declared fields. Writing `n_r+a n=b` gives
`a>=0`, `0<=b<=b_max`, and hence
`0<=n(r)<=integral_0^r b(s)ds<=r b_max`. Therefore

$$
\boxed{A(r)\ge A_0(R)-2\alpha b_{\max}\quad(0\le r\le R).}
$$

Arb evaluates this sufficient bound at 160 bits from exact binary nodal
inputs. The result includes an exact dyadic lower endpoint and a downward
rounded display value. An unresolved or nonpositive bound stops the
reference; failure of this sufficient test does not disprove the geometry.
The certificate concerns the exact radial solution for the fixed declared
fields. It is not an enclosure of the numerical radial integration or a
time-dependent validity interval.

Let `g_chi,g_k` be the mass gradients, `M` the same free-node element mass
matrix and `c=exp(delta_W)`. Define

$$
V=2cM^{-1}g_k,\qquad
\chi_T=V,\qquad
k_T=M^{-1}\{-2cg_\chi+e(B_WRw_+-V_W)\}.
$$

Symmetry of `M` cancels the interior pairings, leaving

$$
\boxed{N_x=\frac{T_x}{2c}V_W(B_WRw_+-V_W)
\le\frac{T_x}{8c}(B_WRw_+)^2.}
$$

This is the physical **weak** wall flux with
`p_hat=V_W/B_W`, `q_hat=f_W+R w_plus-p_hat`. The weak form follows by
integrating the continuum mass-gradient equations: the reaction tail
contributes the volume reaction, and the wall term contributes the same
incoming condition. The factor `c` depends on the state; this skew-gradient
identity alone is not a claim about a canonical Poisson bracket.

The component returns both weak traces and the reconstructed `p,q`, with
their discrepancies evaluated without avoidable `q-f` cancellation. In the
single nonconstant manufactured case at `alpha=.002`, the certified lower
metric value is `0.7837849351318801`. The gradient routes differ by
`3.75699e-13`, and the mass/weak-flux algebra by `7.80626e-18`. Nevertheless,
the wall momentum discrepancy is `-2.06194e-7` and the normal-variable
discrepancy is `-0.320591157262`. Those are properties of the prescribed
manufactured field and incoming data, not measurements of the retained
physical cap. They show why a small mass-balance residual cannot replace
the separate boundary-consistency check.

All [19 controls](../tests/test_nsc_atlas_cap_variational_reference.py) pass.
They include recovery of the exact zero-alpha element functional, a mass
direction checked from perturbed inputs, two fixed ODE tolerances, independent
gradient quadrature, the inherited centre clock and an independent lapse
integral, weak flux conversion, quiet data, exact binary bound
handling, rejected domains/budgets and the corrupted-gradient negative
control. Cumulative development uses **1,704,332 coefficient/derivative
values / 3.859375 active CPU seconds** in 71 executions. Requested tolerances,
observed gradient agreement and paired-run changes are not uniform error
enclosures. This supplies one regular nonlinear weak RHS on a fixed cap;
moving-radius evolution, strong trace compatibility and a coupled physical
cap-wall-exterior trajectory remain open.

Independent Grok review confirms the sensitivity equation, weak form and
all-radius bound. It caught an initial zero-centre-lapse initialization that
the mass checks alone could not detect. Restoring the original clock and
checking the wall lapse and rate normalization independently corrects that
prototype. The cumulative control receipt includes those earlier executions;
no physical evolution was performed with the incorrect clock.

#### Exact wall-momentum projection bound

The [projection-bound component](../src/recursive_horizons/nsc_atlas_cap_projection_bound.py)
bounds one of the preceding trace differences without integrating a radial
ODE. Inputs define exact binary nodal values and the same continuous P1
fields. The result concerns the exact constrained radial solution and exact
mass gradient. Numerical integration, quadrature and linear-solve errors
must be added separately.

Write `beta=B/A0`, `g=beta k=Bp`. The momentum mass gradient gives
`M V=load(g)`, so `V` is the centre-zero L2 projection of `g`. If `I_h g`
is its nodal interpolant, the free-node residual satisfies

$$
M(V-I_hg)=\operatorname{load}(g-I_hg).
$$

For every free hat function let `ell_i=integral phi_i dz`. The mass matrix
has `M_ii-sum_(j!=i)|M_ij| >= ell_i/3`; the eliminated centre makes the
first row stricter. At a maximum component of `V-I_hg`, this gives the
factor-three infinity-norm bound. Piecewise interpolation then gives

$$
|V_W-B_Wp_W|\le\frac38\max_e h_e^2\sup_e|g_{zz}|,
\qquad
|\widehat p_W-p_W|\le
\frac{3\max_e h_e^2\sup_e|g_{zz}|}{8\eta e^{\delta_0}},
$$

where `eta>0` is the preceding metric certificate and
`B_W>=eta exp(delta0)`. Element derivatives are one-sided at internal nodes;
interpolation uses each open element, not derivatives across the jumps.
The centre limit is regular.

The radial constraints simplify the required derivatives:

$$
\begin{aligned}
\beta_z&=\alpha e^\delta C,
& C&=2n/r^2-rP-4n/A_0,\\
\beta_{zz}&=\alpha A_0e^\delta(aC+C_r),\\
C_r&=2\rho/r^2-4n/r^3-P-rU\Phi-4\rho/A_0-8rn/A_0^2.
\end{aligned}
$$

Since `k` is affine on each element,
`g_zz=beta_zz k+2 beta_z k_z`. Also
`|Pi|<=L_k Z/(R A0(R))`. For affine `chi=a_e z+b_e`,

$$
\Phi=a_e K(r)-b_e/r^2,\qquad
K(r)=\frac{r/A_0-\operatorname{atanh}r}{r^2}
=r\int_0^1\frac{2t^2}{(1-r^2t^2)^2}\,dt.
$$

This kernel is positive and increasing on `0<r<1`. Each element has
`|Phi|<=|a_e| K(r_hi)+|b_e|/r_lo^2`; omit the second term on the first
element, where the intercept is exactly zero. Take the smaller of this
elementwise maximum and the earlier global bound. This sharpens the
estimate without changing the mesh or fields.

Let `Q=Pi_max^2+Phi_max^2`, `C_b=P(F)/2+Q/4`, `nu=C_b/3`.
The independent positive-metric certificate implies
`0<=rho<=C_b r^2`, hence `0<=n/r^3<=nu`. With `A_R=A0(R)`, define

$$
\begin{aligned}
C_*&=R(2\nu+P(F)+4\nu R^2/A_R),\\
C_{r,*}&=2C_b+4\nu+P(F)+RU(F)\Phi_{\max}
             +4R^2C_b/A_R+8\nu R^4/A_R^2,\\
a_*&=\alpha RQ/2,\qquad \delta_* =\delta_0+\alpha R^2Q/4,\\
\beta_1&=\alpha e^{\delta_*}C_*,\qquad
\beta_2=\alpha e^{\delta_*}(a_*C_*+C_{r,*}),\\
g_2&=\beta_2 L_k Z+2\beta_1L_k.
\end{aligned}
$$

These bound `|C|`, `|C_r|`, `a`, `delta`, the two beta derivatives and
`|g_zz|`, respectively. Arb at 160 bits supplies exact dyadic outward
endpoints and directed floating-point display values. The bound vanishes
exactly at alpha zero or zero momentum. It fails explicitly if the metric
floor or finite bounds cannot be resolved.

For the same manufactured state the momentum projection bound is
`0.0005850596031365239`; the previous global Phi estimate gave
`0.021594667452231934`. Both exceed the observed momentum difference
`2.06194e-7`. This is a conservative exact-model bound, not an enclosure of
that numerical measurement or of the independent normal-trace discrepancy.
Thirteen [controls](../tests/test_nsc_atlas_cap_projection_bound.py) pass.
Grok independently checked the projection constant, constrained derivative
identities, centre limit and global norm estimate. Its source review also
confirmed the elementwise sharpening and directed endpoints. No physical trajectory,
new raw dataset or observational prediction is generated.

#### Moving-mesh projection remainder

The fixed-cap identity does not automatically extend to a moving wall.
Consider its exact zero-alpha reference with `z_i=xi_i Z(T)`,
`nu=Z_T/Z`, `w(z)=nu z` and centre-zero P1 coefficients `chi,k`. Then
`M=Z M_0`, `K=K_0/Z`, `b=A0(R)/R=2/sinh(2Z)` and `S=K-b ee^T`.
Define the coordinate-advection matrix and its projected primitive rate by

$$
C_{ij}=\int_0^Z\varphi_i w\,\partial_z\varphi_j\,dz,
\qquad a=M^{-1}C\chi.
$$

A natural proposed extension of the fixed-cap equations is

$$
\dot\chi=k+a,\qquad
M\dot k=-S\chi+e(d-k_W)+Ck,\qquad d=A_0(R)Rw_+.
$$

Here dots use the cap time `T`; this is a candidate numerical rule only.
Integration by parts and differentiation of the geometric wall coefficient
give

$$
C+C^T=Z_Tee^T-\nu M,\quad
\dot M=\nu M,\quad \dot K=-\nu K,\quad
\dot b=-Z_T(2A_0(R)+b^2).
$$

Let `s_W=chi_z|W` be the last reconstructed slope. Direct differentiation
of `N=(k^T M k+chi^T S chi)/4` gives the exact identity

$$
\begin{aligned}
\dot N&=J_{\rm weak}+R_T\rho_W+\mathcal E_{\rm move},\\
J_{\rm weak}&=\tfrac12 k_W(d-k_W),\\
R_T\rho_W&=Z_T\left\{\tfrac12 A_0(R)\chi_W^2
              +\tfrac14[k_W^2+(s_W-b\chi_W)^2]\right\},\\
\mathcal E_{\rm move}
 &=\tfrac12\chi^TKa-\tfrac\nu4\chi^TK\chi
   -\tfrac{Z_T}4s_W^2
   +\tfrac12 b\chi_W(Z_Ts_W-a_W).
\end{aligned}
$$

The kinetic mesh terms cancel to the correct boundary density. The primitive
projection leaves the displayed remainder. It vanishes for a globally
affine primitive because `a=nu chi`; it need not vanish for a kinked P1
primitive. The fixed-cap flux is still the weak flux: a disagreement between
the reconstructed normal derivative and the incoming condition is a separate
term, not part of this identity.

An exact counterexample uses nodes `(0,1/2,1)`, primitive `(0,gamma,0)`,
zero momentum and `d=-2 gamma`. The incoming condition then agrees with the
actual boundary derivative at this instant. Direct element integration gives

$$
M=\begin{pmatrix}1/3&1/12\\1/12&1/6\end{pmatrix},\quad
K=\begin{pmatrix}4&-2\\-2&2\end{pmatrix},\quad
C=Z_T\begin{pmatrix}-1/6&1/3\\-5/12&5/12\end{pmatrix}.
$$

Consequently `a=Z_T gamma (1/7,-18/7)`,
`N_T=13 Z_T gamma^2/7`, while `J_weak=0` and
`R_T rho_W=Z_T gamma^2`. The remainder is therefore
`6 Z_T gamma^2/7`. For the timelike expanding speed `Z_T=1/2`, it is
`3 gamma^2/7`; contraction reverses its sign and a static wall removes it.
The construction permits arbitrarily small nonzero gamma.

Seven [exact controls](../tests/test_nsc_atlas_cap_moving_projection_identity.py)
check the element identity, geometric derivative, full mass derivative,
affine limit and the counterexample at positive, negative and zero speed.
This rejects a direct projected-advection extension as an exactly
mass-compatible moving-cap scheme. The next construction must derive its
mesh terms and physical boundary trace together, or retain and bound their
remainder. Subtracting this defect to enforce the desired mass balance would
not establish the field equations. No production evolution is changed.

Independent Grok derivation confirms the matrix shape identities, kinetic
cancellation, full remainder and hand-computed `6/7` coefficient per unit
wall speed. The timelike examples above use half that speed.

#### Dilation-compatible regular field space

The preceding defect comes from projecting a discontinuous slope transport
back into a continuous P1 space. A different regular representation removes
that particular obstruction. Let `xi=z/Z` and use global odd polynomials

$$
\chi=\sum_j c_j\xi^{p_j},\qquad k=\sum_j d_j\xi^{p_j},
\qquad p_j=1,3,\ldots,2m-1.
$$

This changes the finite-dimensional representation, not the action, fields,
potential or wall law. Odd powers give even regular `f=chi/r` and
`Pi=k/(A0 r)` at the centre, with `f_C=c_1/Z`, `Pi_C=d_1/Z` and
`Phi_C=0`. Write `E=diag(p_j)` and `e=(1,...,1)^T`. Dilation stays in
the same space: `z chi_z=sum_j (E c)_j xi^p_j`, and likewise for `k`.

Consider the exact constrained mass `N(c,d,Z)=n(tanh Z)` at fixed alpha
on a positive-metric cap. Changing only the coordinate scale through
`c_T=(Z_T/Z)E c`, `d_T=(Z_T/Z)E d` leaves both physical functions
unchanged at fixed `z`, hence at fixed `r`. It also preserves `c_1/Z` and
the inherited centre lapse. The radial solution on each common interior
prefix is therefore unchanged; only the outer integration limit moves.
The exact shape identity is

$$
\boxed{N_Z+\frac{g_c^TEc+g_d^TEd}{Z}=A_0(R)\rho_W.}
$$

Here the coefficient gradients and `N_Z` are taken in this polynomial
representation. The density is the original radial mass density evaluated
on the reconstructed fields. No new source term or mass correction is used.
The argument requires differentiability of the constrained radial solution
in its declared positive-metric domain.

The mass matrix is `M_ij=Z/(p_i+p_j+1)`, and the wall evaluation vector
is `e`, not the last coordinate vector used for nodal P1 coefficients.
Retain the inherited clock explicitly:

$$
\begin{aligned}
V&=2e^{\delta_W}M^{-1}g_d,\qquad D_W=B_WRw_+,\\
c_T&=V+(Z_T/Z)Ec,\\
d_T&=M^{-1}\{-2e^{\delta_W}g_c+e(D_W-e^TV)\}+(Z_T/Z)Ed.
\end{aligned}
$$

The common gradient cancels the fixed-domain interior pairings. The exact
shape identity supplies the remaining domain derivative, giving

$$
\boxed{N_T=\frac{(e^TV)(D_W-e^TV)}{2e^{\delta_W}}+R_T\rho_W,
\qquad R_T=A_0(R)Z_T.}
$$

This is a moving-domain weak mass identity for exact gradients. It is not
a Poisson-bracket assertion. State dependence of the clock factor is
retained; no derivative of that factor is required in this first derivative
of the clock-independent mass functional.

At alpha zero the independent quadratic check uses
`K_ij=p_i p_j/[Z(p_i+p_j-1)]` and
`N=(d^T M d+c^T(K-b ee^T)c)/4`. Direct integration gives

$$
ME+EM=Zee^T-M,\qquad EK+KE=K+pp^T/Z,\qquad
C=(Z_T/Z)ME.
$$

These identities and `b_Z=-(2A0+b^2)` reproduce the original wall density
for a non-affine odd primitive. The seven executions of
[six exact controls](../tests/test_nsc_atlas_cap_dilation_identity.py) use
zero field samples and 1.296875 active CPU seconds. They check the regular
centre, fixed physical fields and clock, independently integrated matrices,
quadratic shape derivative, general gradient pairing and separate normal
traces. The first run's only failure was an unexpanded symbolic equality;
the expanded zero difference and focused replay pass.

The remaining wall variables remain distinct: the weak momentum is
`p_hat=(e^T V)/B_W`, while the reconstructed momentum is `k_W/A0(R)`.
The weak normal variable follows the imposed incoming characteristic;
the reconstructed one follows `chi_z`. The wall scalar rate combines the
weak velocity with the exact reconstructed coordinate transport. Equality
to the material-interface scalar rate still requires a trace-consistency
check. The previous P1 projection constant does not apply to this new
space. The numerical reference below supplies radial gradients and
conditioning checks. Field truncation, numerical error enclosures, trace
bounds and a coupled physical evolution remain to be established before
this reference replaces any production update.

Independent Grok review confirms the finite-alpha shape argument, inherited
clock, weak pairing and non-affine quadratic limit. It also identifies the
numerical requirements: the monomial Gram matrix is poorly conditioned as
degree grows, and the actual global field slopes need their own certificate.
An orthogonal basis may improve conditioning, but the dilation operator,
wall evaluation and mass gradient must all be transformed consistently.
Using the P1 nodal-slope certificate unchanged would not bound polynomial
derivatives between nodes.

#### Independent numerical shape derivative

The numerical polynomial reference keeps `xi=z/Z` fixed during radial
integration. This supplies a second route to the shape derivative instead
of setting it from the preceding identity. Define `J=Z A0`, so
`n_xi=J rho` and `delta_xi=J a`. With the coefficient vectors held fixed,

$$
\begin{aligned}
r_Z&=\xi A_0,& (A_0)_Z&=-2r\xi A_0,\\
f_Z&=-f r_Z/r,& p_Z&=2r\xi p,\\
q_Z&=(-1/Z+2r\xi)q,& h_Z&=q_Z-f_Z,\\
A_Z\big|_n&=(A_0)_Z+2\alpha n r_Z/r^2,\qquad
J_Z/J=1/Z-2r\xi.
\end{aligned}
$$

Differentiate the original radial density at fixed running mass:

$$
s_Z=rr_ZP+\tfrac12r^2U f_Z
 +\tfrac14(A_Z|_n)(p^2+h^2)+\tfrac A2(pp_Z+hh_Z).
$$

The forward mass-shape derivative starts at zero at the centre and obeys

$$
\boxed{G_\xi=J[-aG+s_Z]+J_Z\rho,\qquad G(0)=0.}
$$

The `-aG` term differentiates the metric's running mass; the last term
differentiates the radial Jacobian. At `xi=1` this gives the total `N_Z`;
there is no additional moving-endpoint term in this fixed-unit-interval
integration. The endpoint contribution appears separately when comparing
with the dilation identity. The mass does not depend on the centre lapse
constant, so no centre-clock sensitivity is added to this mass equation.
The inherited lapse itself still sets the weak rates.

Near the centre the direct difference `q_Z-f_Z` cancels two finite terms.
The equivalent regular form is

$$
h_Z=(-1/Z+2r\xi)h+\xi(2r-\gamma(z))f,\qquad
\gamma(z)=1/z-A_0/r.
$$

For each odd monomial, the existing positive-gap evaluator also gives
`h_j=xi^(p_j-1) [(p_j-1)+z gamma(z)]/(Z A0)`.
Whole-cap derivative bounds are
`L_chi=sum_j p_j |c_j|/Z` and `L_k=sum_j p_j |d_j|/Z`.
Applying the earlier metric inequality to these actual global bounds,
with exact binary coefficients and outward Arb arithmetic, supplies the
polynomial-domain certificate.

The [new adapter](../src/recursive_horizons/nsc_atlas_cap_polynomial_reference.py)
uses the same radial gradient engine as the P1 reference.
It checks the independently integrated shape against dilation, the forward
gradients against attenuation-weighted quadrature, and the Gram solves by
their conditioning and residual. Its monomial one-RHS domain is one to three
odd modes with a positive metric certificate and a prescribed timelike wall
speed. Requested ODE tolerances and these acceptance checks are not uniform
numerical error enclosures. The scalar tangent is compared with both the
weak characteristic pair and the actual hybrid-interface pair.

The hybrid interface combines the supplied `w_plus` with the reconstructed
`w_minus=(p_trace-h_trace)/R`. Set
`dp=p_trace-p_weak`, `dh=h_trace-h_weak`. Its transmitted pair is
`p_interface=p_weak+(dp-dh)/2`,
`h_interface=h_weak+(-dp+dh)/2`. Consequently

$$
f_x^{\rm cap}-f_x^{\rm interface}
=\frac{(R_x-T_xB)\,dp+(R_x+T_xB)\,dh}{2R}.
$$

The smaller expression `R_x dh/R` compares only with the weak pair; it is
reported separately and is not called material transmission. Likewise the
actual interface mass-rate discrepancy is the weak balance defect plus
`T_x (J_weak-J_interface)+R_x (rho_trace-rho_interface)`, with
`J_interface=A B p_interface h_interface/2` and
`rho_trace-rho_interface=A(p_trace^2+h_trace^2-p_interface^2-h_interface^2)/4`.
The potential cancels in this last difference because the field is common.

On the admitted three-mode manufactured state at `alpha=.002`, `Z=.5`,
the mass is `0.00023306432649816832` and the exact-model metric lower bound
is `0.7861214361263889`. The independently integrated shape differs from
the dilation identity by `1.3746870395486832e-15`; forward and adjoint
gradients differ by `2.917457941897794e-14`. The mass/weak-flux-plus-motion
residual is `7.05544563744942e-17`. The Gram condition number is about
`849.576`; a separate 75-digit solve checks its effect on the rates.

The reconstructed-minus-weak momentum discrepancy is `-4.18419e-10`,
the normal discrepancy is `-0.08767909754045838`, and the actual interface
scalar-tangent discrepancy is `-0.05601999753599185`. The interface mass-rate
discrepancy is `-0.0008803230560492166`, decomposed into the displayed trace
contributions. The weak-pair scalar difference is separately `-0.00758934`.
These prescribed fields and
incoming data are a manufactured calculation, not the retained physical
cap. In particular, the tiny mass-balance residual does not erase the
boundary discrepancies or supply material transmission.

All [26 monomial-reference controls](../tests/test_nsc_atlas_cap_polynomial_reference.py)
pass in 54 development executions with **1,252,570 charged values /
2.84375 active CPU seconds**. They
include the independent radius/coefficient direction, quadratic limit,
centre and lapse integral, global metric bound, high-precision Gram solve,
timelike motion of both signs, quiet data and rejection of jointly corrupted
forward/adjoint shape values before the weak solve. The one-mode finite-alpha
case agrees with the existing P1 reference. The common-engine extraction also
passes all 19 P1 controls with unchanged reported values and all 13 P1
projection controls. Independent Grok derivation confirms the shape ODE and
the use of global polynomial slope bounds. Completed Grok source review of
`a36f16a` confirms the common engine and two-trace hybrid comparison. Its
additional checks now exercise the actual `T_x B` timelike threshold, both
Gram solves at high precision, and rejection of a corrupted solve. The
output names `Z_rate=dZ/dx` explicitly and presents both gradients with the
same component names. No numerical mass result changed in that follow-up.
The interface check adds an exact two-trace scalar identity and checks the
mass-rate discrepancy against its independent flux/density decomposition.
A separate local-jet construction calls the existing `hybrid_wall_boundary()`
and recovers the reported incoming/outgoing pair, radius/clock rates and
transmitted scalar tangent. It is a local interface check, not an exterior
source history or a solved coupled evolution.

#### Orthogonal regular cap coordinates

The [basis component](../src/recursive_horizons/nsc_atlas_cap_polynomial_basis.py)
uses the complete odd Legendre tower `P_1(xi),P_3(xi),...,P_(2m-1)(xi)`,
normalized by `P_p(1)=1`. It spans the same regular polynomials as the
monomial reference. No physical field or boundary law changes. Orthogonality
on `[-1,1]` and evenness of the products give

$$
M_{ij}=\frac{Z\delta_{ij}}{2p_i+1},\qquad
K_{ij}=\frac{\min(p_i,p_j)[\min(p_i,p_j)+1]}{2Z}.
$$

The wall value vector is still `e=ones`. Define wall and centre derivatives
`t_j=P'_pj(1)=p_j(p_j+1)/2`, `s_j=P'_pj(0)`. The coefficient dilation
operator is upper triangular:

$$
E_{ij}=\begin{cases}p_i&i=j,\\2p_i+1&i<j,\\0&i>j.\end{cases}
$$

Indeed `xi P'_p=p P_p+sum_(lower odd j)(2j+1)P_j`. In this basis

$$
E^TM+ME=Zee^T-M,\qquad E^TK+KE=K+tt^T/Z,
\qquad e^TE=t^T,\quad s^TE=s^T.
$$

The transposes matter: `E` is not diagonal or symmetric. The centre field
and momentum are `(s dot c)/Z` and `(s dot d)/Z`, respectively. Consequently
pure coordinate dilation preserves the inherited centre clock, while the
same shape and weak-flux identities remain valid. A skipped mode would
break closure under dilation; the API uses the complete tower.

For stable evaluation at the centre let `t=xi^2`, `U_j=P_(2j)(xi)` and
`Q_j=P_(2j+1)(xi)/xi`. Starting with `U_0=Q_0=1`, use

$$
U_j=\frac{(4j-1)tQ_{j-1}-(2j-1)U_{j-1}}{2j},\qquad
Q_j=\frac{(4j+1)U_j-2jQ_{j-1}}{2j+1}.
$$

Differentiate this recurrence in `t`, retaining the product term
`(t Q)_t=Q+t Q_t`. Then `P=xi Q`, `P'_xi=Q+2t Q_t`, and

$$
f_j=\frac{zQ_j}{Zr},\quad p_j^{\rm basis}=\frac{\xi Q_j}{A_0},\quad
h_j=\frac{2t(Q_j)_t+z\gamma(z)Q_j}{ZA_0}.
$$

This avoids subtracting the finite centre values of `q` and `f`.
The global slope certificate follows directly from the Legendre equation:

$$
\frac{d}{dx}\{(1-x^2)(P'_n)^2+n(n+1)P_n^2\}=2x(P'_n)^2\ge0.
$$

At `x=1` the expression is `n(n+1)`, so `|P_n|<=1` on `[0,1]` for
`n>=1`; `P_0=1` separately. The derivative expansion in lower opposite-parity
Legendre polynomials gives `|P'_n|<=n(n+1)/2`. Thus the correct bound here
is `L_chi=sum_j |c_j| p_j(p_j+1)/(2Z)` and likewise for `k`. The same
radial metric inequality applies to these global bounds, with outward
arithmetic. This is not a field-truncation or wall-projection error bound.

`polynomial_cap_reference(..., basis="legendre")` now supports one to
sixteen odd modes. The default remains the one-to-three-mode monomial
reference. The Legendre Gram condition is exactly `(4m-1)/3`, hence 21 at
sixteen modes. The radial engine, independent shape equation, incoming
condition and explicit hybrid trace discrepancies are shared by both bases.
Both forward and adjoint outputs name their `chi`, `canonical_density` and
`Z` components; `Z_rate` is the geometric coordinate rate, while
`radius_rate` remains the areal-radius rate.

Twelve [basis controls](../tests/test_nsc_atlas_cap_polynomial_basis.py) pass
in thirteen executions, with **34,140 charged values / 1.390625 active CPU
seconds**. Exact integrations check the matrices, transpose identities,
centre and wall functionals. Five fixed-point evaluations of sixteen modes
agree with 85-digit arithmetic to a maximum raw difference of `1.95399e-14`.
The normal basis is also checked near the centre using its uncancelled
high-precision expression.

Five [adapter controls](../tests/test_nsc_atlas_cap_legendre_reference.py)
pass with **1,835,038 charged values / 0.6875 active CPU seconds**. The same
three-mode field, represented in the two bases, gives matching mass,
gradients, rates and wall discrepancies after the explicit change of basis.
One fixed sixteen-mode manufactured profile recovers the quadratic mass at
alpha zero. At `alpha=.002`, its mass is `1.6825374175116823e-7`, gradient
discrepancy `1.66228e-13`, shape discrepancy `4.02076e-19` and weak moving
mass residual `-1.28312e-20`. These are operator controls; no prepared pulse
or coupled evolution has been calculated in this representation.

Independent Grok derivation confirms the complete-tower condition,
transpose conventions, centre clock, stable recurrence and derivative bound.
Completed independent Grok source review found no correctness defect in
the extension, including dilation orientation, centre weights, normal-field
cancellation and the global slope certificate. The reviewer did not execute
tests or re-read the shared radial solver body; the controls above supply
execution evidence.

#### Fixed-profile representation diagnosis

The [projection component](../src/recursive_horizons/nsc_atlas_cap_profile_projection.py)
now tests the existing computational profile itself. On a fixed cap slice,

\[
 \chi_0=F(T+Z\xi)-F(T-Z\xi),\qquad
 k_0=F'(T+Z\xi)-F'(T-Z\xi),
\]

and the coefficients of each odd Legendre mode are

\[
 c_j=(2p_j+1)\int_0^1\chi_0 P_{p_j}\,d\xi,\qquad
 d_j=(2p_j+1)\int_0^1k_0 P_{p_j}\,d\xi.
\]

Split at every mapped profile knot, including the quiet-extension edge.
The coefficient integrands have degree at most 36, so 19-point Gauss is
exact in exact arithmetic on each piece. The squared residuals have degree
at most 62, covered by 32-point Gauss. Direct residual integration avoids
subtracting nearly equal energies. The physical norms are

\[
 \|e\|_{L^2(dz)}^2=Z\int_0^1e^2d\xi,\qquad
 \|e_z\|_{L^2(dz)}^2=Z^{-1}\int_0^1e_\xi^2d\xi.
\]

Report keys `chi_xi` and `k_xi` identify derivative integrands; their returned
norms include this `1/Z` factor and are norms of the **physical z derivative**.
The H1 residual is not an energy difference: derivative projection is not
orthogonal in this pairing. The wall comparison retains the value term,

\[
 \Delta w_+=\frac{\Delta k+\Delta\chi_z}{A_0R}
                 -\frac{\Delta\chi}{R^2}.
\]

The last term cannot be dropped when the projection changes the wall value.
Centre values are `2 F'(T)` and `2 F''(T)`, compared with the projected
`s.c/Z` and `s.d/Z`. The certified quiet extension below the first profile
label is admissible; upper extrapolation remains excluded.

Eight [independent controls](../tests/test_nsc_atlas_cap_profile_projection.py)
pass with **234,742 charged values / 0.90625 CPU seconds**. They include
quintic reproduction and exact piecewise integrals of a truncated cubic
profile, its derivative residuals, quiet limits and domain/budget rejection.
Grok confirms the quadrature orders, reflection derivatives and physical
norm factors. Completed independent Grok source review found no correctness
blocker in reflection, projection, split quadrature or trace normalization.
It did not execute the controls or re-read the basis implementation. Its
requested nonzero sixteen-mode coefficient control now also passes; the
control receipt includes that focused addition.

The [one fixed input diagnosis](../scripts/check_nsc_atlas_cap_profile_projection.py)
authenticates the existing operator-comparison payload and 104 bindings,
then uses its **zero-case terminal cut 128**, `x=0.4815042082576364`,
`T=0.19784721605423244`, `R=0.55`, `Z=0.6183813135744636`.
All 65 existing profile intervals/quiet pieces are retained; no source is
refitted. The input payload SHA256 is
`312bc770f7baf63b2aff3d02c26642b32db212febdd34a7e3e6f40fded539d09`;
the selected profile's canonical-array SHA256 is
`04e10695afac696c4b0e2c178ccf17a19a04ce20706b60bdf86f6f44af38e0fe`.

| Quantity | Residual norm | Relative to original norm |
|---|---:|---:|
| Primitive chi | 0.01494870769471052 | 0.06791943818548013 |
| Canonical density k | 0.8700432643748897 | 0.1460260575723819 |
| Primitive derivative | 2.435828693493262 | 0.408823877612646 |
| Density derivative | 151.03382952730925 | 0.5645312602617536 |

The original centre is still quiet: `f=Pi=0`. The projected centre gives
`f=-0.6742457405504162`, `Pi=-38.79888310735645`. The wall errors are
`Delta f=-0.13703417265604292`, `Delta Pi=11.347534097957219`,
`Delta Phi=-161.02998990160847`, and
`Delta w_plus=-149.68245580365124` against original
`w_plus=0.03192911946666193`. These are projection errors relative to the
specified interpolant, not a newly evolved physical signal. The interpolant's
wall value also differs slightly from the saved finite-volume field;
this diagnosis does not identify the two.

The single diagnosis uses **2,328,407 charged values / 2.296875 CPU seconds**,
with one projection, one profile construction and zero radial solves,
source preparations, time steps or trajectory replays. The command prints
its compact result and input/producer hashes; no raw payload is added.

**Outcome:** the sixteen-mode global representation is inadequate for this
narrow background pulse, despite its good matrix condition and mass identities.
Do not promote it to physical evolution or mask the trace error with a
boundary adjustment. The next connection is an enriched representation that
retains the existing characteristic background directly and represents its
departure. Its explicit time and moving-domain dependence must enter the
same mass-gradient calculation; simply subtracting the background from the
current polynomial RHS does not establish that balance. This is a specific
numerical representation limitation, not a change to the action or source.

#### Mass balance with a retained characteristic background

Retaining the narrow pulse directly requires an affine field space,

\[
 \chi=\chi_0(T,z)+\phi(z/Z)^T c,\qquad
 k=k_0(T,z)+\phi(z/Z)^T d,
 \qquad k_0=\partial_T\chi_0.
\]

The same constrained radial mass now has explicit dependence
`N(T,Z,c,d)`. Let `g_c`, `g_d`, `N_T` and `N_Z` be its actual partial
derivatives. At fixed physical `z`, the background does not depend on `Z`.
The coefficient dilation therefore still cancels the fixed-z change of the
**departure**, giving

\[
 N_Z+\frac{g_c^TEc+g_d^TEd}{Z}=A_0(R)\rho_W.
\]

The time dependence adds a separate connection. Put

\[
 a=\partial_T\chi_0,\quad b=\partial_T k_0,\quad
 a_h=M^{-1}\int_0^Z\phi a\,dz,\quad
 b_h=M^{-1}\int_0^Z\phi b\,dz,
 \quad \mu=e^{\delta_W},\quad C=B_WRw_{\rm in}.
\]

A direct background subtraction from the polynomial weak equations gives
fixed-Z departure rates

\[
 v=2\mu M^{-1}g_d-a_h,\quad V_W=a_W+e^Tv,
 \qquad
 w=M^{-1}\{-2\mu g_c+e(C-V_W)\}-b_h.
\]

Adding coordinate motion gives `c_x=T_x v+Z_x E c/Z` and
`d_x=T_x w+Z_x E d/Z`. Differentiating the actual mass, rather than
assuming that subtraction cancels its time dependence, yields

\[
 N_x=T_x\left\{\frac{V_W(C-V_W)}{2\mu}+\mathcal R\right\}
       +R_x\rho_W,
\]

\[
 \boxed{\mathcal R=N_T-g_c^Ta_h-g_d^Tb_h
       -\frac{(a_W-e^Ta_h)(C-V_W)}{2\mu}.}
\]

This is an exact finite-dimensional contraction of the stated rates and
mass derivatives. It supplies a remainder to evaluate or control; it is
not a conservative evolved solver. The flux in this identity uses the weak
pair `(V_W/B_W, R w_in-V_W/B_W)`. Comparison to the actual hybrid trace pair
still requires the two-trace correction derived above.

At alpha zero, `mu=1` and

\[
 N=\frac14\left\{\int_0^Z(k^2+\chi_z^2)dz
                         -\beta\chi_W^2\right\},\qquad
 \beta=A_0(R)/R.
\]

For the reflected background, `b=chi_0,zz`. Orthogonality of the L2
projection and closure of the complete odd polynomial space under the
second derivative give `v=d`, `V_W=k_W`, and

\[
 N_T-g_c^Ta_h-g_d^Tb_h
       =\tfrac12(a_W-e^Ta_h)(\chi_{z,W}-\beta\chi_W).
\]

For example, set `eta=a-Pa`, `theta=b-Pb`. Integrating
`chi_z eta_z+k theta` by parts leaves the boundary term: the remaining
integral vanishes because the departure second derivative and momentum are
in the trial space and `P` is self-adjoint. Consequently

\[
 \mathcal R=\tfrac12\Delta a\,A_0(R)\,dh,\qquad
 \Delta a=a_W-e^Ta_h,\quad
 dh=\frac{k_W+\chi_{z,W}-\beta\chi_W-C}{A_0(R)}.
\]

A compatible incoming trace gives zero remainder in this reference limit.
An arbitrary imposed incoming value does not. This distinction matters
when using the background to carry a pulse whose global polynomial
projection has a large endpoint error.

An independently integrated counterexample uses `F(u)=u^4`, one departure
basis function `z/Z`, and `T=c=d=C=0`. Then `chi_0=0`, `k_0=a=8z^3`,
`b=0`, `a_h=24Z^3/5`, `g_d=4Z^4/5` and `N_T=g_c=0`.
The weak departure rates are `v=0`, `w=-24Z^2`. The direct mass rate is
`-96Z^6/5`, while the weak wall flux is `-32Z^6`, leaving the exact
remainder **`64Z^6/5`**. Choosing the compatible `C=8Z^3` instead gives
zero departure rates and zero remainder. This intentionally incompatible
instantaneous boundary input is an algebraic counterexample, not an
evolved physical failure.

Three [exact controls](../tests/test_nsc_atlas_cap_affine_background_identity.py)
pass, using **6,000 charged algebra values / 0.21875 CPU seconds**, with
zero field samples. They differentiate the polynomial background's
integrated mass, check the general incoming remainder and independently
check the moving shape identity. Independent Grok derivation confirms the
full mass split, moving shape term, geometric-limit remainder and both
incoming choices in the counterexample. A finite-alpha enrichment still needs the
explicit `N_T` and a controlled remainder before a physical time step.

#### Retained-background field jets and metric domain

The [retained-data adapter](../src/recursive_horizons/nsc_atlas_cap_background_jets.py)
keeps the original Hermite values and first two derivatives alongside the
existing numerical `CharacteristicProfile`. The exact C2 model means the
Hermite interpolant of those original binary values. Rounded `PPoly`
coefficients can leave tiny jumps between stored pieces; their conversion
and evaluation errors remain numerical errors. They are not additional
surface sources and are not silently part of the exact-C2 certificate.

Let `I1`, `I2` and `J2` be the existing knot-split integrals of `F'`, `F''`
and `t F''` over `-1<=t<=1`, evaluated at `T+zt`. Let `J0`, `J1` integrate
`F'''` and `t F'''`. With `gamma=1/z-A0/r`, the regular background fields are

\[
 f_0=\frac{z I_1}{r},\quad p_0=\frac{z I_2}{A_0},\quad
 h_0=\frac{z(\gamma I_1+J_2)}{A_0}.
\]

Their time derivatives follow from the same moments:

\[
 f_{0,T}=\frac{z I_2}{r},\quad p_{0,T}=\frac{z J_0}{A_0},\quad
 h_{0,T}=\frac{z(\gamma I_2+J_1)}{A_0}.
\]

The spatial derivatives can avoid a second centre cancellation,

\[
 f_{0,z}=A_0h_0/r,\quad
 p_{0,z}=(I_2+zJ_1)/A_0+2rp_0,\quad
 h_{0,z}=zJ_0/A_0+2r(f_0+h_0)-f_{0,z}.
\]

At the axis, `f_0=2F'(T)`, `f_0,T=p_0,z=2F''(T)` and the other reduced
normal/time limits used here vanish. A merely C2 profile at an axis knot
does not imply an even differentiable expansion of `Pi`; no such claim is
needed for these finite reduced jets. For `z>.01`, the adapter uses the
existing endpoint values and derivatives. Independent Grok derivation and
source review confirm the formulas. Numerical coefficient/evaluation error,
tiny stored knot jumps and the existing gamma-series arithmetic remain
outside an interval error enclosure.

Fourteen jet controls pass in fifteen executions, including the first
roundoff-sensitive assertion corrected from exact equality to a declared
`2e-14` tolerance. Total development uses **9,430 charged values / 0.34375
CPU seconds**. Controls compare against 85-digit polynomial derivatives,
exercise an interior C2 knot on both sides of the evaluation switch and at
a tiny axis span, and check quiet limits and copied original inputs.

For a local metric certificate, quintic Bernstein controls are formed in
outward Arb arithmetic from the original Hermite values. On a label interval
of width `h`, the first three controls are

\[
 b_0=F_0,\quad b_1=F_0+hF'_0/5,\quad
 b_2=F_0+2hF'_0/5+h^2F''_0/20;
\]

the other three use the corresponding right-end expressions. First and
second derivative controls are scaled successive differences. Their convex
hulls bound the exact C2 profile and jets on every intersected input piece,
including the certified quiet extension. Domain endpoints are checked in
exact-binary arithmetic; an inward-rounded slice must be declared by the
caller if its nominal endpoint rounded outside the supplied profile.

On each existing cap interval, bound

\[
 b(r)=r^2P(f)/2+A_0(p^2+h^2)/4\le b_i.
\]

The positive radial solution obeys `n_r+a n=b`, `a>=0`, `n(0)=0`, hence
`0<=n(r)<=integral b`. With `U_i=sum_(j<i) b_j Delta r_j`, the affine
cumulative envelope gives

\[
 \sup_{[r_i,r_{i+1}]}\frac{n}{r}\le
 \max\left(\frac{U_i}{r_i},
           \frac{U_i+b_i(r_{i+1}-r_i)}{r_{i+1}}\right).
\]

For the first interval use `n/r<=b_0`. The comparison ratio is monotone,
so its endpoint values suffice; there is no new sampling grid. The local
metric floor is `A0(r_(i+1))-2 alpha sup(n/r)`. Every floor must be positive
before radial integration. Failure of this sufficient bound does not prove
a physical zero of `A`. Numerical interpolation and ODE errors are separate.
Independent Grok derivation confirms this recurrence. Seven exact and
manufactured certificate controls pass with **7,213 charged values /
0.203125 CPU seconds**. Independent Grok source review is complete and
confirms the outward certificate for the exact C2 input model.
The implementation uses the original Bernstein controls, not rounded
PPoly coefficients; it leaves the existing profile owner unchanged.

#### Explicit-time radial sensitivity

At fixed radius and departure coefficients the retained background supplies
`f_T,p_T,h_T`. Linearizing the same radial mass equation gives

\[
 (n_T)_r=-a n_T+\frac{r^2P_f}{2}f_T
                      +\frac A2(pp_T+hh_T),\qquad n_T(0)=0.
\]

There is no lapse derivative in this source. The lapse still enters the
wall clock and Hamiltonian factor. At fixed `xi`, add `xi` times each
background z derivative to the departure shape partials. Existing knot
splits and the small-z evaluation switch are retained in the shared radial
solver. Exact C2 matching cancels moving-knot boundary terms; the actual
numerical interpolation error is not thereby certified away.

The existing [polynomial adapter](../src/recursive_horizons/nsc_atlas_cap_polynomial_reference.py)
now has the optional `background` and `background_time` arguments for
Legendre departures. The new sensitivity is a further column in the same
forward/adjoint engine. Background time jets use 19-point split Gauss for
their projection. `balance_defect` remains the actual mass-minus-weak-flux
residual. `explicit_time_remainder` records the derived additional term;
`remainder_identity_defect` checks their agreement. No correction has been
installed to force either physical discrepancy to zero.

The 26 monomial controls still pass after this change, using another
578,573 charged values / 1.375 CPU seconds within their original allowance.
Their cumulative record is now **1,831,143 / 4.21875**, eighty executions.
The separately admitted unchanged sixteen-mode regression passes with
**837,921 / 0.25**. Seven retained-background reference controls pass with **378,959 charged
values / 1.0 CPU seconds**. They compare the same full polynomial field in
both representations, including time and radius chain rules, and reproduce
the exact quartic remainder and compatible limit. Quiet and underflow
branches are checked separately. Completed Grok source review found no
blocker in the time source, shape, indexing, projection subtraction or
residual reporting; its scope excluded the separate metric certificate.
The small-z split is shared as one constant. No physical evolution has
been run.

#### Fixed retained-pulse mass-remainder diagnosis

Both reviewed commands have now run once at source commit `0b8d0af`.
They authenticate 104 existing input bindings and retain the zero-case
terminal cut 128, `T=0.19784721605423244`, nominal `R=0.55` and
`Z=0.6183813135744636`. The exact binary profile endpoint requires one
inward ulp, giving `Z=0.6183813135744635`; its floating radius is still
`0.55`. Sixteen departure coefficients in each field are zero. The incoming
trace is reconstructed from this same background:
`w_plus=0.03192911946666193`, with `T_x=1` and `R_x=0`.
This is a fixed background reference, not the evolved finite-alpha state.

The standalone metric certificate passes on the 65 existing intervals:
**`A >= 0.6843674647287019`**, downward rounded from the Arb lower point
`1000204170236996926071558525507717990167727106073 * 2^-160`.
The outward mass upper bound is `27.578152340162468`. This use consumes
**78,388 charged values / 2.15625 CPU seconds** of its 100,000 / 10 allowance,
with no radial solve. The bound concerns the exact C2 Hermite model;
numerical field and ODE errors still require separate control.

The two instantaneous radial references give:

| Quantity | Geometric limit, alpha = 0 | Recorded alpha = 0.00012983127142551443 |
| --- | ---: | ---: |
| Mass N | 17.749674162944892 | 17.633757131127027 |
| Wall lapse delta | 0 | 0.013422275586898612 |
| Explicit time derivative N_T | 3.2756321491130675e-6 | -0.11724584063509795 |
| Weak wall flux | 3.275828645292163e-6 | -8.288200105688344e-6 |
| Mass-minus-weak-flux discrepancy | -1.3992838745002594e-10 | -0.02337491766289285 |
| Derived explicit-time remainder | -1.3992495356829338e-10 | -0.02337491766289628 |
| Difference from remainder identity | -3.4338817325643006e-15 | 3.4312830354821244e-15 |
| Maximum absolute departure rate | 5.495621735462919e-10 | 0.5179424833694384 |
| Forward/adjoint discrepancy | 5.014653083783901e-10 | 4.453792540933671e-10 |
| Independent shape discrepancy | 4.736468941395674e-10 | 4.1815267664191327e-10 |
| Actual hybrid interface mass discrepancy | -1.3994665943215884e-10 | -0.023386447174947225 |
| Wall normal discrepancy | -4.708303191769403e-12 | -0.002584016790169539 |
| Wall scalar-tangent discrepancy | -5.970983256442519e-12 | -0.0032816426664835363 |

The zero-coupling reference recovers the retained pulse and its flux within
the declared numerical acceptance thresholds. At finite coupling, the
explicit-time remainder accounts for the weak mass discrepancy to
`3.44e-15`. This identifies a compatibility failure of the stated affine
background subtraction. Passing the radial gradient and shape comparisons
does not remove it. The separately reported hybrid trace gaps also remain.
The next construction must derive a compatible projection or variational
coupling from the field equations; subtracting the measured discrepancy
from an output would not supply that construction.

The two references consume **15,264,426 charged values / 7.984375 CPU
seconds** of 20,000,000 / 120, with zero time steps or trajectory replays.
All asserted checks pass. The internal coefficient-evaluation counts are
7,620,464 and 7,643,936, with 6,787 and 6,811 field sites respectively.
These are numerical consistency results, not uniform error enclosures or
a finite evolving cap-wall solution. No observational prediction is added.

The compact commands print their receipts and create no raw files. Their
inputs and producers are bound by SHA256:

| Input or producer | SHA256 |
| --- | --- |
| Existing operator payload | `312bc770f7baf63b2aff3d02c26642b32db212febdd34a7e3e6f40fded539d09` |
| Input manifest | `9bdf673c83cc0992f9484a3c0f014337e65cdeb0470d953bc0204028676ec1cd` |
| Original profile arrays | `04e10695afac696c4b0e2c178ccf17a19a04ce20706b60bdf86f6f44af38e0fe` |
| Domain command | `0e62b62704a07fecdc5cd21f56977eb1d1cce348e7046ed1703a2529ee712a5b` |
| Two-reference command | `b4500e686114c905dd9360224807d02c9e1881bdb4534b042dd7f5116699e7d8` |
| Background bound | `cd3013306825187e2c7b9336a56990183eb6dbc1383fad2146f0bb4fe49f3948` |
| Background jets | `f50ec70168f482dc8efaac8cde4a185650170390c26237fb1421fa1c94f39852` |
| Polynomial reference | `35d7273d3ffda6ecb441b0280c07dbd794992950ba609674bd2be1e280b3fb0e` |
| Shared radial engine | `0779ff9783932b2335b9dc09bb64a218daad6b694393d82cfc97f62ce8888f93` |

#### Tangent enrichment and its moving-knot obstruction

The retained-background remainder has a precise fixed-slice interpretation.
Let `P` be the L2 projector onto the common field space, `eta=a-Pa` and
`theta=b-Pb`. For admissible field variations,

\[
 N_T-g_c^Ta_h-g_d^Tb_h=D\mathcal N(\eta,\theta).
\]

If both background tangents are exactly in the trial space, this term and
`a_W-(Pa)_W` vanish. This holds for the nonlinear constrained mass, without
adding an artificial flux. New basis functions have their actual wall
values, not the unit wall values of odd Legendre modes. Affine departure
rates still subtract the represented background rates; total field rates
include them. This establishes a frozen-slice identity only. The independent
review confirms this identity and the need to account for basis motion.

A stronger obstruction occurs for merely C2 input profiles if `b` is also
used as a primitive-field direction. Consider the exact manufactured profile

\[
 F(u)=(u-\kappa)_+^3,\qquad 0<T<\kappa,\quad
 \ell=\kappa-T\in(0,Z),\qquad x=z-\ell.
\]

On the positive cap, `chi0=x_+^3`, `a=k0=3x_+^2` and `b=6x_+`.
Set `chi=chi0+q b`, `k=k0`, with constant nonzero `q`. The primitive is
continuous, but its z derivative jumps by `6q` at the moving knot. At fixed
wall, put `L=Z-ell` and `beta=A0(R)/R`. The actual geometric-limit cap mass is

\[
 N=\frac{9}{10}L^5+3qL^3+9q^2L
             -\frac\beta4(L^3+6qL)^2.
\]

Ordinary pointwise time differentiation away from the knot, including the
wall term, misses **`9 q^2`**. Reynolds differentiation of the piecewise
integral supplies `-ell_T [rho_z]=[rho_z]=9q^2`, where
`rho_z=(k^2+chi_z^2)/4` and brackets mean right minus left. The unmodified
background `q=0` has no such jump. This obstruction is introduced by the
new trial direction; it does not invalidate the retained C2 background
reference or its original smooth time source.

Moreover, `b_T=6 H(z-ell)` is not an H1 variation. For the pure added field
`q b`, translating the knot by epsilon gives the exact spatial-gradient
quotient norm squared `36 q^2/epsilon`. It diverges as epsilon tends to
zero. Thus the expression `D N(psi_T c,psi_T d)` alone is not justified for
this moving basis. The ordinary field derivative must be supplemented by
the moving-interface contribution. With `Z=1`, `ell=1/2` and the exact
space `span{z,a,b}`, the squared L2 projection error of `b_T` is **`3/2`**.
This also directly rules out closure of that particular moving space.

For the nonlinear radial equation `n_r=rho`, continuity of `n` at a moving
internal radius `r_j(T)` gives the necessary sensitivity jump

\[
 [n_T]_j=-r_{j,T}[\rho]_j.
\]

It is zero for the original exact C2 background, whose radial source is
continuous. A new basis that makes `rho` discontinuous needs this jump in
addition to the existing bulk sensitivity. The weak evolution and its
internal fluxes must also be derived; installing this sensitivity jump
alone would not establish a conservative cap-wall evolution.

Three [exact controls](../tests/test_nsc_atlas_cap_background_enrichment.py)
pass across four executions, using **24,000 charged algebra values /
0.765625 CPU seconds**, with zero field samples or time steps. The first
run's difference-quotient assertion compared unsimplified expressions;
its corrected exact algebraic equality passes on the focused repeat.
The mass integration and projection controls passed on their first run.
Independent Grok review confirms the full cap expression, interface sign
and failure of the smooth-only pairing for this moving enrichment.
The next implementation choice must preserve the required regularity or
include derived internal-interface terms and error control. No enriched
physical solver is installed by this diagnosis.

#### Internal wavefront compatibility and evolution

A moving reconstruction kink is not automatically a solution of the scalar
equations. On a finite patch take finitely many isolated internal fronts,
one-sided regular field limits, continuous `f,n,delta`, a single smooth
potential, `A>0`, and the common centre clock. Require `Phi=f_r`, the scalar
equations on every open piece, and no singular scalar source on a front.
Collisions, axis arrival and material-wall intersection end this patch.
These are internal scalar fronts, not additional material walls.

For a front `r_j(T)` with speed `s`, distributional differentiation gives

\[
 E_\Pi=s[\Pi]+B[\Phi]=0,\qquad
 E_\Phi=s[\Phi]+B[\Pi]=0,\qquad B=Ae^\delta.
\]

Brackets mean right minus left. The second condition also follows from
continuity of `f` along the front and `f_T=B Pi`. A nonzero derivative jump
requires `s=+B` with `[Phi]=-[Pi]`, or `s=-B` with `[Phi]=[Pi]`.
The background knot speed `+/-A0` is generally different at finite coupling.
The original C2 background has no such derivative jump and remains valid
as a represented background. Adding a lower-regularity primitive direction
and leaving its kink at the background speed need not obey these equations.

With `P=Pdiff`, the actual density and mass flux are

\[
 \rho=\frac{r^2P}{2}+\frac{Ar^2}{4}(\Pi^2+\Phi^2),\qquad
 J=\frac{ABr^2}{2}\Pi\Phi.
\]

Continuity of `f` makes the potential contribution continuous. Direct
algebra, using arithmetic averages of the two traces, gives

\[
 \boxed{s[\rho]+[J]=\frac{Ar^2}{2}
       \{\operatorname{avg}(\Pi)E_\Pi+
          \operatorname{avg}(\Phi)E_\Phi\}.}
\]

Thus the two scalar jump equations cancel the Reynolds mass contribution.
The converse is false: the stationary traces `Pi_L=-u`, `Pi_R=u`,
`Phi_L=Phi_R=0` have zero density and flux jumps but `E_Phi=2 B u`.
Mass conservation alone does not specify a permissible field update.
In a label x use `v=r_j,x`, `T_x>0` and residuals
`v[Pi]+T_x B[Phi]`, `v[Phi]+T_x B[Pi]`; the left side is then
`v[rho]+T_x[J]`. This is precisely the interface term in the existing
weighted cap-transport identity.

For piecewise fields satisfying the bulk equations and these jump laws,
the existing product identity gives `(n_T-J)_r+a(n_T-J)=0` on open pieces
and `[n_T-J]=0` at fronts. Regular-centre data therefore imply `n_T=J`
throughout the patch. This is a conditional continuum conservation result;
it does not prove that a chosen finite field approximation meets its premises.

The interface amplitude also has an evolution law. Use the already defined
reduced fields `p=r Pi`, `q=chi_r=f+r Phi` and source
`S=B_r f+r exp(delta) U`. The radial constraints cancel their kinetic terms:

\[
 B_r=e^\delta\left(-2r+\frac{2\alpha n}{r^2}-\alpha rP\right).
\]

This derivative and `S` are continuous across the front under the stated
hypotheses. Set `w_sigma=p-sigma q`, where `sigma=+1` denotes the outgoing
family and `sigma=-1` the incoming family. On each side,

\[
 (\partial_T+\sigma B\partial_r)w_\sigma
       =-\sigma B_r w_\sigma-S.
\]

Subtract the traces along the same front. For `j=[w_sigma]`,

\[
 \boxed{r_{j,T}=\sigma B,\qquad j_T=-\sigma B_rj,\qquad
 [\Pi]=\frac{j}{2r_j},\quad [\Phi]=-\sigma\frac{j}{2r_j}.}
\]

No separate lapse or time normalization enters. In the exterior label,
both rates are multiplied by the inherited `T_x`. The common bulk solution
still has to supply `f,n,delta` at the moving front and maintain zero
jump in the complementary characteristic family. The amplitude law contains
no extra `B_T` term. It does not imply `j B` is constant when `B` depends
on time; the geometric closed form below uses the static geometric metric.

The [wavefront stage](../src/recursive_horizons/nsc_atlas_cap_wavefront.py)
implements these rates using the existing normalized potential and the
supplied constrained mass/lapse. It rejects an axis or wall event, a
nonpositive metric, invalid orientation/clock and nonfinite arithmetic.
Its pointwise metric check is not a domain certificate. Five independent
[stage controls](../tests/test_nsc_atlas_cap_wavefront.py) pass across six
executions with **404 charged values / 0.015625 CPU seconds**. The initial
five-control run used 340 values with CPU time below timer resolution;
a focused guard repeat used 64 / 0.015625 after enforcing scalar-real
orientation inputs, including complex and nonscalar rejection. They compare
the exact geometric characteristic
and amplitude, independently reconstruct `B_r` from both radial constraints,
and check orientation, clock scaling, zero jump and domain rejection.
Independent Grok mathematical review confirms the characteristic
orientation, amplitude law, constraint cancellation and geometric curves.
Independent source review finds no implementation blocker. Its suggested
`j B` invariant is not used: the mathematical review and equation above
require the explicit time-dependence qualification.

For the geometric limit, the independent exact curves are

\[
 r_j(T)=\tanh(z_0+\sigma t),\qquad
 j(T)=j_0\frac{\cosh^2(z_0+\sigma t)}{\cosh^2z_0},\quad t=T-T_0.
\]

There is also a complete manufactured finite-interval scalar control.
The previous ramp initial datum `chi=6q0(z-ell0)_+`, `k=0` evolves as

\[
 \chi(t,z)=3q_0\{(z-\ell_0-t)_++(z-\ell_0+t)_+\},\qquad
 k=\chi_t,\qquad 0<t<\min(\ell_0,Z-\ell_0).
\]

It splits into two compatible fronts, with speeds `+/-1` in z. On the three
open pieces the scalar equations hold; both jump equations hold at each
front. With `L=Z-ell0` and `beta=A0(R)/R`, direct piecewise integration gives

\[
 N=9q_0^2L-9\beta q_0^2L^2,
\]

independent of t. Its wall velocity field `k_W` is zero and the compatible
incoming drive is `C=6q0(1-beta L)`, so the wall flux vanishes. If the wall
moves while remaining outside the fronts, using `beta_Z=-2A0-beta^2`
gives `N_Z=A0(R)rho_W` and hence `N_x=R_x rho_W`. This specifies the
manufactured boundary input; it does not replace the retained source or
establish its material-wall transmission.

Five [exact controls](../tests/test_nsc_atlas_cap_wavefront_identity.py)
pass across six executions, using **36,000 charged algebra values /
0.390625 CPU seconds**, with zero field samples or numerical time steps.
One initial assertion compared factored and expanded forms structurally;
the corrected difference simplifies exactly to zero on the focused repeat.
The next coupled construction must evolve the bulk fields and supply their
metric traces to this interface law, with a controlled residual/error account.
A front RHS alone is not an evolved finite-alpha cap-wall configuration.

#### Fixed physical-coordinate bulk space

The route review identifies the restriction that made the preceding
construction difficult: the leading pulse was prescribed as a translating
affine offset, while the finite correction lived in a small global space.
Adding moving time tangents created extra regularity/interface work. The
original C2 pulse has no derivative jumps, and the source-free front law
preserves `j=0`; further front machinery alone does not repair its bulk
projection error.

Three paths were compared: moving tangent enrichment, a fixed local space
for the total fields, and a direct nonlinear characteristic integral solver.
The next bulk candidate is the **fixed physical-coordinate local space**.
It reuses the mass functional, radial metric, common clock and wall source.
The characteristic solver remains an alternative if that representation
fails. The internal-front result above remains a compatibility reference,
not an instruction to add artificial fronts to the smooth pulse.

Take basis functions `phi_i(z)` independent of time and wall radius, with
fixed physical knots. Represent total fields `chi=phi^T c`, `k=phi^T d`.
Restrict the first interval to the required axis parity/regularity subspace;
`phi_i(0)=0` alone is insufficient to preserve higher regular axis jets.
For the first patch, keep `Z=atanh(R)` strictly inside a single last cell
and retain a fixed active set with a positive, controlled Gram matrix.
No basis birth, knot crossing or extrapolated physical data is part of this
patch. A polynomial extension used to define basis support beyond the wall
is not an additional physical exterior state.

With `M(Z)=integral_0^Z phi phi^T dz`, `e=phi(Z)` and the same constrained
mass `N(Z,c,d)`, freezing coefficients freezes the fields at each physical z.
Consequently

\[
 N_Z=A_0(R)\rho_W,\qquad M_Z=ee^T.
\]

There is no explicit background-time derivative or basis dilation term.
For `mu=exp(delta_W)`, `C=B_W R w_in`, define

\[
 V=2\mu M^{-1}g_d,\quad V_W=e^TV,\qquad
 c_x=T_xV,\quad d_x=T_xM^{-1}\{-2\mu g_c+e(C-V_W)\}.
\]

The chain rule and symmetric Gram matrix give

\[
 \boxed{N_x=T_x\frac{V_W(C-V_W)}{2\mu}+R_x\rho_W.}
\]

These are field coefficient rates. If the unknowns were moments `M d`,
their derivative would instead include `M_Z d Z_x`. Inserting that moment
term into the coefficient rates changes the update. The exact quadratic
control shows the spurious mass contribution `-Z_x k_W^2/2` from adding
`-M^-1 M_Z d Z_x` to `d_x`.

At one initial T, the union of mapped original profile knots partitions
`chi0=F(T+z)-F(T-z)` into quintics and `k0=F'(T+z)-F'(T-z)` into quartics.
They are respectively C2 and C1, and thus belong to the common C1 quintic
local space. A cubic endpoint Hermite interpolant plus
`xi^2(1-xi)^2*(a+b xi)` represents every quintic exactly. This retains the
local initial profile without its inaccurate global sixteen-mode projection.
Exact inclusion in an axis-parity subspace additionally requires compatible
initial axis jets, for example a quiet axis or an axis away from a profile
knot. The merely C2 axis-knot case needs its separate regularity treatment.

Four [exact controls](../tests/test_nsc_atlas_cap_fixed_space_identity.py)
pass with **20,000 charged algebra values / 0.265625 CPU seconds**. They
independently integrate the cap mass for a deliberately non-dilation-closed
fixed basis, check its shape and weak flux, distinguish coefficients from
moments, and establish the local quintic representation. Independent Grok
review confirms the fixed-basis cancellation and initial C1 inclusion.

The [local basis](../src/recursive_horizons/nsc_atlas_cap_fixed_basis.py)
and [initial transfer](../src/recursive_horizons/nsc_atlas_cap_fixed_transfer.py)
now implement this representation. There are `4N-1` degrees of freedom per
field for N cells: an odd quintic at the axis, shared node values/scaled
slopes, and two bubbles on each subsequent cell. The zero-based axis cell
has no bubbles. Exact rational mapped knots and shapes define the space;
floating coefficient conversion has a separate value/first-derivative bound.
The final support endpoint extends the last polynomial algebraically, and
the physical wall initially cuts that cell in half. Source evaluation is
restricted to the physical cap. No projected smoothing of the initial
profile is performed.

Six-point Gauss integration forms the degree-ten Gram integrals, with
half-bandwidth five. Diagonal normalization, an eigenvalue condition check
and banded Cholesky define the numerical solve. The last-cell wall fraction
is restricted to `[1/4,3/4]`; numerical SPD/conditioning is checked for each
requested matrix, not certified over an evolving interval. Thirteen
[controls](../tests/test_nsc_atlas_cap_fixed_basis.py) pass, including exact
independent Gram integration, piecewise initial transfer, original-input
domain and axis-parity rejection, and a separate Bernstein/de Casteljau
comparison, now including the physical axis/wall and a banded solve residual.
Their combined receipt across fifteen executions is **197,250 charged values /
0.28125 active CPU seconds**. Numerical evaluation error is not included in the
coefficient-conversion bounds.

The [single retained-input check](../scripts/check_nsc_atlas_cap_fixed_transfer.py)
passes at the previously declared terminal slice `T=0.19784721605423244`,
`Z=0.6183813135744635` (one inward ulp from the nominal wall coordinate).
The mapped space has **65 cells / 259 coefficients per field**. Independent
Bernstein/de Casteljau evaluation checks 1,568 exact equalities at 392
locations; six rational interior points per cell establish **260 polynomial
identities**, with separate axis and physical-wall points. The axis is quiet.
At the wall `chi=-0.008774564904026498`, `k=chi_z=0.0005605323908308457`
and `k_z=0.00011529693710923827`.

| Channel | Largest sampled floating error | Coefficient-conversion bound in the exact basis |
|---|---:|---:|
| chi | `2.220446049250313e-16` | `3.2464582858306955e-17` |
| k | `7.105427357601002e-15` | `1.6585454964236477e-15` |
| chi_z | `7.016609515630989e-14` | `1.3035828315671948e-14` |
| k_z | `1.2505552149377763e-12` | `7.105849423871667e-13` |

The differing columns measure differing errors: floating shape/evaluation
roundoff is excluded from the coefficient-conversion bound. Sampled errors
meet the predetermined `5e-12 * max(1, reference_peak)` channel tolerances;
they are not uniform enclosures. The diagonal-scaled Gram condition is
`30121278.38782701`, giving `condition*eps=6.688267359461932e-9`, below the
`1e-8` limit. This relatively large condition requires continued solve/error
control; it does not certify a moving interval. One transfer and Gram factor
use **936,018 charged values / 2.484375 CPU seconds** of 3,000,000 / 30.
The compact [receipt](../results/summaries/nsc-atlas-cap-fixed-transfer-v1.json)
binds the input, original profile and producer sources. Independent Grok
review confirms the transfer signs, original-input domain and exact
comparison. There is no new raw dataset or numerical time evolution.

The next step after the retained initial transfer check is a numerical
mass/gradient adapter on the fixed cut cell. Its state is the total
coefficient pair `(c,d)` and the existing wall
and clock state; no translating affine field or new wavefront state is
required. Gram conditioning and the endpoint-event window must be checked.
The weak wall value `V_W` still needs comparison with the actual field trace
`B_W k_W/A0(R)`, and field residuals need an error estimate. This identity
removes two representation remainders; it does not yet establish an evolved
finite-alpha solution, its material-wall junctions, or its truncation error.

#### Metric domain for fixed local total fields

The [local certificate](../src/recursive_horizons/nsc_atlas_cap_fixed_bound.py)
interprets the supplied floating coefficient pair as exact binary numbers
in the rational fixed basis. It reconstructs the degree-five local powers
exactly. For `P(xi)=sum a_k xi^k`, the Bernstein controls on `[0,c]` are

\[
 b_j=\sum_{k=0}^j a_k c^k\frac{\binom jk}{\binom nk},\qquad n\ge\deg P.
\]

The maximum absolute control encloses the polynomial. Derivative controls
include the physical cell width. Only the final cell is clipped; no extra
partition or point scan enters the estimate. In the axis cell, chi and k
vanish, so bounds `L_chi,L_k` on their derivatives imply
`|f|<=L_chi*z_hi/tanh(z_hi)` and `|p|<=L_k*z_hi/A0(z_hi)`. Elsewhere the
local value hulls and positive left radius bound these quotients. In every
cell `|h|<=L_chi/A0(z_hi)+|f|`. The potential is bounded at the resulting
nonnegative upper bound for `|f|`.

Write the constrained radial equation as

\[
 n_r+a n=b,\quad
 a=\frac{\alpha(p^2+h^2)}{2r}\ge0,\quad
 b=\frac{r^2P(f)}2+\frac{A_0(p^2+h^2)}4\ge0.
\]

Starting with zero axis mass gives `0<=n<=U=integral b_bound dr`. On each
cell, the cumulative envelope is affine in r, so its ratio U/r takes its
largest value at an endpoint (or equals the source bound on the first
cell). Outward Arb arithmetic then gives `A>=A0-2*alpha*U/r` throughout
the physical interval. Failure to find a positive bound is inconclusive
about the actual metric sign.

Five [controls](../tests/test_nsc_atlas_cap_fixed_bound.py) pass with
**17,600 charged values**; active CPU time is below the reported timer
granularity. They independently expand the Bernstein polynomial, recover
the quiet geometric floor, compare the integrated alpha=0 linear-field
mass with the source envelope, check clipping and reject invalid input.
This certifies the exact declared field/constraint model only. Floating
evaluation, numerical radial error and evolved states remain outside it.
Independent Grok source review confirms the exact cut, hull, axis and
cumulative-envelope arguments. One retained application at the recorded
`alpha=0.00012983127142551443` certifies

\[
 A(r)\ge0.6870371535642128\quad(0\le r\le R),\qquad
 n(R)\le21.971840781627296,
 \quad R=\tanh Z\simeq0.55.
\]

Here Z is the stored exact-binary wall coordinate of the retained check.
This is a bound for the rounded coefficient pair in the exact rational
basis, not a reuse of the original-input certificate. One authentication,
regenerated transfer and certificate use **278,272 charged values /
2.296875 CPU seconds** of 1,000,000 / 15. The certificate itself charges
83,136 values at 160-bit Arb precision. There is no radial solve, Gram
factor, new partition or time step. The compact [receipt](../results/summaries/nsc-atlas-cap-fixed-domain-v1.json)
preserves the dyadic lower bound, source hashes and scope.

For the numerical mass adapter, the existing forward/adjoint radial engine
can retain the fixed local coefficient pair. On a cell write
`z=z_left+h*xi`; an element-index coordinate gives a Jacobian `h*A0`
to the radial equation. The last element ends at its physical cut fraction.
For each supported basis function the coefficient partials use
`f_i=phi_i/r`, `p_i=phi_i/A0`, and
`h_i=(phi_i,z-A0*phi_i/r)/A0`. Near the axis the stable equivalent is

\[
 h_i=\frac{\phi_{i,z}-\phi_i/z+
                  (1/z-2/\sinh(2z))\phi_i}{A_0}.
\]

In the odd axis cell `phi_i=a1*xi+a3*xi^3+a5*xi^5`, the first difference
is `(2*a3*xi^2+4*a5*xi^4)/h`, avoiding subtraction of equal centre limits.
At fixed running mass, the local derivative sources are

\[
 s_{c_i}=\frac{r^2P_f f_i}{2}+\frac{A h h_i}{2},\qquad
 s_{d_i}=\frac{A p p_i}{2}.
\]

These sources enter the existing `G_r=s-a*G` engine and its independently
weighted adjoint quadrature. Only six basis functions contribute on an
ordinary cell, although the global sensitivity vector has `2*(4N-1)`
entries. Retain that existing engine first; a sparse replacement needs its
own demonstrated necessity. The physical shape derivative remains the
boundary value `N_Z=A0*rho_W`, and the already derived weak rates use the
banded Gram factor. No background-time or dilation term is introduced.
The [numerical adapter](../src/recursive_horizons/nsc_atlas_cap_fixed_reference.py)
now implements that connection using the shared engine without changing it.
For the total field and each basis function it first forms the exact local
polynomial `Delta=z*chi_z-chi`, then rounds once for evaluation. The
element-index integration coordinate supplies `h*A0`; it does not move the
physical polynomial basis. The output records the rounded wall-coordinate
offset and the difference between local polynomial and basis-row trace
evaluation. The initial metric certificate excludes those numerical errors.

The first alpha=0 comparison recovered independently integrated mass and
gradients but missed `c_x=T_x*d` by `3.40552204e-9` after the conditioned
Gram solve. One tighter requested radial accuracy, `rtol=2e-13` and
`atol=2e-15`, passed the unchanged rate and wall assertions. These are now
the adapter defaults. The finite coefficient directional derivative and
wall-shape mass differences also pass their original checks.

Gradient differences must be propagated through the operator. For forward
minus adjoint differences `Delta g_c, Delta g_d`, the adapter computes

\[
 \Delta V=2\mu M^{-1}\Delta g_d,\qquad
 \Delta d_x=T_x M^{-1}\{-2\mu\Delta g_c-e(e^T\Delta V)\}.
\]

It reports coefficient infinity norms, physical L2 norms and the signed
wall-velocity difference. With `M=S L L^T S`, the stable numerical norm is
`||L^T S v||^2`; the factor has diagonal scaling and lower bandwidth five.
This is a norm of the represented field using the numerical Gram, not an
interval enclosure. In the manufactured finite case, the raw gradient gap
is `7.70217e-16`, but the canonical coefficient-rate gap is `7.44041e-8`
against coefficients reaching `93205.5378`. Its physical L2 gap is only
`1.01493e-12`, and its wall-velocity gap is `2.31570e-15`. Thus a bare
coefficient absolute threshold would depend on the chosen coordinates.
The dedicated norm control compares independently integrated energy and
requires physical L2/wall discrepancies below `2e-10`; it passes.

The [six manufactured controls](../tests/test_nsc_atlas_cap_fixed_reference.py)
cover these identities, stable axis limits, quiet fields, common centre
lapse, coefficient preservation and domain rejection. The recorded
development work, including diagnosed failures, ten control executions and
one standalone base-state diagnostic, is **5,226,742 charged values /
4.640625 active CPU seconds** under the original allowance plus its
specifically admitted physical-norm follow-up. Independent Grok source and
mathematical reviews confirm the coordinate map, delta algebra, rates and
stable norm. The first retained reference reaches the shared engine's
adjoint quadrature, then stops at **25,000,008 charged values / 7.515625
CPU seconds**. Its forward integrations are complete, but the independent
gradient check and weak RHS are not. The [work diagnosis](prior-art-and-open-claim.md#one-retained-fixed-space-mass-reference)
bounds the full remaining adjoint/operator tail below three million values
and specifies one completion replay at unchanged numerical settings.
Numerical error enclosures, scalar residuals
and finite wall coupling remain outstanding.

The diagnosed completion replay now passes on the original **65 cells /
259 coefficients per field**, with the same supplied background incoming
trace, `T_x=1`, `R_x=0` and `alpha=0.00012983127142551443`:

| Instantaneous quantity | Calculated value |
|---|---:|
| Rescaled constrained mass N | `17.633757131085453` |
| Inherited wall lapse delta | `0.01342227558686664` |
| Mass rate from the coefficient gradients | `3.2413119398946692e-6` |
| Weak wall mass flux | `3.241311948576424e-6` |
| Mass-balance discrepancy | `-8.681754825408308e-15` |
| Actual minus weak wall momentum | `3.358186125013152e-14` |
| Actual minus weak wall normal field | `-2.4560214972879635e-14` |
| Forward/adjoint chi-rate physical L2 gap | `3.524258464086733e-13` |
| Forward/adjoint k-rate physical L2 gap | `8.165981512954437e-10` |
| Forward/adjoint wall-velocity gap | `-2.351542830749427e-14` |

The raw gradient discrepancy is `1.61968e-12`; the coefficient-infinity
k-rate gap is `1.80997e-7`, illustrating why physical norms and wall traces
are checked separately. The strong/weak trace differences are measured,
not forced to zero by an acceptance condition. The exact-field metric
floor remains `0.6870371535642128`; the numerical sampled minimum is
`0.6891748614242271`. The physical endpoint evaluation rounds by
`-1.1102230246251565e-16` in z, recorded separately from the exact input Z.

The completed call has 17,741 field sites and uses **25,463,040 values /
7.734375 CPU seconds** of the admitted 30,000,000 / 120. Together with the
first work-limited attempt preserved in `e5d7cb5`, the retained calculation
uses **50,463,048 values / 15.25 CPU seconds**. The [compact receipt](../results/summaries/nsc-atlas-cap-fixed-reference-v1.json)
binds inputs and producers; its arrays are regenerated by the command rather
than stored as a new field history.

This is the first retained total-field reference where the mass law and
actual wall traces agree at the reported numerical level. It removes the
earlier translating-affine representation remainder at this slice. The
incoming trace and frozen wall rates are declared diagnostic inputs;
admissible material-wall initial data and common wall/exterior evolution
must still be constructed. A time interval, scalar residual estimate and
numerical error enclosure do not follow from this instantaneous result.

#### Connecting the fixed cap to a moving material wall

The fixed cap's canonical primitive is `chi=r*f`; the exterior invariant
field is `c=r*f*sinhc(q_f)`. Their common value is a junction condition.
Consider an instant where that value and the induced metric agree. Write
the actual reduced cap traces as `p=r*Pi`, `h=r*Phi`, and its mass-gradient
pair as

\[
 P=V_W/B,\qquad Q=R w_+-P,\qquad
 d_p=p-P,\quad d_h=h-Q,\quad B=Ae^\delta.
\]

The existing hybrid supplies `w_+` from the outgoing exterior data and
uses the cap's actual `w_-=(p-h)/R`. Its matched derivatives are therefore

\[
 p_m=P+\frac{d_p-d_h}{2},\qquad
 h_m=Q-\frac{d_p-d_h}{2}.
\]

With the hybrid's actual velocity and boost, `H=beta+v>0`,
`beta-v=A/H`, and `T_s=beta/B`, differentiating the cap trace along the
moving wall gives `f_s^cap=(beta*P+v*h)/R`. The hybrid instead transports
`f_s^hybrid=(beta*p_m+v*h_m)/R`. Thus

\[
 \boxed{f_s^{cap}-f_s^{hybrid}
       =\frac{-A d_p/H+H d_h}{2R}.}
\]

Since `du/df=sqrt(D)` for `u=f*sinhc(q_f)`, the cap-minus-exterior
invariant-field rate discrepancy at a common value is
`sqrt(D)*(-A*d_p/H+H*d_h)/2`. A prior value mismatch needs its own record;
this identity does not assume it away after evolution starts.

Subtract the same moving-domain density from the mass rate and divide by
`T_s`. The weak temporal law uses `A^2*exp(delta)*P*Q/2`, while the actual
field product uses `A^2*exp(delta)*p*h/2`. Their residual is

\[
 \boxed{E_N=-\frac{A^2e^\delta}{2}
                   (P d_h+Q d_p+d_p d_h).}
\]

Relative to the hybrid's matched product, put `a=(d_p-d_h)/2`; the residual
is `A^2*exp(delta)*(a*(P-Q)+a^2)/2`. These are distinct comparisons.
Here `N` is the rescaled mass. The metric identity
`1-2m/r=1-r^2-2alpha*N/r` gives `m=r^3/2+alpha*N`, so the physical temporal
Einstein residual is `E_m=alpha*E_N` at fixed radius. On a moving wall,
`m_s=(3R^2/2+alpha*rho)*v+alpha*T_s*N_T`; removing only `alpha*rho*v`
would retain a spurious vacuum contribution. At zero amplitude the
normalized scalar still propagates while its mass backreaction vanishes.
The cap's own weak balance can vanish while either physical discrepancy
is nonzero. Sending the weak outgoing combination `(P-Q)/R` instead makes
the hybrid pair `(P,Q)`, but the moving primitive still has discrepancy
`v*d_h/R`. That substitution alone does not close the common scalar value.

Six [exact controls](../tests/test_nsc_atlas_cap_hybrid_compatibility.py)
pass with **40,000 charged algebra values / 0.484375 active CPU seconds**.
They differentiate the moving trace, apply the nonlinear frame derivative,
separate actual and matched mass products, and exhibit a counterexample
with zero weak balance but nonzero scalar/mass residuals. This identifies
the quantities that a coupled stage must preserve or bound. It does not
replace the actual traces, add a physical source or prescribe wall motion.
Grok independently confirms the trace residuals. Its normalization question
is resolved by the metric identity above. Admissible initial junction data
are checked separately before any common update.

The [initial-event audit](../results/summaries/nsc-atlas-cap-initial-junction-v1.json)
finds that the retained terminal cap and exterior are not a shared event.
Its invariant scalar differs from the stored zero-case exterior by
`-1.339906988e-5`; the finite-case exterior also has a different radius and
canonical time. Combining these records would change the common data.
The quiet corner already present in the source/gauge cache does pass:
`x=.4,R=.45,T=0`, zero scalar and all quiet departures, with material rates
`R_x=1.084575413,T_x=1.840737097`. The largest algebraic junction residual
is `9.325873407e-15`. This one audit costs 2,180 charged values / 2.1875 CPU
seconds, with no radial or prefix solve and no time step.

The current coupled implementation obtains the exterior endpoint scalar
from the cap at every stage; it does not evolve an independent wall scalar.
The fixed-space connection should retain that common-value elimination and
report the derived tangent residual, together with the temporal mass
residual. Its finite-dimensional field evolution remains an approximation
until the residual and error over a nonzero interval are bounded.

The [instantaneous adapter](../src/recursive_horizons/nsc_atlas_cap_fixed_hybrid.py)
now makes this connection. The radial solve passes an immutable scalar
packet to one wall callback; the exterior callback receives its common
`r,c` and returns three outgoing first jets. The hybrid supplies `T_x,R_x`
and the incoming cap field. Independently prescribed versions of those
three inputs are rejected on the coupled route. The physical wall lapse
rate follows from the existing fixed-coordinate boost law. A second cheap
hybrid evaluation fills in the incoming lapse jet with that derived rate;
it does not change the scalar port or trigger another radial/exterior solve.

The returned diagnostics compare the actual coefficient tangent with the
hybrid tangent and its exact factored expression. They retain both actual
and matched temporal mass products. The implemented mass derivative is
used, so its small numerical weak-balance discrepancy is included rather
than silently set to zero. Physical mass rates include the moving vacuum
contribution. Exterior bulk constraints and all time-error estimates remain
the caller's work.

Four [focused controls](../tests/test_nsc_atlas_cap_fixed_hybrid.py) pass
with **375,431 charged values / 0.375 active CPU seconds**. They include
an independently integrated zero-coupling energy/rate oracle, callback/direct
agreement, quiet vacuum wall clocks, a nonzero manufactured state, a moving
field finite difference, immutable inputs and invalid-port rejection. No
retained nonzero state or time trajectory is evaluated by these controls.
The common stage below supplies the actual exterior prefix. A finite state
update and its residual bounds remain required. The original retained reference remains bound to its producer
at `8f6c73c`; extending the adapter does not relabel that earlier evidence.

#### Common exterior, wall and fixed-cap stage

The [common stage](../src/recursive_horizons/nsc_atlas_fixed_cap_stage.py)
joins the fixed cap to the existing `departure_slice`. Its state contains
three exterior departure arrays, four wall/clock departures `r,L,T,s`,
and two total-field coefficient vectors `chi,canonical_density`. The wall
scalar remains the image of the cap coefficients. The cap's actual floating
radius is used at the exterior endpoint; its conversion difference from the
evolved radius is recorded and guarded by four machine epsilons times the
radius. No mismatch is removed from a field derivative.

The exterior restores source and wall values, evaluates one transverse
departure derivative, and uses it in both ordered prefixes and coordinate
transport. The extracted orchestration agrees bit for bit with the frozen
stage's exterior outputs. Original field equations and historical producers
remain unchanged.

Write `Q_X` for the exterior fixed-Ubar longitudinal derivative, `P_X` for
its numerical incoming derivative, and `P_X^req` for the hybrid's required
derivative. For `X=r,L,c`, the hybrid kinematic defect remains

\[
 K_X=X_x^{hybrid}-Q_X-W_xP_X=-W_x(P_X-P_X^{req}).
\]

The reconstructed cap scalar has the separate tangent discrepancy `E_c`
derived above. The actual shared-endpoint update defect is `K_c+E_c`.
The stage reports these separately as `wall_kinematic_defect`, the
connection's `conformal_wall_flux`, and `shared_wall_update_defect`.
Numerical incoming jets are retained, and no independent wall scalar is
advanced. When the source supplies its own `r_xx`, its pure constraint is

\[
 E_{x,source}=r_{xx}-L_xr_x+
 \frac{\alpha r}{2D}\left(\frac{c_x-(c/r)r_x}{r}\right)^2,
 \qquad D=1+\alpha(c/r)^2/6.
\]

Missing curvature is reported as unavailable, rather than inferred from a
mixed equation. Four [controls](../tests/test_nsc_atlas_fixed_cap_stage.py)
pass, plus one focused repeat after output naming was aligned with the
existing stage: **250,248 charged values / 0.265625 active CPU seconds**
in total. They cover vacuum recovery, frozen-prefix equality, one shared
derivative, nonzero scalar/jet residuals, source normalization, immutability
and domain rejection. Grok source review finds no actionable defect. This
is a common rate map, not a finite evolution or pure-bulk-constraint bound.

#### Time-operator bound in the frozen linear cap limit

The old face-based Courant criterion does not control clipped quintic
coefficients. Start from the independently checked limit `alpha=0`, fixed
`Z>0`, `T_x=1`, `R_x=0`, and zero incoming field. Define

\[
 M=\int_0^Z\phi\phi^Tdz,\quad K=\int_0^Z\phi_z\phi_z^Tdz,
 \quad e=\phi(Z),\quad b=\frac{1-R^2}{R}=\frac2{\sinh(2Z)},
 \quad S=K-bee^T.
\]

Every represented field vanishes at the axis. Cauchy--Schwarz gives
`ee^T <= Z K`; also `bZ<1`, since `sinh(2Z)>2Z`. Hence
`S >= (1-bZ)K > 0`. For coefficient vector `u=(c,d)`,

\[
 N=\tfrac14u^THu,\quad H=\operatorname{diag}(S,M),\quad
 u_x=Ju,\quad
 J=\begin{pmatrix}0&I\\-M^{-1}S&-M^{-1}ee^T\end{pmatrix}.
\]
\[
 J^TH+HJ=\operatorname{diag}(0,-2ee^T),\qquad
 N_x=-\tfrac12(e^Td)^2\le0.
\]

The exact constant-operator flow is contractive in the H norm without
requiring J to be normal. Its induced norm is bounded by its weighted
Frobenius norm:

\[
 \|J\|_H^2\le\operatorname{tr}(H^{-1}J^THJ)
 =2\operatorname{tr}(M^{-1}S)+(e^TM^{-1}e)^2\le L_*^2,
\]
\[
 L_*^2=2\operatorname{tr}(M^{-1}K)+(e^TM^{-1}e)^2.
\]

For rational knots and cut, these defining matrix entries are rational.
This is a dynamical bound; the scaled Gram condition alone is not one.
For the cubic Runge--Kutta polynomial, the exponential series gives

\[
 \|P_3(hJ)-e^{hJ}\|_H\le
 \delta=\frac{(hL_*)^4e^{hL_*}}{24},\qquad
 \|P_3(hJ)\|_H\le1+\delta.
\]

Telescoping n identical steps bounds the flow discrepancy by
`(1+delta)^n-1 <= n*delta*(1+delta)^(n-1)`. This does not assert that the
polynomial itself is contractive. For `hL_*<1`, replacing the exponential
by `1/(1-hL_*)` supplies a rational upper bound for an exact control.

Four [exact matrix controls](../tests/test_nsc_atlas_cap_time_operator.py)
pass using **60,000 charged algebra values / 6.5 active CPU seconds**.
Independent monomials span the same one-cell odd-quintic space. At `Z=1/2`
the exact upper bound is `L_*^2=3444`. The controls check positivity,
dissipation, the weighted norm, the cubic polynomial, and a positive-definite
rational energy-amplification bound for one manufactured Robin value and
one predetermined step. They run no numerical trajectory.

This estimate is restricted to constant geometry and zero coupling. A moving
wall changes M and H; finite coupling and exterior feedback change the rate
map and source terms. Those effects and the spatial residuals must be bounded
before this result can support a finite physical interval.

#### Moving-cut energy and continuous time residual

The zero-coupling estimate extends to a prescribed moving cut inside the
same last cell, without changing basis functions or coefficient meanings.
Put `l=phi_z(Z)`, `v=l-b*e`, `A0=sech^2 Z`. Since `b_Z=-b^2-2A0`,

\[
 M_Z=ee^T,\qquad S_Z=vv^T+2A_0ee^T,\qquad
 H_Z=\operatorname{diag}(S_Z,M_Z)\succeq0.
\]

For supplied `T_x>0` and incoming drive `C=A0*R*w_plus`, the same coefficient
equations give the exact moving-domain identity

\[
 N_x=\frac{T_x}{2}(VC-V^2)+\frac{Z_x}{4}u^TH_Zu,
 \qquad V=e^Td.
\]

The second term is `Z_x*A0*rho_W=R_x*rho_W`. It accounts for the moving
domain. Adding `-M^-1*M_x*d` to the coefficient rate would introduce the
extra work `-Z_x*(e^T*d)^2/2`. In fact,

\[
 \frac14u^TH_Zu=
 \frac14\{(v^Tc)^2+(e^Td)^2\}+\frac{A_0}{2}(e^Tc)^2.
\]

Define the nonnegative geometry factor

\[
 \gamma(Z)=\operatorname{tr}(H^{-1}H_Z)
 =v^TS^{-1}v+2A_0e^TS^{-1}e+e^TM^{-1}e,
 \qquad g(x)=\max(Z_x,0)\gamma(Z).
\]

Positive semidefiniteness and completing the flux square imply
`N_x <= g*N + T_x*C^2/8`. Writing `G(x)=integral_(x0)^x g(t)dt`,

\[
 N(x)\le e^{G(x)}\left[N(x_0)+
 \int_{x_0}^{x}e^{-G(t)}\frac{T_x(t) C(t)^2}{8}\,dt\right].
\]

This is conditional on the stated geometry and input bounds. For the same
prescribed `Z(x),T_x(x),C(x)`, let a continuous coefficient approximation
have residual `r=u_h'-F(x,u_h)` and error `epsilon=u_h-u`. Dissipation of
the homogeneous generator gives, without assuming normality,

\[
 \|\epsilon(x)\|_{H(x)}\le
 e^{G(x)/2}\|\epsilon(x_0)\|_{H(x_0)}+
 \int_{x_0}^{x}e^{[G(x)-G(t)]/2}\|r(t)\|_{H(t)}\,dt.
\]

Thus an independently bounded continuous time defect can control a cap
coefficient error. The estimate is not obtained from a saved endpoint
difference alone. It does not bound spatial truncation or replace the
exterior feedback with prescribed physical data.

One additional exact control verifies the matrix derivatives, moving energy
identity and single domain contribution using **16,000 algebra values /
1.484375 active CPU seconds**. All five time-operator controls together use
**76,000 / 7.984375** of the admitted 80,000 / 20. The common finite-amplitude
problem still requires source/exterior closure and a nonlinear error estimate.
Independent Grok review confirms the displayed moving energy and residual
estimates. The error estimate uses the dissipative pairing directly; it
does not need a separate Lipschitz or generator-norm bound for this linear,
prescribed-geometry problem.

### Fixed-coefficient time step and continuous reconstruction

The [time adapter](../src/recursive_horizons/nsc_atlas_fixed_cap_time.py)
uses the existing SSPRK(3,3) coefficients and exact-time/combination helpers.
Stages run at left, right predictor, and midpoint. A fourth evaluation at
the accepted right state supplies its own slope. Historical `Ip/Iq`
drivers remain unchanged. There is no retry, hidden source change, or cap
Courant assertion inferred from the exterior speed.

For accepted endpoints and slopes, put `h=x_R-x_L`, `t=(x-x_L)/h`,
`d=y_R-y_L`. The continuous path and its residual are

\[
y_h=y_L+htF_L+t^2a+t^3b,\qquad
a=3d-h(2F_L+F_R),\quad b=-2d+h(F_L+F_R),
\]
\[
y_h'=F_L+(t/h)(2a+3tb),\qquad
\mathcal R(x)=y_h'(x)-F(x,y_h(x)).
\]

`fixed_time_residual` makes an independent stage call on this path. It
returns each component in its own derivative units, with physical stage
diagnostics. A selected-point residual is not a uniform bound; the preceding
linear estimate requires control throughout the interval with identical
prescribed inputs. The coupled finite-amplitude estimate remains open.

The derived wall scalar needs a separate derivative. With `a=sqrt(alpha/6)`,
`u(f)=sinh(a*f)/a` and `chi_W=R*f`, its moving-trace chain rule is

\[
c_W'=\cosh(af)\chi_W'+[u-f\cosh(af)]R',\qquad
\chi_W'=e^T c_{\rm coeff}'+\phi_z(Z)^Tc_{\rm coeff}\,R'/(1-R^2).
\]

Zero-amplitude limits are continuous. This distinguishes the dense path's
derivative from the stage's actual scalar derivative and the hybrid target.

Four [controls](../tests/test_nsc_atlas_fixed_cap_time.py) pass, plus one
focused extension: **531,005 charged values / 0.53125 active CPU seconds**
over five executions. Known cubic and exponential solutions check time
integration and the continuous derivative. Failure/mutation controls preserve
the accepted input and count attempted calls. One manufactured cubic-source,
finite-amplitude step produces nonzero cap and exterior states and a finite
midpoint residual. It checks `K=-W_x J`, `E_c=R sqrt(D) E_f`, physical temporal
mass factors, and the derived scalar tangent against an independent centred
difference. This source is not the retained physical preparation. Grok review
finds no adapter defect; its request for derived-scalar and physical-residual
coverage motivated the focused extension.

### Retained source at the first scalar front

The [quiet-front audit](../scripts/check_nsc_atlas_quiet_source_front.py)
authenticates the retained source and reconstructs its original **129**
Hermite knots from `exact_rows[::2]`. It never refits the 1025-row evaluation
cache or integrates another source. Its [compact receipt](../results/summaries/nsc-atlas-quiet-source-front-v1.json)
records exact binary-input polynomial coefficients and three source queries.

The first 35 radial-correction polynomials are identically zero. Interval
35, from `x=0.42228630694544744` to `0.42292305857246026`, straddles the
scalar front `x=0.42269179059454887`. Its nonzero Hermite polynomial starts
before that front. At the front, `c=c_x=0` exactly, but
`h=1.433684912015276e-14`, `h_x=-8.464169193503768e-13`, and the actual
source constraint is `-6.90103555838494e-11`. A source radius rounded to
its vacuum value does not erase this curvature residual.

At the last pre-front knot, `h=h_x=0` and the reference wall has
`R=0.4748482346624701`, `T=0.04377128443765652`. At the scalar front its
reference values would be `R=0.4753138162777697`, `T=0.044623011279717725`.
Vacuum transport to the scalar front therefore does not reproduce the full
retained source. The earlier knot supplies a common value/first-jet candidate
for a numerical start, but its right-sided source curvature already has
residual `6.846851463133793e-11`. The source interpolation defect belongs
in the subsequent error account. This is not a constraint-compatible
continuum evolution. Grok review confirms the polynomial, source residual,
and distinction between the source label and wall radius. Transverse source
splines were not audited here; they are not inputs to `departure_slice`.

Completion uses **13,260 values / 2.34375 CPU seconds**. The first attempt
stopped at an unsuitable result-relative polynomial comparison; its recorded
upper-bound cost is 12,423 values / 20 CPU seconds (actual CPU was not printed).
The corrected comparison uses the absolute coefficient scale under
cancellation. Source coefficients and queries are unchanged. Neither
attempt integrated a source or evolved a physical state.

#### Whole-cell source representation bound

The [source-cell component](../src/recursive_horizons/nsc_atlas_source_cell_bound.py)
converts the interpolation residual into errors in the six quantities used
by the actual exterior stage: `r, log_C, c, r_x, log_C_x, c_x`.
`departure_slice` reconstructs its transverse jets from the evolved exterior
state; the separate source transverse splines are not inputs to this map.

Let `a,b` be the original cell endpoints, `l=b-a`, and let h-tilde be the
actual cubic defined by the stored binary `SourceSnapshot` coefficients.
The exact source obeys

\[
 h''=Q(x,h,h'),\qquad Q=\frac{gg''-2g'^2}{6s}+\frac{gg'v}{3s^2},
 \qquad s=x+\alpha h,\quad v=1+\alpha h'.
\]

Its left h and h' vanish because the analytic pulse is quiet there. The
saved polynomial has exactly that left jet. On a common tube with
`s>=s_min`, `|v|<=v_max`, take uniform pulse bounds G0,G1,G2. Then

\[
 L_h=\frac{\alpha(G_0G_2+2G_1^2)}{6s_{\min}^2}
       +\frac{2\alpha G_0G_1v_{\max}}{3s_{\min}^3},\qquad
 L_w=\frac{\alpha G_0G_1}{3s_{\min}^2}
\]

bound the two partial derivatives of Q. Enclose the actual polynomial
residual `D=ess sup |h_tilde''-Q(x,h_tilde,h_tilde')|`. For the error e,
integration of `|e''|<=D+Lh|e|+Lw|e'|` twice gives

\[
 \theta=L_hl^2/2+L_wl<1\quad\Longrightarrow\quad
 \boxed{\|e\|_\infty\le\frac{l^2D}{2(1-\theta)},\qquad
 \|e'\|_\infty\le\frac{lD}{1-\theta}.}
\]

These are source representation bounds, not a coupled stability estimate.
A source continuation argument
uses the corresponding bound Qmax on |Q|: the inequalities
`a-alpha*l^2*Qmax/2>s_min` and `1+alpha*l*Qmax<v_max`, together with
`1-alpha*G0^2/(6*s_min^2)>F_min`, keep the exact solution in the tube.
The polynomial is separately enclosed in that same tube.

The next PPoly piece owns the final knot in the actual source interface.
Its binary value and first jet can differ from this cubic's exact left limit
by coefficient-rounding jumps J0,J1. The endpoint-inclusive errors are
therefore `H+|J0|` and `W+|J1|`. The endpoint is checked in the common tube
and included in the frame box; the source coefficients themselves stay
unchanged. This does not use the next piece's second derivative.

The frame map is

\[
 F=1-\frac{\alpha g^2}{6s^2},\quad r=\sqrt{s^2-\alpha g^2/6},\quad
 L=\log F,\quad r_x=\frac{sv-\alpha gg'/6}{r},\quad
 L_x=-\frac{\alpha}{3F}\left(\frac{gg'}{s^2}-\frac{g^2v}{s^3}\right).
\]

Its interval Jacobian in `(s,v,g,g')` converts the radial error bounds to
the used frame quantities. The declared front can differ slightly from
the exact analytic support root. The calculation encloses the compact
pulse on the intersection of its coordinate interval with `[-1,1]`, adds
zero for the quiet branch, and separately bounds the sliver omitted by
the snapshot's additional declared-front mask. It does not set the exact
analytic field to zero at a rounded front.

Arb at 160 bits uses just the two original panels separated by that front.
Endpoint powers preserve positive denominator bounds; repeated wide-ball
multiplication must not erase a known positive floor. Six independent
controls initially passed in eleven development executions, including one diagnosed
interval-dependency failure, plus a small arithmetic diagnostic. Recorded
charges are 120,144 values; allowing a further 8,192 for the two initially
uncharged failing frame evaluations gives 128,336 within the admitted
150,000-value ceiling. Active CPU totals 0.40625 seconds. Charges now precede
the interval work, including failures. The controls include
symbolic Q derivatives, a closed-form forced ODE, 85-digit independent frame
derivatives and compact-source jets, zero-coupling recovery and rejected
domains. No physical source or trajectory has been replayed for them.

The [first cell certificate](../results/summaries/nsc-atlas-source-cell-bound-v1.json)
uses 10,372 values / 2.4375 CPU seconds and binds 59 producer files at
`ac98a98`. On the declared closed cubic, it bounds radius error by
`7.86156057e-16` and radius first-jet error by `2.48452710e-12`.
Its exact-record audit verifies theta and both integral inequalities, and
finds right-owner minus left-limit differences `-2.08603e-30` and
`-2.85265e-28`. That rational audit costs 1,024 values / 0.15625 CPU seconds.
The optional endpoint-owner correction now passes its manufactured control
with another 11,240 values and CPU below the timer resolution. Seven distinct
controls pass; the development value ceiling including failed-work allowance
is 139,576. The original v1 receipt is preserved for its declared cubic.

The [endpoint-inclusive certificate](../results/summaries/nsc-atlas-source-cell-bound-v2.json)
at frozen `191a05b` completes with 10,372 values / 2.5625 CPU seconds.
It retains the same residual, pulse envelopes, two panels and theta as v1;
the exact endpoint jumps are included in the radial and frame bounds.
Grok's separate review confirms the continuation, interval map and final
owner correction. All 59 source bindings authenticate. In the existing
dimensionless source coordinates, the whole-cell upper bounds are:

| Source input | Representation error bound |
| --- | ---: |
| Einstein-frame radius r | `7.861560568382867e-16` |
| First radius jet r_x | `2.4845270927685856e-12` |
| Logarithmic lapse log C | `8.186704230948388e-32` |
| First lapse jet (log C)_x | `3.0996481306140273e-27` |
| Reduced scalar c | `9.419625771084104e-59` |
| First scalar jet c_x | `1.1117684512743183e-41` |

The radial inverse has `theta<=4.184327547902973e-15`, with
`|h_tilde-h|<=2.863505873123005e-12` and
`|h_tilde'-h'|<=8.99410618409098e-9`, including the endpoint owner.
The scalar entries bound only the declared-front mask difference from the
analytic pulse. They are bounds on the prescribed value/jet representation,
not a claim of smoothness across that rounded mask or a floating-evaluation
error estimate. The possible mask jump is covered by the same entries.

The certificate performs no source readout calls or integrations. Its 10,372
charged values include 132 existing-knot/geometry validation values and
10,240 interval-work charges. It writes 14,751 bytes with exact dyadic upper
bounds and directed floating displays; no raw sampled arrays or time history are added.
This supplies one input to the coupled estimate. Uniform implementation,
exterior residual and inverse-evolution bounds remain required.

### First retained-source fixed-cap interval

The [exact small-matrix helper](../src/recursive_horizons/nsc_atlas_cap_linear_bound.py)
integrates the fixed rational basis to form M,K,e at a declared wall cut.
It computes the initial frozen bound `L_*^2=2tr(M^-1 K)+(e^T M^-1 e)^2`
and its integer square-root ceiling. Two [independent controls](../tests/test_nsc_atlas_cap_linear_bound.py)
recover the one-cell invariant `3444` and the exact mass, stiffness and
wall value of a global linear field on two cells. Both pass with
**38,144 charged values / 0.03125 active CPU seconds**. Grok mathematical
review confirms the integrals, bound and integer ceiling.

The [onset runner](../scripts/run_nsc_atlas_fixed_onset.py) uses the same
retained source interval 35 and splits time at the analytic scalar front.
It starts with quiet total cap fields and exterior/wall departures at the
last checked quiet value/first-jet knot. The original finite amplitude and
129 Hermite knots stay fixed. The exterior retains its existing 64-cell
quadratic derivative.

Let `dT,dZ` be the reference endpoint changes, `Z0` the initial wall
coordinate, and `a=max(4*dZ,2*dT)`. The two fixed cells have knots
`[0,Z0-a,Z0+a]`. The initial cut is one half and the reference end cut is
at most five eighths. The quiet reference characteristic starting at the
initial wall ends at `z=Z0-dT`; therefore `a>=dT` contains this reference
incoming layer inside the last cell. These are geometric sizing rules;
two cells are not asserted to resolve the pulse, and `T+z` is not the
exact finite-amplitude characteristic label.

Each time piece uses `ceil(16*L*T_max*piece_length)` equal steps, where L
is the initial integer cap bound and `T_max` is the integer ceiling of
the reference endpoint clock rates. Actual positive `T_x<=T_max` and the
cap cut domain are checked at every stage. This selects a time scale from
the initial frozen operator. It does not extend the frozen energy theorem
to the moving, state-dependent cap/exterior coupling.

Every accepted step includes a fresh midpoint Hermite residual. Cap
components use separate physical `L2(0,Z)` norms through the current Gram;
wall and exterior components retain their own units. The derivative of
the derived wall scalar is checked separately. Scalar, temporal mass,
source-curvature and wall-jet defects remain visible. The independent
pure-U constraint algebra is evaluated on the original centre and wall
patches, using the quadratic radial curvature of the departures. It is
a selected spatial diagnostic, not a continuum constraint enclosure.

The runner keeps only aggregate maxima, physical endpoint values, source
and producer identities, an accepted-state history hash, and execution
counts. A stopped attempt records its phase and reason without claiming
the unfinished endpoint. Source/code review found no actionable defect in
the stage physics or time schedule. The selected pure-U diagnostic and
nullable initial metric minimum were added afterward and checked by the
parent; they do not change the state update.

The [first attempt](../results/summaries/nsc-atlas-fixed-onset-v1.json)
has `dZ=0.0009446994071695691`, `dT=0.0013384741798292757`, `L=12106`
and `T_max=3`, giving **236+135 steps**. It stopped at the CPU check after
**341/371 accepted steps**, at `x=0.42287166568847995`, using
**235,384,613 values / 180.015625 CPU seconds** and 1,708 attempted common
stage calls. The initial maximum rate was `9.73353e-14`; scalar zero
recovery held through the source front. No right-endpoint result or
completed interval is claimed. The summary includes extrema from attempted
intermediate stages as well as accepted states.

The minimum stage metric certificate was `0.773878364523775`. Selected
exterior pure-U residuals reached `2.02810e-11`, and source residuals
`9.78749e-11`. Weak cap mass balance stayed within `1.02438e-20`, but the
conformal wall-flux gap reached `7.07677e-4` (the scalar-tangent gap
`E_f` was `1.48821e-3`). The sampled continuous
time residual in the derived wall scalar was `1.17754e-9`; cap density
and chi residuals had maximum L2 norms `9.26385e-8` and `1.46171e-11`.
These are different residuals, not interchangeable error enclosures.
They motivate a boundary-representation diagnosis before another interval
attempt. Completing the remaining steps would not by itself explain the gap.

### Finite-space onset lag at the scalar boundary

Consider the same finite-space equations in the zero-coupling limit, with
smooth prescribed geometry, nonsingular `M(0)`, zero initial scalar
coefficients, and a drive `C(t)=c_m t^m+O(t^(m+1))`, `m>=1`, `c_m!=0`.
Here t measures time from onset and `ell=phi_z(Z)`, `b=A0/R`. Writing the
coefficient vectors as q,d avoids confusing q with the wall scalar c:

\[
q'=T_xd,\qquad
d'=T_xM^{-1}\{-Sq+e(C-e^Td)\}.
\]

Leading-order integration gives

\[
d=\frac{T_{x0}M_0^{-1}e_0c_m}{m+1}t^{m+1}+O(t^{m+2}),
\qquad
q=\frac{T_{x0}^2M_0^{-1}e_0c_m}{(m+1)(m+2)}t^{m+2}+O(t^{m+3}).
\]

Consequently the strong incoming boundary residual satisfies

\[
r_B=(\ell-be)^Tq+e^Td-C=-c_m t^m+O(t^{m+1}).
\]

The trace of the finite-dimensional solution responds one time power later
than its boundary drive. Reducing time-integration error at fixed space
cannot remove this leading residual. In this limit the weak and strong
momentum traces coincide, while `dh=r_B/A0`. The previously derived
moving-wall identity therefore becomes

\[
E_f=\frac{(R_x+T_xA_0)r_B}{2RA_0},\qquad E_c=R E_f.
\]

Weak mass balance can hold while this strong boundary relation fails.
Two [exact controls](../tests/test_nsc_atlas_cap_onset_lag.py) verify the
cubic-drive Taylor recurrence and this scalar-tangent conversion using
**24,000 algebra values / 0.40625 active CPU seconds**. Independent Grok
review confirms the leading orders and factors. The scalar coefficients
are exactly zero before the front in the numerical attempt; its nonzero
pre-front gravitational interpolation residual is a separate matter.

This result establishes a mechanism in the zero-coupling finite-space
equations. The following calculation identifies the actual hybrid drive's
leading coefficient. Neither result attributes the entire measured
finite-amplitude gap. No time or mode scan, weak-flux replacement, or
interval retry follows from this diagnosis.

### Actual hybrid onset from the retained source

Let x_star be the mathematical root of `zeta(x_star)=-1`, not the rounded
stored `support_x` label, and put `t=x-x_star>0`. In the original source
coordinate,

\[
 V(x)=\frac{(1-x/a)e^{x/a}}{U_0},\quad a=0.8,\qquad
 \zeta(x)=\frac{4[V(x)-(V_a+V_b)/2]}{V_b-V_a},
 \quad \kappa=\zeta'(x_*)=
 -\frac{4x_*e^{x_*/a}}{a^2U_0(V_b-V_a)}>0.
\]

The original C3 source therefore has the one-sided expansion

\[
 g(x)=16\kappa^4t^4+O(t^5),\qquad
 g_x(x)=64\kappa^4t^3+O(t^4).
\]

To carry that coefficient through the exterior, use its actual scalar
prefix rather than equating the wall jet with g_x. Consider the
exact-arithmetic finite equations with the constrained mass and lapse
evaluated as mathematical functionals, zero scalar onset, fixed finite
basis and bounded spatial differentiation, on a regular local interval:
positive radii, lapses and H, nonsingular M, and a wall cut away from the
representation's domain edges. Geometric coefficients are right-smooth;
the exterior geometry need not be exactly vacuum. For the cap coefficients
q,d and exterior scalar state e, write y=(q,d,e). The scalar right-hand
side is locally Lipschitz, vanishes with y,g,g_x, and obeys

\[
 |y'(t)|\le K(|y(t)|+|g(t)|+|g_x(t)|).
\]

The integral inequality gives `y=O(t^4)`. The cap equation
`q'=T_x V(q,d)` has no direct source-derivative forcing, with `V(0,0)=0`,
so `q=O(t^5)`. This also makes the physical exterior scalar and its fixed
spatial first jets O(t^4), including the supplied source node. These are
local orders of the finite equations, with constants depending on the
chosen space; they are not uniform spatial-convergence estimates.

Two implemented trace laws are essential to that argument. The exterior
U-jet is the bounded derivative of its scalar state and endpoint values,
not an independently supplied source derivative. In the cap,

\[
 \chi_W=\phi_W^Tq,\quad k_W=\phi_W^Td,\quad
 \chi_{z,W}=\phi_{z,W}^Tq,\qquad
 w_-=\frac{k_W-\chi_{z,W}}{A_0R}+\frac{\chi_W}{R^2}.
\]

The actual common scalar is `c_W=R*u(f_W)`, with `f_W=chi_W/R`,
`u(f)=sinh(sqrt(alpha/6)*f)/sqrt(alpha/6)` and its continuous alpha-zero
limit. Its coefficient tangent is

\[
 (c_x)_{\rm cap}=\sqrt{D_W}
 \left(\phi_W^Tq'+\chi_{z,W}\frac{R_x}{A_0}\right)
 +(u_W-f_W\sqrt{D_W})R_x=O(t^4).
\]

It does not use d' or the hybrid-requested c_x. Together these identities
give `w_-=O(t^4)` and close the state-to-jet argument. Independent Grok
review identifies these as necessary hypotheses: the coefficient would
not follow from the source order alone if either trace law were changed.

Write `u=c/r`, `D=1+alpha*u^2/6`. The existing prefix gives, at the wall,

\[
 c_V=u_Wr_{x,W}+\sqrt{D_W}\left[
 \frac{g_x-u_0r_{x,0}}{\sqrt{D_0}}
 -\int_0^W r_x f_{\bar U}\,d\bar U\right],\qquad
 f_{\bar U}=\frac{c_{\bar U}-u r_{\bar U}}{r\sqrt D}.
\]

Here c_V is the hybrid's outgoing jet in its relabelled x coordinate,
as supplied by `departure_slice`; its name does not change that coordinate.
The fixed quadrature preserves the same orders. Since the integral and
u_0 are O(t^4), while u_W=O(t^5), this proves
`c_V=g_x+O(t^4)` under the stated local assumptions.

At zero cap state the constrained metric has A=A0=1-R^2 and mu=1.
For fixed finite alpha, its first changes are quadratic in the scalar
state: `N=O(t^8)`, `delta=O(alpha*t^8)`, and
`B=A0+O(alpha*t^8)`. The exact hybrid identities are

\[
 C_- =\frac A H(R\sqrt{D_W}w_- -u_W),\qquad
 w_+=\frac{2V_+c_V/H-u_W}{R\sqrt{D_W}},\qquad
 (c_x)_{\rm hybrid}=c_V+\frac{C_-}{2V_+}.
\]

The strong outgoing cap trace is O(t^4), so C_minus is O(t^4). Consequently

\[
 C=B Rw_+=\frac{2A_{0*}V_{+*}}{H_*}\,64\kappa^4t^3+O(t^4),
 \qquad
 \boxed{E_c=(c_x)_{\rm cap}-(c_x)_{\rm hybrid}
       =-64\kappa^4t^3+O(t^4).}
\]

The geometry cancels from the leading scalar-tangent defect. At alpha zero
this also follows from `T_x A0+R_x=H/V_plus` and the previous r_B formula.
This is a finite-space onset obstruction, including fixed finite alpha
under regularity, rather than a failure of the analytic cap-mass identity.
It supplies neither the remainder constant nor a bound over the attempted
window. Floating masking at a rounded support edge shifts the analytic
onset; the existing source-jet roundoff controls address numerical replay,
not equality of that binary edge with x_star.

### Incoming cap reference that carries the onset

Reuse the [reduced canonical incoming solution](#reduced-canonical-field-and-an-exact-incoming-wave-control).
In the alpha-zero cap, z=atanh(r) and the bulk primitive obeys
`chi_TT=chi_zz`. Before the front reaches the centre, its incoming family is

\[
 \chi(T,z)=F(T+z-\lambda_*),\qquad k=F',\qquad
 F(\lambda)=0\quad(\lambda\le0).
\]

This is the early part of the existing reflected characteristic solution;
it introduces no independently adjustable incident waveform. At the wall,

\[
 p=\frac{F'}{A_0},\quad h=\frac{F'}{A_0}-\frac F R,
 \quad w_-=\frac F{R^2},\quad C_-=0,\quad c_{\bar U,+}=0.
\]

Put `lambda_W=T+Z-lambda_star` and `b=A0/R`. The hybrid gives
`lambda_W,x=H/(V_plus*A0)>0` and

\[
 C=\frac{2c_V}{\lambda_{W,x}}-bF=2F'-bF,
 \qquad \frac{dF(\lambda_W(x))}{dx}=c_V.
\]

Thus the same exterior supplies the incoming primitive. It has a quartic
wall onset rather than the finite fixed-space quintic onset. For the
local quartic example `F(lambda)=a_F*lambda_+^4`, let
`d=T+Z-lambda_star` with `0<d<Z`. Direct integration gives

\[
 N=\frac14\left[\int_0^Z(k^2+\chi_z^2)dz-b\chi_W^2\right]
  =\frac{8a_F^2d^7}{7}-\frac{b a_F^2d^8}{4},\qquad
 N_T=\frac{k_W(C-k_W)}2.
\]

The moving derivative adds `Z_x N_Z`, retaining the derivative of b(Z).
The leading physical incoming layer's mass is seventh order, whereas the
fixed-space scalar density is eighth order. This comparison concerns their
different onset representations, not a claim that either exponent bounds
the numerical error over the full interval.

The reusable input is [PreparedFieldSolution.wall](../src/recursive_horizons/nsc_atlas_prepared_response.py),
which returns `wall_chi_over_epsilon` and `wall_dchi_dV_over_epsilon` from
the original source and exterior equation. Its argument is
`(V-Va)/(Vb-Va)`, and its derivative uses original V. Conversion to x
therefore multiplies that derivative by V_x. The wall-arrival map supplies
the incoming label. These are bounded numerical outputs, not an exact
continuum field or a finite-alpha characteristic. No prepared solve or
terminal-profile fit is run for this derivation.

Three [new controls](../tests/test_nsc_atlas_hybrid_onset.py) pass once using
**38,000 charged values / 1.28125 active CPU seconds**. They independently
expand the source and prefix, check the hybrid coefficient, integrate the
quartic layer's mass, and exercise the production hybrid on a compatible
incoming trace. The next representation must retain the incoming history
and account for its explicit time dependence in the mass derivative. The
earlier affine-background remainder and moving-knot obstruction still
apply; this alpha-zero reference alone does not authorize a finite-alpha
lift or another trajectory.

### Finite-coupling characteristic transport and common wall port

The [reduced canonical equations](#reduced-canonical-field-and-an-exact-incoming-wave-control)
already specify the finite-alpha dynamics. Their useful transport variables
in z=atanh(r) are

\[
 E^+=A_0(p+q),\quad E^-=A_0(p-q),\quad
 q=\chi_r,\quad \beta=B/A_0,\quad A_0=1-r^2.
\]

Here q is a field derivative, not the finite-space coefficient vector used
above. Put `mu=exp(delta)`, `P=3+alpha*P1`, `m=r^3/2+alpha*N`. The common
reaction is

\[
 S=B_r f+r\mu U_\alpha
  =\mu\left[r(U_\alpha-2f)
       +\alpha\left(\frac{2N}{r^2}-rP_1\right)f\right].
\]

Adding and subtracting the polar equations gives

\[
 \boxed{E^+_T-\partial_z(\beta E^+)=-A_0S,\qquad
 E^-_T+\partial_z(\beta E^-)=-A_0S,}
 \qquad \chi_T=\frac\beta2(E^++E^-).
\]

The plus family moves inward and the minus family outward. Their identical
sources cancel in the primitive constraint:

\[
 K=\chi_z-\frac{E^+-E^-}{2},\qquad K_T=0.
\]

Regularity at the axis gives `chi=0`, `E_plus=-E_minus`, hence zero primitive
and mass flux there. At alpha zero, beta=1 and S=0; the previous incoming
cap has `E_plus=2F'`, `E_minus=0`. Finite coupling generally generates the
minus family and changes the characteristic speed. Keeping `T+z` as an
exact finite-alpha characteristic or forcing the outgoing family to zero
would change these equations.

The [smooth-point operator](../src/recursive_horizons/nsc_atlas_cap_characteristics.py)
`characteristic_point_rhs()` implements these equations at at most 513
supplied points with `0<r<1`. It reuses `reference_potential()` and the
stable quotient `(U_alpha-2f)/alpha` from the existing reduced cap. The
metric derivative is evaluated with its analytic vacuum cancellation,

\[
 \beta_z=\alpha\mu\left(\frac{2N}{r^2}-rP_1-\frac{4N}{A_0}\right).
\]

Rates include `T_x partial_T + z_x partial_z` at moving query points.
An independently supplied chi_z is compared with `(E_plus-E_minus)/2`;
the operator reports a discrepancy rather than replacing that derivative.
The axis, internal interfaces and the history reconstruction remain the
caller's responsibility. This is a local PDE operator, not a radial solve
or a claim that a discretization preserves K.

At the material wall, use the actual chi and outgoing E_minus. With
`f=chi/R`, `u=sinh(a*f)/a`, `a=sqrt(alpha/6)`, the existing hybrid receives

\[
 c=Ru,\qquad w_-=\frac{E^-/A_0+f}{R},\qquad
 E^+=A_0\left[\frac{2V_+c_V/H-u}{\sqrt D}+f\right].
\]

The last expression is the one incoming boundary condition. It supplies
neither an independent wall-value waveform nor a second scalar condition.
The hybrid identities imply

\[
 \chi_x=\frac{H E^+}{2A_0V_+}+\frac{A E^-}{2A_0HV_+},\qquad
 c_x=\sqrt D\chi_x+(u-f\sqrt D)R_x
     =c_V+\frac{C_-}{2V_+}.
\]

`characteristic_wall_port()` makes the shared-value exterior callback,
calls the original hybrid and derives the incoming trace and wall rates.
It uses the same boost-gauge conversion as the fixed-cap adapter to obtain
the exterior lapse rate. Its scalar tangent residual checks this local
conversion; an independently evolved bulk trace must still be compared
with the returned incoming value.

The small frame term is evaluated without subtracting equal leading terms.
For `v=a*f`,

\[
 u-f\sqrt D=-f\left(\frac{v^2}{3}+\frac{v^4}{30}
 +\frac{v^6}{840}+\frac{v^8}{45360}+\frac{v^{10}}{3991680}+\cdots\right).
\]

On `|v|<=.1`, the displayed truncation has a real-arithmetic remainder
bounded by `|f|*v^12/[518918400*(1-v^2/180)]`. Floating function and
arithmetic errors remain separate. The original hybrid itself is unchanged;
the reported tangent residual includes any difference in its arithmetic.

The existing constrained mass identity, in these variables, has

\[
 p=\frac{E^++E^-}{2A_0},\quad h=\frac{E^+-E^-}{2A_0}-f,\quad
 \rho=\frac{r^2P_1}{2}+\frac A4(p^2+h^2),\quad
 N_T=\frac{AB}{2}ph,\quad N_x=T_xN_T+R_x\rho.
\]

The implementation returns rho and the temporal/moving mass **targets**.
It does not calculate an independent radial mass derivative, so their
presence is not a verified numerical mass balance. Physical mass remains
`m=r^3/2+alpha*N`, including its moving vacuum contribution.

All six [controls](../tests/test_nsc_atlas_cap_characteristics.py) pass in
seven executions using **27,948 charged values / 0.234375 active CPU
seconds**. The first five pass with 25,680 values; the stable-frame control
and one affected wall-port repeat add 2,268 (active CPU below the timer's
reported resolution). They check the exact conservative transformation,
primitive identity, an incoming alpha-zero profile, independent 70-digit
polar-equation arithmetic at finite alpha, moving query rates, reported
primitive inconsistency, the production hybrid's common value and both
traces, 80-digit frame-gap arithmetic, quiet recovery and invalid
metric/domain/budget rejection. No source, radial ODE or physical time
evolution is run. Independent Grok derivation confirms the continuum
signs, characteristic directions and nonlinear wall conversion.

The following section supplies the incoming-history field representation
and its actual constrained mass derivative. The existing characteristic profile
and regular background jets supply the reference history; finite-alpha
departures must retain the generated outgoing field and explicit-time
mass contribution. A common numerical history and boundary flux are still
needed before invoking these local identities in a trajectory.

In particular, [CapLeadingField.profile](../src/recursive_horizons/nsc_atlas_cap_field.py)
already maps the prepared wall solution through the arrival clock and
returns F,F' and optional F''. Its second jet uses the differentiated
exterior integral, not saved-waveform differences; it is not the exact
second derivative of the numerical field interpolant. This is the existing
source-to-cap route to reuse. `RetainedCapBackground` instead supplies
self-consistent derivatives of its declared C2 Hermite representation.
The original background/profile adapter made its second jets by a declared
gradient of supplied first jets. Those two numerical derivative meanings
must remain distinct when forming the history and its error account.

### Incoming history with fixed departures and an actual mass derivative

The [fixed-space radial adapter](../src/recursive_horizons/nsc_atlas_cap_fixed_reference.py)
now accepts the existing `RetainedCapBackground`, with fixed physical C1
quintic departures:

\[
 \chi=\chi_0(T,z)+\varphi(z)^Tc,\qquad
 k=k_0(T,z)+\varphi(z)^Td,\qquad z=\operatorname{atanh}r.
\]

Here \(k_0=\partial_T\chi_0\). The incoming front stays in the background
representation. Only its time derivatives are projected when forming the
departure rates. The radial forward sensitivity and independent adjoint
use the same total fields and inherited centre lapse. Alongside \(g_c,g_d\),
they calculate \(N_T\) at fixed coefficients and fixed physical radius.
This derivative is calculated from the radial constraint, not assigned
the wall-flux target.

Let \(a_h=M^{-1}\int\varphi k_0\,dz\),
\(b_h=M^{-1}\int\varphi\partial_Tk_0\,dz\),
\(e=\varphi(Z)\), \(a_W=k_0(T,Z)\), \(\mu=e^{\delta_W}\), and
\(C=BRw_+\). The implemented rates and total wall velocity are

\[
 v=2\mu M^{-1}g_d-a_h,\quad V_W=a_W+e^Tv,\quad
 w=M^{-1}[-2\mu g_c+e(C-V_W)]-b_h,
 \qquad c_x=T_xv,\quad d_x=T_xw.
\]

The actual derivative and its retained remainder obey

\[
 \begin{aligned}
 N_x&=T_xN_T+g_c^Tc_x+g_d^Td_x+R_x\rho_W\\
 &=T_x\left[\frac{V_W(C-V_W)}{2\mu}+\mathcal R\right]+R_x\rho_W,\\
 \mathcal R&=N_T-g_c^Ta_h-g_d^Tb_h
 -\frac{(a_W-e^Ta_h)(C-V_W)}{2\mu}.
 \end{aligned}
\]

Symmetry of the same Gram solve cancels the two bulk products. Neither
completeness of the departure space nor vanishing remainder is assumed.
The inherited central lapse affects the wall rate through mu; it does not
add a term to the constrained mass sensitivity because only the integrated
radial lapse difference enters the zero-centre-mass solution. Moving C2
history knots have continuous radial source, so their Leibniz jump is zero.
The wall endpoint contributes \(R_x\rho_W\). These are first variations;
the coefficient rate itself is not differentiated a second time.

The [metric certificate](../src/recursive_horizons/nsc_atlas_cap_fixed_bound.py)
combines original Hermite control hulls and exact fixed-cell polynomial
hulls over their existing rational cuts. Centre slope bounds remove the
apparent inverse-radius singularities. The radial source envelope supplies
an all-radius lower bound for A. Numerical interpolation, ODE and time
evolution errors remain outside that certificate. Radial integration and
six-point polynomial projection additionally split at the existing
small-z formula switch; no resolution scan is used. A solver endpoint
rounded at most one ULP beyond its panel is restored to that endpoint;
larger excursions fail.

Four [new controls](../tests/test_nsc_atlas_fixed_background_reference.py)
and six existing fixed-reference controls pass in 15 development executions,
using **4,569,376 charged values / 5.28125 active CPU seconds**. They include
identical total-polynomial versus background/departure calculations, an
independent clock difference for N_T, a nonrepresentable incoming quartic
front, a nonzero finite-coupling remainder, quiet/cancelling fields and
domain rejection. Grok independently verifies the affine mass identity
and the combined certificate's hull/envelope construction. Its axis
conditions are enforced by `FixedCapBasis`: the first knot is zero and
the first cell is odd.

For the quartic control \(F=(\lambda-1/2)_+^4/64\), \(T=1/8,Z=1/2\),
the computed mass is `1.268548382311671e-10` and explicit time derivative
`7.05434095864415e-9`, agreeing with the preceding independent formulas.
The common hybrid's conformal scalar-flux discrepancy is `1.35278e-13`;
its **actual** temporal mass discrepancy is `-8.81777e-18`. This is one
manufactured zero-coupling event. Departure rates are checked in the
physical Gram norm and at the wall, since raw coefficients amplify
quadrature error in a clipped basis. A separately rounded cancelling field
may leave mass below numerical resolution; negative/unresolved mass is
rejected rather than declared a physical zero.

The [hybrid adapter](../src/recursive_horizons/nsc_atlas_cap_fixed_hybrid.py)
passes these total fields through its existing common-value wall callback.
No prepared physical source or finite-time trajectory has been replayed.
The next connection is to supply the incoming history through the existing
prepared-field interface, distinguish its continuum second jet from its
numerical interpolation derivative, and keep the finite-coupling departure
and time-remainder checks in that common evolution.

### Preserving the source onset when retaining its history

The successful quartic control above does not establish that ordinary
quintic interpolation preserves a general physical onset. Let the first
active interval have length h, use \(s=x/h\), and supply
\(y=F(h),v=hF'(h),w=h^2F''(h)\), with the first three jets zero at x=0.
The usual quintic Hermite coefficient of \(s^3\) is

\[
 C_3=10y-4v+\tfrac12 w.
\]

For the independent profile \(F(x)=Ax^4+Bx^6\), this equals \(Bh^6\).
The interpolant therefore has a lower-order cubic term even though all
supplied endpoint jets are correct. Its small coefficient does not restore
the quartic asymptotic order. Merely setting the third front derivative
to zero and using a sextic is insufficient to preserve the coefficient:
for \(F=Ax^4+Bx^7\), the resulting \(s^4\) coefficient is
\(Ah^4+Bh^7\).

There is a unique degree-seven first piece preserving the known quartic
coefficient as well as the next value and first two jets. Set

\[
 \begin{gathered}
 C_4=Ah^4,\quad Y=y-C_4,\quad V=v-4C_4,\quad W=w-12C_4,\\
 P(s)=C_4s^4+C_5s^5+C_6s^6+C_7s^7,\\
 C_5=21Y-6V+W/2,\qquad C_6=11V-W-35Y,\qquad
 C_7=(W-10V+30Y)/2.
 \end{gathered}
\]

The eight conditions are the four zero derivatives of orders 0..3 at the
front, the specified fourth derivative \(24A\), and the three right jets.
The coefficient matrix is nonsingular. This piece joins the quiet zero
branch with C3 continuity and ordinary following quintics with C2
continuity. The fourth derivative at the front is **one-sided**; imposing
it on the preceding quiet interval would destroy that quiet interval.

For the existing compact preparation, write
\(\zeta=-1+4(V_{\rm null}-V_f)/\Delta V\) and
\(q-q_f=q_V(V_{\rm null}-V_f)+O((V_{\rm null}-V_f)^2)\).
Since \(g=(1-\zeta^2)^4\) and the initially quiet wall correction starts
at fifth order, the inherited reference profile has

\[
 A=16\left(\frac{4}{\Delta V\,q_V}\right)^4.
\]

This fixes A from the source and reference arrival clock. It is not a
parameter fitted to a later wall signal. The floating source mask and
inverse-arrival bracket need their own arithmetic account; a numerical
sample close to a mathematical support edge is not an exact zero jet.

Three [exact controls](../tests/test_nsc_atlas_history_onset.py) pass with
**zero field samples / 0.21875 active CPU seconds**. They derive the
quintic and sextic counterexamples, independently solve all eight
conditions, verify reproduction through degree seven, and recover A under
the arrival change of variable. This establishes onset order and
coefficient, not an interpolation error bound for the rest of the interval.
Independent Grok algebra confirms the three coefficient systems, the
source-to-arrival factor and the one-sided front restriction.

Consequently the next source adapter must first extend the existing
profile's **first active piece**, its centre moment sums and original-input
Bernstein certificate together. The present centre sums have quintic
degree limits; changing only interpolation would silently omit terms.
Keep the quiet interval zero and subsequent C2 pieces unchanged. Use
existing prepared-field cuts and arrival-geometry knots, with the source
edges explicit, rather than a new knot scan. Preserve the continuum
integral meaning of `CapLeadingField.profile(..., second=True)` and record
its difference from the resulting interpolant. No such prepared-history
transfer or physical evolution has yet been evaluated.

The [profile](../src/recursive_horizons/nsc_atlas_characteristic_background.py)
now implements the optional onset as one existing interval, using Bernstein
controls

\[
 (b_0,\ldots,b_7)=\left(0,0,0,0,\frac{Ah^4}{35},
 y-\frac{2hd}{7}+\frac{h^2e}{42},y-\frac{hd}{7},y\right),
\]

where d and e are the next physical first and second derivatives.
The fourth derivative is imposed only from the active side. Original
binary A and its interval index are copied and retained alongside the
Hermite inputs. Both metric certificates use the same original-input
controls, differentiating them with factors 7/h and then 6/h. The regular
centre moment routines include every term of the actual degree. The
older exact transfer into a quintic **total-field** space rejects this
history; it must remain in the retained-background route instead.

Eight [new controls](../tests/test_nsc_atlas_retained_onset.py) and 65
affected existing controls pass in 75 development executions, using
**787,030 charged values / 2.0625 active CPU seconds**. The first run's two
failed assertions were corrected: a C2 join does not match third jets,
and the underflow fixture initially triggered an earlier nonfinite guard.
The corrections add seven values, with CPU below the reported timer
resolution. The checks include an independent degree-seven polynomial,
one-sided quiet/front jets, C2 joins, degree-six moments, 80-digit centre
and formula-switch values, original-input Bernstein reconstruction,
analytic incoming mass/time derivatives and a finite-coupling independent
time difference. Grok independently confirms the Bernstein and moment
formulas. These are manufactured controls; the prepared history is next.
Grok's subsequent source review finds no concrete defect in the onset,
centre sums, original-input hull selection or guarded total-field transfer.

`build_prepared_background()` in the existing
[profile adapter](../src/recursive_horizons/nsc_atlas_background_profile.py)
now supplies the prepared-field connection. It uses the union of existing
field cuts, arrival-spline knots and source edges, rejecting more than 513
labels. One `CapLeadingField.profile(..., second=True)` call supplies the
value and two jets. The source and arrival derivative fix A; no waveform
differences or fitted amplitude enter this route. The original supplied
quiet-jet discrepancies must lie within explicitly supplied numerical
tolerances, zero by default, before the analytic quiet anchor is used.
Those discrepancies, the front coordinate roundoff, inverse bracket,
physical input domain and knot interpolation discrepancies remain visible.

Five new [adapter controls](../tests/test_nsc_atlas_prepared_background.py)
and 28 existing wall-profile controls pass once, using **8,596 charged
values / 0.09375 CPU seconds**. The public prepared-field budget exposes
its same count to cap callers; it is not reset between preparation and
later evaluation. The [one-source producer](../scripts/check_nsc_atlas_prepared_history.py)
is prepared to evaluate the 257/257 source and compare its native-interval
midpoints under the [separate finite allowance](prior-art-and-open-claim.md#one-prepared-physical-history-transfer-17-october-2026).
The [first physical transfer](../results/summaries/nsc-atlas-prepared-history-v1.json)
completes the source, then stops at the strict quiet check before any
interpolation comparison: the quiet jets are
`[0,5.37477524314984e-40,1.5572687669723015e-39]`. It uses 1,052,512
values / 2.03125 CPU seconds. The source shape and tolerances are unchanged.
The resulting C2 interpolation is still distinct from the continuum-integral
second jet sampled to construct it.

The [causal field interface](../src/recursive_horizons/nsc_atlas_cap_field.py)
previously returned exact zero only for q no later than the preparation
corner. Its actual fixed source has a longer known quiet interval. On
\(V_a\le V\le V_f\), the source is zero and its Volterra correction obeys

\[
 z_V(U,V)=-\int_{U_0}^{\min(U,U_W(V))}K(U',V)z(U',V)\,dU',
 \qquad z(U,V_a)=0.
\]

The finite-domain bounded kernel and Gronwall uniqueness give z=0.
The compact source is C3 at its front, so F,F',F'' are zero there too.
The interface now uses this known solution through the declared front
arrival. It leaves subsequent active-field evaluations in place. This
avoids converting an exactly quiet label through several floating-point
coordinates: a single V step near the nominal front changes the computed
source mask from zeta `-1.0000000000000016` to `-0.9999999999999956`.
The input front's `1.11022e-16` coordinate roundoff remains recorded; the
numerical arrival is not claimed to be the exact continuum clock.

Three targeted fake-prepared controls pass. The [second fixed-input transfer](../results/summaries/nsc-atlas-prepared-history-v2.json)
completes with zero quiet tolerances and 393 retained knots. Its 392
native-interval midpoint comparisons give maximum absolute discrepancies
`5.96343e-12`, `1.82137e-10`, `9.60804e-5` for F,F',F'', respectively.
Their global peak-relative values are `5.99031e-12`, `7.69983e-12`,
`8.03301e-8`. The onset coefficient is `429711.59455406887`, fixed by
the source and arrival derivative, and the largest F'' discrepancy lies
in its first active interval. Knot jet discrepancies are also retained.

The transfer and comparison use 1,362,327 values / 2.3125 CPU seconds.
Both physical attempts total 2,414,839 values / 4.34375 CPU seconds, with
no time steps, grid changes or saved full histories. These are finite
comparisons against the same prepared numerical field, not continuum
error bounds. Independent Grok review confirms the causal quiet argument;
an evolved solver's own quiet residual remains a separate diagnostic.
The common-stage connection is implemented below.

### Retained history at the evolving cap clock

The existing [common stage](../src/recursive_horizons/nsc_atlas_fixed_cap_stage.py)
accepts the retained background alongside the fixed cap departures. Each
stage evaluates it at `T_history=gauge.T+state.wall.T+background_clock_offset`.
The offset is a declared constant coordinate origin. The existing hybrid supplies
the physical T_x; the background time derivative is subtracted once from
the departure rates. The constrained mass rate remains

\[
 N_x=T_xN_T+N_c\cdot c_x+N_d\cdot d_x+R_x\rho(R).
\]

The total scalar trace enters the exterior callback before its outgoing
jets and wall rates are calculated. No second background is added at the
wall, and the explicit-time affine remainder remains visible. The existing
fixed-coefficient SSPRK3 adapter evaluates the same background at the
physical clock of each of its predictor, midpoint and accepted states.

Two [manufactured controls](../tests/test_nsc_atlas_fixed_background_stage.py)
check an independently represented total cubic history, a direct mass
directional difference, and one common time step with nonzero clock departure.
The mass-direction values are `2.2976403795706864e-5` (difference) and
`2.2976377177617093e-5` (equations); omitting the `5.14562e-6` clock term
would fail. Both new and eight existing controls pass with the
[development accounting and test corrections](prior-art-and-open-claim.md#retained-history-in-the-common-stage-17-october-2026)
recorded. This connects the representation and evolution interface; a
physical prepared-source departure trajectory and its error control remain.

The radial partition now omits only interior labels preceding a declared
onset whose endpoint jets are all exactly zero. Those pieces are the same
zero polynomial, so their internal cuts do not affect the degree bound
for the projection quadrature. The original first label, active onset,
later labels, basis cuts and centre formula switch remain. The original
inputs and all metric-certificate cuts remain unchanged: its closed
Bernstein hulls can see an active piece even at point contact.

A [single independent control](../tests/test_nsc_atlas_quiet_history_partition.py)
reduces five radial panels to four while preserving the certificate
exactly, and checks both projected time derivatives against rational
polynomial integrals. It uses 386,961 values / 0.609375 CPU seconds.
This is an exact zero-piece integration simplification; physical-window
cost and accuracy have not yet been measured with it.

The prepared field and wall gauge have different time origins. The gauge
defines T=0 at its quiet x=.4 anchor, while the retained field's absolute
cap time there is `Ta=-1.0664592301893958`. That inherited datum is the
prepared history's `background_clock_offset`; it is not selected by matching
an observed signal. At the source front, the shifted label differs from
the retained label by `1.01030e-14`, and the independent arrival derivative
differs by `-6.38496e-11`. These remain floating comparisons, not exact
continuum clock bounds. A constant translation of the history and its
evaluation clock preserves the common fields, rates and mass in the added
control. The clock rate T_x and explicit-time mass term acquire no extra
scale factor.

#### Derived time residuals and the full retained-history interval

Let `delta` denote Hermite derivative minus common-stage rate. Both paths
use the same gauge, so the wall departure residuals are physical `delta R`
and `delta T`. For departure coefficients `a,b`, total wall field
`f=chi/R`, and background `k0=partial_T chi0`, the chain rules are

\[
 \delta\chi=e\cdot\delta a+
 \frac{\chi_z}{1-R^2}\delta R+k_0\delta T,
 \qquad
 \delta c_W=\cosh(\kappa f)\delta\chi+
 (u-f\cosh(\kappa f))\delta R,
\]
\[
 u=\frac{\sinh(\kappa f)}{\kappa},\quad
 \kappa=\sqrt{\alpha/6},\qquad
 \delta N=N_T\delta T+N_a\cdot\delta a+N_b\cdot\delta b+\rho\delta R,
 \quad \delta m=\alpha\delta N+\frac32R^2\delta R.
\]

The zero-coupling limit is `u=f`. The radial trace is the total history
plus departure trace; `N_T` is the actual explicit-time constrained mass
derivative. Its affine remainder is not added a second time. The
[time adapter](../src/recursive_horizons/nsc_atlas_fixed_cap_time.py) now
provides these contractions. Two independent shifted-slice controls,
with and without the cubic background, agree within `6e-13` for the wall
scalar and physical mass, and `5e-15` for normalized mass. They use
630,439 values / 0.796875 CPU seconds and discriminate omission of the
clock terms.

The [admitted physical diagnostic](prior-art-and-open-claim.md#one-full-fixed-interval-evolution-with-retained-history-17-october-2026)
uses the original full onset interval and 371-step schedule with this
representation. A three-event source check found a pre-front metric
interpolation constraint of about `6.9e-11`; dropping the quiet prefix
would not be a justified replacement for evolving it. The source values
and first derivatives enter the common RHS, while its curvature and
transverse-jet interpolation defects remain separate diagnostics. The
same runner retains all of them, plus the affine mass remainder and the
clock-aware residuals, in a compact receipt. Its new creation-only mode
binds committed sources before computing and charges one field preparation
into the original overall clock and count. The complete coupled error
estimate and physical prediction remain outstanding.

The [completed receipt](../results/summaries/nsc-atlas-fixed-history-onset-v1.json)
uses the same full source cell and 236+135 time steps, with 1,855 common
stages. It costs 1,176,614,645 values / 619.4375 CPU seconds, including one
prepared field and transfer. The retained input hash exactly matches the
v2 history, and all 70 source bindings authenticate at `2551c9f`. The
minimum stage metric certificate is `0.773824049621412`.

| Sampled quantity | Maximum absolute value |
| --- | ---: |
| Conformal scalar wall-flux residual `E_c` | `5.40982898e-8` |
| Canonical scalar-tangent residual `E_f` | `1.13752336e-7` |
| Derived wall-scalar time residual | `5.38806476e-13` |
| Actual physical temporal-mass residual | `1.54459677e-15` |
| Affine explicit-time mass remainder | `2.50884930e-12` |
| Weak-balance minus explicit-time remainder identity | `2.39405933e-20` |
| Selected exterior pure-U constraint | `2.02810185e-11` |

The previous, shorter 341-step attempt had `E_c=7.07677230e-4` and
`E_f=1.48820981e-3`. The new complete interval covers its domain and has
recorded maxima about 13,000 times smaller. This supports the onset-lag
diagnosis; it does not supply a convergence order or a uniform error bound.
The new weak-balance defect is `5.28415752e-12`: its nonzero affine remainder
is recorded, not erased to reproduce the polynomial-only balance.

The final source knot has two metric curvatures. Writing the actual local
Hermite polynomial as `a3 t^3+a2 t^2+a1 t+a0`, the left endpoint curvature
is `6 a3 Delta+2 a2`, while the next piece gives `2 a2_next`. At source
knot 36, these are `-1.13548391e-6` and `0.00272394368`; the common source
ODE target is `-5.82457777e-6`. Their actual null constraints are
`6.08790298e-10` and `3.54409283e-7`. The run's maximum uses the second,
as specified by the snapshot's right-interval convention. The direct frame
constraints agree with `alpha*(h_xx-Q)/sqrt(F)` to below `1e-15` on both
sides. This independent 134-value check changes neither the saved result
nor the input. Curvature is not a first-jet RHS forcing, and different
constraint and scalar residual units preclude assigning causation by
comparing their magnitudes. The remaining junction gap needs its own
incoming-characteristic/weak-projection decomposition before another run.

#### Separating the remaining scalar boundary channels

With the same total cap traces and the common hybrid input, define

\[
 d_v=v_W-Bp_W,\qquad
 d_+=\frac{p_W+h_W}{R}-w_+.
\]

The existing matched-pair equations give the exact decomposition

\[
 \boxed{E_f=\frac{T_x}{R}d_v+
 \frac{R_x+T_xB}{2}d_+,\qquad E_c=R\sqrt D\,E_f.}
\]

This separates the weak-velocity trace from incoming-characteristic
compatibility. The latter involves both cap traces and the common incoming
data; it is not automatically an exterior-only error. Independently
differentiating the moving field proves the split, including cases with
only one nonzero channel. The existing hybrid adapter reports both terms
and their identity defect without changing any state rate or wall law.

The velocity term has a further exact description. Hold chi, T and R fixed
while varying a canonical-density coefficient `d_i`. Set `k=A0 p=k0+k_h`,
`A0=1-r^2`, `z=atanh(r)`, and let `phi_i` be the fixed basis. The mass and
lapse constraints imply

\[
 a=\frac{\alpha(p^2+h^2)}{2r},\quad \delta_r=a,\qquad
 (\partial_{d_i}N)_r+a\,\partial_{d_i}N
   =\frac{A p\,\phi_i}{2A_0}.
\]

The central mass derivative is zero. Its integrating factor, `dr=A0 dz`,
and `mu=exp(delta(R))` therefore give

\[
 2\mu\,\partial_{d_i}N
  =\int_0^{Z}\phi_i(z)\frac{B(z)}{A_0(z)}k(z)\,dz.
\]

Here mu multiplies the already evaluated mass gradient; no derivative of
the product mu N is taken. The chosen central lapse depends on the fixed
field chi, so a density variation adds no central-lapse term. Let `P_h` be
the exact L2 projection on the fixed basis and set `b=B/A0`. Subtracting
the background's projected time derivative in the departure equation yields

\[
 v_h=k_h+P_h[(b-1)k],\qquad
 \boxed{d_v=\{P_h[(b-1)k]\}(Z)-(b(Z)-1)k(Z).}
\]

At alpha zero the anchored lapse is zero and `b=1`, so this channel vanishes
for every admitted retained history. At finite coupling it is the endpoint
projection error of the metric correction times the total momentum, with
numerical constraint-gradient, projection and arithmetic errors entering
separately in an implementation. The existing forward/adjoint wall-velocity
comparison is useful diagnostic evidence, not an error enclosure.

Two exact controls pass with 14,000 algebra values / 0.390625 CPU seconds.
The second differentiates the mass constraint, checks the integrating
factor, and evaluates an independent nonconstant polynomial projection.
The quiet and nonzero numerical channels pass under the original hybrid
budget; their development record includes one corrected size assumption.
The completed v1 summary does not retain the coefficient state needed to
numerically divide its physical gap into these channels. Future receipts
retain one small final state, wall rates, trace and gradient diagnostics,
and its existing state hash. The v1 receipt remains unchanged, and no
trajectory has been repeated for this algebraic result.

The [single endpoint limit control](../results/summaries/nsc-atlas-fixed-history-endpoint-velocity-v1.json)
uses the original prepared history, endpoint radius/time and fixed basis,
with alpha zero and zero departures. It gives `d_v=2.10085855e-15`.
The independently contracted forward and adjoint wall discrepancies are
`2.10098798e-15` and `-3.23501807e-18`, with a contraction identity defect
`1.29427e-19`. The two projected L2 discrepancies are `3.76349e-17` and
`3.75948e-19`. It completes in 1,696,469 values / 2.546875 CPU seconds.
These are errors measured against the exact zero-coupling identity at this
event, not bounds on the earlier finite-alpha trajectory. No change to the
mass-gradient solver follows from this check. The next local comparison
traces the leading prepared field through its clock, history and common port.

At alpha zero the same mass functional also has the exact form

\[
 N=\frac14\int_0^Z(k^2+\chi_z^2)\,dz
       -\frac{A_0(R)\chi(Z)^2}{4R}.
\]

Expanding the original radial density gives this identity by one integration
by parts; regularity removes the central boundary term. Its first variation
is `2g_chi=-load(chi_zz)+e*A0*h_W`. For a retained linear-wave history
with zero departures, `chi_zz=k_T`, `v_W=k_W` and the projected interior
forcing cancels. The remaining canonical-density rate is therefore

\[
 \boxed{d_x=-T_x A_0(R)R\,d_+\,M^{-1}e.}
\]

Here `e` is the wall evaluation vector and M the dz Gram matrix. This is
an instantaneous statement at a moving-wall event, without a time derivative
of M. An independent exact density/variation control passes with 6,000
algebra values / 0.203125 CPU seconds. Numerical constraint-gradient and
projection errors remain separately measurable against this boundary force.

The [leading reference endpoint](../results/summaries/nsc-atlas-prepared-endpoint-connection-v1.json)
identifies a history-transfer contribution: retained-minus-direct prepared
material rate `6.76460751e-9`, versus direct-minus-exterior-prefix
`-3.22726177e-12`. The common conformal gap is `6.76138468e-9`, and its
weak-velocity channel is at roundoff. The two incoming labels differ by
`2.54241e-14` and their rates by `2.14395e-11`; no clock was adjusted.
The computed canonical-density rate has L2 norm `6.52894299e-7`, with
only `4.03972e-13` remaining after subtracting the predicted boundary force.
It completes with 1,708,414 values / 2.5625 CPU seconds and no time steps.

This is a leading reference state, not the earlier evolved finite-alpha
endpoint. A further local split must distinguish differentiation of the
retained polynomial from differentiation of the prepared dense output and
the prepared equation-derived first jet. The earlier midpoint comparison
does not bound discrepancies everywhere within an interval.

The [local dense-polynomial audit](../results/summaries/nsc-atlas-prepared-endpoint-connection-v2.json)
isolates that difference without changing any equation. At history interval
88 (the special onset piece is 87), retained-minus-dense F' is
`1.88328538e-9`, while dense-minus-equation F' is `-3.06454081e-14`.
The two differences sum exactly to the observed derivative gap. The
analytic nested-polynomial value agrees exactly with the installed dense
evaluator. This audit uses 1,710,662 values / 2.46875 CPU seconds and records
only the relevant small dense/history polynomials, their source identity
and existing endpoint summaries.

For an ordinary quintic Hermite piece the leading degree-six remainder is
proportional to `(q-a)^3(q-b)^3`. Its derivative is zero at the midpoint,
which explains why midpoint-only derivative comparisons can miss it.
A degree-seven correction can retain the original endpoint jets while
matching an additional interior value and first derivative. With
`xi=(q-a)/h`, `psi=xi^3(1-xi)^3` and interior coordinate eta, use

\[
 H_7=H_5+\psi[A+B(\xi-\eta)],\quad
 A=\frac{\Delta F}{\psi(\eta)},\quad
 B=\frac{h\Delta F'-\psi'(\eta)A}{\psi(\eta)}.
\]

It leaves both endpoint jets through second order unchanged and reproduces
degree-six and degree-seven remainder terms exactly. Two independent
symbolic controls pass with zero physical field samples / 0.125 CPU seconds.
The actual binary interior label must be retained, rather than silently
replaced by an unrepresentable exact midpoint. Implementation must carry
these inputs into every certificate; the special quartic-onset piece stays
separate. No physical evolution with this enrichment has been run.

Optional interior jets are now implemented in the existing profile and
retained-background owners. The original binary labels and jets enter the
runtime polynomial and both Arb metric certificates through one Bernstein
formula. The prepared adapter queries ordinary active pieces only, leaving
the quiet prefix, first active quartic piece and knot grid intact. A
quintic total-field transfer rejects enriched histories explicitly. Default
constructors retain their existing inputs and behavior.

The [independent controls](../tests/test_nsc_atlas_interior_history.py)
check exact rational reconstruction, C2 joins, stable moments, mass,
immutable original inputs and a bubble invisible to endpoint-only bounds.
Six new and 105 existing controls pass in 112 development executions;
the [admission record](prior-art-and-open-claim.md#degree-seven-interior-jet-enrichment-17-october-2026)
retains the corrected numerical test assumption and compute costs. The next
one-event comparison has fixed numerical decision limits; no enriched
trajectory or continuum error estimate is implied by these controls.

The [first enriched reference](../results/summaries/nsc-atlas-prepared-endpoint-enriched-v1.json)
passes all predeclared local limits: its conformal gap is `-4.70422e-13`,
down from `6.76138e-9` in the unchanged leading reference with quintic
ordinary pieces. History material-rate difference is `3.25911e-13` and
the prepared/prefix rate difference is `-8.02150e-13`. The density-rate
difference from the incoming boundary force is `5.81967e-13` in L2.
One preparation and event use 1,866,626 values / 2.75 CPU seconds.

The source, geometry, 393 knots, original nodal inputs, quartic onset and
direct prepared field are unchanged; 304 ordinary pieces acquire the
additional interior jets. Both old and enriched hashes are recorded.
This resolves the diagnosed transfer contribution at one reference event.
The [full enriched interval](../results/summaries/nsc-atlas-fixed-enriched-history-onset-v1.json)
now completes the same finite-alpha 371 steps / 1,855 stages and passes the
predeclared `1e-10` scalar decision. All 70 source bindings authenticate at
`42c7741`; original input, geometry and schedule match the earlier full
interval exactly. Both prepared-history hashes are unchanged from the
checked enriched reference.

| Sampled quantity | Maximum absolute value |
| --- | ---: |
| Conformal scalar wall-flux residual `E_c` | `2.75429853e-13` |
| Canonical scalar-tangent residual `E_f` | `5.79146106e-13` |
| Incoming scalar channel | `5.69213831e-13` |
| Weak-velocity scalar channel | `2.00642153e-14` |
| Forward/adjoint wall-velocity discrepancy | `4.53418108e-15` |
| Actual physical temporal-mass residual | `7.75251382e-21` |
| Derived wall-scalar time residual | `3.27464746e-15` |
| Affine explicit-time mass remainder | `1.32922763e-17` |
| Weak-balance minus explicit-time remainder identity | `2.75795871e-20` |
| Selected exterior pure-U constraint | `2.02810185e-11` |

These maxima include the four SSPRK/endpoint states and the continuous-path
midpoint check per step, not a supremum over time. The channel maxima need
not occur at the same event. The scalar channel identity defect is below
`2.33e-18`, and the metric certificate stays above `0.7738240496214122`.
Compared with ordinary quintic history, the conformal gap is about 196,000
times smaller. The radial mass and wall laws are unchanged. The final cap
departure L2 norms are `1.30859e-19` for chi and `1.11268e-15` for density;
their smallness is a result in this short onset domain, not an all-pulse
or finite-amplitude remainder bound.

One preparation and evolution cost 1,244,922,621 values / 660.703125 CPU
seconds, with no retries. The 28,312-byte receipt retains the small final
coefficient/exterior state, its SHA256, and summary diagnostics, not the
trajectory. Its source and selected exterior constraint residuals persist;
the already identified source-knot one-sided curvature distinction is
unchanged. The numerical scalar decision therefore does not supply the
continuous source/exterior error account or the coupled validity estimate.
Those are the next required connections; a further history/grid scan is
not indicated by this result.

### Central clock in general null labels

The observer conversion can be explicit. Parameterize the regular cap
axis by `lambda`, with positive rates
$d_-=dQ_-/d\lambda$ and $d_+=dQ_+/d\lambda$. Regularity gives

$$
d_-r_{Q_-}+d_+r_{Q_+}=0,\qquad
-4r_{Q_-}r_{Q_+}/C=1,\qquad
d_-c_{Q_-}+d_+c_{Q_+}=0.
$$

Consequently the Jordan central clock and scalar readout obey

$$
\boxed{\frac{dT_J}{d\lambda}
=\sqrt{C d_-d_++\frac{2\alpha}{3}(d_+c_{Q_+})^2},\qquad
p_{J,C}=\frac{2d_+c_{Q_+}}{dT_J/d\lambda}.}
$$
$$
\mathcal F_C=\frac{C d_-d_+}{(dT_J/d\lambda)^2}.
$$

These follow from `u_C=2 d_+ c_Qplus/sqrt(C d_- d_+)` and
`dT_J/dlambda=sqrt(C d_- d_+/F_C)`. Both physical outputs are invariant
under constant rescalings of either null coordinate. For a diagonal axis
`d_-=d_+=1`, the clock rate is
$\sqrt{C+2\alpha c_{Q_+}^2/3}$; the earlier anchored-gauge readout
`p_J,C=2 c_Qplus` follows when this rate is one. In another gauge, fix
`T_J` at the same quiet event and integrate its derived rate rather than
imposing a second clock normalization. This conversion assumes a regular
axis and does not by itself construct one.

## Connected scalar characteristic inverse

The source certificate can enter a first-order continuum comparison without
differentiating its six error bounds. The following estimate connects the
exterior, cap and wall on a **prescribed compatible geometry**. It does not
yet close the changes in that geometry caused by scalar errors. It uses
characteristic integration, rather than a finite cap mode matrix.

### One canonical reduced field on both sides

Write `X=r*f=r*asinh(a*c/r)/a`, `a=sqrt(alpha/6)`, with its continuous
zero-alpha limit. X equals the cap's existing chi. On the exterior the
canonical scalar is massless, so its exact reduced equation is

\[
 X_{\bar U x}=-K_E X,\qquad
 K_E=\frac{r_{\bar U}r_x+C/4}{r^2}.
\]

This follows by differentiating `X=r*f` and using both the radius and
canonical scalar equations. Unlike the invariant-c equation, this scalar
equation has no scalar first-jet product in its potential. It is an exact
change of variables, not a new source, action or discretization.

The certified invariant source inputs need this same change of variables.
Put `u=c/r`, `h=(1+a^2 u^2)^(-1/2)` and `j=asinh(a*u)/a-u*h`. At a source
point the complete differential, including the variation of h, is

\[
 dX=h\,dc+j\,dr,\qquad Q=h c_x+j r_x,\qquad
 dQ=h\,dc_x+j\,dr_x+
 \frac{h_u(c_x-u r_x)}r(dc-u\,dr),\quad
 h_u=-a^2u h^3.
\]

On the segment joining two source inputs, assume `r>=r_*`, `abs(u)<=u_*`,
`abs(c_x)<=c_*`, `abs(r_x)<=v_*`. Since `h<=1` and
`j_u=a^2*u^2*h^3`, one has `abs(j)<=alpha*u_*^3/18`. Therefore

\[
 e_X\le e_c+\frac{\alpha u_*^3}{18}e_r,\qquad
 e_Q\le e_{c_x}+\frac{\alpha u_*^3}{18}e_{r_x}
 +\frac{\alpha u_*(c_*+u_*v_*)}{6r_*}(e_c+u_*e_r).
\]

`canonical_source_error_bound()` implements these rational upper bounds.
They supply canonical d_X and d_q contributions from four of the six source
errors once the common tube is bounded; the lapse inputs belong to the metric
comparison. Neither `c_x` nor its error can simply be relabelled Q at finite
alpha. No second source derivative is introduced, and runtime floating
evaluation error remains additional to the existing representation bound.

The [original-cell application](../results/summaries/nsc-atlas-canonical-source-bound-v1.json)
now supplies that conversion. The authenticated source certificate has
`F_min=1/2`, `rJ_min=1/5`, `abs(rJ_x)<=3/2`; hence a conservative Einstein
radius floor is `r_min=F_min*rJ_min=1/10`. Its exact g and g_x envelopes give
`u_max=g_max/r_min`, `cx_max=gx_max`, and
`rx_max=(3/2)/F_min+alpha*g_max*gx_max/(6*r_min)` from the source frame
derivative. These bounds cover the exact source, the binary cubic and its
right-owner jet. The mask can only remove scalar values or first jets.
Positive radius makes c/r monotone along the straight input segment, so the
endpoint tube covers the required Jacobian path too.

Frozen `d365a13` yields the outward decimal bounds

\[
 e_X\le4.658\,10^{-38},\qquad e_Q\le2.576\,10^{-33}.
\]

The exact rational results are retained in the 6,437-byte summary, SHA256
`6817ecccf46d2eb0620ea081d1ffd1cc67a4d319431a9b1e9cdd80b41dfe3b9a`.
All 13 producer bindings authenticate. The calculation uses 1,024 charged
values, with active CPU below timer resolution; a separate 128-value rational
audit verifies the displayed outward bounds. The first command had stopped
before calculation because it passed the JSON input to a code-only source
graph; the corrected runner authenticates that input separately. No field
readout, integration or time step occurs. These small source representation
errors are inputs to the inverse, not yet bounds on the final wall signal.

At the wall put `U_+=d Ubar/ds`, `V_+=dx/ds`, `W_x=U_+/V_+>0`,
`A0=1-R^2`, `H=2 V_+ r_x+S_w R`, `S_w=2 cosh(a f)^3`. Define

\[
 j(f)=f-\frac{\tanh(af)}a,\qquad J=S_wRj(f),\qquad
 P=X_{\bar U},\quad Q=X_x.
\]

The existing hybrid law and the chain rule give both transmission directions:

\[
 \boxed{E^+=\frac{A_0}{H}(2V_+Q+J),\qquad
 P_W=\frac1{2U_+}\left(\frac{A}{HA_0}E^-+J\right).}
\]

The common value is `X_W=R*f`, with
`(X_W)_x=Q_W+W_x P_W`. Canonical normal derivatives have jump `-J`,
where the jump is exterior minus cap using the same outward normal. Setting
that jump to zero would alter the finite-alpha wall law. The invariant c
still has its homogeneous normal transmission.

At fixed R the full J derivative is

\[
 J_X=2\cosh y\sinh y\,(3y\cosh y-2\sinh y),\qquad y=af.
\]

Thus J and its derivative vanish at alpha zero. For a prescribed coefficient
comparison, an enclosure L_J of `abs(J_X)` gives value gains
`j_e>=sup(A0/H)*L_J`, `j_p>=sup(1/(2U_+))*L_J`.
Changes in H, null rates, radius and the remaining metric coefficients add
their actual algebraic differences to the boundary forcing. H itself depends
on the wall and scalar through the hybrid law; holding H fixed here is a
comparison decomposition, not a claim that it is independent physical data.

The principal cross gains are `a_e=2 V_+ A0/H` and
`b_e=A/(2 U_+ H A0)`. Their product is `A/(H^2 W_x)`. In proper-time
normalization it is `A/H^2=(beta_minus-v)/(beta_minus+v)`. An expanding
timelike wall makes the latter less than one, but this product alone is
not an inverse estimate: it omits exterior propagation and the cap source.

### Continuum inequalities and the instantaneous inverse

Use `x0<=x<=x1`, `0<=Ubar<=W(x)`, and `z_*<=z<=Z(x)` in the cap,
where `z=atanh r`, `Z=atanh R`, and the wall expands. The inner edge is
strictly positive and lies inside the quiet region. Only the incoming cap
family E-minus is supplied at that edge. The outgoing E-plus is zero there
as a consequence of no arrival, not a second boundary condition. A sufficient
quiet-core condition is that the wall's earliest disturbed radius in z,
minus `integral T_x beta_max dx`, remains greater than z_*. This needs a
uniform coefficient bound, not a sampled travel time.

All suprema below are running suprema through the indicated x, on the common
physical domains. For a scalar difference or linearized error, order them as

\[
 u=(M,w,p,X,q,E,Y)
  =(\|e_{E^-}\|,\|e_{X_W}\|,\|e_P\|,\|e_{X,+}\|,
       \|e_Q\|,\|e_{E^+}\|,\|e_{X,-}\|).
\]

Let `L>=sup W`, `k>=sup abs(K_E)`, `v_c>=sup T_x beta`,
`gamma>=sup abs(T_x beta_z)`, `v_w>=sup W_x`. On a common scalar tube let
`s_c>=sup T_x A0 abs(dS/dX)` for the cap reaction
`S=B_r f+r exp(delta) U_alpha(f)`, so
`dS/dX=B_r/r+exp(delta)*U_alpha'(f)` at fixed metric.
Choose uniform nonnegative upper bounds `a,b,j_e,j_p` for the wall gains.
All coefficient variations, moving-domain pullback terms and implementation
residuals must be included in the forcing. They are not assigned zero.
Piecewise absolutely continuous fields suffice when the characteristic
integrals exist; jumps in an approximate field contribute their accumulated
absolute jump sizes to the corresponding residual terms.

With nonnegative d_i bounding initial/source/boundary errors and the time or
space integrals of these residuals over the declared domain, characteristic
integration and the two primitives give

\[
\begin{aligned}
 M&\le d_M+\int(\gamma M+s_cY)\,dx,\\
 w&\le d_w+\int(q+v_wp)\,dx,\\
 p&\le d_p+bM+j_pw+\int kX\,dx,\\
 X&\le d_X+Lp,\qquad q\le d_q+LkX,\\
 E&\le d_E+aq+j_ew+\int(\gamma E+s_cY)\,dx,\\
 Y&\le d_Y+w+\int\tfrac12v_c(E+M)\,dx.
\end{aligned}
\]

For example, P is integrated along fixed-Ubar characteristics from the
initial cut or the advancing wall. Q is integrated from the source at fixed
x. The exterior value follows `X=X_src+integral P dUbar`. For the cap value,
a fixed-z line starts at the initial slice or enters through the expanding
wall; integrating `X_x=T_x beta(E^++E^-)/2` gives the last inequality.
Initial and boundary terms have been added rather than maximized, which is
a conservative bound. The d terms need finite uniform or integrated bounds;
stored point maxima do not supply them.

Write these inequalities as `u<=d+B*u+integral A*u`. B has only

\[
 B_{pM}=b,\ B_{pw}=j_p,\ B_{Xp}=L,\ B_{qX}=Lk,
 \ B_{Eq}=a,\ B_{Ew}=j_e,\ B_{Yw}=1.
\]

Its longest dependency chain is `M -> p -> X -> q -> E`, with four edges;
the chain starting at w is no longer. Therefore

\[
 B^5=0,\qquad G=(I-B)^{-1}=I+B+B^2+B^3+B^4\ge0,
 \qquad u\le Gd+\int GAu\,dx.
\]

This closes the instantaneous exterior/wall/cap scalar feedback without a
smallness assumption on `a*b` or `L*k`. The remaining feedback passes through
time integration. With positive component scales h_i, define

\[
 \lambda=\max_i\frac{\sum_j(GA)_{ij}h_j}{h_i},\quad
 \theta=(x_1-x_0)\lambda<1,\quad
 D_* =\max_i\frac{(Gd)_i}{h_i}.
\]

Taking running suprema gives the explicit conditional inverse

\[
 \boxed{u_i\le h_i\,D_*/(1-\theta).}
\]

The [exact rational comparison helper](../src/recursive_horizons/nsc_atlas_transmission_bound.py)
forms B, A and G and evaluates this sufficient short-time bound. Its inputs
are **already established upper bounds**; it does not turn sampled physical
values into certificates. Failure of theta<1 is inconclusive about stability.
Only first-jet source errors enter this uncommuted estimate. A comparison on
different geometries also needs bounded first spatial derivatives of the
reference characteristics for the transport-coefficient forcing, and a
compatible fixed-domain identification. This is distinct from demanding a
time derivative of every source error bound.

For the subsequent geometry comparison, one explicit fixed-domain map is
`Ubar=eta*W(x)` outside and `z=z_*+xi*(Z(x)-z_*)` inside. The exterior P
equation has transport coefficient `-eta*W_x/W`. The cap equations have
coefficients `-(T_x*beta+xi*Z_x)/(Z-z_*)` for E-plus and
`(T_x*beta-xi*Z_x)/(Z-z_*)` for E-minus, when the transport term is placed
on the left. Their difference multiplies the reference spatial first jet.
The map needs positive lower bounds for W and `Z-z_*`; wall timelikeness
fixes the incoming/outgoing signs at xi=1. These are the additional map
coefficients to bound, not an assumed fixed physical wall.

### The remaining metric connection has a radial bound

On the quiet-core collar, let r_* be its fixed inner areal radius. At a common
physical radius the two exact cap constraints are

\[
 N_r+a_cN=b_c,\quad \delta_r=a_c,\qquad
 a_c=\frac{\alpha(p_c^2+h_c^2)}{2r},\quad
 b_c=\frac{r^2P_1(f)}2+\frac{A_0(p_c^2+h_c^2)}4.
\]

Here a_c and b_c are radial coefficients; p_c and h_c are the
cap variables, not the error norms above. The quiet core anchors `N=delta=0`.
For two states, subtraction and `a_c>=0` give the exact variation equation
and its integrating-factor estimate

\[
 (\Delta N)_r+a_{c,1}\Delta N=\Delta b_c-N_2\Delta a_c,
 \qquad |\Delta N(r)|\le\int_{r_*}^r
       (|\Delta b_c|+|N_2|\,|\Delta a_c|)\,d\rho,
 \qquad |\Delta\delta|\le\int_{r_*}^r|\Delta a_c|\,d\rho.
\]

Uniform state bounds `abs(p_c)<=P_*`, `abs(h_c)<=H_*`,
`abs(P1')<=L_P` on the segment joining the two states imply

\[
 |\Delta a_c|\le\frac\alpha r(P_*|\Delta p_c|+H_*|\Delta h_c|),
 \qquad
 |\Delta b_c|\le\frac{r^2L_P}2|\Delta f|
       +\frac{A_0}2(P_*|\Delta p_c|+H_*|\Delta h_c|).
\]

At the same r, the characteristic variables supply
`Delta f=Delta X/r`, `Delta p_c=(Delta Eplus+Delta Eminus)/(2 A0)` and
`Delta h_c=(Delta Eplus-Delta Eminus)/(2 A0)-Delta f`.
Thus the radial metric difference uses these same scalar norms without
another scalar derivative. Also `Delta A=-2 alpha Delta N/r`, and

\[
 |\Delta\beta|\le\frac{\mu_*}{A_{0,*}}
 \left(\frac{2\alpha}{r_*}|\Delta N|+A_*|\Delta\delta|\right),
 \qquad \mu_*=\sup e^\delta.
\]

Here `A_*` is an upper bound and `A0,*` a positive lower bound. At zero alpha,
the finite normalized N may vary but its physical metric contribution and
the lapse difference vanish. Nonzero inner data or radial constraint
residuals add their anchor and integrated forcing terms. Different radii
or clocks add endpoint/pullback terms; the displayed same-radius estimate
does not silently drop them.

Seven [independent controls](../tests/test_nsc_atlas_transmission_bound.py)
pass in eight executions, using **99,024 charged values / 0.796875 active
CPU seconds**. They check the canonical exterior equation, exact J derivative,
both production-hybrid transmission directions and signed nonlinear jump,
the rational inverse against direct elimination with large cross gains, an
exact incoming quartic on an expanding cap, and the damped radial variation.
The source conversion additionally checks its full symbolic Jacobian,
an independent 75-digit finite difference against the bound, and zero-alpha
recovery. Grok's review confirms the transmission and inverse proof, and
identified the need to convert the certified invariant source data explicitly.
Its follow-up confirms the complete source Jacobian and rational bound.
The first run exposed a proper-time/source-time conversion omitted in the
test's wall tangent comparison; dividing by V-plus corrects it. The production
hybrid and derived transmission equations did not change. There is no
physical source readout, evolution, sampled coefficient enclosure or new raw
data in these controls.

This makes the next full-system obstruction specific: cap metric errors
depend instantaneously on E-plus, while the incoming wall law also depends
on cap A and the exterior metric prefix. Those extra edges can create a
cycle in the comparison matrix. Their uniform gains, exterior metric
transport, clock/wall variations and continuous residuals must be bounded
before inverting that larger system. The scalar nilpotence result does not
assert nilpotence, or an inverse, for the complete Einstein-wall system.

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
