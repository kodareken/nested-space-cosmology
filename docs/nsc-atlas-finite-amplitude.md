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

This mixed-chart prescription has a noncircular constraint implication.
The actual cap normal acceleration, computed from its metric and the
wall curve, is

$$
K_{\tau,-}^{\rm actual}
=\frac{v_s+A_{-,r}/2}{\beta_-}
+\beta_-\delta_{-,r}+\alpha r v\Pi\Phi.
$$

Let $d=K_{\tau,-}^{\rm actual}-K_{\tau,-}^{\rm predicted}$, where
the prediction is the normal-force expression used above. The exterior
rapidity equation imposes its predicted acceleration, while the cap
constraints follow from radial reconstruction. With the exterior outgoing
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
