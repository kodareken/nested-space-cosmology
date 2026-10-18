# A shared variational form of the original cap--wall action

The [original action and invariant transmission](nsc-atlas-finite-amplitude.md#exact-material-interface)
already require one field with continuous value and normal derivative at
the wall. P001's separate numerical representations do not preserve that
strong trace: the [actual-state diagnosis](nsc-atlas-finite-amplitude.md#pointwise-equations-on-an-archived-evolved-state)
finds the natural weak boundary lift. This note supplies an equivalent
action/interface and local Legendre map for a coherent next approximation.
It does not supply a new evolved solution or close P001's bulk error.

Use the same Jordan action, both outward `2fK` terms and material tension.
Here `R` is the **Jordan** areal radius, `chi=R*phi_J` is the unnormalized
invariant field, and `kappa=kappa_4`. Define the unit-sphere metric and
curvature coefficient

$$
g_J=R^2(\bar g_2+d\Omega^2),\qquad
X=\frac{R^2}{2\kappa}-\frac{\chi^2}{12},\qquad
\mathcal F=1-\frac{\kappa\chi^2}{6R^2}>0.
$$

The common radial metric is frame invariant:
`g_J/R_J^2=g_E/R_E^2`. This does not identify the two areal radii or masses;
existing Einstein-frame numerical inputs still require their exact frame
conversion. In fact `2*kappa*X=R_E^2`, connecting this curvature coefficient
directly to the existing Einstein area variable. No old spectral/Maxwell
candidate is imported here.

## Exact bulk and boundary reduction

After angular integration the original action is exactly

$$
\begin{split}
\frac{S}{4\pi}={}&\sum_\pm\left\{
\int_{\mathcal M_\pm}\sqrt{-\bar g_2}
\left[X(\bar R_2+2)+\frac3\kappa(\bar\nabla R)^2
-\frac12(\bar\nabla\chi)^2-V_{0,\pm}R^4\right]
+\int_{\partial\mathcal M_\pm}\sqrt{-\bar\gamma_1}\,2X\bar K_{\rm out}
\right\}\\
&-\sigma\int_\Sigma\sqrt{-\bar\gamma_1}\,R^3.
\end{split}
$$

The region vacuum densities remain `V0_minus=3*H^2/kappa`, `V0_plus=0`.
The direct conformal substitutions leave the bulk total divergence

$$
\bar\nabla_a\left[
\left(-\frac{3R}{\kappa}+\frac{\chi^2}{2R}\right)
\bar\nabla^aR\right].
$$

The transformed outward boundary term contains precisely its negative,
in addition to `2X*Kbar_out`; discarding it would change the interface.
The [component](../src/recursive_horizons/nsc_atlas_conformal_action.py) and
[controls](../tests/test_nsc_atlas_conformal_action.py) retain this identity.

In a local common normal collar `d nbar^2-exp(2*omega)*d taubar^2`, define
`A=[nbar(R)]`, `B=[nbar(chi)]`, `C=[Kbar2]`. Variation of the shared induced
metric, radius and field gives

$$
RA-\frac{\kappa\chi}{6}B=-k_{\rm wall}R^3,\qquad
3A+RC=-3k_{\rm wall}R^2,\qquad
B+\frac\chi3 C=0,
$$

with `k_wall=kappa*sigma/2`. The coefficient determinant is `-R^2*F`.
On the admitted branch its unique solution is

$$
[\bar nR]=-k_{\rm wall}R^2,\qquad
[\bar n\chi]=0,\qquad [\bar K_2]=0.
$$

These recover the existing Jordan jumps after undoing the rescaling. A
shared `C1` invariant field and regular conformal metric therefore retain
the correct interface, while the radius must carry its physical kink.
Forcing the radius into that same `C1` space would erase the tension.
The region potential and wall motion also enter the normal force/Bianchi
balance; these boundary variations alone are not an evolution scheme.

## Actual velocities from a common time action

For `d sbar2=-N^2 dt^2+a^2(dy+beta dt)^2`, set

$$
v_a=a_t-(a\beta)_y,\qquad
v_R=R_t-\beta R_y,\qquad
v_\chi=\chi_t-\beta\chi_y.
$$

The reduced time kinetic density is

$$
L_{\rm kin}=\frac1N\left[
-2\left(\frac R\kappa v_R-\frac\chi6v_\chi\right)v_a
-\frac{3a}\kappa v_R^2+\frac a2v_\chi^2\right].
$$

Canonical momenta per `4*pi` obey `p=M*v`, where

$$
M=\frac1N\begin{pmatrix}
0&-2R/\kappa&\chi/3\\
-2R/\kappa&-6a/\kappa&0\\
\chi/3&0&a
\end{pmatrix},\qquad
\det M=-\frac{4aR^2\mathcal F}{\kappa^2N^3}.
$$

The local Legendre map is invertible for `R,a,N,kappa>0,F>0`. Lapse and
shift constraints remain essential; this gravitational kinetic matrix is
not a positive matter-energy matrix. The axis is a separate regularity
limit. Physical field derivatives must be reconstructed from these
velocities and the represented field, rather than treating canonical
coefficients as pointwise velocities.
For a finite basis, this means assembling and differentiating the same
finite action before solving its Legendre map. Pointwise inversion followed
by an independent projection would not preserve that identity automatically.

Eleven controls check the bulk/boundary cancellation, direct normal
variation, determinant, shifted two-metric curvature/ADM sign, manufactured
velocity recovery, smooth-wall limit and branch restrictions. Their only
generated results are test output; no trajectory or extra JSON is retained.

## What a finite method still has to establish

A common invariant-field representation removes a particular two-sided
trace ambiguity. It does not make a finite bulk projection exact. The next
numerical method must derive its total energy and interface power, retain
the radius kink, handle the physical centre, and resolve the active wall
region. Moving-wall differentiation also includes the change of normal and
evaluation point. In covariant form,

$$
\nabla_{u_w}(\bar n^a\nabla_a\chi)
=(\nabla_{u_w}\bar n)^a\nabla_a\chi
+\bar n^a u_w^b\nabla_b\nabla_a\chi.
$$

An additional scalar phase-space term `lambda*C(q,p,t)` is insufficient on
its own: canonical evolution gives `Cdot={C,H}+C_t`, because `{C,C}=0`.
For `H=(p^2+k*q^2)/2`, `C=p+a*q-g(t)`, its rate is
`a*p-k*q-gdot`, independent of that multiplier. A replacement needs a
complete variational/port and constraint derivation, not a force chosen
to cancel the archived residual. The algebraic reduction alone does not
establish accuracy or convergence of a new P001 evolution.

## A shared finite field action on a quiet wall-normal collar

The [finite component](../src/recursive_horizons/nsc_atlas_collar_field.py)
uses the same quiet sources, `H=1,k_wall=2,M=.4`, on
`.45<=Rwall<=.55`. The physical Gaussian coordinates are wall proper time
`tau` and normal distance `|n|<=eta=.02`, with
`ds2=-N^2 dtau^2+dn^2`. This is a prescribed-background scalar control;
the metric momenta and backreaction are not evolved here.

On the cap,

$$
r=R\cos n+\beta_-\sin n,\qquad
N=\cos n+\beta_{-,R}\sin n.
$$

On the exterior the signed normal geodesic and Jacobi equations are

$$
r_{nn}=M/r^2,\quad N_{nn}=-2MN/r^3,\quad
w_{nn}=-2Mw/r^3+6MN^2/r^4,\qquad w=N_R.
$$

Initial wall data are `(r,r_n,N,N_n,w,w_n)=
(R,beta_plus,1,beta_plus_R,0,beta_RR)`. They imply `r_R=N`,
`r_tau=v*N`, `N_tau=v*w`. The checks retain both
`r_n^2-v^2=f(r)` and `r_n*N_n-r_nn*N=v*v_R`, including the signed
turning point `beta_plus=0`.

These coordinates have a direct positive bound over the declared collar.
For `rlo=Rmin+min(beta_plus(Rmax),0)*eta`, the exterior satisfies
`r>=rlo` and, before any first zero,
`N>=1+min(beta_plus_R)*eta-M*eta^2/rlo^3>0`. The cap has
`N>=cos(eta)` and `r>=Rmin*cos(eta)-max(beta_minus)*sin(eta)>0`.
The bound prevents a first lapse zero. It is a local chart with its own
slicing, not a repair or continuation of the failed wide coordinate map.

The invariant real scalar has the following action per `4*pi`:

$$
L_w=\frac12\int\!dn\left(\frac{\chi_\tau^2}N
-N\chi_n^2-V\chi^2\right)
-\frac12\sum_s b_s\chi_s^2,
\quad V_-=0,\quad V_+=\frac{2MN}{r^3},
$$

$$
b_s=\frac{\epsilon_s}{3}\left(N_n-\frac{Nr_n}{r}\right),
\qquad \epsilon_s=(-1,+1).
$$

The last terms are the scalar part of the original timelike GHY terms
at `n=+-eta`. At the material wall the two outward contributions cancel,
since `[N_n]=-k_wall` and `[r_n]=-k_wall*R`. They are retained at the
outer ends, including their radius derivatives in the work calculation.

Two cubic Hermite elements share value and slope at the wall. The six
coordinates are `(chi_left,eta*chi_n_left,chi_wall,eta*chi_n_wall,
chi_right,eta*chi_n_right)`. If `B,D` evaluate value and derivative,

$$
M_h=\int B^TN^{-1}B,\qquad
K_h=\int(D^TND+B^TVB)+\sum_s b_s e_se_s^T.
$$

The field is `C1`; the independently represented radius retains its
physical kink. The restricted original ADM momentum is not `M_h*qdot`:

$$
\pi_h=M_h\dot q+S_hq,\qquad
S_h=-\frac v3\int B^T r^{-1}B.
$$

The time-end generating function is `G=q^T*S_h*q/2`. The restricted
ADM action differs from `L_w` by `dG/dtau`; its scalar density momentum
is `chi_tau/N-v*chi/(3r)`, independently checked against the full
three-variable kinetic matrix above. This scalar canonical transformation
applies to the prescribed background, not a removal of the coupled metric
momenta.

With `u=solve(M_h,pi_h-S_h*q)` and `v=dR/dtau`, the finite equations are

$$
q_R=u/v,\qquad
(\pi_h)_R=(S_h)_Rq+S_hq_R+(F_h-K_hq)/v.
$$

The wave-action energy and its work identity are

$$
E_w=\tfrac12u^TM_hu+\tfrac12q^TK_hq,\qquad
(E_w)_R=-\tfrac12u^T(M_h)_Ru+\tfrac12q^T(K_h)_Rq+u^TF_h/v.
$$

This is not the total Einstein mass. Shared field traces cancel the
internal field power; a finite bulk mass constraint still needs its own
check in any coupled evolution.

Five controls check the exact geometry and independently differentiated
parameter dependence, cubic reconstruction and original transmission,
the full original ADM momentum, independently projected continuum forcing,
energy/work and the admitted domains. The single frozen application uses
`chi=.01*sin(pi*(R-.45)/.1)^2*(1+.3*n+.2*n^2+.1*n^3)`.
Its volume force comes from differentiating that field in the continuum
equation; its endpoint forces retain both normal flux and GHY. The field
lies in the finite space, so recovery tests the action and time integrator,
not spatial convergence or physical predictions. This manufactured force
is never supplied to the P001 model.

The [first frozen application](../results/summaries/nsc-atlas-collar-field-v1.json),
from commit `8cf250dd229621c2870b3be40ef81fabff9fa8e0`, completes the full
strip in 82 accepted steps, 1232 field RHS calls, 59549 spatial-background
RHS calls and 1.984375 process-CPU seconds. The pre-sampling lapse bound is
`.8854464933`; the radius bound is `.4234086834*ell`.

| Check | Absolute maximum | Normalized maximum |
| --- | ---: | ---: |
| Field coordinates | 6.7521e-14 | 6.7113e-12 |
| Actual proper-time velocity | 1.9712e-11 | 7.1079e-11 |
| Energy minus metric and input work | 2.6902e-18 | 1.7452e-15 |

The normalized denominators are the maximum exact field, exact velocity,
and energy/work scale over the declared 33 diagnostic cuts, respectively.
The original scalar transmission residual is `1.041e-17`; the shared
internal field-power difference is zero in reference arithmetic. The
momentum shift reaches `.0001074389`, so omitting it would fail the
original ADM momentum comparison. The preserved radius kink reaches `1.1`.
The single 16/32 assembly comparison at `Rwall=.5*ell` differs by at most
`2.558e-13` across all six matrices. These are numerical controls on the
specified manufactured field, not an error range for a physical P001
observable. The control is closed; it supplies no reason to rerun P001 grids
before the coupled metric/action replacement exists.
