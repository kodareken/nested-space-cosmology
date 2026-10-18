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
+\int_{\partial\mathcal M_\pm}\sqrt{-\bar\gamma_1},2X\bar K_{\rm out}
\right\}\\
&-\sigma\int_\Sigma\sqrt{-\bar\gamma_1},R^3.
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
