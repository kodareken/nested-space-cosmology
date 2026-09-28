# The compact warp changes the regulated Dirac source

The finite-cutoff determinant of the declared five-dimensional carrier now
has an evaluated compact spectral weight. It composes with the existing
four-dimensional Dirac operator, including a curved g4. This supplies a
previously missing part of the common source functional without adding a
scalar, a vacuum offset or a second gravitational weight.

The numerical application evaluates the homogeneous spherical potential
and its radius, compact-size, cutoff and warp variations. This is the
Dirac determinant magnitude for a specified Euclidean prescription. The
full compensating functional, finite matching terms, Lorentzian state and
transmitting interaction still have to be combined before solving a
self-sourced throat. The [preceding reduced source equation](nsc-horizon-source.md)
cannot identify its full U(r) with this potential alone.

## Reused ingredients and the one new connection

Reuse the [classical compact reduction](nsc-compact-mass-map.md),
[Euclidean adjoint domain](nsc-compact-boundary-action.md), and
[common coefficient normalization](nsc-vacuum-charge-matching.md).
The new question is how that same compact warp changes the finite quantum
weight, given that it cancels from the canonical classical evolution.
The monopole spectrum is imported from
[Borokhov–Kapustin–Wu, section4.1 and the appendix](https://arxiv.org/html/hep-th/0206054v2).
The proper-time and conformal variation methods are established in
[Vassilevich, sections2,3.3 and7.1](https://arxiv.org/html/hep-th/0306138v3).
None of these spectra, heat coefficients or generic horizon results is
claimed as new mathematics here.

## Keep the quantum measure and the adjoint domain

Use two free massless complex 5D Dirac copies with complementary compact
chiralities, unit gauge charge, no link mass and no torsion. In natural units,

\[
g_{5,s}=e^{2s\sigma(Y)}(g_4+dY^2),\quad
\sigma=-a(2Y/\ell)^2,\quad a=18/1015,
\quad Y\in[-\ell/2,\ell/2],\quad 0\leq s\leq1.
\]

The full point is s=1 and ell=2 L_star. This compact reflecting interval
is distinct from a radial transmission boundary. The permitted endpoint
values are P=(1-tau3 gamma5)/2; its complementary choice gives the same
free paired result.

The geometric Dirac operator is e^(-3s sigma) D0 e^(2s sigma). Transport
to fixed L2 measure gives instead

\[
D_s=e^{-s\sigma/2}D_0e^{-s\sigma/2},\qquad
f=e^{5s\sigma/2}\Psi=e^{s\sigma/2}\chi.
\]

D_s has the chiral domain; its adjoint has the complementary domain. The
positive determinant operator is L_s=D_s^dagger D_s. Its form in chi is

\[
\boxed{\int e^{-s\sigma}|D_0\chi|^2,dY
   =\epsilon\int e^{s\sigma}|\chi|^2,dY.}
\]

Only the chiral value constraint is imposed on this form. The compatible
Robin condition follows from it, as in the existing boundary calculation.
Discarding either weight would compute a different regulated operator.
This is a physical warp deformation at fixed cutoff, not a simultaneous
change of units.

In a paired nonzero D4 sector choose D4=lambda sigma1 and gamma5=sigma3.
The compact first-order expression is lambda sigma1-i sigma3 d/dY.
Cosine and sine functions provide the allowed and Dirichlet components.
Their weighted quadratic and Gram matrices give a symmetric generalized
eigenproblem. The implementation retains their mixing; it does not insert
the unwarped compact masses into a warped heat exponent by hand.

The real matrix uses the unitary phase basis chi=(f,i g), not a restriction
to physically real spinors. Direct component multiplication gives
|D0 chi|^2=(lambda g-f')^2+(lambda f-g')^2. Its off-diagonal matrix element
is -lambda integral w(f'g+fg'), w=exp(-s sigma). This is lambda integral
w'fg because g vanishes at the endpoints. It is zero for a constant weight
but generally nonzero with this warp. A direct complex-spinor norm and
independent scalar quadrature check this term: dropping it changes the
selected trial-field norm by about 0.0321555. The numerical record retains
that control to prevent a real-component restriction from changing the
Dirac operator.

## An effective spectral weight for the same D4

For y=lambda^2>0 let epsilon_j(y;s) be the compact singular values squared.
With E1(z)=integral_z^infinity exp(-t) dt/t, define

\[
h_{\Lambda,s}(y)=\sum_j E_1(\epsilon_j(y;s)/\Lambda^2),
\qquad
\boxed{\Gamma_{5,\Lambda}[g_4,s]=\tfrac12
                 \operatorname{Tr}_4 h_{\Lambda,s}(D_4^2).}
\]

This already includes the two complementary 5D copies. At s=0 the weight
is E1(y/Lambda^2)+2 sum_{p>=1} E1((y+(p pi/ell)^2)/Lambda^2), consistently
with one 4D Dirac zero level and two fields at each nonzero level. Adding
another factor of two would double-count them.

For a variation of g4, the existing spectral functional calculus gives
delta Gamma=Tr4[D4 h'(D4^2) delta D4], after transporting the varying
Hilbert space/domain consistently. Thus the weight and its derivative
can feed the existing covariant metric-response owner. A discrete D4
zero mode requires its separate IR/state prescription. The numerical
application uses a density on noncompact R2, where the large-proper-time
heat density vanishes; it does not discard a finite-box zero mode.

An independent warp derivative checks this matching. Put R_s=D_s D_s^dagger.
The closed-operator identity relating the L and R heat kernels gives

\[
\partial_s\Gamma_{\Lambda}
=\tfrac12\operatorname{Tr}\sigma
  [e^{-L_s/\Lambda^2}+e^{-R_s/\Lambda^2}-P_{\ker L_s}-P_{\ker R_s}].
\]

For the identical paired copies, tau1 exchanges the two adjoint domains
and commutes with sigma, so the smeared L/R traces agree. The kernel
density vanishes in the noncompact application. Differentiating the
generalized eigenproblem and evaluating this dual-heat trace are distinct
checks; their finite Galerkin difference decreases with compact resolution.
The finite trace identity does not require a claim about cancellation of
all boundary a5 invariants. A renormalized anomaly statement is a separate
limit and cannot remove this finite-cutoff contribution.

## Evaluated potential and homogeneous source

Set g4=R2 x S2(r), with magnetic flux integral F=2 pi q, q an integer.
The imported angular eigenvalues squared and degeneracies are

\[
\alpha_n^2=n(n+|q|)/r^2,\quad d_0=|q|,
\quad d_n=2(|q|+2n)\ (n\geq1).
\]

For q=0 the zero sector is absent. Nonzero d_n includes both angular
signs. The spin trace on R2 and the determinant's one-half produce

\[
V_{2,\Lambda}(r,\ell,s)=\frac1{4\pi}
 \sum_n d_n\int_0^\infty dz\;h_{\Lambda,s}(z+\alpha_n^2),
\qquad
\partial_r V_2=\sum_{n\geq1}
  \frac{d_n\alpha_n^2}{2\pi r}h_{\Lambda,s}(\alpha_n^2).
\]

V2 has units length^-2 per reference R2 area. The radius derivative uses
the integral's spectral endpoint; it is checked against a metric finite
difference. The compact-size and warp derivatives use both matrices in
the generalized eigenproblem. Dimensional covariance requires
r V_r+ell V_ell-Lambda V_Lambda=-2 V at fixed a,s,q.

At the unfitted development point q=1, r=1, ell=2, Lambda=2, s=1:

| Projection | Value |
|---|---:|
| Flat-product V2 | 2.21888250509 |
| Warped V2, 26 compact modes | 2.15498055956 |
| Radius derivative | 4.37111539727 |
| Compact-length derivative | 1.03018040018 |
| Warp-path derivative | -0.06242437994 |

The displayed potential changes by about -2.88% with the physical cutoff
held fixed. This does not contradict the unchanged classical KK masses:
the covariant quantum norm and regulated weight differ from canonical
free evolution.

For the homogeneous variational projection, rho=V2/(4 pi r^2),
p_parallel=-rho and p_sphere=-V_r/(8 pi r). The parallel null contraction
is zero. These projections are integrated over compact Y; they are neither
pointwise 5D stress nor the complete in-in stress on a varying throat.

## Bounds, uncertainty and source bookkeeping

Because -a<=sigma<=0, the quadratic forms imply, by the min-max principle,

\[
\epsilon_j(y;0)\leq\epsilon_j(y;s)
 \leq e^{2sa}\epsilon_j(y;0),
\quad
V_{2,0}(\Lambda e^{-sa})\leq V_{2,s}(\Lambda)\leq V_{2,0}(\Lambda).
\]

These are comparison inequalities for the stated free modulus, not for
the complete action. The application checks q=0,1,2 without selecting a
physical flux. Positive flat heat tails bound omitted angular, compact
and momentum sectors. Those bounds do not bound the errors in retained
Galerkin eigenvalues. The record separately varies compact basis size,
compact quadrature, angular cutoff, momentum extent and momentum tolerance.
From 18 to 26 compact modes the potential changes by about 3.93e-8;
the continuum warp-variation residual is about 3.82e-8 at 26 modes.
These are observed convergence measures, not a proved remainder estimate.

The same computation fixes a source contribution, not permission to add
another arbitrary u0. If canonical massless LLL modes are retained at a
4D matching cutoff nu, their homogeneous potential is |q| nu^2/(4 pi).
Define the complementary contribution by subtracting that same low action
from Gamma5. The sum is independent of that partition. This is a mode
matching prescription, not a sharp energy projection of the warped 5D
eigenvalues.

Gravity also must be counted once. For the already induced coefficient
A4, the spherical Einstein term inside this potential is -8 pi A4.
With G=(16 pi A4)^-1, its contribution 2G V_Einstein equals -1.
Including it in U(r) while also retaining that same Einstein functional
on the geometric side would duplicate the spherical curvature term.
The magnetic term is 2 pi C4 q^2/r^2 and obeys the existing gauge/Newton
normalization. Move or retain each full invariant functional consistently;
a constant-radius potential cannot determine every kinetic or curvature
variation by itself.

The next assembly therefore uses h(D4^2) in the common covariant source
and matches the remaining compensator, finite terms and state in the
same domain. It does not restart a collapse simulation or choose U(r)
from the desired neck. The full self-sourcing and observational goals
remain open.

## Reproduction

`scripts/check_nsc_warped_source.py --check` authenticates its source and
input closure, then compares every field of
`results/development/warped-source.json`, using the existing comparator
(absolute tolerance 3e-9, relative tolerance 3e-8; exact structures,
strings and hashes). It does not run prior scientific generators.
New output paths refuse overwrite. The primary sources, hypotheses,
numerical uncertainty and remaining physical inputs are retained together.
