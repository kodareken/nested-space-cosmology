# Upstream metric directions span one diagonal preparation response

For the existing raw KS Dirac vertices at a selected energy `E*`,

\[
V(\rho)=
\begin{pmatrix}
-m&0&0\\
\ell/r&0&-\ell/r^2\\
-E_*/a&E_*/a^2&0
\end{pmatrix},\qquad
\boxed{\det V=-\frac{m\ell E_*}{a^2r^2}\ne0}
\]

when `m,ell,E*` are nonzero. Consequently a sufficiently narrow common
positive parent bump produces three independent **actual metric-response**
columns. Real amplitudes of those lapse, axial-scale and radius directions
can cancel one specified diagonal traceless K while leaving every incoming
metric jet unchanged. This is an analytic first-order construction map;
no amplitudes, width, physical profile or field solve is selected here.

## Reused operator and missing connection

Reuse the [owned raw KS vertices](../src/recursive_horizons/nsc_incoming_cauchy_jets.py),
the [fixed reference chart](nsc-pg-ks-metric-pullback.md), the
[short-shell response](nsc-retarded-radial-response.md), and the
[Fourier matching equations](nsc-incoming-fourier-matching.md). The
[diagonal sign certificate](nsc-retarded-diagonal-sign-bound.md) concerns
the original compact radius direction. The new connection is whether an
upstream change of the same metric can compensate that preparation response
without changing its incoming normal data or the source law.

At the reference `N_K=1,beta_K=0`, the three raw variations give

\[
\delta H=
\left(-m\sigma_1+\frac\ell r\sigma_2-\frac{E_*}{a}\sigma_3\right)\delta N_K
+\frac{E_*}{a^2}\sigma_3\delta a_K
-\frac\ell{r^2}\sigma_2\delta r.
\]

These are derivatives of the existing Dirac Hamiltonian, not independently
adjustable matrix forces. The spatial kinetic operator uses the owned Weyl
symmetrization. Its matrix element has midpoint momentum
`k_mid=-(Ei+Eo)/2`; at zero transfer it gives exactly the displayed terms.
Compact axial-gradient terms are included by this rule, not discarded from
the operator. Here a_K is a metric field, not a refitted action coupling or
the magnetic flux.

## Common bump and transported rank

Consider real metric directions
`(delta N_K,delta a_K,delta r)=b(rho) w(z) (alpha_N,alpha_a,alpha_r)`.
The common smooth radial bump b is nonnegative, not identically zero, and
supported strictly upstream of rho1 and away from the horizon. The real
compact axial function has `w_hat(0)>0`. No particular bump is chosen.

Let `U_E(1,rho)` be reference canonical transport and let its adjoint action
on Pauli matrices be `O_E(rho) in SO(3)`. With the downstream Duhamel
orientation, the response is

\[
K_{\rm comp}(E_*)=-i\mu_b(\overline V_b\alpha)\cdot\sigma,
\qquad
\mu_b=\widehat w(0)\int\frac{b(\rho)}{a_0(\rho)}\,d\rho>0,
\]
\[
\overline V_b=
\frac{\displaystyle\int\frac{b}{a_0}\,O_E(\rho)V(\rho)\,d\rho}
{\displaystyle\int\frac{b}{a_0}\,d\rho}.
\]

The same positive scalar weight multiplies all three columns. Choose an
abstract interior center rho_c and a support halfwidth h lying in a compact
positive-metric neighborhood. In the inherited raw-field coordinates,

\[
\|\overline V_b-O_E(\rho_c)V(\rho_c)\|_2
\le h L,\qquad
L=\sup_{\rm support}
\left(2\|H_{E_*}/a_0\|_2\|V\|_2+\|V'\|_2\right).
\]

This follows from positive weighted averaging and the unitary-conjugation
bound `norm(O_E')<=2 norm(H_E*/a0)`. Smooth fixed geometry makes L finite,
while the local determinant makes `sigma_min(V(rho_c))` strictly positive.
Therefore

\[
hL<\sigma_{\min}(V(\rho_c))
\]

is a sufficient rank certificate. Such a width exists without knowing the
common transport rotation at the center: rotations preserve singular
values. Fixing any qualifying nontrivial positive bump then gives an
invertible real response matrix. A zero-integral bump or zero axial mean
would not satisfy this argument.

For `K_original(E*)=-i k_original.sigma`, the real linear equation
`mu_b*Vbar_b*alpha=-k_original` has a unique solution. This establishes the
available compensation directions; it does not assign that solution as
physical metric data. Narrower support can reduce mu_b and require larger
relative amplitudes, so rank alone gives no finite-amplitude admissibility
claim.

## Fixed clock, boundary data and implementation

The chart `T=T0(rho), z=tau+S0(rho)` remains fixed during variation.
In `G_rho=i H_K/a0`, a0 is the **reference coordinate factor**: varying raw
a_K does not differentiate this clock factor. The canonical half-density
Hamiltonian already includes its field/measure transformation; no separate
measure force is added.

Direct linear KS Fourier vertices suffice for this zero-transfer response.
A PG implementation instead requires the full rho-dependent KS-to-PG
Jacobian and the existing principal-operator variation. Raw KS lapse/axial
directions cannot simply be inserted into PG logarithmic slots or handled
by the radius-only adapter. Their interior characteristic velocities may
change. For a fixed compact direction and sufficiently small overall family
amplitude, positive N,a,r and the connected PG chart condition `F>0` persist.
This is a local qualification, not a selected finite history.

Because the parent bump vanishes on an open neighborhood of Sigma, all its
incoming metric derivatives vanish. It leaves the original radius
direction's incoming jets and the fixed restriction map unchanged. Compact
support also allows the unchanged reference past and original source law
to supply initial data. Apply the existing retarded initial-time and right
inflow support guards to the **union** of supports; extending support does
not automatically preserve a previous numerical time window or grid edge.
Initial and incident tangents vanish only where that unchanged-past/upstream
construction proves they do.

## Single-pair development condition and next test

The amplitudes and metric functions are common to **all** energies and
angular/compact channels. They cannot be chosen separately by spectral pair
or channel. Cancelling K at E* does not cancel other diagonal or off-diagonal
pairs. It establishes neither full C0 preservation nor preservation of the
incoming stress source, nonlinear preparation matching, parent constraints
or extended stationarity. The original source law and past data stay fixed;
the resulting complete incoming covariance remains an equation to solve.

This is legitimate off-shell history selection through the existing metric
operator, not a new force, action term, source-state adjustment or coupling
fit. All other equations of the declared action must still be imposed.

The smallest next test is a directed rank/width certificate for one parent
control interval, using the existing profile and derivative bounds in the
inequality above. It needs no chosen amplitudes or field propagation. Only
after that certificate should actual response columns be evaluated for a
candidate history. No metric timestep, physical initial tuple or full-goal
completion follows from this lemma.
