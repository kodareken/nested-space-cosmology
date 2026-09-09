# Which inherited Dirac state contributions survive child expansion?

For the fixed free compact Dirac sector, a defined class of homogeneous
state differences has vanishing stress density in the expanding child.
The bound is O(r^-4) for massless modes and O(r^-3) when finite massive
occupation moments are present. It uses unitary evolution and the actual
proper volume, so it does not require a simulation of every state.

There is also a direct application to the existing parent-state control:
the free lowest-Landau-level (LLL) Unruh minus Hartle-Hawking stress
difference decays as r^-4 in the interior. Its absolute common vacuum
stress is not thereby fixed. These statements do not establish the state
of the interacting recursive theory or identify its full source with
the round-sphere Euclidean calculation.

## Existing results and the specific new map

Reuse the [free conformal/compact reduction](nsc-compact-mass-map.md), the
[normalized covariance evolution](nsc-influence.md), the
[interior geometry and clocks](nsc-clock-horizon.md), and the
[computed null-state constants](nsc-horizon-source.md). Their spectra,
evolution algorithms, horizons and result generators are not repeated.

The renormalized difference is governed by the standard Dirac stress
construction: on the same geometry and in the same prescription, the
local subtraction terms cancel between Hadamard states. The remaining
stress is the point-split differential operator applied to their smooth
two-point-function difference. This is explicit in
[Dappiaggi–Hack–Pinamonti, Definition 4.1 and Theorem 4.1](https://arxiv.org/html/0904.0612).
It fixes neither the absolute reference stress nor its finite coefficients.

The new calculation binds those results to the actual anisotropic
interior volume and momentum redshift. Trace-norm conservation and
Holder's trace inequality are imported operator facts; no new universal
cosmic-no-hair theorem is claimed.

## The actual future conformal geometry

Use the existing interior coordinates:

\[
ds^2=dT^2-a_\parallel(T)^2dz^2-r(T)^2d\Omega_2^2,
\quad a_\parallel=\sqrt{-A},\quad dT=-d\rho/\sqrt{-A}.
\]

The parent static coordinate is the spatial coordinate z here. Spatial
sections are unwrapped R x S2. Define s=a_parallel/r. On rho<=0,

\[
s^2=3(\pi/2-\arctan\rho)-\frac{1+3\rho}{1+\rho^2},
\qquad (s^2)'=\frac{2(\rho-3)}{(1+\rho^2)^2}<0.
\]

Thus s_min=sqrt(3 pi/2-1)<=s<sqrt(3 pi), with increasing s toward the
future. The angularly integrated proper volume per coordinate dz is
4 pi s r^3. This is a local volume factor, not a finite global volume.
The displayed numbers use hbar=c=L_throat=1. Restoring a profile length
L gives s_min=sqrt(3 pi/2-1)/L; k and m have inverse-length units while
the angular eigenvalue kappa is dimensionless. The mode measure includes
longitudinal momentum density, so the bound has energy-density units.

At u=-1/rho approaching zero from above,

\[
W(u)=s^2=3\pi-u^2-2u^3+O(u^4),\qquad
\frac{g}{r^2}=\frac{du^2}{(1+u^2)^2W(u)}-W(u)dz^2-d\Omega_2^2.
\]

The conformal metric is nondegenerate at u=0. The remaining conformal
time obeys 0<eta_future-eta(u)<=arctan(u)/s_min. No new reflecting
condition is imposed on this future boundary.

For the massless Dirac field, conformal covariance gives the spinor
weight r^-3/2. Between states in the same prescription, local anomaly
terms cancel; orthonormal stress differences transform with r^-4.
Bounded stress of the conformally transported difference therefore gives
physical O(r^-4) decay. Smoothness at each interior point alone does not
prove a uniform bound at conformal infinity. The mode estimate below
states a sufficient condition explicitly and also covers fixed masses.

## A state bound that retains all free angular and compact modes

For the fixed free tower, m_n=n pi/ell. The prior chiral-domain choice
and physical angular sectors are retained. The interior Dirac expansion
connection is H_parallel/2+H_sphere; the canonical field is
chi=(a_parallel r^2)^(1/2) psi. Each channel then has the standard form

\[
H_j=m_j\beta+\frac{k_j}{a_\parallel}\alpha_z
                 +\frac{\kappa_j}{r}\alpha_\Omega,\qquad
\|H_j\|=\sqrt{m_j^2+k_j^2/a_\parallel^2+\kappa_j^2/r^2}.
\]

The three Hermitian matrices anticommute. This notation is for the
actual spinor channel decomposition; physical degeneracies belong in
the mode measure once, and a channel label is not itself a chirality.

Compare two covariances evolving under the same closed Hamiltonian:
Delta C_j(T)=U_j(T) Delta C_j(T0) U_j(T)^dagger. Their trace norm
N_j=||Delta C_j||_1 is constant, even with noncommuting time-dependent H.
The hypothesis is finite weighted moments

\[
M=\sum_j\!\int d\nu_j\,|m_j|N_j<\infty,\qquad
K=\sum_j\!\int d\nu_j\,
 (|k_j|/s_{\min}+|\kappa_j|)N_j<\infty.
\]

The notation includes continuous longitudinal momentum and the compact
and angular multiplicities. It refers to a homogeneous density or the
specified spatially averaged density. Arbitrary inhomogeneous pointwise
stress requires additional spatial estimates.

Holder's inequality and a_parallel>=s_min r give

\[
\boxed{|\Delta\rho(T)|\le
 \frac{M+K/r(T)}{4\pi s_{\min}r(T)^3}.}
\]

The corresponding pressure bound is
|Delta p_parallel|+2|Delta p_sphere|<=K/(4 pi s_min r^4).
If M=0, the energy bound is r^-4. For finite M it has an r^-3 envelope.
No claim that every state saturates these rates is made. The finite
moments justify taking the limit through this state-difference sum;
they do not bound a vacuum zero-point sum.

This is compatible with unitarity: N_j need not decay. Expansion reduces
the contribution to local stress through momentum redshift and proper
volume growth; it does not erase the full state information.

## The source and geometric work obey the same energy ledger

For each fixed-mass channel,

\[
\dot H_j=-H_\parallel\frac{k_j}{a_\parallel}\alpha_z
         -H_\perp\frac{\kappa_j}{r}\alpha_\Omega.
\]

The commutator term in d tr(H Delta C)/dT vanishes by cyclicity.
The same metric derivatives define the two pressure components, giving

\[
\Delta\dot\rho+(H_\parallel+2H_\perp)\Delta\rho
 +H_\parallel\Delta p_\parallel+2H_\perp\Delta p_\perp=0.
\]

If compact size or a link changes the mass, m_dot tr(beta Delta C) is
additional work and must be retained. An open reservoir, changing domain,
interacting source or evolving compact geometry is outside this fixed
free bound. A vanishing free state-difference flux does not imply that
acceleration ends: a persistent common vacuum stress remains possible.

## A direct comparison of the already specified parent states

For the existing null coordinates u=z-x_star, v=z+x_star in the trapped
region, partial_T u=-1/a_parallel and partial_T v=1/a_parallel.
Only the constants t_u,t_v change between the two free LLL states.
Their common geometric term cancels. In the interior orthonormal frame,

\[
\Delta\rho=\Delta p_\parallel
 =\frac{\Delta t_u+\Delta t_v}{4\pi a_\parallel^2r^2},\quad
\Delta p_\perp=0,\quad
\Delta T_{\hat T\hat z}
 =\frac{\Delta t_v-\Delta t_u}{4\pi a_\parallel^2r^2}.
\]

For Unruh minus Hartle-Hawking, per unit |q|,
Delta t_u=0 and Delta t_v=-kappa_h^2/(48 pi), using the recorded horizon
surface gravity. Both density and this covariant flux component decay
as r^-4. The asymptotic r^4 density is
-kappa_h^2/(576 pi^3). This directly resolves the relative contribution
of those two existing state choices in the free projected sector. It
does not identify either one's complete four-dimensional vacuum stress.

The same stored angular spectrum gives the first nonzero angular gap
kappa_1^2=(1+|q|)/r^2. At fixed finite q its ratio to H_child^2 tends
to zero. A lowest-mode-only approximation therefore cannot remain a
uniformly separated description of the late child. Full vacuum matching
needs the angular sum; taking the late limit of each vacuum channel
before summing is not justified by the relative-state bound.

## What remains to identify the absolute endpoint state

The [round-sphere response](nsc-spectral-endpoint.md) remains an
authenticated input. To use it as the actual late source requires a
Hadamard reference state on the actual domain whose absolute stress
approaches that same-regulator invariant endpoint, and control of the
parent/reference difference. Neither a future instantaneous ground
projector nor Hadamard regularity alone supplies that construction.

Grok checked the scope of the closest primary state constructions.
[Dappiaggi–Hack–Pinamonti's spinor boundary construction](https://arxiv.org/abs/1009.5179)
addresses spatially flat FRW, while the previously identified
[Allen–Lutken reference](https://pure.mpg.de/view/item_153497) is maximally
symmetric. Their applicability to this anisotropic R x S2 interior is
an additional mathematical step. Scalar-only results cannot be substituted.

The canonical free-tower bound does not by itself establish the Lorentzian
continuation of the finite proper-time functional, the global transmission
state, the interacting gauge/metric vacuum, or the completed quantum
measure. It supplies an explicit class in which inherited excitations
cannot create a persistent stress correction, and identifies the angular
and reference-state work still required for the common source.

Evidence: [record](../results/development/child-state.json), reproduced by
scripts/check_nsc_child_state.py --check. The symbolic chart, canonical
work and null-state maps and all numerical bound fields are compared.
No previous scientific generator or endpoint spectral sum is run.
A bounded read-only Grok review independently checked the volume
normalization, pressure split, null-state signs and stated limit conditions;
its assessment is additional review evidence, not a substitute for the
explicit derivation or the unresolved absolute reference construction.
