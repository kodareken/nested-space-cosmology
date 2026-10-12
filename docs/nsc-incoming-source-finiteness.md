# Finite retained source and the local constraint relation

For the exact source defined by the current retained channel prescription,
the incoming stress moments are finite. This concerns the mathematical
source; it does not certify the errors of the archived finite-collar modes,
quadratures or their floating evaluation.

The distinction matters for the established local constraint plane:

\[
\mathcal E_N=S_N+c_u u+c_{v^2}v^2,\qquad
\mathcal E_\beta=S_\beta+c_vv.
\]

Both response coefficients are now positively enclosed. Within the
established retained incoming bulk functional, finite real source constants
therefore give a unique finite point-germ pair

\[
v_*=-S_\beta/c_v,\qquad
u_*=-\bigl(S_N+c_{v^2}v_*^2\bigr)/c_u.
\]

Knowing a numerical approximation to each source constant to `3e-11` is not
a logical prerequisite for this algebraic statement. Such accuracy, with
coefficient and evaluation errors propagated, is required to certify a
reported numerical root to that tolerance. Neither statement supplies a
Cauchy-surface solution or parent-compatible physical initial data.

## Reuse and decision

Reuse the exact fermionic source fibers and current-normalized physical
sewing, canonical unitary interior transport, the existing incoming vertices,
gapped fourth-order reference, analytic LLL allocation and vacuum-tail bound.
The missing connection is whether unfinished numerical source errors also
leave the source's mathematical finiteness unestablished. Combine these
existing bounds, without new modes or quadrature. Stop after establishing
finite real `S_N,S_beta,c_v2` in the declared retained model and identifying
the separate coefficient and surface-solution requirements.

## Finite energy intervals

For the positive-energy folded fibers used by the incoming source, the owned
[source covariance](../src/recursive_horizons/nsc_paired_horizon_preparation.py)
has a two-dimensional horizon block with eigenvalues0 and1, since
`s^2=f(1-f)`, and an independent incoming occupation `0<=n<=1`.
On a closed channel the incoming occupation vanishes and its source
projection removes that column. Thus the exact source fiber obeys
`0<=C_source<=P_source` in both open and closed sectors.

The [physical sewing](../src/recursive_horizons/nsc_pg_massive_modes.py)
uses the exact flux identity `abs(R)^2+T=1` and the appropriate open/closed
projection. On the common canonical incoming frame it is a coisometry;
the exact interior Dirac evolution is unitary. Consequently

\[
0\le C_g(E)\le I_2.
\]

This is a property of the defined exact physical mode prescription. A
sampled coisometry residual or ODE tolerance does not prove that its numerical
approximation has a specified error. No numerical columns are normalized
to obtain this inequality.

At the fixed incoming surface, each non-LLL group has gap

\[
\omega_g(E)^2=m_g^2+(\ell_g/r)^2+(E/a)^2
\ge M_g^2=m_g^2+(\ell_g/r)^2>0.
\]

The finite retained inventory contains32 such groups. The existing formal
reference through order4 is a smooth matrix function of real energy on any
compact interval including zero, with gap denominators bounded away from
zero. Let `B_g=max_[0,L] norm(P_ad4,g)`; this maximum is finite. This does
not treat that formal subtraction as a positive covariance.

The [incoming vertices and measure](../src/recursive_horizons/nsc_incoming_state_moments.py)
are

\[
H=-m\sigma_1+(\ell/r)\sigma_2-(E/a)\sigma_3,
\quad V_\parallel=-(E/a)\sigma_3,
\quad V_{01}=(E/a)I,
\quad V_\perp=\ell\sigma_2/(2r),
\]

with the fixed finite multiplicity factor per actual angular sign. Hence
each vertex has a finite maximum `V_g,j` on `[0,L]`, and

\[
\int_0^L
\left|\operatorname{tr}[(C_g-P_{\rm ad4,g})V_j]\right|\,dE
\le 2L(1+B_g)V_{g,j}<\infty.
\]

The source fibers and scattering modes are measurable; threshold changes
and closed-channel reflection phases cannot violate this norm estimate.
The local measure is the owned constant multiple of `dE`, not an extra
inverse asymptotic group velocity. No small-energy numerical refinement is
needed to prove this finite-interval statement. Those refinements remain
necessary to determine the integral accurately.

## Infinite tails and the LLL

The [existing leading matching identity](../src/recursive_horizons/nsc_incoming_source_tail.py)
and [analytic tail certificate](nsc-incoming-tail-quadrature-bound.md)
give `source(x)=x^2*g(x)`, with analytic `g(0)=0` and `x=1/E`, for the
owned order16-minus-ad4 approximation. Its kernel is therefore `O(E^-3)`.
The exact-state difference is integrable by the already certified
[vacuum-mode remainder](nsc-incoming-vacuum-tail-bound.md) and the positive
exponentially decreasing thermal bound. Their combination proves absolute
tail integrability; a numerical tail-integral error certificate is a
different requirement. The matching identity is algebraic in the general
channel parameters, so this argument is not limited to the group22
quadrature pilot.

The massless LLL is excluded from the gapped argument. Its finite physical
state moments use the separate
[analytic allocation](../src/recursive_horizons/nsc_incoming_source_assembly.py):
with `t_u=abs(q)*kappa^2/(48*pi)` and `t_v=Omega^2*t_u`. The state vector is
`(t_u+t_v,t_u+t_v,t_v-t_u,0)/(4*pi*r^2*a^2)`; all its parameters are finite
and `a,r>0`. It is already included once by
[baseline_matter_source](../src/recursive_horizons/nsc_incoming_joint_constraints.py).
Its geometric terms remain with the local action. Summing the32 gapped
groups and that LLL allocation is finite; this is not a convergence theorem
for an undeclared infinite angular or compact complement.

## Consequence and remaining boundary

The local action density and Euler coefficients are finite at the fixed
positive intrinsic `N,a,r`. On the
[restricted plane](nsc-incoming-constraint-germ.md), the changed reference
has a finite gap and an absolutely integrable momentum tail. Its exact
polynomial dependence on `u,v` therefore has finite real coefficients,
including `c_v2`. Thus `S_N,S_beta,c_v2` are finite real quantities for the
exact retained prescription, despite unresolved errors in their numerical
approximations.

The [lapse coefficient](nsc-incoming-lapse-coefficient.md) lies in
`[0.05739927758627631,0.057399277586276436]`, and the
[shift coefficient](nsc-incoming-shift-coefficient.md) lies in
`[16.575975622058497,16.575975622058518]`. Both are strictly positive, so the
Jacobian determinant `c_u*c_v` never vanishes on this plane. The finite
source constants and finite `c_v2` establish the unique algebraic zero above.
No numerical value of that pair is asserted: source and `c_v2` error bounds
are still needed to certify a reported numerical solution. This result
uses the explicit fixed-C0 experiment and included retained operators;
it does not assign any remaining `Gamma_rest` term or establish full
functional stationarity.

A surface solution additionally needs spatial normal-data functions and
their compatible mixed derivatives, a stated spatial domain and boundary
or asymptotic conditions, and constraints between evaluation points. The
constant-jet point plane cannot simply be repeated at every spatial point:
for example, `r_Tz` is the spatial derivative of the function `r_T(z)`.
Parent-state compatibility, transmitting endpoint variation and extended
stationarity remain independent requirements. No source value, physical
initial datum, action term, elapsed duration or metric timestep is selected
by this finiteness argument.
