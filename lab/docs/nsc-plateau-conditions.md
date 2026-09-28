# A dark-sector plateau as a dynamical requirement

The [reproducible calculation](../scripts/check_nsc_plateau_conditions.py) and
[record](../results/nsc-5-plateau-conditions.json) turn the proposed late-time
plateau into equations that the common operator must satisfy. They are
requirements and analytic controls, not a derived NSC cosmological history.

The baseline is flat FLRW with constant Newton coupling, positive densities,
negligible late-time radiation and a conserved total effective stress. Visible
means fixed-mass baryonic dust. Dark means the combined unresolved response,
including both the effects ordinarily assigned to dark matter and dark energy.
This effective stress decomposition does not postulate a new dark substance.
Asymptotic flattening here means curvature invariants tend to zero.

## Conservation determines what a plateau needs

Write `N=ln(a)`, with prime denoting differentiation in N. Adopt the transfer
sign convention

\[
\dot\rho_b+3H\rho_b=Q,\qquad
\dot\rho_d+3H(\rho_d+p_d)=-Q.
\]

Thus positive Q transfers energy into visible matter. Set

\[
r=\rho_d/\rho_b,\quad f=r/(1+r),\quad
w_d=p_d/\rho_d,\quad \epsilon=Q/[H(\rho_b+\rho_d)].
\]

Directly differentiating the two densities gives

\[
\boxed{r'=-3w_dr-(1+r)^2\epsilon},\qquad
\boxed{f'=-3w_df(1-f)-\epsilon
=(1-f)(1-2q_{\rm dec})-\epsilon},
\]

where `q_dec=-1-d(ln H)/dN=(1+3w_d f)/2`. The last equality uses the flat
Einstein background equations, constant G and the stated matter accounting.
It cannot be carried unchanged into a different gravity or transfer convention.

For a finite positive ratio fixed point,

\[
\epsilon_*=-\frac{3w_{d*}r_*}{(1+r_*)^2}
=(1-f_*)(1-2q_{{\rm dec}*}).
\]

A scalar autonomous equation `r'=F(r)` is locally attracting when `F'(r*)<0`.
If the stress depends on geometry, boundary state or other evolving variables,
the complete coupled physical Jacobian must instead be tested. A trajectory
approaching a limit is not by itself a demonstration of attraction.

For separately conserved baryons, Q=0, the exact solution is

\[
r(N)=r(0)\exp\!\left[-3\int_0^N w_d(s)\,ds\right].
\]

A finite positive limiting ratio requires the integral to converge to a finite
value. In a regular asymptotic dynamical state this implies `w_d -> 0`, hence
`w_total -> 0` and **`q_dec -> 1/2`**. A flat expanding Einstein solution then
has `H^2 proportional to a^(-3)`, `a(t) proportional to t^(2/3)` and curvature
tending to zero. Current acceleration can therefore be followed by a finite
ratio plateau and asymptotic flattening, provided acceleration ends.

For autonomous `w_d=w_d(r)` with Q=0, the positive root obeys `w_d(r*)=0` and
is locally stable when `w_d'(r*)>0`. Smooth asymptotic evolution is an explicit
assumption: convergence of a function alone does not guarantee convergence of
its derivative.

Maintaining acceleration at a finite fraction instead requires a derived
positive transfer `epsilon*=(1-f*)(1-2q_dec*)`. This is a necessary condition,
not a prescription to add a transfer term. The common action must provide its
physical origin and conservation law. With nonzero transfer an accelerating
power-law expansion can also have H and curvature tending to zero.

## One present snapshot permits different futures

The rounded values `rho_b0=0.05`, `rho_d0=0.95` and `w_d0=-0.7/0.95=-14/19`
are illustrative inputs in units `H0=1`, `8*pi*G/3=1`. They are not a new
measurement, an inferred NSC coefficient, or a predicted ratio of 19.

Use the future control histories

\[
w_d(N)=-\frac{14}{19}e^{-kN},\quad Q=0,
\qquad k\in\{1,4,20\},\quad N\ge0.
\]

Conservation fixes

\[
\rho_d(N)=\frac{19}{20}
\exp\!\left[-3N+\frac{42}{19k}(1-e^{-kN})\right],
\]

\[
\boxed{r_*=19\exp\!\left(\frac{42}{19k}\right)},\qquad
\boxed{f_*=\frac{r_*}{1+r_*}}.
\]

All three have identical present H, dark fraction and `q_dec=-0.55`. They
approach different plateau values at different rates. The runner determines
when acceleration ends, checks the analytic densities against independent
quadrature, and reconstructs deceleration by differentiating H independently
of the pressure expression. These prescribed time histories are controls;
they do not supply the perturbation law needed to establish an attractor.

The tail is controlled exactly:

\[
\log\frac{r_*}{r(N)}=\frac{42}{19k}e^{-kN},\qquad
0\le f_*-f(N)\le\frac14\log\frac{r_*}{r(N)}.
\]

The bound follows from the logistic derivative `f(1-f)<=1/4`. It illustrates
how an indefinitely continuing expansion can yield a finite observable with
an explicit remainder. An infinite decimal expansion or an infinite number of
recursive steps is not itself a convergence proof; the particular sum,
integral or response needs a limiting argument and an error bound.

If Q=0 and `w_d<=0` throughout the future, f cannot decrease. Current
acceleration over a nonzero interval therefore places any final finite plateau
above today's fraction. A plateau exactly equal to today's rounded 0.95 would
need a later decrease, transfer or changed assumptions. Neither the current
fraction nor an observational fit selects a topological resistance value.

For comparison, `w_d=-1/(1+N)` tends to zero but has
`integral w_d dN=-ln(1+N)`, so `r=19(1+N)^3` diverges and `f -> 1`. Its
deceleration tends to 1/2 and curvature still flattens. Decaying negative
pressure alone does not establish a finite ratio.

## Recursion, fraction divergence and curvature are separate tests

At `Omega=1`, `b=1`, `K(i)=i`, the scalar recursion is

\[
\Gamma=i-\Gamma^{-1},\quad
\Gamma_*=i\varphi,\quad \varphi=(1+\sqrt5)/2.
\]

Its upper-half-plane branch is locally contracting under room-depth iteration,
with derivative `-1/varphi^2`. The same equation contains no N or cosmic-time
stress law. It can coexist as an algebraic equation with the independent
history `rho_b=exp(-3N)`, `rho_d=1`, whose ratio diverges. This logical
countermodel disproves an implication from the Schur fixed point alone; it is
not a proposed solution of the complete self-sourcing NSC action.

The standard flat dust-plus-positive-Lambda control further distinguishes
ratio growth from a singularity. With rounded current fractions `(0.05, 0.25,
0.70)` for baryons, cold dark matter and vacuum,

\[
r=5+14e^{3N}\longrightarrow\infty,\qquad f\longrightarrow1.
\]

Using signature `(-,+,+,+)` and `R=6(dot H+2H^2)`, its future limits are

\[
H^2/H_0^2\to0.7,\quad R/H_0^2\to42/5,\quad
R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}/H_0^4\to294/25.
\]

They are finite. The density ratio diverges because the baryon denominator
dilutes. Its fraction is bounded, and this asymptote has no curvature blow-up.

## Matching test for the eventual common operator

The next physical calculation must supply `H(a)`, the full state-defined
metric-variation stress and any Q from the same action. The spectral link b,
off-diagonal Phi and recursion eigenvalues cannot simply be labeled as
`rho_dark`, pressure or energy transfer.

For the conserved-baryon baseline, independently measured baryon abundance and
the predicted expansion provide the audit

\[
r_H=\frac{E(a)^2a^3}{\Omega_{b0}}-1,\quad
f_H=1-\frac{\Omega_{b0}a^{-3}}{E(a)^2},\quad E=H/H_0,
\]

\[
w_{{\rm total},H}=-1-\frac23\frac{d\log H}{dN},\qquad
w_{d,H}=w_{{\rm total},H}/f_H.
\]

A finite positive plateau requires
`E^2*a^3 -> Omega_b0/(1-f*)` and `d(ln H)/dN -> -3/2` on the regular branch.
Compare these reconstructions with the stress outputs, both conservation
equations and the Friedmann/Raychaudhuri residuals. The same frozen operator
must match current expansion and the combined fraction, determine future
curvature, and pass physical perturbation tests. This audit can be specified
before the solution exists; it does not supply the missing solution.

Reproduce with `python3 scripts/check_nsc_plateau_conditions.py --check`.
Every record field and both source-file hashes are compared. The script uses
exclusive creation for `--output`; no historical result is overwritten.
