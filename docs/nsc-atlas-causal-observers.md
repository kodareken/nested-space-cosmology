# Centre, moving wall and parent: event-specific causal access

This successor asks which named events can exchange light in the isolated
reference associated with the [explicit connection slice](nsc-atlas-connection-slice.md).
Its answer must retain the observers' different domains and the time of
emission. Initial-data attachment is not by itself a communication channel.

## Admission before calculation

- **Borrowed:** the signed vacuum-wall law from
  [0711.2198v1](https://arxiv.org/html/0711.2198v1), Section 2 equation (2),
  and the patch distinctions in
  [gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2), Section II.
  Spradlin, Strominger and Volovich,
  [hep-th/0110007v2](https://arxiv.org/html/hep-th/0110007v2), Section 2.1,
  equations (15)--(21), supplies the conformal and planar charts. Their
  displayed planar metric contracts; here its time reverse is chosen
  explicitly. This coordinate geometry is not a new horizon mechanism.
- **Missing:** the v4 source-compatible slice did not give emission- and
  reception-specific access for its centre, moving wall and right parent.
- **Decision/unlock:** identify allowed historical light paths and distinguish
  them from future inaccessible paths. This supplies observer labels for the
  atlas and prevents a spatial bridge from being promoted to a signal link.
- **Fixed theory and events:** H ell=1, k ell=2, M/ell=0.4, expanding wall.
  The initial wall event W0 is R/ell=1.5 with flat de Sitter time t=0 and
  scale factor one. C0 is the centre at the same flat time. The right parent
  observer is the separate Schwarzschild R=4M worldline.
- **Finite computation:** R/ell in [0.05,3], with finitely many samples and
  paired 32/64-point quadratures, less than 20,000 scalar nodes and 60 CPU
  seconds. A horizon is crossed with regular boost variables, never by
  integrating a singular static time coordinate. No mass sweep or PDE run.
- **Stop/check:** stop on non-real or incorrectly oriented motion, failed
  boost normalization, inconsistent null propagation, or a historical
  emission preceding the retained wall interval. Future no-access statements
  require an analytic inequality valid beyond the plotted interval; a grid
  alone cannot certify a complete future light cone.
- **Attribution and independent checks:** inherited wall tests, two regular
  boost formulas, direct null-coordinate equalities, monotonic characteristic
  intercepts, quadrature-order agreement, unit rescaling and domain exits.
  This is an observer map for borrowed solutions, not a priority claim.

## Flat de Sitter chart and regular wall clock

Use the expanding flat chart

$$
ds^2=-dt^2+e^{2Ht}(dr^2+r^2d\Omega^2)
      =e^{2Ht}(-d\eta^2+dr^2+r^2d\Omega^2),
\qquad \eta=-e^{-Ht}/H<0,\qquad R=e^{Ht}r.
$$

The centre r=0 has proper time t. This chart's future conformal boundary is
eta=0. The wall's proper time is a separate clock. Let s=dot R>0 and retain
the signed beta_- from the imported junction. Its regular boost relative to
the flat slice can be evaluated without division by the static horizon:

$$
z={\beta_-+s\over 1+HR},\qquad
\gamma={z+z^{-1}\over2},\qquad
\gamma v={z-z^{-1}\over2},\qquad
{dt\over dR}={\gamma\over s},\quad
{d\tau_W\over dR}={1\over s}.
$$

An independent rationalized expression is
gamma=(1+s^2)/(beta_-+HR s). Both agree with the Lorentz relations
s=gamma(HR+v) and beta_-=gamma(1+HR v). Thus the de Sitter static horizon
R=1/H is a regular event in this chart.

## Null intercepts distinguish historical and future access

At an emission from the wall, an ingoing radial light ray reaches the centre
only if its intercept eta_W+r_W is negative. If it is zero, arrival would
lie at infinite centre proper time; a positive value lies beyond the chart's
future boundary. At W0, eta_W=-ell and r_W=1.5 ell, so the intercept is
0.5 ell: no finite future centre reception is possible in this expanding
de Sitter reference.

Conversely, a centre emission reaching W0 has

$$
\eta_C=\eta_W-r_W=-2.5\ell,\qquad
t_C=-\ell\log(2.5).
$$

This is a historical emission, not C0. It counts as an actual reference path
only after its entire segment is shown to lie within the retained matched
child. For a timelike future wall, d(eta_W-r_W)/dt=(1-v)e^{-Ht}>0. Once the
record verifies that t_C lies after the retained earliest wall event, this
strict monotonicity proves that the outgoing ray stays inside the wall until
its unique reception at W0. No extrapolation through the wall's birth is
needed for that finite path.

An additional earlier wall event at R/ell=0.5 is below the de Sitter horizon.
Its ingoing ray reaches the centre at finite proper time. An even earlier
centre emission reaches this same wall event on its outgoing null ray. Both
paths are checked against the retained beginning at R/ell=0.05. This earlier
cut is needed to cover the latter emission; the calculation never continues
to the wall's zero-radius past boundary.

There is also a distinguished outgoing centre ray that meets the wall at
R=2M=0.8 ell. Its constant eta-r equals eta_H-r_H, and its centre emission
time is t_H-ell log(1.8). The adapter checks that this entire segment is
within the retained history. The companion Schwarzschild reference identifies
the wall crossing with U_K=0: earlier wall crossings have U_K<0 and can reach
the right static parent, while U_K>=0 cannot reach that observer at finite
time. Consequently the centre emission threshold is strict. The equality
ray continues along the reference event horizon, rather than predicting a
finite-time parent reception. This combined interpretation applies only to
the isolated completion and uses its separately checked exterior signs.

For R>=1.5 ell on this fixed expanding branch, beta_- is increasing and is
already greater than one. Together with HR>1 this gives gamma v>0. Hence the
wall's comoving radius never falls below its value 1.5 ell at W0. A light ray
emitted from C0 has at most one ell of remaining comoving radius before
eta=0, and cannot catch the wall. Future wall emissions also have
eta_W+r_W=(HR-1)e^{-Ht}/H>0, so cannot reach the centre. These inequalities
hold on the entire future branch, independently of the finite plotting
cutoff. The branch has s>0 for all R>=1.5 ell and s<2R/ell there, so its
proper time to unbounded radius diverges. The absence of a future arrival
is not caused by stopping the numerical trajectory early.

The right-parent observer is handled in the companion reference's regular
Schwarzschild null coordinates. The de Sitter calculation alone does not
decide access from or to that observer. In particular past-parent emission
and future-parent emission are distinct questions, and an isolated reference
completion is not automatically the development of a glued many-body parent.

The [adapter](../src/recursive_horizons/nsc_atlas_causal_observers.py) and
[tests](../tests/test_nsc_atlas_causal_observers.py) own these finite data and
the explicitly stated analytic extension.

## Finite result and remaining boundary

The retained wall history has $t/\ell$ in
$[-1.5439442061,\,0.6185152317]$. At the early wall event $R/\ell=0.5$,
$t_W/\ell=-1.1187517160$ and $\tau_W/\ell=-1.1138558046$. Its outgoing
centre witness was emitted at $t_C/\ell=-1.5242168241$, which is
$0.0197273820\ell$ after the retained beginning. Its inward ray reaches
the centre at $t_C/\ell=-0.4256045355$. These are explicit events on one
reference history, with separately named local clocks.

At $R/\ell=0.8$, the wall has $t_W/\ell=-0.6435077876$,
$\eta_W/\ell=-1.9031450127$ and $r_W/\ell=1.5225160102$. The outgoing
parent escape horizon therefore starts from the centre at
$t_C/\ell=-1.2312944525$, with constant
$(\eta-r)/\ell=-3.4256610229$. Its central emission is
$0.3126497537\ell$ after the retained beginning. At the distinct de Sitter
observer horizon crossing $R/\ell=1$, the wall has
$t_W/\ell=-0.3946760303$ and $\eta_W=-r_W=-1.4839033733\ell$.

The paired clock integrations use 11,616 quadrature evaluations and agree
to floating-point precision. This comparison is an integration diagnostic,
not a rigorous interval error bound. Direct Lorentz normalization has
maximum residual $2.28\times10^{-13}$; the null intercept identities have
maximum residual $4.45\times10^{-16}$. All 28 focused tests pass in 0.53
seconds, including the finite-history margins, horizons, signed junction,
future analytic bounds, clock derivatives and unit rescaling. The Windows
test runner reports that its POSIX CPU guard is inactive; the explicit
bounded node count and measured short run are the execution evidence.

Reading coverage for this addition is the cited de Sitter chart section,
the previously read signed wall equations, and the companion
[classification review](nsc-atlas-causal-classification.md). The parent
reference and classification own the continuation of the escape ray into
Schwarzschild coordinates. No new horizon theorem or priority is claimed.
Field transmission through the wall, emission amplitudes, regular formation,
perturbative stability and the future development of the glued many-body
data remain outside this calculation. No result here supplies a physical
signal between siblings or a common age for different universes.
