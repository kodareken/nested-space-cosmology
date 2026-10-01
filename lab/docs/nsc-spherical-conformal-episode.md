# Same-source conformal episode at $T=0.05$

The optional same-action formulation in
[the conformal-gauge proof](nsc-spherical-conformal-gauge.md) reaches
$T=0.05$ on the unchanged saved $T=0$ preparations at $n_f=256,512$ and
$dt=0.001,0.0005$. Every stage evolves the same geometry and Gaussian columns
through the existing Galerkin coupling. The gauge is $L=Q$, $\beta=0$,
with the derived $-C$ chain-rule term in $\dot p_Q$ and the complete
coordinate work $F_Q\dot Q+F_L\dot L$. The two sealed result owners are
[episode-v1](../results/development/nsc-spherical-conformal-episode-v1.json)
and [frames-v1](../results/development/nsc-spherical-conformal-frames-v1.json),
with their corresponding NPZ payloads.

The finite measured branch remains localized while carrying resolved internal
energy and probability currents. Exchange through the four outer regional
boundaries is unresolved over this short window. That distinction determines
what this episode demonstrates. It does not establish renewal, a continuum
Cauchy solution, or a total error bound for the physical observables.

## Measurement domain

The same initial $Q,r,\chi$, momenta, occupations and columns are loaded
bitwise from the v5 payload. The original $T=0$ region-0 pair remains the
fixed mode observer. The frozen control uses the exact constant-initial-$Q$
Dirac group on the same Gaussian columns and retained band. It is a control
for that initial geometry; it does not inject a different source.

Four width-2 windows are $[0,2],[2,4],[4,6],[6,8]$. Their contents are exact
integrals of each periodic nodal trigonometric interpolant, using the
spectral interval integration weights. Normal shell energy is $F_L/r$.
The regional ledger uses the active $L=Q$ system, the actual projected
lifted rate, normal boundary flux, pressure work and momentum-lapse-gradient
work. The numerical directional derivative varies both the Cauchy state and
$L=Q$. Its residual is retained. The source $F_L$ is independent of $L$, so
the direct variation of this normal shell has no omitted $F_L\dot L$ term;
the coordinate field energy does contain that term.

The normal clock is anchored at $x=1$ and obeys $d\tau=rQ\,dt$. Its frozen
control has the same initial clock rate. The equal-proper-clock comparison
interpolates both curves on their own per-step clocks over their common
clock interval. It is not a comparison at equal coordinate endpoints, and
its interpolation error is an indicator rather than a certified enclosure.
The original mode observer has finite spatial support; this clock names its
reference protocol rather than identifying the mode with a single worldline.

## Measured effects and refinement

The fine $dt=0.0005$ case gives the following. Relative indicators use the
measured change or control separation, and retain its own physical domain.

| Quantity | Fine measurement | Time indicator | Space indicator |
|---|---:|---:|---:|
| Field energy change | $-1.41447\times10^{-7}$ | $6.7445\times10^{-11}$, about $0.0477\%$ | $1.42\times10^{-14}$ |
| Original observer occupation change, each mode | $-0.0842690387$ | $2.74\times10^{-11}$ | $2.80\times10^{-9}$ |
| Coupled minus frozen occupation at coordinate $T=0.05$ | about $4.1447\times10^{-9}$ | about $2.74\times10^{-11}$, $0.661\%$ | below $10^{-15}$ for this paired comparison |
| Coupled minus frozen occupation at equal normal clock | about $-2.16786\times10^{-6}$ | $8.21\times10^{-9}$, $0.379\%$ | about $1.14\times10^{-11}$ |
| Region-0 normal shell change | $+0.0027172381$ | $3.11\times10^{-11}$ | $3.18\times10^{-11}$ |
| Maximum internal normal energy flux in region 0 | $3.6359135$ | $2.15\times10^{-9}$ | $2.64\times10^{-4}$, $0.00727\%$ |
| Maximum internal canonical probability flux in region 0 | $1.1084188$ | about $1.96\times10^{-10}$ | $2.35\times10^{-5}$, $0.00212\%$ |

The field and gravity energy changes cancel to about $2.20\times10^{-12}$.
The field-energy work integral and its actual energy change agree within the
reported temporal quadrature indicators. Most of the original mode-occupation
change is shared with the frozen control; the coupling separation is a
smaller measured effect. Neither a large drift nor a small control separation
is by itself renewal.

The region-0 normal-shell leader remains the same on every saved frame. Its
share stays above about $0.66103$. The original region-0 probability content
stays at about $2$ of the total $3$, so its width-2 share remains $2/3$.
The circular concentration changes from about $0.76970948$ to $0.76911609$;
the circular location stays near $1.635703405$. These are finite-window
measurements, with their sampled refinement indicators. No minimum leader
share or reversal was used as a stop requirement.

Pressure and lapse-gradient work account for the shell rise. On the fine
case their frame-Simpson integrals are about $0.0016122863$ and
$0.0011049530$, totaling $0.00271723937$. The difference from the shell
change is about $1.26\times10^{-9}$. Their combined Simpson-versus-trapezoid
indicator is about $7.94\times10^{-7}$, or $0.0292\%$ of the shell change.
The normal continuity residual remains reported, including projection,
quadrature and the numerical directional derivative; it is not zeroed.

The **outer** boundary flux has a different outcome. The coarse region-0
maximum net boundary flux is about $2.53\times10^{-9}$, while the fine value
falls to about $1.57\times10^{-13}$, comparable with the boundary identity
roundoff. Its integrated fine net flux is about $-2.82\times10^{-16}$.
This does not resolve transfer between those width-2 regions. The robust
internal-current maxima refer to the interior profile and cannot replace
that missing boundary exchange.

## Constraint forcing domain

For $h_c=QC$ and $D$, the continuum same-action formulation has
$\dot h_c=D_x$, $\dot D=(h_c)_x$. The actual Galerkin ODE has the explicitly
differentiated forcing $f=(\dot h_c-\mathsf D_xD,\dot D-\mathsf D_xh_c)$.
The skew-adjoint periodic derivative gives

$$
R(t)\le R(0)+\int_0^t\|f(s)\|_2\,ds,
\qquad R=\|(h_c,D)\|_2.
$$

The fine endpoint has $R\simeq2.09543\times10^{-6}$, with full $C,D$ norms
about $7.94571\times10^{-6}$ and $6.35330\times10^{-7}$. The sampled actual
forcing integral is about $7.30987\times10^{-5}$, leaving a positive sampled
budget margin of about $7.30987\times10^{-5}$. The unprojected quadrature
forcing integral is about $3.55\times10^{-9}$; the dominant contribution is
the actual rate-projection defect. The fine timestep movement of the forcing
integral is about $6.80\times10^{-9}$, while its frame-free trapezoid/Simpson
indicator is about $1.17\times10^{-10}$.

That forcing integral increases from the coarse value about
$1.68498\times10^{-5}$, while the actual constraint norm falls substantially
with refinement. It is a valid finite ODE control primitive with a measured
forcing term, not a claim that all constraint indicators converge at one
percent. Temporal quadrature, a continuum norm enclosure and the map from
constraint control to a physical observable error remain distinct needs.
The source records explicitly set both continuum and temporal certification
flags to false.

All runs keep $r,Q>0$. Fine endpoint proper radial motion ranges from about
$-0.0370994$ to $+0.00663986$. The fine Gram gap stays below
$2.8\times10^{-14}$. Neither the historical $10^{-8}$ initial line nor the
prescribed-gauge coordinate endpoint was used to veto this formulation.

## Producers and budget

The sealed producers are
[`derive_nsc_spherical_conformal_episode.py`](../scripts/derive_nsc_spherical_conformal_episode.py)
and
[`derive_nsc_spherical_conformal_frames.py`](../scripts/derive_nsc_spherical_conformal_frames.py).
The episode uses $120.66$ CPU seconds. Episode plus the finalized physical
postprocessing uses about $126.62$ CPU seconds; an earlier postprocessing
pass adds about $5.90$ CPU seconds of preparation. The combined saved payload
is about $8.44$ MB, within $64$ MiB. The original v5, prescribed episode,
source and producer hashes are retained; neither trajectory nor initial
radius is regenerated for physical frame postprocessing.

The independent owners
[`test_nsc_spherical_conformal_gauge.py`](../tests/test_nsc_spherical_conformal_gauge.py),
[`test_nsc_spherical_conformal_episode.py`](../tests/test_nsc_spherical_conformal_episode.py)
and
[`test_nsc_spherical_conformal_frames.py`](../tests/test_nsc_spherical_conformal_frames.py)
check the action, lapse chain rule, dynamic work, actual constraint derivative,
source and initial-data equality, exact frozen group, interval quadrature and
saved forcing/physical ledgers. They do not produce a new trajectory or
certify the remaining continuum and observable-error gaps.
