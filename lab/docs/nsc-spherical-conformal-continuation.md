# Finite conformal continuation and local response

The same-action conformal realization connects regional source differences,
energy transfer, geometric response, maintained localization, and a quantitative
local response on one saved trajectory. The primary runs reach $T=0.3$ from
the original v5 preparation. The accepted finite precision criterion is movement
below one percent of each reported effect; refinement and quadrature comparisons
are numerical indicators unless separately enclosed. This result demonstrates
the maintained branch. Renewal is false.

## Preparation, action, and trajectory

The [initial episode](nsc-spherical-conformal-episode.md) runs from $T=0$ to
$0.05$. Each of its four final Cauchy arrays hands off bitwise to the matching
[continuation record](../results/development/nsc-spherical-conformal-episode-v2.json)
and [payload](../results/development/nsc-spherical-conformal-episode-v2.npz),
which reach $T=0.3$ at $n_f=256,512$ and step caps $0.001,0.0005$.
The source columns and weights remain unchanged. There is no midtrajectory
reset or imposed radius law.

The [conformal formulation](nsc-spherical-conformal-gauge.md) varies the same
action before setting $L=Q$, $\beta=0$. The actual Galerkin rate retains the
gauge chain rule, source force, and coordinate work $F_Q\dot Q+F_L\dot L$.
Positive $r,Q$ and the column Gram remain admissible. For the fine confirmation
continuation, field/gravity exchange is about $2.36029\times10^{-5}$ and the
total-energy change is about $-1.09\times10^{-11}$.

## Regional transfer and maintained structure

Four fixed width-2 windows use exact integrals of the periodic nodal
trigonometric interpolants. Normal shell energy is $F_L/r$; surface flux uses
the owned flux_nodal / dx normalization. The
[surface record](../results/development/nsc-spherical-conformal-transport-v2.json)
keeps entry and exit separately from their net contribution.

| Surface | Absolute normal-energy-flow integral, $[0.05,0.3]$ | Largest frame indicator |
|---|---:|---:|
| $x=0$ | $0.03937795$ | $0.0950\%$ of that flow |
| $x=2$ | $0.07025487$ | $0.0427\%$ |
| $x=4$ | $0.06351828$ | $0.0314\%$ |

Their spatial movements are below $0.0064\%$, and timestep movements are
smaller. Flow at $x=6$ remains unresolved in its own domain. A zero net flux
does not veto nonzero entry and exit. Internal current maxima are separate
quantities and are not substituted for surface exchange.

Region 0 remains both the shell and probability leader, with shares at least
$0.65098$ and $0.65585$. Its shell content changes by $-0.0413133882$.
The Simpson ledger combines net boundary contribution $-0.10956538$,
pressure work $+0.04884147$, and lapse-gradient work $+0.01941054$.
Its gap is about $1.99\times10^{-8}$; the trapezoid indicator is $0.1641\%$
of the shell change. Conservation, work, probability, and occupation retain
their different measures.

## Source dependence, metric, and clocks

The [source controls](nsc-spherical-conformal-source-controls.md) change only
the regional pair occupations, preserving the same six columns and total
occupation 3. Each initial radius and periodic shift momentum is prepared
from its own source; current means are retained. Reversing the imbalance
reverses signed flow at $x=2$ and moves the maintained shell leader to region 1.
Perturbing imbalance by $\pm5\%$ changes that flow by about $\pm5.04\%$.
The paired frame indicator is about $0.043\%$ of the perturbation difference.
These are exploratory coarse controls with their recorded Newton floor and
frame indicators, not new fine trajectories for every source.

The [metric curvature](nsc-spherical-conformal-curvature.md) uses the
directional Jacobian of the actual projected $\dot Q$, with $\dot L=\dot Q$.
Independent differences and a Christoffel contraction check the time jet and
scalar. The auxiliary $\chi$ is compared afterward. At $T=0.3$ the fine
maximum Weyl invariant is about $0.11331608$, with spatial movement about
$0.00110\%$. Dense/FFT spatial-derivative conditioning is reported separately.

The continuous normal clock is at $x=1$, $d\tau=rQ\,dt$. The
[clock record](../results/development/nsc-spherical-conformal-clock-v2.json)
evaluates the analytic frozen group at the matched clock, retaining the
original observer. Its occupation separation is approximately
$(-1.19109,-1.06846)\times10^{-5}$, with time and space indicators below
$0.00079\%$ and $0.00014\%$. It supersedes the sealed producer's interpolated
frozen-control comparison without replacing that historical record.

## Same-realization local response and finite scope

The [local successor](nsc-conformal-local-response-successor.md) consumes the
same generated geometry and actual $\Phi(0.2)$ on $[0.2,0.3]$, measured by
the original $T=0$ observer. The actual initial cross norm is $0.33863658$.
Occupation changes by $0.07342952$; streamed/full error is $6.3508\times10^{-6}$
or $0.00865\%$ of the effect. Memory, exterior drive, and actual cross
omission move occupation by $0.01681854$, $0.04854120$, and $0.09206143$,
each resolved against its own error or refinement. The representative
generator has no isotropic source-copy factor, and no time-indexed exterior
propagator is stored.

The independent memory-control supplement uses the same state and schedule
with its triangular omission generator. Its occupation discrepancy from the
streamed memory-off series is about $2.2011\times10^{-6}$, or $0.01309\%$
of the memory effect. This is a negative control, not a replacement physical
Hamiltonian.

The fine initial-chart certificate and source-derived controls remain in
their declared domains. Actual quadrature and rate-projection forcing are
tracked through the finite norm relation $R'\le\|f\|_2$; its sampled time
integral is not a validated integral. A propagated bound for an unknown
geometric path, continuum certification, and stress are stronger claims not
established here. They are not additional prerequisites for the accepted
finite measured realization. The incoming-gate campaign remains paused.

## Read-only replay

The episode, frame, continuation, surface and clock producers refuse existing
output before compute. Their scientific source bytes are authenticated at explicit checkpoint
c2d5fb6f02722d9f69ec70a7145ea3f2287266df; numerical input bytes remain
resident and hash-bound. Current-request checks remain strict.

~~~sh
python scripts/lab.py scripts/derive_nsc_spherical_conformal_continuation.py --check
python scripts/lab.py scripts/derive_nsc_spherical_conformal_transport.py --check
python scripts/lab.py scripts/derive_nsc_spherical_conformal_clock.py --check
python scripts/lab.py scripts/derive_nsc_conformal_local_response_successor.py --check
~~~
