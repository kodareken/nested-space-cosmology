# Same-action source controls on the conformal realization

The four coarse controls in
[the source-control record](../results/development/nsc-spherical-conformal-source-controls-v1.json)
reach $T=0.3$ at $n_f=256$, $n_q=1024$, with step cap $0.0005$. They use the
same conformal action and the unchanged six prepared $T=0$ columns. The
source parameter is

$$
c(\alpha)=(0.5+0.25\alpha,0.5+0.25\alpha,0.5,0.5,
0.5-0.25\alpha,0.5-0.25\alpha).
$$

The tested values are $0,-1,0.95,1.05$. The total occupation is analytically
$3$, with only floating representation roundoff in the stored sum. The
original $\alpha=1$ trajectory, preparation and producers remain sealed.
These are changes of initial source state under the same law. No source is
injected later, no radius is imposed, and no evolving state is reset.

## Constraint preparation and domain

The existing `nsc_regeneration_controls.prepare_population` recomputes the
source from the same columns with each new occupation vector. Its positive
radius Newton and projected periodic shift antiderivative prepare each
initial state. $Q$ remains the same constant initial ratio; $\chi,p_r,p_\chi$
start at zero. The current array is unchanged by the shift solve. A periodic
momentum cannot remove its mean, so that mean remains in the residual.
The four means are approximately $-1.41,-0.58,-2.20,-2.28$ times $10^{-14}$.
No current mean was deleted.

All four preparations are source-derived and positive. They retain a
Newton solver floor, with full Hamilton maxima about $0.0013742,0.0016135,
0.0017485,0.0017878$, respectively. The projected, full and held-out domains
are reported separately. The constructor's older tolerance flags remain
reported; those flags were not turned into a new scientific veto.
The controls follow the actual Galerkin rates of the
[conformal formulation](nsc-spherical-conformal-gauge.md). Coordinate work
includes $F_Q\dot Q+F_L\dot L$, with $L=Q$. Normal-observer regional work
keeps pressure, momentum-lapse-gradient work and the actual balance defect.

## Observed source dependence

The following are signed time integrals of the outward normal energy flux
at the physical surfaces. For $x=0$ the sign is outward from $[0,2]$ on its
left; for $x=2,4$ it is outward to the right from the preceding region.
The baseline is the joined sealed coarse series from $T=0$ to $0.3$.

| Source | $x=0$ | $x=2$ | $x=4$ | Initial/final normal-shell leader |
|---|---:|---:|---:|---|
| Sealed $\alpha=1$ | $+0.0393771$ | $+0.0702504$ | $+0.0635174$ | $0\to0$ |
| Uniform $\alpha=0$ | $+0.0266388$ | $+0.00000701$ | $+0.1241505$ | $1\to0$ |
| Reversed $\alpha=-1$ | $+0.0136283$ | $-0.0702352$ | $+0.1835064$ | $1\to1$ |
| $\alpha=0.95$ | $+0.0387420$ | $+0.0667135$ | $+0.0666025$ | $0\to0$ |
| $\alpha=1.05$ | $+0.0400123$ | $+0.0737913$ | $+0.0604253$ | $0\to0$ |

Reversal changes the sign at $x=2$ and moves the maintained normal-shell
leader to region 1. Both five-percent imbalance perturbations retain
region 0 as the leader throughout. Their $x=2$ signed-flow changes relative
to the baseline are approximately $-5.035\%$ and $+5.040\%$; their $x=4$
changes are about $+4.857\%$ and $-4.868\%$. At $x=0$ the changes are
approximately $\mp1.613\%$. These are measured source-dependent responses,
rather than an assumed uniform-source disappearance rule.

Uniform occupations on the six columns are half of their rank-six
projector. They are not $\tfrac12 I$ on the full retained fermion band.
The empty complement still permits dynamics and outer regional flow.
The uniform run has two comparably occupied width-2 regions and changes
its normal-shell leader, while its probability leader remains the same.
Its minimum maximum-window normal-shell and probability shares are about
$0.49732$ and $0.49750$. This outcome is recorded directly; no mandatory
minimum share or reversal was imposed on the controls.

The reversed source keeps its leader with minimum normal-shell and
probability shares about $0.63676$ and $0.64212$. The two perturbations
have minimum normal-shell shares about $0.64323$ and $0.65874$.
All controls maintain a positive evolving chart, and their maximum column
Gram gaps stay below $1.7\times10^{-13}$. Their finite constraint-forcing
budgets have positive endpoint margins about $8.81$ to $9.83$ times
$10^{-5}$.

The field/gravity coordinate exchange is approximately $2.41,2.44,2.38,2.37$
times $10^{-5}$ for the four sources, with total coordinate-energy changes
about $-1.3\times10^{-11}$. The regional shell ledger reports pressure,
lapse-gradient and boundary-flow integrals separately. The original fixed
$T=0$ observer is used for the occupation measurements; its initial
occupation changes with the declared source vector.

## Coarse uncertainty and what is established

The controls are exploratory coarse trajectories. They do not have a new
independent fine trajectory or a control-specific continuum error enclosure.
The sealed baseline's spatial and time indicators provide reference scales,
not certificates for every changed source.

For the measured five-percent differences at $x=2$, the paired
Simpson-versus-trapezoid change indicator is about $1.51\times10^{-6}$,
against a flow difference about $0.00354$: approximately $0.043\%$ of that
observed difference. At $x=4$ the corresponding indicator is about
$9.71\times10^{-7}$ against $0.00309$, or $0.031\%$. At $x=0$ it is about
$6.04\times10^{-7}$ against $0.000635$, or $0.095\%$.
These are numerical comparisons of the same frame protocol and signed
quantities, not rigorous temporal integration bounds. The baseline spatial
indicators for the underlying flows are below $0.007\%$, while its time
indicators are substantially smaller. The observed perturbation effects
are therefore separated from those coarse reference scales. This does not
close a total physical-observable error bound.

The controls remove an assumption that maintained surface flow would be
source-independent or must vanish at uniform occupation. They do not prove
renewal, a continuum limit, an infinite continuation or an absolute vacuum
stress. Those claims are not used to describe these finite observations.

## Reproduction and preserved evidence

The producer is
[`derive_nsc_spherical_conformal_source_controls.py`](../scripts/derive_nsc_spherical_conformal_source_controls.py).
Preparation took about $14.61$ CPU seconds; the measured forecast admitted
all four $T=0.3$ controls, and the complete run used $387.78$ CPU seconds of
the $600$-second exploratory budget. The compact payload is about $925$ kB.
It stores full initial and final Cauchy checkpoints, the original columns
and observer, every global numerical sample and every physical regional
frame. It does not duplicate the earlier trajectories.

[`test_nsc_spherical_conformal_source_controls.py`](../tests/test_nsc_spherical_conformal_source_controls.py)
independently checks the source family, unchanged inventory/columns,
constraint-preparation records, current-mean retention, actual-time forcing
integrals and signed surface accounting. The original baseline file hashes
are bound before and after the run.
