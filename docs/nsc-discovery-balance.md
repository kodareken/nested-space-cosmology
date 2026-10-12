# Disjoint RK4 normal-energy accounting

The first-batch observation ledger summed pressure and lapse over child
$[1,3]$ and parent $[0,4]$, double-counting their overlap, and omitted an
integrated finite-projection defect. This successor accounts for disjoint
windows on the same full-state ODE:

| Window | Coordinate interval |
|---|---|
| left | $[0,1]$ |
| child | $[1,3]$ |
| right | $[3,4]$ |
| ambient | $[4,8]$ |
| parent | sum of left, child, right |
| global | sum of all four |

The owner is [nsc_discovery_balance.py](../src/recursive_horizons/nsc_discovery_balance.py),
with [CLI](../scripts/derive_nsc_discovery_balance.py) and
[tests](../tests/test_nsc_discovery_balance.py). Existing trajectories and
source files are preserved.

At each of the four owned RK4 stages, one `model.rates(...,
return_bundle=True)` call supplies both the exact state rate and its full
source bundle. Stresses reuse that source. The state arithmetic is identical
to `model.rk4_step`; coupled and frozen controls are checked bitwise.
Production applies the existing episode step restriction at each accepted
state. A separate endpoint source evaluation measures the final stock.

The regions owner's analytic polarization of $F_L/r$ supplies the shell
slope. For a window, the instantaneous identity is

$$
\dot E=F(a)/\Delta x-F(b)/\Delta x+P+L+D_{\rm finite}+R.
$$

$D_{\rm finite}$ retains both the pointwise projection residual and the
gap between interpolated endpoint flux and the zero-Nyquist derivative
window. $R=\eta(b)\dot b-\eta(a)\dot a$ is zero for these fixed cuts.
There is no extra piston work on a measurement cut. Multiplicity $M=4\kappa$
is included once in the source, and nodal shell summands are divided by
$\Delta x$ once.

All channels integrate with $h(f_1+2f_2+2f_3+f_4)/6$, at the actual
RK4 stages. The endpoint stock difference minus their sum is recorded as
`quadrature_closure_error`; it is not filled with zero. Parent/global
channels derive from the disjoint pieces, and internal fluxes cancel.
Coordinate metric work $F_Q\dot Q+F_L\dot L$ remains a separate ledger.
Frozen metric work evaluates to zero while spatial lapse work remains.

Probability current uses the weighted $\psi^\dagger\sigma_2\psi$.
Each endpoint's positive inward/outward current is integrated at the
stages and reported as a gross **lower bound**: unresolved opposite-moving
traffic can exceed that net-current split. Their difference equals the
signed boundary integral. The projected probability slope and its defect,
endpoint probability change and RK quadrature remainder are also retained.

Normal clocks at $x=1,2,3$ keep the episode's trapezoid protocol per accepted
step. RK-stage clock increments are an additional comparison, not a change
of observer protocol. Common coordinate stations are not common proper clocks.

## Production and replay

No production trajectory was run during implementation. After freezing
producer code, the following commands create exclusive scalar records under
the new `lab/results/development/nsc-discovery-balance-v1` directory:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_balance.py --check --target-time 3
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_balance.py --run --target-time 1 --spatial \
  --write results/development/nsc-discovery-balance-v1/coupled-T1.json
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_balance.py --run --target-time 3 --spatial \
  --write results/development/nsc-discovery-balance-v1/coupled-T3.json
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_balance.py --check \
  --record results/development/nsc-discovery-balance-v1/coupled-T3.json
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_balance.py -q
```

`--run` advances nf256 at caps $0.001$ and $0.0005$; `--spatial` adds the
saved nf512 handoff at $0.0005$. `--control frozen_geometry` selects that
control. Default endpoint is $T=1$; another endpoint is explicit.
`--write PATH` alone implies a run; `--check --write PATH` saves a forecast.
Existing destinations and paths outside the new balance directory are
refused before numerical work. Individual records are capped at 64 MiB.

Input and producer hashes are checked before and after production.
`--check --record` authenticates those hashes and replays the recorded
channel sums, endpoint identities, defect splits and partition sums without
re-evolving a trajectory. It does not independently regenerate numerical
endpoint states. Source changes cause strict replay authentication to fail.

The forecast uses historical accepted-step CPU and nf512 pilot CPU;
cap-based counts are lower bounds when the owned restriction requires
smaller steps. Added ledger overhead is unmeasured. An aggregate six
CPU-hour limit covers all cases within one run and is checked before each
step. A production pilot remains necessary to refine the historical estimate.
Tests use manufactured states only, independently compare the analytic
slope with endpoint differences, verify bit-identical state updates and
source reuse, and check that halving the step reduces energy quadrature
error and global probability drift.

The resulting ledger measures finite projected accounting, with a stated
quadrature remainder. Continuum certification, renewal and a stability
theorem remain unclaimed; tidal growth is assessed by its separate owner.
