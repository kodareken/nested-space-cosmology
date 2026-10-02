# Finite inherited grandchild handoff

The [new consumer](../src/recursive_horizons/nsc_discovery_grandchild.py) loads
full saved analytic state/tangent and the independently rooted $\alpha=+0.01$
arm from `nsc-discovery-prediction-v1`. It extends that finite realization
without reinitializing the source, geometry or constraints. The
[CLI](../scripts/derive_nsc_discovery_grandchild.py) previews by default;
[tests](../tests/test_nsc_discovery_grandchild.py) use a short future interval.
The future increment is declared explicitly on the existing normal clock at
$x=2$. The original $\Delta\tau=0.05$ record remains sealed. Its primary
predicted change is $-1.59545\times10^{-11}$ and measured change
$-1.68151\times10^{-11}$, with error about $5.12\%$ of that tiny effect.
This is a readout/timing blocker: regional content is insensitive before
the exterior front arrives. It is not prediction success and does not call
for more precision on that early tail.

A new successor declares $\Delta\tau=0.75$, with absolute target
$\tau_*=1.1078554632$ near baseline coordinate time $T=1$, where exterior
transport can enter the grandchild region. The baseline forecast must again
be locked before measuring the new held-out future. The interval change
does not change the operator, full tangent, source or opening geometry.

## Observer and physical content

Declare $I_G=(1.5,2.5)\subset I_C=(1,3)\subset I_P=(0,4)$ on the original
period-8 carrier. Half-density dilation about $x=2$ of the original middle
profiles gives two grandchild seeds. These are evaluated on the fine AP
carrier and projected into the retained band. A single ordered QR of those
seeds, the original middle pair, and the four original outer modes supplies
fixed nested ranks $2,4,8$. Ambient fields retain the entire original band.
No observer is selected again during evolution.

The record reports the dilated seed's unresolved projection tail and each
observer's actual interpolant mass outside its assigned spatial interval.
Finite Fourier observers have tails; exact compact support is unclaimed.
The large nf128 dilated-seed projection tail leaves modal-basis comparisons
unresolved there. A coarse spatial $N_G$ readout can still be reported;
neither the basis nor the readout receives an acceptance label from QR.
The primary readout is

$$N_G=\int_{I_G}\sum_k c_k(|\phi_{0k}|^2+|\phi_{1k}|^2)\,dx.$$

Its stored half-density samples include the existing canonical normalization.
The rank-two modal covariance and its trace are secondary measurements.
For example, nf256 opening spatial content is approximately $0.5000000013$,
while its modal trace is approximately $0.1213977463$.

## Correct inherited clock

The saved NPZ has the analytic tangent but omits its accumulated $\delta\tau$.
The consumer recovers it from the saved matched child-content derivative:

$$\delta\tau=(\delta N_C|_t-\delta N_C|_\tau)\dot\tau/\dot N_C.$$

It independently cross-checks that value using the saved proper-mean-radius
derivative. An ill-conditioned recovery or inconsistent cross-check refuses
continuation and names unchanged-prefix replay as the required fallback;
it does not insert zero or borrow the nonlinear arm's Taylor clock.
On the opening nf256 data the two recoveries agree within $2.5\times10^{-16}$.
The inherited tangent is aligned once by

$$\delta y_\tau=\delta y_t-F\delta\tau/\dot\tau,$$

preserving the occupation direction $(1,1,0,0,-1,-1)$. The new segment's
clock tangent starts at zero. All future clock increments and their tangents
use the existing RK4 stage clock routines and analytic full-state Jv.

The future absolute target is `case.tau_target + proper_increment`. The held-out arm's
local increment subtracts its **actual saved rooted clock**, including its
small opening residual. A fresh bracketed decreasing-step event produces
the nonlinear endpoint; a Taylor clock correction is not the measurement.

## Driven causal composition

Let $V_G$ be the first two columns, $V_D$ the next two, and $V_F$ the next
four. Direct elimination retains $V_G$ against its full complement.
Sequential elimination retains $V_G$ through the child and parent layers,
with separate initially present exterior drive and generated memory.
Every route reconstructs full $\Phi$ before the common-action source forces.
All geometry uses the same generated full-state stage.

In particular, the detail-memory equation includes the exterior:

$$i\dot a_{D,m}=V_D^\dagger H(V_Da_{D,m}+V_Ga_G+
V_F(a_{F,d}+a_{F,m})+e_{P,d}+e_{P,m}).$$

The parent detail has the corresponding $V_F^\dagger H(e_{P,d}+e_{P,m})$
drive. Dropping these terms breaks the composition identity. Streaming
columns supply the causal drive/memory equations; no dense exterior
propagator or ambient covariance is stored.

The record compares reconstructed fields, actual full source forces at
stages, grandchild covariance, spatial content, clocks and source CAR.
It records force differences from omitting off-diagonal source cross terms.
Signed normal-energy stocks are read over six disjoint spatial pieces;
normal-energy closure is not asserted from those stocks.

## Lock, measure and replay

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py \
  --prepare --nf 256 --proper-increment .75 \
  --output results/development/nsc-discovery-grandchild-v1/nf256-dtau0p75
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py \
  --run --output results/development/nsc-discovery-grandchild-v1/nf256-dtau0p75
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py \
  --check --output results/development/nsc-discovery-grandchild-v1/nf256-dtau0p75
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_grandchild.py -q
```

`--prepare` probes baseline steps, admits a CPU forecast, evolves only the
baseline plus analytic tangent, and creates immutable `forecast.json/npz`.
It locks predicted $N_G$ before generating any new held-out future state.
`--proper-increment` accepts finite values in $(0,1]$; the default remains
$0.05$. Set the new increment during preview or preparation, never during
measurement or replay of an existing locked forecast. The primary cap is
$0.0005$. After that forecast has actually finished, its measured step count
and elapsed CPU determine the reserved held-out cost and whether a
$0.00025$ baseline/measurement pair still fits. Otherwise the record declares
that later step confirmation is required. `--no-matched-step` requests only
the primary cap. The sealed short-run probes give a historical planning
estimate of $184.23$ CPU seconds for the new primary interval; fresh measured
probes and the aggregate runtime guard own actual admission.

`--run` first replays the locked forecast, then advances the independently
saved nonlinear arm and both causal routes to the proper-clock event.
It creates immutable `measurement.json/npz`, preserving the forecast bytes.
A half-step measurement can be deferred if the remaining aggregate budget
is insufficient. Probes, preparation and measurement share a 300 CPU-second
budget; each JSON/NPZ pair is capped at 64 MiB. Existing prefixes are refused.

`--check` binds producer bytes and old input hashes, replays the
forecast from saved full state/tangent, recomputes measured $N_G$ from its
endpoint, and reconstructs both causal-route endpoints. It takes no steps
and never heals historical JSON. New measurement requires its current
premeasurement producer hashes. Historical replay after a source edit
requires `observed-run-binding.json`: both immutable artifact hashes and
every recorded producer are authenticated against that record's local Git
commit. It reports the historical producer commit and explicitly states
that readout replay uses current code without recomputing old trajectories.
Preparation records the producing commit only when all declared source
hashes match that commit's blobs; tests or uncommitted source changes retain
explicit working-tree hashes without assigning a false producer commit.

The exact composition identity and the new held-out prediction error are
separate evidence. A successful physical prediction is conditional on this
finite inherited state, observer and clock. It establishes neither universal
$\Omega$ scaling nor continuum certification, recurrence or renewal.
