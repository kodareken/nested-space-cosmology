# Finite inherited grandchild handoff

The [new consumer](../src/recursive_horizons/nsc_discovery_grandchild.py) loads
full saved analytic state/tangent and the independently rooted $\alpha=+0.01$
arm from `nsc-discovery-prediction-v1`. It extends that finite realization
without reinitializing the source, geometry or constraints. The
[CLI](../scripts/derive_nsc_discovery_grandchild.py) previews by default;
[tests](../tests/test_nsc_discovery_grandchild.py) use a short future interval.
The production experiment is a newly measured grandchild content forecast
at another $\Delta\tau=0.05$ of the existing normal clock at $x=2$.

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

The future absolute target is `case.tau_target + 0.05`. The held-out arm's
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
  --prepare --nf 256 --output results/development/nsc-discovery-grandchild-v1/nf256
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py \
  --run --output results/development/nsc-discovery-grandchild-v1/nf256
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_grandchild.py \
  --check --output results/development/nsc-discovery-grandchild-v1/nf256
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_grandchild.py -q
```

`--prepare` probes baseline steps, admits a CPU forecast, evolves only the
baseline plus analytic tangent, and creates immutable `forecast.json/npz`.
It locks predicted $N_G$ before generating any new held-out future state.
The primary cap is $0.0005$; cap $0.00025$ is included if the aggregate
forecast fits. `--no-matched-step` requests only the primary cap.

`--run` first replays the locked forecast, then advances the independently
saved nonlinear arm and both causal routes to the proper-clock event.
It creates immutable `measurement.json/npz`, preserving the forecast bytes.
A half-step measurement can be deferred if the remaining aggregate budget
is insufficient. Probes, preparation and measurement share a 300 CPU-second
budget; each JSON/NPZ pair is capped at 64 MiB. Existing prefixes are refused.

`--check` binds current producer bytes and old input hashes, replays the
forecast from saved full state/tangent, recomputes measured $N_G$ from its
endpoint, and reconstructs both causal-route endpoints. It takes no steps
and never heals historical JSON. Source changes fail strict replay.

The exact composition identity and the new held-out prediction error are
separate evidence. A successful physical prediction is conditional on this
finite inherited state, observer and clock. It establishes neither universal
$\Omega$ scaling nor continuum certification, recurrence or renewal.
