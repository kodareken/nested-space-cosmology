# Assessment of discovery episode v1

The [independent consumer](../scripts/assess_nsc_discovery_episode.py) reads
finished first-batch inputs in `lab/results/development/nsc-discovery-episode-v1`:
six cases, 18 immutable chunks, stations $T=0.3,1,3$, and summed
child-process CPU $277.185634$ seconds. The numerical producer is
`1c9e7705c53e123497ffd40f94fb365e250e8e79`. Its binding envelope and
`authenticated-closure.json` authenticate the consumed episode artifacts.
The assessment binds every consumed input and its own script by SHA-256.
It takes no evolution steps and creates no production assessment until
`--write` is explicitly invoked after the consumer is frozen.

Spatial comparisons use `nf128_coupled_dt0.0005` and
`nf256_coupled_dt0.0005`; time comparisons use the latter and
`nf256_coupled_dt0.001`. Relative gaps are $|a-b|/\max(|a|,|b|)$.
There is no universal one-percent gate or continuous evolution error bound.

## Physical movements and control subtraction

At fine coupled stations $T=0.3,1,3$, child proper length is
$2.375218143$, $1.488905010$, $0.001020437674$. Parent proper length
moves from $4.709426854$ to $0.002156434760$. The fine coupled minimum
areal radius at $T=3$ is $20.17060796$, versus the frozen geometry's
$4.305851583$. Child length and the areal-radius floor have spatial gaps
$0.046149\%$ and $0.0097647\%$ at $T=3$.

The independent consumer reconstructs the saved $W$ geometry map and
periodic geometry interpolant on the declared $4n_f$ quadrature carrier.
It integrates $rQ$ for proper lengths and evaluates the three fixed
coordinate clock rates. It separately prolongs the two antiperiodic
matter blocks with their canonical normalization, integrates their weighted
probability, and checks these recomputed values against the saved rows.
This preserves the observer's interpolant quadrature; the additional nodal
box sum has its own cut error and is reported separately.

The weighted source occupation trace remains $3$. The fine coupled child
fractions of that trace are $0.4255971504$, $0.4999799261$,
$0.00100306588$. The frozen child fraction at $T=3$ is $0.00414463483$.
The fine coupled localization peak moves from $x=0.53125$ through
$x=1.671875$ to $x=3.25$. Child probability leaves both controls; the
shrinking proper length is a coupled geometry movement.

Each control effect is explicitly
$[O_c(T)-O_c(0.3)]-[O_f(T)-O_f(0.3)]$. The initial values remain distinct:
for example, the coupled and frozen initial $R_h$ maxima are
$15.52337813$ and $0.06877947517$, despite their common spatial state.

| Initial-subtracted fine control effect | $T=1$ | $T=3$ | Spatial gap of the $T=3$ effect |
|---|---:|---:|---:|
| Child proper length | $-0.8863131330$ | $-2.374197706$ | $0.0043464\%$ |
| Child probability / trace | $3.344150958\times10^{-5}$ | $-0.003141568955$ | $0.0537085\%$ |
| Child probability / proper length | $0.3759579809$ | $2.943693543$ | $0.362889\%$ |
| Child normal energy | $3.040300476$ | $-1.128567471$ | $9.29279\%$ |

The last two quantities have distinct meanings: probability per proper
length rises because the window contracts, and normal energy is a signed
matter measurement with an incomplete balance ledger.

## Curvature and clocks

At $T=1$, the $R_h$ maximum spatial gap is $0.0232510\%$. At $T=3$,
$R_h$ maxima are $34240.39162$ and $29030.98879$, a $15.2142\%$ gap;
Weyl $C^2$ maxima are $1780.692331$ and $1279.442408$, a $28.1492\%$
gap. The fine time-step comparison for $R_h$ max at $T=3$ is
$9.08345\times10^{-7}$. The initial-subtracted control effects have
spatial gaps $15.2211\%$ and $28.1509\%$ for these two maxima.
These scalar comparisons do not classify the radial and angular tidal
response; a separate native-Jv assessment owns those computations.

The first sampled increase of the relative $R_h$ spatial gap is $T=0.35$
(from $1.16449\times10^{-4}$ to $1.23851\times10^{-4}$). The Weyl gap
first increases at $T=0.4$. Neither gap is monotone across the entire
cadence. A separate scale comparison finds the spatial gap first exceeding
the fine case's previous cadence movement at $T=2.90$ for $R_h$ and
$T=2.65$ for Weyl. These are cadence diagnostics, not acceptance thresholds.
The sealed coupled chunks contain no $\ddot Q$ summary; the consumer does
not substitute the atlas acceleration or the frozen zero into late curvature.

All station/control comparisons use coordinate time. From $T=1$ to $3$,
the fine coupled proper-time increments at $x=1,2,3$ are
$0.2429581484$, $0.2494421047$, $0.2644711378$; frozen increments are
$2.417449034$, $2.379804788$, $2.316074043$. Thus these control effects
are not measurements at matched proper clocks.

## Work and missing energy closure

The independent coordinate-time trapezoids reproduce every saved work
ledger with absolute discrepancy zero. Fine coupled integrals through
$T=3$ are coordinate fieldwork $-0.01024838243$, pressure
$9342.171467$, lapse $-360.8553489$, net child boundary
$-1613.835707$, and net parent boundary $-546.1244030$.
Pressure and lapse channels sum the child $[1,3]$ and parent $[0,4]$
windows, which overlap. They are neither per-child summands nor a global
disjoint energy ledger. Agreement with the stored trapezoid checks sampled
bookkeeping, not an integrated normal-energy identity.

The finite-projection normal-energy defect is unrecorded and remains null.
Net boundary samples do not separate incoming and outgoing contributions.
At fine coupled $T=3$, field and gravity energies are $11.9897660831$
and $-11.9897660881$; their small total does not establish regional closure.
Separate $Q,r,\chi$ and source-column energy accounts are absent.
Constraint projection forcing is a different ledger. The atlas's final
child balance interval has residual $0.0115388171$ without its unrecorded
defect. Normal-energy closure, continuum curvature certification,
regeneration and a stability theorem remain unclaimed.

The next assessment should combine these coordinate comparisons with the
separate native-Jv radial/angular response, choose matched proper-clock
comparisons, and provide a finite-projection-defect normal-energy ledger
before extending trajectories.

## Reproduction

```sh
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_episode.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_episode_assessment.py -q
```

`--write PATH` exclusively creates one new JSON file after authentication.
An existing destination, any destination inside the sealed episode, and
atlas JSON/payload destinations are refused. `--check` prints deterministic
JSON and exits unsuccessfully when binding or consistency checks fail.
