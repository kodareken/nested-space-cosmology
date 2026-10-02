# General-Q, chi=0 dynamic preparation

The [owner](../src/recursive_horizons/nsc_discovery_dynamic_preparation.py)
prepares initial data for the existing common-action conformal evolution.
It imposes no static force condition. The fixed Gaussian has six positive
real standing eigenmodes of the actual retained Dirac operator, unchanged
occupations $(.75,.75,.5,.5,.25,.25)$, trace 3 and multiplicity $4\kappa$
once. Its covariance commutes with that initial operator. Its source is
recomputed when selecting a different initial Q, never reset during dynamics.
Mean current and the full shift residual are measured without deletion.

The balanced n=8 [saved geometry](../results/development/nsc-discovery-stationary-balance-seed-v2.json)
supplies only the declared mean and cosine shape of Q. The mean is lifted by
1.1. A finite seven-cosine polynomial avoids exponential interpolation aliases.
The same-mean uniform control removes that pattern. Both preparations set
chi and all geometric momenta to zero. Q is an initial input, not a claimed
dynamically selected or holding geometry. The stationary
[critical-v3 result](../results/development/nsc-discovery-stationary-critical-v3.json)
retains its noncompact tendency and unsatisfied local critical equations;
this constructor does not continue its optimizer.

## Positive continuum seed

At chi=0 and zero momenta, the smooth periodic lapse constraint with
$y=\sqrt r>0$ is

$$
-4\partial_x(y_x/Q)+Qy-
\frac{\rho+2\pi C_F\mathrm{flux}^2 Q}{8\pi A y^3}=0.
$$

Write $B=(\rho+2\pi C_F\mathrm{flux}^2 Q)/(8\pi A)$.
For positive Q, A and B, its energy

$$
\mathcal E[y]=\int\left[2y_x^2/Q+Qy^2/2+B/(2y^2)\right]dx
$$

is strictly convex on the positive domain. This supplies a positive seed,
not a claim that a spectral product rule holds. The implementation solves
the pulled-back finite variation of this energy, maintaining fine positivity,
then maps its y to an initial radius.

## Exact finite correction and scope

The next solve calls the existing sampled `hamilton_constraint` on its
actual general Q and chi=0 fields. Its Jacobian is the adjoint pullback of

$$
\delta C=2Z(Dr)D\delta r/Q-QV_r\delta r
-2D\bigl[D(F_r\delta r)/Q\bigr].
$$

The source rho is independent of r in canonical conformal variables.
Neither the old constant-Q radius residual nor its Jacobian is reused outside
that domain. Positivity-preserving Newton has finite iteration, backtracking
and CPU limits. A stalled correction is a named preparation blocker, not
nonexistence. The record preserves actual represented and unresolved lapse
residuals, source/eigenmode tails, current, CAR, geometry positivity and
nonzero dynamic forces. No continuum Cauchy certificate or holding claim is
made. No trajectory is executed.

## Checkpoint and replay

The builder is pure and preflights by default. `--prepare` explicitly constructs
both cases under a shared 30 CPU-second cap. `--write` is separate and creates
new files only. The checkpoint exposes each full nodal state, its nested
canonical encoding, geometry frame, eigenvalues and fine source/constraints.
`load_case` returns the existing pair and encoded state for the existing
episode interface; it does not call the old initializer. `--check` authenticates
hashes and replays the owned source and constraints without a solve or write.

```sh
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py --prepare --producer-commit HEAD \
  --write results/development/nsc-discovery-dynamic-preparation-v1
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --check results/development/nsc-discovery-dynamic-preparation-v1
```

The focused [tests](../tests/test_nsc_discovery_dynamic_preparation.py) use
explicit smaller-resolution/harmonic test preparations and temporary
checkpoints only. Root owns frozen NF128 execution and any later episode.

## Frozen controls and one coherent successor

The immutable v1 checkpoint was produced at `2dab719`. Its original JSON has
eight source hashes but no `producing_commit`. Extending this owner does not
change that checkpoint or repair its header. A separately supplied observed
binding names the exact JSON and payload SHA256 plus `producing_commit`.
Authentication checks the old hashes against immutable Git blobs and labels
this as a post-run external observation. Numerical replay additionally
requires the live physical owners to match their historical bytes.

The successor copies both uniform and patterned controls by value. Its one
additional source rotates the **stored** real eigenframe columns 0 and 4 by
$\theta=\pi/6$ before solving its own chi=0 constraint. The old frame's order,
signs and hash are retained. The rotation preserves Gram, CAR eigenvalues,
occupations and trace 3, while changing the physical covariance because the
two occupations differ. The covariance Frobenius difference is
$\sqrt2/4$. Its nonzero Hamiltonian commutator is expected physical coherence,
not a rejected spectral preparation. It is not charge conjugation.

Source density differences are measured at the same Q. The new radius solves
the source of the rotated frame; no old radius is rebound as a solution.
Both timestep caps fork this same stored state. No independent degenerate
eigenframe is chosen for a timestep comparator. Spatial convergence is not
claimed; any later spatial confirmation must preserve a compatible physical
source frame and phase, then solve that source's radius.

## Thin diagnostic episode

`prepare_episode` freezes six cases: the three preparations at NF128, each at
caps 0.001 and 0.0005. Stations are coordinate times 0.3, 1 and 3. The saved
native geometry frame, canonical momenta, full field columns, actual clocks,
source columns and fixed reference are retained in the existing episode
checkpoint format. Initial owned lapse and shift constraints must be solved
within the explicitly declared resolved error, default $10^{-6}$ in their
owned density units. Static critical forces are not an admission condition.

`run_episode` delegates the existing episode pool worker, coupled RK4,
principal step restriction, checkpoint publication and locked aggregate CPU
ledger. At most four workers share one pool. Default budget is 300 CPU
seconds, including preparation and observations; an explicitly accepted
six-hour budget is the upper limit. Last admissible states, finite budget
stops and immutable chunks of at most 64 MiB are preserved. Resume uses the
stored state and cumulative ledger, without source reselection, constraint
reset, or another initial solve.

The scoped adapter adds station 0.3 and redirects diagnostics to the frozen
period-aware extent consumer (`real-periodic-rfft-v2`, source SHA256
`0b9adea69e33f4e1683d50f44d6961f76ab71cea4533d64c025d5c8aa7080170`,
committed by root at `6562313`). Real clock samples still implement the same
$d\tau=rQ\,dT$ and owned endpoint trapezoid per RK4 step. Rows use the actual
returned runner clocks. No old atlas implementation is selected dynamically.

Primary localization uses **all columns**: positive probability per proper
length, its contrast and peak (flat profiles have no assigned peak), the
effective proper participation width and fixed spatial contents. Columns 2
and 3 are only selected spectral labels, with no prepared-child ancestry.
The observer retains actual metric tides, the constraint pair $h_c=QC$ and
$D$ separately from raw C, actual rates and source forces, fixed-window
normal energies, disjoint pressure/lapse/flux accounting and trapezoid work.
No moving boundary or physical wall is introduced.

A first sampled sign change of an initially positive local lapse-source
density is saved as an event and evolution continues. CAR positivity does
not imply local Dirac energy-density positivity later. That event is neither
an instability nor a sector failure. Stops are chart failure, nonfinite or
diagnostic admissibility failure, numerical restriction and budget stops.
Startup force alone is not evidence of holding. Each returned state retains
the finite observed effect or its named blocker; no continuum or propagated
observable error certificate is claimed.

The returned assessment reports changes in child proper length, areal radius,
all-column probability, normal energy, source contrast and proper participation
width, together with actual tides and the last constraint pair. Timestep
comparisons join equal timestamps and verify their identical initial-state
pins. They remain finite timestep indicators, with no spatial convergence or
holding inference. A case with only a startup row returns its recorded stop or
the explicit blocker `no evolved physical observation`.

Every new production checkpoint requires a frozen producing commit and the
complete observed local dependency closure. The parent freezes source bytes
before these explicit commands; the preserved v1 binding remains external:

```sh
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --check results/development/nsc-discovery-dynamic-preparation-v1 \
  --observed-binding /path/to/v1-observed-binding.json
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --coherent results/development/nsc-discovery-dynamic-preparation-v1 \
  --observed-binding /path/to/v1-observed-binding.json \
  --execute --producer-commit HEAD \
  --write results/development/nsc-discovery-dynamic-preparation-v2
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --episode results/development/nsc-discovery-dynamic-preparation-v2 \
  --execute --producer-commit HEAD \
  --write results/development/nsc-discovery-dynamic-episode-v1
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --run-episode results/development/nsc-discovery-dynamic-episode-v1 \
  --workers 4 --episode-cpu 300
```

The observed-binding JSON has `producing_commit`, `record_sha256` and
`payload_sha256`. Its assertion is checked, not copied into the old record.
Tests use temporary checkpoints and one owned RK4 step; scientific production
remains root-owned.
