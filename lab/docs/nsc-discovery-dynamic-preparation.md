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
  scripts/derive_nsc_discovery_dynamic_preparation.py --prepare \
  --write results/development/nsc-discovery-dynamic-preparation-v1
.venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_dynamic_preparation.py \
  --check results/development/nsc-discovery-dynamic-preparation-v1
```

The focused [tests](../tests/test_nsc_discovery_dynamic_preparation.py) use
explicit smaller-resolution/harmonic test preparations and temporary
checkpoints only. Root owns frozen NF128 execution and any later episode.
