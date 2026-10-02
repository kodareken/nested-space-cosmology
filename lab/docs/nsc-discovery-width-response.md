# Physical source-width response

The owner is [nsc_discovery_width_response.py](../src/recursive_horizons/nsc_discovery_width_response.py),
with [driver](../scripts/derive_nsc_discovery_width_response.py) and
[focused tests](../tests/test_nsc_discovery_width_response.py). This is a finite
source-scale response of the existing one-action conformal system. It is not
an affine coordinate pullback or a universal scale exponent. The existing
family and prediction producers and their evidence are preserved.

## Declared experiment

The baseline source half-width is `1`. The held-out width is **1.05**, declared
before its source is constructed or measured. The two source-frame derivative
increments are `0.001` and `0.0005`; neither is the held-out displacement.
Occupations stay `(0.75,0.75,0.5,0.5,0.25,0.25)`, with trace 3. The angular
channel, action coefficients, carrier, Fourier operators, fixed observer and
frozen complete geometry basis W stay fixed. NF128 is available as an optional
spatial control; NF256 is the production resolution. The default invocation
prepares a short NF128 preview at coordinate duration 0.01.

The continuous source adapter uses the owning family bump/lobe primitives,
normalizes sampled half-density lobes, preserves source ordering and phase,
then projects the outer columns off the declared child and applies the same
outer Löwdin operation. The nominal child is exactly the owned reference pair.
At other widths any needed child orthogonalization is explicitly reported.
There is no alignment rotation between unequal occupation blocks. Source
supports are `[0,2-w]`, `[2-w,2+w]`, `[2+w,4]`; Fourier interpolation tails are
measured rather than called exact support. Full-precision float hex values,
frame hashes and a canonical numerical-input digest avoid catalogue rounding
or case-ID collisions.

## Initial direction and propagation

The initial state is the baseline source's own dense radius construction.
For each h the *complete* frame direction, including the outer orthogonalization,
is

\[
V_h=\frac{S(1+h)-S(1-h)}{2h}.
\]

The record measures source Gram defects, outer Gram positivity, tangency
`S^H V + V^H S`, and the Frobenius norm of

\[
\partial_w C=V\operatorname{diag}(c)S^\dagger+
 S\operatorname{diag}(c)V^\dagger.
\]

A nonzero covariance derivative is evidence of a change of the finite Gaussian
state on the fixed chart, rather than an invisible source-frame rotation. The
frame difference is a numerical width derivative; it is not mislabelled as an
analytic preparation derivative. Its later full-state Jv transport is analytic.

`prepared_nested_radius_tangent` already accepts nonzero column directions.
It computes the source derivative and solves the *same* projected initial
constraint Jacobian,

\[
J_r\partial_w r_0=-P_g\partial_w\rho,\qquad
\rho=\frac{M}{\Delta x_q}(K/Q_0+\kappa S_{\rm spin}).
\]

Both the initial column direction and implicit radius direction are retained.
The preparation family keeps `Q0=b0/a0` fixed and starts with zero chi and
canonical momenta. The occupation direction is zero. This restriction is an
initial-data choice, not suppression of subsequent geometric response. The
full coupled state and direction use the existing `coupled_rk4_step` and
`nested_rate_jacobian_vector`. There is one shared baseline and two width
directions, with no new dynamics or force.

The primary observable is **total six-column probability in the fixed child
window `(1,3)`**, not a tagged-source fraction, a moving preparation window, or
a child-observer kernel trace. At baseline coordinate time 0.3, stage-integrated
proper time supplies the matched-clock coefficient

\[
D_w O|_\tau=D_w O|_t-\dot O\,D_w\tau/\dot\tau,
\qquad \dot\tau=(rQ)(2).
\]

Both h coefficients and their difference are retained. The frozen held-out
forecast is `O(1,tau*) + 0.05 D_w O|tau*`.

Actual parent, fixed-child and source-support proper lengths and their ratios
are recorded. The width derivative of the source-support length includes its
moving-endpoint contribution,

\[
D_w L_{\rm source}=\int_{2-w}^{2+w}D_w(rQ)\,dx
 +(rQ)(2+w)+(rQ)(2-w).
\]

Ratios are computed from the actual metric, not from coordinate width alone.
Reported length derivatives are labelled at fixed coordinate time.

## Immutable stages and independent measurement

`prepare.json/.npz` save the initial state, basis, observer, four derivative
frames, two initial directions and solver diagnostics. `prediction.json/.npz`
seal the baseline terminal state, propagated directions, clock data and held-out
forecast **before** `--measure` can construct the width-1.05 source.

The measurement uses a new own-source dense radius constructor with unchanged
observer and W, then the existing nonlinear RK4 and bracketed decreasing-step
proper-clock root. It does not copy the baseline radius or reset the field.
The target is the baseline's actual proper time at the declared coordinate
endpoint. A crossing can require coordinate time beyond 0.3; the existing
25% allowance is declared as an event ceiling of 0.375. Preview measurement
never exceeds 0.01. A target not bracketed by the ceiling remains unavailable.

Preparation failures, budget stops and chart/numerical stops retain their
status and finite last admissible arrays. An initial-constructor blocker does
not imply non-existence. All files use exclusive creation, become read-only,
and are authenticated by payload and array hashes. Each JSON or NPZ chunk is
at most 64 MiB. A stopped immutable stage requires a new successor, not repair
or replacement. Producer pins authenticate committed historical bytes
independently of the current verifier; uncommitted source pins are explicit.

The aggregate budget is at most **300 CPU seconds**, with a forecast multiplier
of **1.5**. One root executor owns the sequential preparation, two-direction
forecast and one held-out arm. Budget and remaining-cost forecasts are checked
between steps; an individual dense solve or bounded clock-root call cannot be
interrupted internally. The saved total includes CPU spent in earlier stages.
No campaign process pool is opened by this owner.

The measurement retains observed and predicted effects and the signed nonlinear
remainder. Available indicators include the h-forecast spread, proper-clock
root residual and bracket interpolation difference. Effect size is compared
with these named indicators, without a mandatory 1% gate. Those indicators are
not an evolution error bound. Spatial and timestep refinement remain separate
controls; until measured, the record explicitly leaves total evolution
uncertainty unresolved. A positive covariance derivative or held-out agreement
is not a continuum certificate, a stationary state, or a universal power law.

## Commands

After root freezes and commits the producer, run from the repository root.
The launcher changes the scientific working directory to `lab`, so evidence
paths begin with `results/`:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_width_response.py \
  --prepare --output results/development/nsc-discovery-width-response-v1 \
  --nf 256 --production --duration 0.3 --cpu-budget 300
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_width_response.py \
  --predict --output results/development/nsc-discovery-width-response-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_width_response.py \
  --measure --output results/development/nsc-discovery-width-response-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_width_response.py \
  --check --output results/development/nsc-discovery-width-response-v1
```

`--run` combines those stages sequentially and still seals the forecast before
measurement. A separate successor with `--nf 128` provides a spatial control;
its own target proper time is reported and must be accounted for when comparing
resolutions. Focused tests use short preview evolutions only:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_width_response.py -q
```
