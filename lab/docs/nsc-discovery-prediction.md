# Stage-5 prediction: matched child content

This note specifies the short-time response of retained child content under
the existing full-state tangent. It does not publish a production history
and it does not form the future effective stress.

The executable owner is
[nsc_discovery_prediction.py](../src/recursive_horizons/nsc_discovery_prediction.py).
Checks live in
[test_nsc_discovery_prediction.py](../tests/test_nsc_discovery_prediction.py).
The directional primitive remains
[nsc_discovery_response.py](../src/recursive_horizons/nsc_discovery_response.py).
Carriers remain
[nsc_discovery_backend.py](../src/recursive_horizons/nsc_discovery_backend.py).
The nested state, its RK4 step and the dense source initializer remain
[nsc_nested_parent_child.py](../src/recursive_horizons/nsc_nested_parent_child.py).

## What is propagated

The baseline seed is the saved `nf128` or `nf256` Cauchy state at coordinate
time 0, together with the frozen replay basis `W`. The four saved cases are
`nf128_baseline_dt0.001`, `nf128_baseline_dt0.0005`,
`nf256_baseline_dt0.001` and `nf256_baseline_dt0.0005`, in
`results/development/nsc-nested-parent-child-v1` and
`nsc-nested-parent-child-replay-basis-v1`. Loading either endpoint checks
the recorded state hash. The `T=0.3` frame is a comparison target. Its
occupations are not replaced and its radius constraint is not solved again.

The occupation direction is the ordered contrast

\[
(+,+,0,0,-,-) = (1,1,0,0,-1,-1)
\]

on the separated packets (left parent, child, right parent). The sum is
zero, so the total occupation is unchanged. Child packet weights stay at
the saved values. The radius tangent is one linear solve,

\[
J\,\delta r = -\mathrm{pull}(\delta\rho),
\]

through `prepared_nested_radius_tangent`. Column, lapse, shift and momentum
tangents start at zero. No Newton loop is used on this tangent.

The coupled step evaluates the nested rate and
`nested_rate_jacobian_vector` at the four RK4 stages. The state combiner is
the nested model's combiner, so one unperturbed step matches `rk4_step`.
Occupation increments have no rate. Evolution uses the FFT carrier with one
FFT worker. The dense grid is reserved for the initializer and for the
radius Jacobian.

## Clock

The child clock is the existing sample \(\dot\tau=(rQ)(x=2)\). Over one
step the proper time and its tangent are the RK4 stage quadrature

\[
\int_{t}^{t+\Delta t} s = \frac{\Delta t}{6}(s_1+2s_2+2s_3+s_4)
\]

applied to \(\dot\tau\) and \(\delta\dot\tau\). Summing those increments is
the propagated clock. The frozen-slice product \(\delta\dot\tau(t_0)\,\Delta t\)
stays available as a comparison and is not the reported increment.

At the final slice the retained response is

\[
\delta O\big|_{\tau} = \delta O\big|_{t} - \dot O\,\frac{\delta\tau}{\dot\tau}.
\]

The primary observable is the child regional content on \((1,3)\). The
secondary observable is the child proper mean of \(r\). Both use this
formula. \(\dot O\) and \(\dot\tau\) are the final-slice rates. \(\delta\tau\)
is the stage sum from the seed.

## Independent control

Derivative steps are \(h=0.001\) and \(h=0.0005\). They are occupation
amplitudes, not the saved time-step caps, even though the numerals coincide
with those caps. For each sign the dense `initial_state` constructor builds
a new source-consistent radius at \(c_0+h\,v\). The saved baseline is not
passed through that constructor. Those states evolve on the FFT carrier.

Two comparisons are stored side by side.

The Taylor indicator evaluates both trajectories at the same coordinate
time and applies

\[
\delta O\big|_{\tau}\approx\delta O\big|_{t}-\dot O\,\frac{\delta\tau}{\dot\tau}.
\]

That is the same formula as the analytic prediction. It is labeled
`linear_clock_correction_not_the_equal_tau_measurement` and is not the
control that judges the prediction.

The equal-proper-time control stops each finite trajectory where its
integrated child clock equals the baseline proper time at the requested
coordinate duration. The crossing step is bracketed and the step is
decreased by bisection. A linear interpolation of the stored
\((t,\tau,O)\) samples is recorded beside the rooted event. The rooted
child content, minus the baseline content, is the measured effect. Its
difference from the predicted effect is the space error. The proper-time
residual and the coordinate-time offset are the time error. Neither number
is compared with an acceptance threshold.

The held-out amplitude is \(\alpha=0.01\). It is not a derivative step.
Its equal-time change and its Taylor indicator are both reported.

Source weights of every prepared arm stay in \((0,1]\). Covariance
eigenvalues stay in \([0,1]\).

## Run record

`--create` writes `prediction-manifest.json` and does not evolve.
`--run` creates `prediction-run.json` and `prediction-run.npz` by
exclusive write. The JSON holds the settings, source pins, producer file
hashes, CPU forecast, Taylor indicators, and equal-time errors. The NPZ
holds the final equal-time states, the analytic state and clock tangent,
and the selected scalar histories. `--check` on a run directory recomputes
the child content and the proper mean of \(r\) from those states. It does
not take a step and it does not call the initializer.

The run is one process. It does not start a process pool. The forecast
uses either one probed FFT step or the recorded `0.037 s` for three
`nf256` RK4 steps, scaled by resolution and by the rate count, including
the root iterations. Jacobian-vector products are counted as one rate, so
the forecast can sit low. Dense initial solves are named separately. The
forecast does not stop the run.

Accepted durations are a short test of at most `0.01`, the first target
`0.3`, and optional `1`. Late `T=3` is not a target. Default `--run` is
`0.3` for `nf` 128 and 256, step cap `0.0005`. Passing `--duration 0.3`
or `--duration 1` is the executor's production call.

When the saved step cap binds, one case of duration \(T\) and cap \(\Delta t\)
takes \(N=\lceil T/\Delta t\rceil\) steps. The analytic advance costs \(4N\)
rate evaluations and \(4N\) Jacobian-vector products. The nonlinear baseline,
which is also the saved-frame comparison, costs another \(4N\) rate
evaluations. Each of the four derivative arms and the held-out arm costs
\(4N\) rate evaluations. Those five arms also cost five dense initial
solves. Solves depend on `nf` and the occupation shift, so the two caps of
one `nf` can share them. The future effective stress is not evaluated. A
recorded one-thread sample in the FFT backend is \(0.037\,\mathrm{s}\) for
three `nf256` RK4 steps. Jacobian-vector products and the dense solves are
additional. The count is not a wall-clock certificate.

For all four saved baselines through \(T=0.3\), with caps prepared
separately, the cap-binding totals are \(1800\) steps, \(7200\) analytic
rate evaluations, \(7200\) Jacobian-vector products, \(7200\) nonlinear
baseline rate evaluations, \(28800\) finite-difference rate evaluations,
\(7200\) held-out rate evaluations and \(16\) dense derivative solves.
Sharing solves across the two caps of each `nf` reduces the derivative
solves to \(8\) and the held-out solves to \(2\).

## What this wave does not do

- It does not form `ctp_future_effective_stress`.
- It does not re-solve or edit the saved `T=0.3` source.
- It does not treat the Taylor correction as the equal-proper-time test.
- It does not install an acceptance threshold on the time or space error.
- It does not rederive the induced spectral coefficient.
- The \(2\times 2\) child kernel remains a modal covariance, not spatial
  exterior content.
- A short FFT prefix agreeing with a saved frame is a carrier check, not a
  continuum error certificate.
- Late \(T=3\) curvature stays unresolved. Optional coordinate time `1` is
  the later target after `0.3`.

## Reproduction

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_prediction.py -q
```

The tests stay at or below coordinate time `0.01`. The executor writes the
`T=0.3` record with

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_prediction.py \
  --run --output <successor> --duration 0.3 --nf 128 --nf 256
```
