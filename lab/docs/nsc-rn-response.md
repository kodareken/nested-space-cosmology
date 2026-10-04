# Coupled response of the radial RN step

`nsc_rn_response` differentiates the current rank-one RK4 step. The discrete
owners remain `nsc_rn_pg`, `nsc_rn_source`, `nsc_rn_observables`, and
`nsc_rn_reference`. The parent march remains
`derive_nsc_rn_parent`. This adapter adds no evolution equation, pump, or
reset, and it writes no result record.

## Baseline

The active numerical point is \(\mu=M/r_m=1.04\), \(\varepsilon=0.001\),
\(\kappa=1\), shell multiplicity 4 once, nine frequency samples, and a
rank-one column. The outer boundary is \(32 r_m\). The march step is the cap

\[
\Delta t=0.001\,r_m,
\]

so the matched half step is \(0.0005\,r_m\). On the archived 869-point
initial state, half the full CFL is about \(0.00875\,r_m\). That larger
number is not the production step.

The truncated window is

\[
E_j=\omega(1+0.9\,\ell_j),\qquad \ell_j=\mathrm{linspace}(-1,1,9).
\]

Width changes the Gaussian weights at those fixed nodes. Frequency moves the
nodes and the scattering modes. At \(12 r_m\) the vacuum ingoing speed is
\(-1.4079079416840141\) and the span is \(1.2784926817352849\). A \(+5\%\)
width or frequency stays on the truncated branch. Killing-frequency metadata
on the PG state is not that spectrum direction.

## Derivatives

The filling, with the factor 4 already inside \(\nu\), is

\[
\nu=\frac{\varepsilon M}{G_N\,4\,E},\qquad
\frac{\delta\nu}{\nu}=\frac{\delta\varepsilon}{\varepsilon}-\frac{\delta E}{E}.
\]

A normalized column obeys \(\operatorname{Re}\langle\phi,\delta\phi\rangle=0\).
Gram, CAR, and probability do not contain the shell factor. The
reconstruction weight in `stage_rates` is 1.

Unit-lapse stress uses the owner bilinears. For one column,

\[
\delta I_\pm=\delta\nu\,H(\chi^\pm,\chi^\pm)
+\nu\bigl(H(\chi^\pm,\delta\chi^\pm)+H(\delta\chi^\pm,\chi^\pm)\bigr),
\]

and the angular bilinear uses the same two-order sum. Multiplicity 4
multiplies the resulting \(F_N\), \(F_\beta\), and \(F_r\) once.

The mass and lapse tangent is `nsc_rn_pg._linearised_rates` at weight 1.
\(F_r\) does not enter that profile. The inner-mass rate keeps the owner's
second \(v F_N\) copy,

\[
\partial_t m_Q=G_N N\bigl[2v F_N+(1+v^2)F_\beta+v r F_r\bigr],
\]

with \(\delta v=\delta m/(r v)\) and \(\delta N(r_{\mathrm{out}})=0\).

The Dirac tangent uses `nsc_rn_pg._split` on each coefficient and on each
spinor, including the outer SAT. The replay tangents of the two base
scalars, which do not move the metric, are

\[
\delta\dot t=0,\qquad
\delta\dot\tau_K=\frac{N\delta N-\beta\delta\beta}{\dot\tau_K},\qquad
\delta\dot N_{\mathrm{excision}}=\delta N(r_{\mathrm{excision}}).
\]

\(\dot\tau_K=\sqrt{N_{\mathrm{out}}^2-\beta_{\mathrm{out}}^2}\) is the
static proper clock at the fixed outer node while \(g_{tt}>0\). The base
scalar called `normal` is \(\dot N_{\mathrm{excision}}=N[0]\), the legacy
integral of the lapse at the fixed excision node. There is no observer
radius on that scalar. Between the horizons its curve is spacelike, so it
is not a transported proper clock.

`rk4_tangent` repeats the four original-state stages of `rk4_step`. The
updated remaining probability is \(2\operatorname{Re}\langle\phi,\delta\phi\rangle\),
not an RK4 sum of a remaining rate. Initial \(\delta m_{\mathrm{inner}}\),
clock tangents, and stock tangents are zero.

The frozen arm holds the initial sourced lapse and shift. Their parameter
derivatives are the initial reconstruction tangents and stay in every Dirac
stage. Inner mass and the mass-flux stocks do not advance on that arm.

A station compared at fixed proper time uses the chain rule

\[
\delta O_\tau=\delta O_t-\dot O\,\frac{\delta\tau}{\dot\tau}
\]

only for a declared proper protocol. The admitted protocol here is
`outer_static`. `matched_proper_derivative` and `declared_proper_clock`
refuse the legacy excision lapse integral and refuse `transported_normal`.
That transported clock is a separate owner. It is not evaluated here.

`observer_metric_stage` hands one field stage's radius, lapse, shift, and
their tangents to that owner. Field RK4 stages start from the original
state with fractions \((0,1/2,1/2,1)\) and weights \(1,2,2,1\). The legacy
`normal` scalar is not a sample of that profile.

Width and target energy are algebraic on one prepared column. Frequency is
the centered difference of `prepare_radial_packet`, including its phase
alignment, at a step and at half that step. The call refuses a grid above
128 points, so it is not a production scattering campaign.

## API

`baseline`, `baseline_window`, `production_step`, `filling_derivative`,
`energy_initial_tangent`, `width_initial_tangent`,
`frequency_preparation_derivative`, `stress_tangent`,
`reconstruction_tangent`, `quadratic_energy_tangent`, `stage_tangent`,
`rk4_tangent`, `frozen_arm_tangent`, `clock_semantics`,
`declared_proper_clock`, `matched_proper_derivative`,
`observer_metric_stage`, `current_car_tangent`, `coordinate_energy_tangent`,
`trapping_horizon_tangent`, `displace_state`.

## Checks

One numerical thread. The focused test is

```text
python scripts/lab.py -m pytest tests/test_nsc_rn_response.py -q
```

It compares the reconstruction, one coupled RK4 step, and the frozen arm
with two centered steps on dense and CSR operators. It also checks
multiplicity, normalization, CAR, the fixed-energy initials, the archived
pilot `t0` at \(\Delta t=0.001\), and a rank-9 frequency preparation on a
17-point grid. The process CPU time of that test stays under 60 seconds.

## Limitation

The missing piece is the integrated forecast: capture, content, and delay
under the finite \(+5\%\) width and \(+5\%\) frequency controls at fixed
initial incoming Killing energy. Preparing that comparison requires the
production scattering campaign. This adapter does not run it, and a passed
tangent check is not that forecast.
