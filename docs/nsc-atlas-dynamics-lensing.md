# Predicting lensing from finite freely falling motion

Status: the precheck, timelike calibration, frozen NSC-P004 forecast and
independent optical validation are complete. The preparation below is
retained as written before the target ray was evaluated.

## Preparation and frozen test design — 26 October 2026

Objective: use motion in the already formed source to predict a previously
unevaluated image direction, without fitting an independent lens mass.
The [formed observer](nsc-atlas-formed-observer.md) supplies the metric,
field derivatives and authenticated short source interval. Formation and
held-clock checks remain closed.

Three routes were compared: an equilibrium rotation curve, an adjustable
lens profile, and finite timelike velocity impulses. The third retains
the breathing field and makes every calibration input clear.
First inspect the force's variation time relative to a dynamical time;
do not assume circular orbits where that condition fails.

Roadmap: add nsc_atlas_dynamics_lensing.py around the existing SourceSlab,
with timelike worldlines (radius, momentum, proper clock), release-impulse
calibration and a data-only Born prediction. A two-stage producer first
writes an immutable forecast, then separately evaluates the full null ray.
The existing prediction register records NSC-P004 before validation. The
manuscript and scientific owners receive the result after the check.
No new physical coupling, source evolution or external library is needed.

Fixed before calculation:

- Same formed 10% source, A*=1.2, inherited action, weak metric and parent.
- New receiver azimuth Theta=pi-.12; it was absent from the previous
  optical fan. Initial source/receiver comoving radii .25. Both freely
  fall from zero peculiar velocity at the existing slab start; clocks
  start there together. Emission proper age is the previous homogeneous
  reference age at A=1.2.
- Calibration consists only of timelike motion and background geometry.
  Radial test probes are released with zero peculiar velocity at the
  retarded reference-chord events; centered velocity changes over
  log-scale half-width5e-5 infer the local force. Controls use half that
  interval, doubled quadrature and the retained source refinements.
- Two Gauss panels with24 nodes each sample the force along the chord;
  the refined calibration uses48 nodes per panel. A release-speed
  error of +/-1e-6 is tested per node and summed as an absolute weighted
  sensitivity envelope. It is a declared state uncertainty, not survey noise.
  Source/receiver trajectories are calibrated inputs in this synthetic
  test; their required measurement precision is reported separately.
- Predict the moving observer's image direction using the finite-endpoint
  Born kernel, calibrated source/receiver trajectories and exact radial
  aberration. No withheld optical number enters this stage.
- Freeze the forecast before evaluating the exact null boundary problem
  in the same metric. The comparison tolerance is the numerical/state
  envelope plus a predeclared0.5% of the lens contribution for the Born
  approximation. This is a test tolerance, not a full-GR confidence interval.

Expected before calculation: an attractive lens contribution of order1e-7
radians, with a motion correction that must be calculated. The freely
falling frequency signal may be dominated by endpoint Doppler terms.
No numerical value for the new target ray has been evaluated or supplied.

Criteria: authenticate the retained field, recover analytic timelike and
boost limits, show calibration convergence, preserve the frozen forecast,
and compare its new image prediction with the independent null solution.
Keep synthetic validation distinct from empirical confirmation. Budget:
300 CPU seconds, four useful workers, no new formation grids.

Failure forecast: (1) force evolution invalidates equilibrium or an overly
long impulse, visible in the timescale and half-interval control; (2)
small differences are unresolved, visible in source and integration
comparisons; (3) release-state uncertainty overwhelms the lens prediction,
visible in the signed velocity controls. Any failed comparison is reported
without retuning the lens. A later real measurement with the same defined
geometry, preparation and error model could falsify the tested relation.


## Forecast checkpoint

Before optical evaluation, the calibration predicts19.43948287 mas with
comparison tolerance .18327870 mas. The source-force timescale check gives
instantaneous change-per-nominal-orbit ratios about3.84 at r=.0073 and2.87
at r=.015. These do not support an adiabatic circular-orbit assumption.
The [forecast](../results/summaries/nsc-p004-dynamics-lensing-forecast-v1.json)
is creation-only, with the new receiver and tolerance fixed. Its full
timelike calibration is preserved; no target optical value was supplied.

## Exact release law and the finite motion measurement

Retain the preceding weak metric. Write alpha=sqrt(1+2Phi),
s=A L sqrt(1-2Phi), B=1-C L^2 r^2. For a unit-mass radial test body with
canonical momentum P, the physical-time Hamiltonian is

\[
\mathcal H_t=\alpha\sqrt{1+BP^2/s^2}.
\]

The implementation integrates its timelike Hamilton equations and proper
clock using N=ln A. With a fixed momentum unit ell, p=P/ell,
G=B ell^2/s^2 and E=sqrt(1+Gp^2), these are

\[
r_N=\frac{\alpha Gp}{H\ell E},\qquad
p_N=-\frac{\Phi_r E/\alpha+\alpha G_rp^2/(2E)}{H\ell},
\qquad \tau_N=\frac{\alpha}{HE}.
\]

The measured local radial velocity relative to the cosmological orthonormal
frame is v=sqrt(B) P/(s E). At a release with v=0, the geodesic equation gives

\[
-\left.\frac{dv}{d\tau}\right|_0
=\frac{\sqrt B\,\Phi_r}{A L\sqrt{1-2\Phi}(1+2\Phi)}.
\]

The body's proper four-acceleration is zero. The quantity here is the
change of its velocity relative to the stated reference frame, not a
nonzero accelerometer reading carried by a freely falling body.

Each calibration probe has zero peculiar velocity at its own central
event. Its before/after velocities and elapsed proper time give a centered
estimate of that derivative. Initial source/receiver worldlines are a
different pair: they are released once at the slab start, then keep falling.
Their positions, velocities and clocks are recorded as calibration inputs.
Spherical symmetry gives the same radial trajectory at the same time,
but emission and reception occur at different times, so their endpoint
radii are not set equal.

To leading weak order the reconstructed force is

\[
\Phi_r^{\rm rec}=-\frac{A L}{\sqrt B}
\frac{v_+-v_-}{\tau_+-\tau_-}.
\]

It uses motion and geometry only. The exact potential-dependent factor is
not silently read from the source to improve this prediction. The finite
impulse bias is checked by halving its interval. The calibration refinement
uses96 probes, split between the incoming and outgoing halves of the
retarded reference chord. It is a motion profile, not a claim that one
velocity measurement determines the entire lens.

## A finite-endpoint light prediction from those motions

Let the background chord run from z_e to z_o with length d, impact b and
r(z)=sqrt(b^2+z^2). For small transverse displacement eta, fixed endpoints
and equal leading metric potentials give

\[
\eta''=-2\partial_\perp\Phi,\qquad
\eta(z_e)=\eta(z_o)=0.
\]

Integrating this boundary problem yields the normal-frame image-angle shift

\[
\Delta\beta_{\rm lens}
=\frac{2}{d}\int_{z_e}^{z_o}(z-z_e)\,
\Phi_r^{\rm rec}(N(z),r(z))\frac{b}{r(z)}\,dz.
\]

The endpoint weight is necessary; the total deflection along a ray is a
different observable. N(z) follows the homogeneous conformal travel time,
not the withheld lensed trajectory. The calibration therefore samples the
right time-dependent force without using the optical answer.

Known source/receiver positions determine the curved-background geometric
angle. A receiver radial speed v_o then gives the exact local aberration

\[
\cos\beta_{\rm obs}=
\frac{\cos\beta_{\rm normal}-v_o}{1-v_o\cos\beta_{\rm normal}}.
\]

The data-only prediction function accepts the recorded impulses, endpoint
kinematics and geometry. It has no source-potential or optical-result
argument and introduces no separate fitted lens mass. The equal-potential
relation is a premise of the retained metric. In a model with distinct
time and spatial potentials, dynamics and light would probe different
combinations; this test does not derive that equality afresh.

## Prospective synthetic result, with its forecast unchanged

The forecast was committed and pushed in
[e626f2f](https://github.com/kodareken/nested-space-cosmology/commit/e626f2f)
before the target null problem was evaluated. That checkpoint contains the
new implementation, frozen calibration and prediction register entry.
The forecast's base source-commit field refers to the preceding23d76bf
source context; its producer file hashes are preserved by the checkpoint.

| Quantity | Milliarcseconds |
| --- | ---: |
| Motion-calibrated forecast | 19.43948287 |
| Frozen comparison tolerance | .18327870 |
| Independently evaluated moving-endpoint null ray | 19.43946116 |
| Null result minus forecast | -.00002171 |
| Null-grid comparison range | .08049 |

The test passes without changing a source parameter, lens normalization,
prediction or tolerance. The discrepancy plus the null-grid comparison
range also stays inside the original window. That small discrepancy is
an internal closure residual; it does not replace the much larger
calibration, numerical or physical uncertainty. Both channels are simulated
in the same prescribed metric. This is not an independent astronomical
confirmation or a unique identification of nesting.

The frozen tolerance is the sum of3.9020e-10 radians of numerical
comparison,2.3549e-11 of probe-release sensitivity and4.7481e-10 of
predeclared Born allowance. The latter is0.5% of the lens contribution.
The allowance is a declared comparison accuracy target, not a rigorous
bound on all omitted GR terms. Tests preserve the forecast's byte identity
and reject changing its value in the register.

The normal-frame lens contribution is9.49624e-8 radians. Calibrated
endpoint geometry and observer aberration add a net -7.17128e-10 radians.
Thus the prediction includes a motion correction rather than applying
the old held-clock value to moving observers.

## What the moving frequency and stationarity checks establish

The source and receiver velocities at their respective events are
-4.27921e-8 and -6.67110e-8 in units of c. Their exact radial Doppler
contribution to log frequency is1.09306e-7. The expansion-removed normal-frame
gravitational term is -2.93250e-12. The motion contribution dominates that
small gravitational term; dropping it would change the problem.

With local photon direction cosine mu at each endpoint,

\[
\log\frac{\nu_o}{\nu_e}
=\log\frac{\nu_{o,\rm normal}}{\nu_{e,\rm normal}}
+\log\frac{\gamma_o(1-v_o\mu_o)}
{\gamma_e(1-v_e\mu_e)}.
\]

The moving-clock arrival identity still holds:
d(tau_o-tau_e)/d tau_e=1/nu_o-1 for unit emitted frequency.
The finite difference at proper-age step1e-4 differs from that prediction
by3.96e-13. Halving the step gives5.52e-12 instead, so no monotonic
convergence claim is made for this cancellation-sensitive auxiliary
derivative. Both are small relative to the1.12e-4 flight-time derivative.

The force-timescale diagnostic is
|partial_t log(g)|/Omega, with g the source acceleration at fixed physical
radius and Omega^2=g/R. Its values are .611 near r=.0073 and .457 near
r=.015. Multiplying by2pi gives order-one changes per nominal orbital
period. This is a local timescale test, not an integrated orbit theorem;
it does not justify an equilibrium rotation curve in this breathing phase.
The finite timelike method does not require that assumption.

For observational use, endpoint kinematics cannot be treated as exact
without measurement. Around beta=.06, an observer velocity error delta v
changes the angle by approximately sin(beta) delta v. Keeping this single
error contribution below the frozen tolerance requires roughly4.4 m/s
velocity precision if everything else is fixed. Position/distance and
initial-state errors require their own propagation. This is a requirement,
not a claim that such an astronomical calibration has been achieved.

## Attribution, limits and reproduction

Combining lensing and dynamics is established practice; see
[Reyes et al.,1003.2185](https://arxiv.org/abs/1003.2185) for the distinct
large-scale statistical approach. Moving-lens frequency/direction effects
are treated in [Korzynski--Kulejewski,2601.13912v2](https://arxiv.org/abs/2601.13912v2).
Their abstracts and bibliographic records were checked. No survey data or
statistical estimator from those papers is imported into this finite test.
The release law, moving clocks, Born kernel and boosts are derived here
for the retained NSC metric.

Grok independently checked the release law, centered impulse and weighted
Born kernel, explicitly without computing the target ray. It also warned
against interpolating the old optical fan. The implementation uses no
optical fan values for its calibration. Independent tests check cosmological
momentum redshifting, the release-law limit, aberration inversion and the
data-only predictor. The original source and its physical assumptions
remain unchanged.

The calibration covers the inner chord, not a direct measurement of the
whole compensating reservoir. An alternative source with the same metric
on these probe paths would give the same result. The next outward test
should investigate whether measurements of the surrounding compensation
and its inherited scale relation break a specified same-core degeneracy.
It should not claim that this universal weak-gravity consistency relation
already distinguishes nesting.

The [forecast](../results/summaries/nsc-p004-dynamics-lensing-forecast-v1.json)
and [validation](../results/summaries/nsc-p004-dynamics-lensing-validation-v1.json)
are distinct records. Preparation used about11 CPU seconds and validation
about15, in four workers; raw formation data were reused. To reproduce:

    python scripts/run_python.py scripts/replay_nsc_p004_calibration.py
    python scripts/run_python.py scripts/check_nsc_atlas_dynamics_lensing.py --validate
    python scripts/run_python.py scripts/check_prediction_register.py
    python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_dynamics_lensing.py tests/test_prediction_register.py

The calibration replay uses a verified temporary workspace directory and
leaves the frozen bytes untouched. Validation requires the forecast to be
committed first. Regenerate the preceding seed-collapse raw fields if
absent; no full formation evolution is hidden in either new stage.
