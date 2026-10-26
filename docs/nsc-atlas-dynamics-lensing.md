# Predicting lensing from finite freely falling motion

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
