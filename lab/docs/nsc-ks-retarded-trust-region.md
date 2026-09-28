# Corrected retarded trust-region primitive

The production local gate uses 129 residual nodes and 64 history directions,
so its Jacobian is rectangular: 258 equations by 64 coefficients. This owner
provides the physical-metric Levenberg step, measured-value acceptance,
rank-one Broyden update and refresh policy independently of any campaign.

The trust metric is the existing geometric principal scaling of
`H^3(w) + H^1(U)`. A candidate is accepted only after a fresh value evolution
reduces the joint maximum norm and preserves positive radius and preparation
support. Predictions are never residuals. The Jacobian is refreshed after
five accepted Broyden updates or two poor prediction ratios.

This module does not launch a scientific run. The single production driver
will bind it to the corrected DOP853/grid1024/degree48 evaluator only after the
certification-method feasibility gate permits the nonlinear search.
