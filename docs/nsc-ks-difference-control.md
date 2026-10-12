# Actual-source reference/difference control

This control reuses the saved `14_1` upstream source, weights, coherent source
matrix, mass and angular label from the authenticated KS/PG comparison. It
reuses the stored 512-node total-field result; no previous field run is repeated.

The changed dependency is the algebraically equivalent
[joint reference/difference evolution](nsc-ks-difference-envelope.md).
Three runs, capped at 30 CPU seconds, compare 256 with 512 envelope nodes
and tighten the 512-node ODE tolerances. The history is the existing
`alpha=.001` Chebyshev `w` control with `U=0`. Both raw matter components
and their full retarded derivatives are compared on the 21 saved nodes of
`I=S(1)+[.12,.18]`. The predeclared value-indicator target is `1e-11`.

The stored NPZ includes the reconstructed fields, axial derivatives,
tangents, homogeneous reference and difference arrays. The JSON binds the
implementation, source owners and upstream record; `--check` reconstructs
the numerical comparisons from saved arrays without evolving a field.

Passing these comparisons establishes numerical consistency for this
sampled source and history. It does not supply continuum, full-source or
between-node bounds, nor a constraint root. The physical local gate remains
OPEN. No measured stress drift is subtracted, no scale is fitted and no
metric timestep is taken.
