# Continuous-trajectory input for the residual enclosure

One 64-node capture reuses the exact upstream columns and source matrix of
the existing nonzero `14_1` control, with amplitude .001 and U=0. Its purpose
is to retain every accepted local reconstruction polynomial; the previous
endpoint-only artifacts cannot bound the continuous equation residual.
The cap is one field run and 30 CPU seconds. No source-mode or old PG run
is repeated and no resolution sweep is scheduled by this producer.

The payload stores the homogeneous reference and difference state at every
accepted endpoint, the six anchored DOP853 correction coefficients per
step, and the final prepared-state arrays. The actual PDE residual and its
first two axial derivatives are sampled at step fractions 1/4,1/2,3/4 on
an axial mesh doubled for inspection. Original source weights are included.
The verifier reconstructs these samples from the saved polynomials without
another history solve.

The sampled norms and rectangle sums are indicators, not upper bounds on
the continuous residual. Their role is to identify where a directed
continuous enclosure is needed and whether it is likely to resolve the
claimed tolerance. No continuum error component is filled with these
numbers. The physical local incoming gate remains OPEN.
