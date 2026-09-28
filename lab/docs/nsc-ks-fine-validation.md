# Continuous validation of the saved fine history

This evaluator validates the saved 2048-node trajectory one accepted rho
step at a time. It checks source, analytic profile, period and data bindings,
then combines the actual field polynomial, time-coefficient enclosures,
certified profile convolution and all omitted remainders. No field or source
equation is evolved during validation.

Each completed segment proof is immutable and hashed. Interrupted work
resumes from completed proofs with exactly the same code/input signature;
it does not restart completed cells. The pilot is segment122, comparable
in rho to the earlier coarse control. The full run covers all246 saved
steps. Proof cells may subdivide in time when the flat cutoff or scalar
Taylor remainder requires it; the original trajectory remains unchanged.

Residual integral bounds add over disjoint rho cells. Only after every cell
is covered is the local characteristic estimate evaluated. Its zero initial
numerical error is relative to the same saved upstream columns, whose
physical preparation accuracy is a separate source requirement. A numerical
field-error bound is not a physical source or incoming-constraint certificate.

`--pilot` tests the fine representation, `--run` continues all missing cells,
and `--check` replays bindings, coverage and error accumulation. The CPU
budget stops dispatching new cells after a completed cell; completed work
is preserved. No partial sum is promoted to a full-history bound.
