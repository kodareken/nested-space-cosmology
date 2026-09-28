# Centered enclosure of the retained middle-band energy error

The current step reuses the [rank-one energy estimate](nsc-incoming-energy-error-budget.md),
the same order16 Riccati projector, and the optimized normalized derivative
jets from the [centered order24 owner](nsc-incoming-centered-order24.md).
Only the interval enclosure changes. The physical state, source values,
adiabatic reference, interval endpoints and geometric scales remain fixed.

## Decision and stopping condition

The first target is group12 on its authenticated `[40,320]` band. Its old
eight-cell energy bound is `1.69e-9` in lapse units;128 natural cells reduce
that to `1.59e-10`, still above `3e-11`. Apply the centered fourth-degree
remainder on those same128 cells, using the unchanged order16 approximation.
Stop the pilot after one verified coefficient record. Extend it to the
remaining retained groups only if that result establishes a useful budget
improvement. Do not regenerate old modes or scattering data.

The enclosure on each radial cell intersects the ordinary interval with a
Taylor polynomial about its midpoint plus the whole-cell derivative
remainder. Derivative jets contain the factorial normalization. Intermediate
derivative arrays are truncated only above the degree needed by later
recurrence stages. Geometry through degree20 supplies the order16 defect
and its fourth coordinate derivative; it does not introduce more terms in
the physical mode approximation.

Risks are interval dependency in the remainder, a genuinely broad absolute
transport-defect estimate, and unchanged low/subgap or quadrature errors.
A bound larger than tolerance remains OPEN rather than being interpreted
as actual source error or physical NON-EXISTENCE.

The [record](../results/development/nsc-incoming-centered-middle.json) keeps
the natural pilot and centered result separately and authenticates their
coefficient data. Replay contracts saved coefficients through the energy
vertex without repeating the radial calculation. Numerical-order24 group14
has a separate certificate and explicit source correction.

```sh
python3 scripts/derive_nsc_incoming_centered_middle.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_centered_middle.py
```
