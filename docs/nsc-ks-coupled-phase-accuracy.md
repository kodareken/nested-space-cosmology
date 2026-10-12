# Directed phase-value accuracy for the coupled history

This record encloses the phase derivative at every solve/verification node
of the [coupled-history pilot](nsc-ks-coupled-history-pilot.md), with the
same two functions, preparation slice and analytic profile identity.
It uses the [directed Taylor owner](nsc-dirac-source-phase-bound.md), the
original full-source angular-square weights and a ball enclosure of the
actual incoming axial scale.

The numerical edge gradient is compared **directly** with the resulting
true-gradient ball. The error includes displacement of the numerical value
from the ball midpoint, the ball radius and contraction arithmetic. A
radius around a different approximation is not transferred to the stored
value. The component target is1e-11 for both N and beta at every node.

The production32-node phase quadrature is checked first. Only if its bound
fails is192-node quadrature evaluated. This changes numerical accuracy,
not the source, history, action or physical parameters. No source or field
evolution is performed.

This is a nodal phase-value component. Source/field errors, the finite
source tail, phase-tangent accuracy and the between-node remainder are
separate. The local physical gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_phase_accuracy.py --check
```
