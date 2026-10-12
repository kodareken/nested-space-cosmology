# Incoming source with explicit certified component refinements

The next incoming-source approximant is the previously assembled matter
source plus three disjoint, authenticated increments:

$$
T_{\rm updated}=T_{\rm previous}
+\Delta T_{14,[16,160]}+\Delta T_{12,[40,320]}
+\Delta T_{22,[32,40]}.
$$

The first two are numerical order24-minus16 vacuum corrections, retaining
their previous thermal insertions. The third replaces the selected LOW32
approximation on one existing cell with a direct vacuum approximation and
an explicit positive thermal error bound. No term is added to
the action. The target affine-horizon source law and physical C0 are unchanged.
Original source arrays and records stay intact.

## Reuse and current decision

Reuse the existing finite source, group13 subgap replacement, approximate
tail, physical LLL state, local actions, and changed-normal reference.
The three new source increments now have matching numerical-order24
physical energy certificates. The missing connection is their application
to the same incoming source and constraint equations, with honest residuals.
Compose them once, verify independent before/after identities and unchanged
light allocation, and replay the two incoming data controls without another
mode, reference or local-action generation. Stop at this updated record;
source accuracy is still OPEN wherever a component is not bounded or its
aggregate does not meet `3e-11`.

Risks are duplicated intervals or light geometry, a bound from the wrong
numerical order, and confusing quadrature indicators with error bounds.
The owner validates group/interval identities, binds the separate certificates,
and keeps physical bounds, numerical indicators and unresolved terms separate.

The [record](../results/development/nsc-incoming-source-update.json) includes
the old and updated matter vectors, explicit increments, updated lapse/shift
coefficients and partial error budgets. Local light geometry remains in its
existing action owner, counted once. The physical group22 thermal current
is bounded rather than declared zero. Low/subgap accuracy, complete spectral
scope, a physical incoming solution, endpoint variation and extended
stationarity remain OPEN. No metric timestep or physical IV is selected.

```sh
python3 scripts/derive_nsc_incoming_source_update.py --write
python3 scripts/derive_nsc_incoming_source_update.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_source_update.py
```
