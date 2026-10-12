# Actual prepared-state response on the local incoming interval

This control uses the resolved7601 grid and unchanged initial columns from
the local reference-accuracy record. It replaces the reference geometry by
one explicit member of the approved local two-function family: w is alpha
times Chebyshev T1 times the flat axial plateau; U=0. The incoming interval
is I=S(1)+[.12,.18], with cutoff support S(1)+[.09,.21] and normal .007/.03.
These are numerical coordinates, not a cosmological duration. The unchanged
past support follows from these exact windows and the owned retarded chart.

The same exact-phase prepared-state owner evolves full fields and retarded
tangents, then restricts F,F_z,dF,dF_z from the actual node PDE. The original
group14_1 covariance, four signed energies, three fibers and multiplicity12
are retained. Canonical source weights are applied once. No incoming C0,
free stress, normalization, or history-dependent incident source is inserted.

Reuse the geometric coefficients and their analytic Jacobian, and form the
partial two-constraint diagnostic as baseline+geometry+G(current)-G(reference).
The stationary computational reference comes from the same initial columns;
it has zero parameter tangent. Its analytic -iE derivative is used only for
that reference. Current columns always use the actual PDE derivative.

Four solves are predeclared: alpha=.001 and its centered pair +/-.0001 at60
midpoint steps on[0,.18], then alpha=.001 at120 steps. Center runs carry one
actual direction, and finite-difference runs have zero tangent count. Total
CPU cap300seconds. Each completed case is checkpointed; no extra case starts
after a failure or cap. The local full constraint-derivative tolerance is3e-8.
The time-change indicator target for raw matter is1e-11; that comparison is
not a continuum error bound. Its result determines whether the approved
higher-order time method is needed.

Every physical gate remains OPEN: one four-energy family cannot establish
full source accuracy, complete changed-history spectral coverage or a uniform
residual between nodes. The supplied alpha is a response control, not a
constraint solution. No measured reference drift is removed from matter.

```sh
python3 scripts/derive_nsc_local_prepared_response.py --run
python3 scripts/derive_nsc_local_prepared_response.py --check
```
