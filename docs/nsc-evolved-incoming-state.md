# Source-fixed evolved incoming state

Douglas's adopted formulation is that the same quantum state is the same
upstream/source preparation, not a covariance pinned to Sigma independently
of the metric. For an allowed history family the incoming kernel is the
evolved restriction of the original source columns:

\[
C_\Sigma[g]=F[g]\,C_{\rm src}\,F[g]^\dagger,
\qquad
\delta C_\Sigma
=\delta F\,C_{\rm src}F^\dagger
+F C_{\rm src}\delta F^\dagger,
\qquad
\delta C_{\rm src}=0.
\]

`F` is the actual rho1 normal/half-density image of those columns. All source
channels, coherences and cross-position kernels are retained. A frozen
incoming `C0` or a legacy matter dictionary is not this state and is rejected
by the new interface. Equality of a genuine evolved kernel with a reference
matrix is not a reason to reject the state.

## Physical entry

[evolve_incoming_state](../src/recursive_horizons/nsc_evolved_incoming_state.py)
calls the owned prepared-field evolution for the supplied
`CompatibleIncomingMetric` family, then restricts the rho1 characteristic
trace with `MODE_TO_CURRENT` and the existing normal/half-density map.
The radius family leaves the intrinsic metric and normal fixed at Sigma, so
the restriction map has zero parameter derivative. An unsupported family that
changes that intrinsic map is rejected; the derivative is not omitted.

`EvolvedIncomingState` exposes `z(nz)`, `columns(nz,2,nsrc)`,
`column_tangents(ndir,nz,2,nsrc)`, `source_covariance(nsrc,nsrc)` and
`column_weights(nsrc)=repeated sqrt(w/(2 pi))`. Properties
`weighted_columns` and `weighted_column_tangents` apply that quadrature once.
`covariance()` and `covariance_tangent()` return the full cross-`z` kernels
with shapes `(nz,2,nz,2)` and `(ndir,nz,2,nz,2)`. Source energies remain
initial labels; they are not substituted as outgoing momenta `k=-E`.
Angular multiplicity is not added. Coverage records sampled-source and
continuum error as OPEN; this is not a boolean proof of a physical history
solution.

Finite arrays are stored read-only and bound to the generating source digest
and metric-family amplitudes, so a changed geometry cannot reuse an arbitrary
old matrix as the physical state.

## Bounded adapter control

The control reuses the authenticated retarded group14 `14_1` 801-node grid,
four signed source energies and compact window, through the existing
`inputs`/`support_guard` helpers. It does not regenerate baseline, scattering,
horizon or source grids. Time coordinates `[0,0.3]` and window `.007/.03` are
the same numerical controls as the retarded pilot, not a physical duration.

The new runs are three 801/64 field solves of the same compact radius family
at amplitude `0.001` and the centered pair `0.001±1e-4`. Initial, incident
and source tangents are explicit zeros by the unchanged past/upstream
support. A chart/radial cache evaluates the actual finite radius at every
time; it is checked against `CompatibleIncomingMetric` at several nodes and
nonzero metric values and is not a zero-radius freeze.

Field, incoming-trace and coherent covariance derivatives are compared with
the centered pair against `3e-8`. Flux algebra is compared against `3e-11`.
Executed residuals live in the
[development record](../results/development/nsc-evolved-incoming-state.json).
Omitting the column tangent or inserting a fake `dC_src` changes the kernel.
An array-only comparison with the saved retarded response checks the same
restriction and `dC_src=0` formula without rerunning that solver.

This is adapter verification, not a physical root. Finite source and
discretization errors remain OPEN. Full CAR, cosmology EXISTENCE, metric
evolution, new stress/Gamma and parameter refits are out of scope.

```sh
python3 scripts/derive_nsc_evolved_incoming_state.py --run
python3 scripts/derive_nsc_evolved_incoming_state.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_evolved_incoming_state.py
```

`--check` replays saved arrays and provenance and does not propagate fields.
The consumer of this state contracts actual `z` derivatives of `F`; it must
not restore a stationary `k=-E` substitution on these columns.
