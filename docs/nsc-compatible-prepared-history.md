# Compatible history with explicit preparation derivatives

The new adapter differentiates a supplied whole-field family using the
existing Dirac operator:

\[
\dot\Phi=L[g]\Phi+f,\qquad
\dot Y_j=L[g]Y_j+(\partial_jL[g])\Phi+\partial_j f,
\qquad Y_j(t_0)=\partial_j\Phi(t_0).
\]

The initial columns and both preparation derivatives are explicit inputs.
There is no default that silently freezes them. The
[development record](../results/development/nsc-compatible-prepared-history.json)
checks the entire supplied family against independent centered changes,
including nonzero initial, incident and source-covariance variations.
This is an executable derivative connection, not a physically prepared
history or a constraint solution.

## Geometry connection

[CompatibleIncomingMetric](../src/recursive_horizons/nsc_compatible_history_geometry.py)
uses the existing reference chart and the completed surface ansatz:

\[
s=T(\rho)-T(1),\qquad z=\tau+S(\rho),\qquad
\delta r=\chi(s)\left[s\,w(z)+\frac{s^3}6U(z)\right].
\]

The explicit caller window is smooth, one near s=0 and zero beyond its
outer radius. Therefore the six incoming derivative slots are exactly
`w,w_z,w_zz,w_zzz,U,U_z`. N, beta and a in the KS chart stay fixed, while
the radius direction in PG logarithmic fields is `delta r/r_actual`.
The window does not force the incoming normal derivatives to vanish.
Nonzero normal jets imply a changed neighborhood on the parent side as
well; the adapter makes no claim of an unchanged reference past.

The window radii and finite PG time interval in the record are numerical
control choices. They do not choose a physical duration, global spatial
boundary, parent profile or surface initial tuple. Supplied functions must be smooth and their derivative-order callbacks
must return coherent derivatives of those same functions. Finite evaluations
cannot certify an arbitrary callable's regularity or derivative coherence.

## Prepared field derivative

[evolve_prepared_field_jets](../src/recursive_horizons/nsc_prepared_history_jets.py)
reuses the canonical generator and its actual-metric direction. It evolves
explicit canonical characteristic columns and tangent columns with midpoint matrix
exponentials. The incident callback supplies the already SAT-weighted
forcing and its derivative on the two right inflow rows only. The metric
and its variation retain the existing exterior boundary values. Thus
boundary speed/SAT derivatives vanish, while incident preparation
variations remain present. Interior forcing is rejected. Tangents refer to
these canonical columns; holding an underlying physical spinor fixed is
not automatically the same convention.

The finite numerical control uses two smooth geometry directions, two
synthetic source columns and four midpoint steps. Its angular label is nonzero so the radius geometry affects the Dirac
operator. Removing the geometry derivative or either initial or incident
preparation derivative changes the answer measurably. The
control tolerance is `3e-8` for finite-difference derivatives; it is not a
stationarity tolerance or a continuum error bound. The comparison tests
the derivative of the implemented midpoint scheme. The underlying discrete
flux algebra is checked against `3e-11`. No full spectral calculation or
old scientific producer is rerun.

For arbitrary supplied restriction columns, the second adapter retains the
complete coherent derivative

\[
\delta C_\Sigma=\delta F\,C_{\rm src}F^\dagger
+F\,\delta C_{\rm src}F^\dagger
+F C_{\rm src}\delta F^\dagger.
\]

The matrix control preserves off-diagonal coherence and cross-position
entries on the incoming surface. The solver records both characteristic
spin components at its exact rho1 node, initially and after each step.
The existing current-frame and normal/half-density maps convert these to
canonical incoming columns, labeled by `z=tau+S(1)`. No interpolation or
stationary `exp(iES)` factor is applied to the time-dependent fields. The
radius family leaves the intrinsic metric and normal fixed there, so this
restriction map has zero parameter derivative. A more general family
would also require that derivative.

These are actual surface traces of synthetic numerical control columns;
they do not become physical horizon/infinity modes. The finite sampled
kernel and its derivative supply neither a physical source law nor the
full continuous spectral kernel.

## Physical connection still required

For the fixed-C0 surface experiment, a physically prepared family must obey

\[
F_\Sigma[g_{\rm past}]C_{\rm src}F_\Sigma[g_{\rm past}]^\dagger=C_0.
\]

Here Sigma is rho1 with its actual normal and half-density restriction;
PG initial-time data are a different surface. The adapter now carries the actual evolved trace and fixed surface frame
map. Physical F_Sigma additionally requires horizon/infinity-prepared
columns for the changed history and the full source-frequency kernel
including mixing and complement correlations. If the source law is fixed, its derivative is
zero by that law, not because the solver omitted the input.

The preparation equation, global compatible incoming data, reference and
local action derivatives on that same family, and full extended stationarity
remain OPEN. The eight legacy diagnostic endpoint coordinates acquire no
new physical embedding. No stress, action term, fitted parameter or metric
timestep is produced here.

## Reproduction and provenance

The [execution map](nsc-compatible-history-plan.md) records the bounded gap
and ownership. Old authenticated owners remain byte-for-byte unchanged.
The new receipt hashes its adapters, reused generator/chart owners, control
producer, integration test and this document. Inspect it without rerunning:

```sh
python3 scripts/derive_nsc_compatible_prepared_history.py --check
python3 -m pytest -q tests/test_nsc_compatible_prepared_history.py
```

A default producer invocation runs the finite control when a changed
implementation or concrete derivative concern justifies doing so. It is
not a source generator or a metric evolution run.
