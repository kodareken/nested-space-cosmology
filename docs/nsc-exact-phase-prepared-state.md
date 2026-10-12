# Exact-phase prepared incoming state

Phase A owner for the approved local incoming sprint. The missing connection
is value accuracy of source-fixed columns under the already owned harmonic
incident law: the generic prepared adapter freezes `f(t_mid)`, while the
source phases are known exactly. This module reuses that exact phase block
on an actual compatible radius history and returns PDE node derivatives
`F_z`, `dF_z` on Sigma. CF4, full-source error, continuum bounds, metric
timesteps and constraint roots are out of scope.

## Reuse

| Input | Owner |
|---|---|
| SBP(4,2)/SAT Dirac generator | `FourthOrderModePropagator` |
| Metric direction `dL` | `generator_direction` |
| Right-inflow amplitude `B` | `incoming_reference_columns` |
| Chart cache, actual `r` | `CachedCompatibleIncomingMetric` |
| `ndir=0` | `AmplitudeOnlyMetric` |
| Restriction `R` | `rho1_restriction_map` |
| Nested state `F`, `dF`, `C_src` | `EvolvedIncomingState` |
| Source fibers/weights | `FixedSourcePreparation` |
| Exact phase diagnostic | `scripts/derive_nsc_exact_source_phase.py` (5cf59e0) |

No old hashed producer, source window, optimizer or 801/3201 field run is
executed here.

## Equation

On the existing characteristic grid, with source-fixed preparation,

\[
\frac{d}{d\tau}
\begin{pmatrix}\Phi\\ Y_j\\ Q\end{pmatrix}
=
\begin{pmatrix}
L[g] & 0 & B\\
\partial_j L[g] & L[g] & 0\\
0 & 0 & -i\,\mathrm{diag}(E)
\end{pmatrix}
\begin{pmatrix}\Phi\\ Y_j\\ Q\end{pmatrix},
\qquad
Q(t_0)=I,\quad Y_j(t_0)=0.
\]

`B` is the SAT-weighted right-boundary amplitude of the initial reference
columns. There is no free inflow callback and no interior counterforce.
`L` and `\partial_j L` are the actual midpoint metric and radius direction
on each step; the harmonic block stays inside the sparse `expm`. This is
not midpoint-held forcing (`Q'=0`, `f=f(t_{\rm mid})`).

The initial metric equals the reference and its metric directions vanish
at `t_0`. The caller supplies the unchanged-past retarded support argument;
the numerical guard checks only the initial slice and cannot establish that
continuum history from samples. The supplied preparation and incident law
are fixed, so their tangents are zero. Unsupported intrinsic families (changing `N,q,\beta,r` at
rho1) are rejected; `R` is not differentiated.

At every time node, including `t_0`,

\[
F=R\Phi,\qquad
F_z=R\bigl(L[g]\Phi+BQ\bigr),\qquad
dF=RY,\qquad
dF_z=R\bigl((\partial_j L[g])\Phi+L[g]Y\bigr).
\]

On rho1, `z=\tau+S(1)`, so these PDE values are the spatial axial
derivatives used by the existing matter contraction. Source energies remain
labels; `F_z=-iEF` is not substituted.

## API and shapes

`prepare_exact_phase_incoming(owner, times, provider, initial_fields, source)`

| Object | Shape |
|---|---|
| `state.columns` | `(nz, 2, nsrc)` |
| `state.column_tangents` | `(ndir, nz, 2, nsrc)` |
| `axial_columns` | `(nz, 2, nsrc)` |
| `axial_tangents` | `(ndir, nz, 2, nsrc)` |
| `weighted_*` | same, times `sqrt(w/(2\pi))` once |
| `state.z` | `(nz,)` equal to `times+S(1)` |
| `sampled_node_metrics` | `(nz, 4, n)` |
| `sampled_midpoint_metrics` | `(nsteps, 4, n)` |
| `sampled_node_directions` | `(nz, nfam, 4, n)` |
| `sampled_midpoint_directions` | `(nsteps, nfam, 4, n)` |
| `diagnostics['phase']` | `(nz, nsrc, nsrc)` |

`ExactPhaseIncoming.validate` / `require_history` bind the generating source,
initial columns, grid, start time, mass/ell, amplitudes, windows, **and**
the actual sampled node/midpoint metric and direction profiles. Equal
amplitudes and window radii with a different evaluated callback cannot reuse
the state. Those finite samples do not identify an arbitrary callable or
prove continuum equality.

Profile binding uses the same chart cache as the evolution, including the
full family directions in zero-tangent primal runs. It does not recompute
the uncached chart at every time node or change the finite radius.

## Status

Physical full spectrum, source-error and continuum field error: **OPEN**.
CF4 is **not implemented**; root assesses actual time error. No metric
timestep. Frozen-C0 and homogeneous exclusions remain regressions.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_exact_phase_prepared_state.py
```
