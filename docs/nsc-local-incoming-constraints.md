# Bounded local incoming constraint assembly

Local source-fixed N,beta residuals on a declared connected interval I,
using the actual [exact-phase owner](nsc-exact-phase-prepared-state.md).
The production local gate is I=S(1)+[.12,.18]. Global matching, extended
stationarity, metric stepping and observations are out of scope.

## Reuse

| Input | Owner |
|---|---|
| Weighted PDE columns F, F_z, dF, dF_z | `ExactPhaseIncoming` |
| Preparation digest from initial fields and source | `_fixed_preparation_digest` |
| Raw N,beta matter and complete source tangent | `source_column_matter` |
| Compatible slots (w,w_z,w_zz,w_zzz,U,U_z) | `compatible_history_slots` |
| Local/reference geometry, D[0]=0 | `surface_geometry_response` |
| Geometry family | `LocalIncomingFamily.metric()` |
| Coherent C_src, weights once, μ | `FixedSourcePreparation` / retained ledger |

No frozen incoming C0, reference-stress dictionary, interpolated F_z, or
authenticated 7601-node generator is executed here.

## Equation

\[
\mathcal E_B[g]
=
\mathcal E_B[g_{\rm ref}]
+\Delta\mathcal E_B^{\rm local+reference}[j_\Sigma g]
+\sum_{\mathrm{families}}
\bigl(G_B[C_\Sigma[g]]-G_B[C_{\rm ref}]\bigr),
\qquad B=N,\beta.
\]

Baseline and geometry are added once. The computational reference of each
family is the stationary phase of the **same** initial fields

\[
F_{\rm ref}(t)=R\,\Phi(t_0)\,e^{-iE(t-t_0)},
\qquad
F_{{\rm ref},z}=-iE\,F_{\rm ref},
\]

with fixed restriction R, times and signed source labels. Its history
tangent is identically zero. Current matter uses the weighted exact-phase
PDE pair (F, F_z) and (dF, dF_z), never k=-E on evolved columns.

Signed-family multiplicity is the ledger identity

\[
\mu=\frac{\texttt{copy\_count}\times\texttt{degeneracy}}{n_{\rm actual\ angular\ signs}}.
\]

Horizon coherences in C_src are retained. Source energies remain labels.

For ell=0 the pure-radius Dirac generator is independent of r, so the
history response may be omitted by that operator identity. If such a family
is included it is documented from ell=0; a small numerical tangent is not
an identity proof. Its baseline contribution is still in E_B[g_ref].

## API

`assemble_local_incoming(prepared_families, provider, baseline_gradient, coefficients, channels, interval, error_budget, coverage=None)`

Each prepared family is an `ExactPhaseIncoming` bound to an actual
group and `angular_sign`. Channel mass and signed ell must rebuild the
stored preparation digest from `diagnostics['node_fields'][0]` and the
actual source; m and ell are not public fields of the prepared state.
Duplicate family IDs, duplicate (group, sign), and mismatched
times/grids/source/group m/ell/intrinsic data are rejected.

`coverage` is the caller's declaration of the currently sampled family IDs.
Missing families and energy ranges stay explicit OPEN. Finite energy
samples are not a complete spectrum: even every retained ID cannot certify
energy coverage.

Returns raw N,beta residuals on the prepared z-nodes, the full history
Jacobian, per-family matter corrections, and an explicit local/error scope.
Physical constraint status is OPEN.

## Error budget and assessment

`local_error_budget` requires every component, each a nonnegative N,beta
pair or None:

- baseline covered regions
- baseline low/subgap
- changed-history finite energy
- changed-history tail
- preparation/field/axial derivative
- coefficient+arithmetic
- between-node remainder

Missing or None cannot become 0.

`assess_local_residual` takes the max sampled |N,beta| on a connected
positive I, adds every declared bound, and requires sample coverage of both
endpoints. Combined value ≤ 3e-11 is a **conditional numeric tolerance**.
It is never a physical EXISTENCE certificate. Zero residual with missing
bounds is not a PASS. Production remains OPEN until an independent
source/error certificate is actually supplied.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_local_incoming_constraints.py
```
