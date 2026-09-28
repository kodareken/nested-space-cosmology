# KS-envelope local incoming constraint assembly

Numerical adapter of the owned [local incoming split](nsc-local-incoming-constraints.md)
to genuine [`KSEnvelopeIncoming`](nsc-ks-source-envelope.md) families. The
local incoming gate remains OPEN; this is residual/Jacobian assembly only.
No archived science is rerun.

For the joint reference/difference prepared state, the adapter uses the
[coherent difference contraction](nsc-ks-matter-difference.md). It retains
the difference between the jointly integrated reference and the cached
reference, preserving the same raw residual while reducing cancellation.

## Reuse

| Input | Owner |
|---|---|
| Weighted columns F, F_z, dF, dF_z | `KSEnvelopeIncoming` |
| Raw N,beta matter and complete source tangent | `source_column_matter` |
| Compatible slots (w,w_z,w_zz,w_zzz,U,U_z) | `compatible_history_slots` |
| Local/reference geometry, D[0]=0 | `surface_geometry_response` |
| Error budget and residual assessment | `local_error_budget`, `assess_local_residual` |
| Homogeneous radial generator | `ks_generator` |
| Incoming intrinsic (a, r) at rho=1 | `reference_chart` / `ks_to_pg` |

No frozen incoming C0, reference-stress dictionary, PG `HistoryBinding`, or
`node_fields` is accepted. A later difference-envelope owner may subclass
`KSEnvelopeIncoming`; this adapter uses `isinstance`.

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

Baseline and geometry are added once. Families share the target incoming `z`
and the metric history. Each family may use its own numerical envelope mesh
and source digest. Signed-family multiplicity is the ledger identity

\[
\mu=\frac{\texttt{copy\_count}\times\texttt{degeneracy}}{n_{\rm actual\ angular\ signs}}.
\]

Horizon coherences in `C_src` are retained. Source energies remain labels;
evolved `F_z` is never replaced by `-iE F`.

## Computational reference

The reference of each family is generated from the same `A_up`, energies,
mass, angular and `rho_up` by the homogeneous radial ODE

\[
\partial_\rho A = G_E(\rho)\,A,\qquad A(\rho_{\rm up})=A_{\rm up},
\]

integrated with the bound DOP853 settings to `rho1=1`. Then

\[
F_{\rm ref}=A_{\rm ref}\,e^{-iEz},\qquad
F_{{\rm ref},z}=-iE\,F_{\rm ref},
\]

with identically zero history tangent. `A_ref` is cached across histories by
the full preparation digest and the bound `rtol`, `atol`, `max_step`. The
numerical `z` mesh is not part of that cache key.

## API

`assemble_local_incoming(prepared_families, provider, baseline_gradient, coefficients, channels, interval, error_budget, coverage=None)`

Each prepared family is a `KSEnvelopeIncoming` (or subclass) bound to an
actual group and `angular_sign`. Channel mass and signed ell must rebuild the
stored preparation digest from `A_up` and the actual source. Optional
`expected_source_digest` is checked per family. Duplicate family IDs,
duplicate `(group, sign)`, mismatched target `z`, and mismatched intrinsic
`(a, r)` are rejected. Distinct computational meshes are allowed.

Returns raw N,beta residuals on the target z-nodes, the full retarded
history Jacobian, per-family matter corrections, and an explicit local/error
scope. Every missing error-budget component stays `None` and is listed under
`scope['unresolved']`. Physical constraint status is OPEN.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_local_constraints.py
```
