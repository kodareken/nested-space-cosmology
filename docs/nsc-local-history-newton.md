# Local incoming history Newton CONTROL

Numerical control for a candidate [`LocalIncomingFamily`](nsc-local-incoming-family.md)
on the fixed incoming interval `I=S(1)+[.12,.18]`. This owner does not run
source evolution, a physical optimizer, the old generator, or a gate-closed
report. Every result is OPEN. A failed or budgeted step is not NON-EXISTENCE.

The production physical gate still needs both
`sup|N,beta| + error <= 3e-11` and a fresh source-fixed evolution. A small
control residual here is not that gate.

## Reuse

| Input | Owner |
|---|---|
| Two Chebyshev `(w,U)` functions, shape `(2,n)`, `n=8,16,32` | `LocalIncomingFamily` |
| Fixed `I=S(1)+[.12,.18]`, axial cutoff `.03/.06`, normal `.007/.03` | `physical_incoming_interval` / family domain |
| Compatible slots `w,w_z,w_zz,w_zzz,U,U_z` | `compatible_history_slots` |
| Local geometric N,beta map, `D[0]=0` | `surface_geometry_response` |
| Raw N,beta, full retarded Jacobian `(ndir,nz,2)` | `KSConstraintAccumulator.finalize` schema |
| Error-budget slots, missing remains `None` | `local_error_budget` |

`C_up`, the source, `A/q/Omega/zeta/V_full` and the seam stay fixed. No
`Gamma_rest`, invented stress, frozen C0 or matter-dict shortcut is accepted.
Geometry Jacobians may precondition the linear step. They do not replace the
retarded state response.

## Mixed-order free data

The included geometric operators are mixed order: `N` sees `w''` and `U`,
`beta` sees `w'''` and `U'`. A unique square collocation solve is not assumed.

The caller must declare the free jet at a chosen point of `I`,

\[
\bigl(w,w',w''\bigr)\big|_{z_\star}
=
\bigl(w_\star,w'_\star,w''_\star\bigr),
\]

or an equally explicit boundary specification. Those three numbers are
solver controls. They are recorded. They are not physical matching data.

Default collocation is the tau layout on the family's Chebyshev-Lobatto
nodes of count `n`:

- `n` N-equations,
- `n-3` beta-equations on Lobatto indices `1:n-2` (drop the left endpoint
  and the two rightmost nodes),
- 3 declared jet conditions,

for `2n` coefficients. Start at `n=8`. `n=16` or `32` is an explicit caller
choice. The solver never raises the resolution by itself.

Verification target nodes are a separate grid (default `2n+1` Lobatto
nodes). All normal and mixed derivatives still come from the same two
functions. Between-node residual is reported on that grid and is not
certified.

## Rectangular Gauss-Newton

If the residual were a local mixed-order ODE, dropping three beta nodes
would be the standard tau replacement by the jet. The actual residual also
contains retarded state response, so it is not that ODE. Rectangular
least-squares on all `n` N-nodes, all `n` beta-nodes and the same 3
declared conditions (`2n+3` equations) keeps the full retarded Jacobian.
That is the better production linearization. Both layouts are implemented
and checked on the synthetic coupled problem. Neither layout is a
between-node certificate.

## Evaluator protocol

`evaluate(family, z)` is the production API. The requested node array is
passed in and the returned `z` must equal it exactly (`np.array_equal`).
Offsets at `3e-12` are rejected; they can exceed the `3e-11` residual
tolerance at large gradient.

The mapping must contain:

- `z`, `action_gradient` of shape `(nz,2)` in N,beta order,
- `history_jacobian` of shape `(ndir,nz,2)` with `ndir=2n`,
- an explicit matching `history_identity` or `profile_identity` (the
  family key or coefficient profile). Geometry slots alone do not prove
state evolution; the identity is not fabricated from the requested family,
- `full_retarded_state_derivative=True`. This is a production-owner
  declaration, not a physics certificate. A scope that admits frozen
  matter or a removed state derivative is rejected,
- upstream/source identity: either an immutable snapshotted
  `source_identity`, or production `batch_records`/`family_records` that
  bind **both** `source_digest` and `preparation_digest` for every sorted
  batch key. Same `C_src` with a changed `A_up`/rho/mass must not pass.
  Partial missing batch bindings are rejected,
- optional `incoming_slots`, slot tangents, and
  `local_reference_gradient_tangent` as a preconditioner.

A typed `HistoryEvaluatorReceipt` copies arrays and snapshots identities
before caching. Reused evaluator buffers cannot change an earlier result.
Missing error-budget components stay `None`. A changed source identity,
mismatched nodes, mismatched or missing history proof, a dropped Jacobian,
or a frozen C0/matter dictionary is rejected. The
[KS history evaluator](nsc-ks-history-evaluator.md) supplies this protocol
by evolving each candidate from the original upstream batches.

## Newton control

Damped Newton on the tau square system, Gauss-Newton on the rectangular
system. The linear step consumes the full retarded Jacobian. Geometric
column/row scales are the preconditioner only. The initial family's
`radius_lower_bound` must be positive before the first evaluator call.
The line search also rejects nonpositive-radius trials. Armijo is in
`[0,1)`; the minimum line-search scale is in `(0,1]`. Step and iteration
budgets are finite.

A tau solve can drive the selected beta rows to tolerance while other
solve or verification nodes remain large. The result records three
separate flags:

- `control_converged` — selected control residual (tau rows + jet),
- `solve_all_nodes_pass` — every N and beta sample on the solve nodes,
- `verification_nodes_pass` — the separate verification grid.

`converged_numerical_control` is returned only when all three pass.
Otherwise a `resolution_or_verification` stop is returned. The solver
does not raise `n`. Physical status stays OPEN. Between-node remainder
is uncertified.

Stop reasons, exactly one of:

| Token | Phrase |
|---|---|
| `converged_numerical_control` | converged numerical control |
| `failed_line_search` | failed line search |
| `rank_or_conditioning` | rank/conditioning |
| `iteration_budget` | iteration budget |
| `resolution_or_verification` | resolution/verification |

None of these is NON-EXISTENCE. Stopping at the numerical tolerance does
not certify the residual between nodes or any missing error.

## API

`DeclaredJetBoundary(point, w, w_z, w_zz)`

`HistoryCollocation(coefficient_count=8, layout='tau'|'rectangular')`

`LocalHistoryNewtonSettings(...)` numerical knobs only.

`bind_history_evaluator_receipt(raw, family, expected_nodes=None, expected_source_identity=None)`

`solve_local_history(evaluator, initial_family, *, boundary, collocation=None, settings=None, verification_nodes=None, expected_source_identity=None)`

## Production adapter

[`KSHistoryEvaluator`](nsc-ks-history-evaluator.md) wraps the owned energy
propagator and `KSConstraintAccumulator.finalize`. Every candidate is
evolved from fixed upstream batches on the shared solve/verification node
set. Source and numerical configuration are snapshotted; operator cache
keys bind the analytic profile, directions, mesh, channel and solver data.
The [original-source integration control](nsc-ks-history-evaluator-control.md)
checks its retarded derivative and saved-operator replay. The physical
gate remains OPEN until the complete error bounds and source-limit bridge
close, even if the numerical control reports convergence.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_local_history_newton.py
```

## Limitations

Synthetic algorithm tests only. No NSC source campaign, PDF, metric step
or physical gate claim. Rank, line-search and budget stops remain OPEN.
The between-node remainder and source/UV error bounds are owned elsewhere.
