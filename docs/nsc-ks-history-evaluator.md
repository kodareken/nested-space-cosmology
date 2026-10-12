# Source-fixed KS history evaluator for local incoming Newton

Production adapter from a candidate [`LocalIncomingFamily`](nsc-local-incoming-family.md)
to the owned KS residual/Jacobian schema used by
[`solve_local_history`](nsc-local-history-newton.md). The incoming law remains

\[
C_\Sigma[g]=U_g C_{\rm up}U_g^\dagger
\]

with the full retarded state variation in both N and beta. `C_up`, `C_src`,
weights, mass, ell, `rho_up`, `A/q/Omega/zeta/V_full` and the seam stay
fixed. This is a finite retained-source CONTROL. Missing error bounds stay
`None`. No physical EXISTENCE or NON-EXISTENCE claim is made.

## Reuse

| Input | Owner |
|---|---|
| Frozen upstream `A_up`, `C_src`, weights, mass, ell, `rho_up` | `KSUpstreamBatch` |
| Two identity columns per energy, original source reconstructed on apply | `evolve_energy_propagator` |
| Opposite-angular negative family, explicit `C_neg`, S3 map, no factor two | `negative_angular_partner` / `map_negative_batch` |
| Baseline and local/reference geometry once | `KSConstraintAccumulator` |
| ell=0 groups 0,13,23 keep their baseline allocation | `declare_radius_invariant_channel` |
| Analytic profile/directions identity | `profile_identity` |
| Receipt contract, `full_retarded_state_derivative=True` | `bind_history_evaluator_receipt` |

No frozen incoming C0, invented stress, `Gamma_rest`, global matching or
metric timestep is accepted. LLL columns are not invented.

## Target grid and cache

The operator is evolved on the strictly increasing union of the solver nodes
and the held-out verification nodes. `evaluate(family, z)` may return an
exact subset of that union, in the requested order. A node offset of `3e-12`
is absent and is rejected.

Every new `g` is evolved with all `2n` retarded history directions. The
operator cache key binds the complete analytic profile and directions, the
computational mesh, interpolant interval/degree, solver options, `rho_up`
and the `(mass, angular)` channel. Amplitudes alone are not a key: two
metrics with the same amplitude vector and different `w`/`U` functions do
not share an operator.

Construction snapshots the already prepared source arrays, baseline,
coefficient arrays, error inputs and channel ledger. This does not repeat
source preparation. Caller-owned buffers and nested returned receipts
cannot change an earlier cached computation or its provenance.

Positive and negative signed families are grouped first. One operator
evolution is reused for a negative partner only when the owned S3 /
opposite-angular identity actually applies to that frozen batch. An
opposite-angular positive-energy family is a different channel.

## API

`KSHistoryEvaluator(batches, *, baseline_gradient, coefficients,
energy_interval, interpolant_degree, z_grid, solve_nodes,
verification_nodes, rtol, atol, max_step, ...)`

`batches` are already bound `KSUpstreamBatch` objects, prepared once.
`evaluate(family, z=None)` returns raw N,beta, the full retarded Jacobian
`(2n, nz, 2)`, the separate geometric preconditioner tangent, immutable
source/preparation identity, and `profile_identity` compatible with
`LocalIncomingFamily`. Mapped negative covariances are the supplied
`C_src`, never `I-C_positive` and never a doubled positive result.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_history_evaluator.py
```

## Limitations

Interpolation, node-evolution, preparation/field and source-cutoff /
coincidence errors remain `None`. The local physical gate stays OPEN until
those bounds close, even if a Newton control residual is small. This owner
does not run an all-family campaign or a metric step.
