# Rank-general physically scaled parent admission

This numerical owner retains the unchanged leading-action rate and analytic
Jv. It generalizes the sealed six-column step owner to a live source width
`R`, including nonorthonormal source columns and `R=0`. It creates no state,
trajectory, pool, checkpoint, force, damping, or runtime replacement of an
old owner. The local system has `4+4R` real components; every field component
and reciprocal source link remains present.

For `a=-8pi A`, use the inherited scales

```text
omega* = max(actual geometry-band k, AP field k, 1)
S_geo = (Q,r,2|a|r² omega*/Q,2|a|r omega*)
psi = Phi/sqrt(dx_q).
```

The helper sums absolute rows of `S^-1 J0 S`, where `J0` is the actual
leading analytic Jv with spatial operators zeroed. Field order is
`phi0.real,phi0.imag,phi1.real,phi1.imag`, each with its actual `R` columns.
The live fine-system weights enter the owned source derivative once;
the representative field generator is not multiplied by `M`.

Moving-scale rates use the **actual projected and lifted** geometry rate:

```text
g_geo = (Qdot/Q,rdot/r,2 rdot/r-Qdot/Q,rdot/r)
reaction = max row_sum(abs(S^-1 J0 S)) + max abs(g_geo)
gradient = 2 max(norm(DQ/Q,inf),norm(Dr/r,inf))
omega = max(quadrature principal k,AP field k) + reaction + gradient
dt = min(cap,1.4/max(omega,1)).
```

The absolute moving-scale term avoids signed cancellation. Growth remains
observable; the formula is a finite admission indicator, **not** a global
stability or continuous error certificate. A frozen-geometry control has
zero geometric scale velocities but retains all declared field variables.

The API is

```python
scaled_step_admission(pair, state, cap, control_mode="coupled")  # metadata dict
step_restriction(pair, state, step_cap, control_mode="coupled")  # (dt, dict)
scaled_step_restriction(...)                                    # tuple alias
stable_timestep(pair, state, step_cap, control_mode="coupled")   # dt only
scaled_local_reaction(grid, fine, system)                       # inherited tuple
```

`pair.grid`, `pair.geometry_map`, and `state` must obey the existing leading
canonical encoding, including `pi=dx_g W.T p`. The source arrays have exactly
their declared width; the helper never adds dummy columns. Fine occupations
must match that width and be finite and nonnegative. The parent preparation
owns global CAR checks: for nonorthonormal columns the nonzero covariance
eigenvalues are those of `sqrt(weights) Gram sqrt(weights)`, not the weights
alone. An orthonormal observer is distinct from the unchanged source columns.

The historical FFT grid factory still requires its original six-weight
metadata. A new preparation can resolve its carrier before replacing
`grid.fine.occupations` with the actual source weights, and before evaluating
any source. No six-column source arrays are needed by that carrier operation.
The old nested-pair constructor and six-column observation/pool wrappers are
not broadened or reused for a rank-general parent.

Validation uses in-memory manufactured nonorthonormal sources of ranks
0, 1, 2, and 6. It checks independently transformed unit-column local
Jacobian row sums, a rank-two finite directional derivative, both dense and
FFT carriers, actual moving-scale rates, frozen controls, immutable state/RHS
and weights, and exact numerical agreement with the sealed rank-six saved
fixture. These fixtures are validation inputs, not physical preparations or
scientific trajectory evidence.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q \
  tests/test_nsc_discovery_parent_step_control.py
```
