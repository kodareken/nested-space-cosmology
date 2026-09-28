# Primal-block DOP853 step control

The n=64 rank-64 evolutions leave a node-constant N shift near `-3.1714e-4`
against iterate6. Damping does not move it, and beta follows the Jacobian.

Saved family solves split that shift. Iterate8 and iterate9 take the same
number of accepted steps on every family. On forty families that count equals
iterate6, and their common-mode N sums to about `-1.9e-7`. On the other
twenty, both n=64 solves take one fewer step than iterate6, and that
common-mode N sums to the whole offset. The large pieces are the high-angular
families. Same-step families do not carry it.

The difference envelope integrates `(A, D, Y)` as one DOP853 state. The RMS
error norm and the initial-step rule divide by the square root of the full
length. `Y` grows with the number of retarded directions, so the declared
`rtol` no longer controls `(A, D)`. A 128-direction solve is looser on the
primal than a 64-direction solve of the same radius. High angular eigenvalues
sit near a step boundary, drop one step, and the matter contraction turns
that primal change into N. Edge and baseline are not the carrier. This is not
a missing constant mode and not NON-EXISTENCE of `(w, U)`.

`evolve_primal_controlled` runs the same difference-envelope code object in a
private namespace whose solver factory accepts steps from the `(A, D)` prefix
only. Tangent count cannot drop a primal step. The global DOP853 binding is
unchanged, so older joint-norm registers stay replayable.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_offset.py --probe
python3 scripts/derive_nsc_ks_n64_primal_norm_offset.py --check
```

The probe evolves the iterate6 coefficients at 64 directions and the zero-padded
64-coefficient metric at 128 directions, for family `(32, -1)` and family
`(1, +1)`, under the primal-block controller. It does not take a coefficient
step. Agreement of those two matters is the cancellation of the representation
offset on those families. The published iterate6 residual remains the best
measured gate point until the same controller is used for a full family sum.
