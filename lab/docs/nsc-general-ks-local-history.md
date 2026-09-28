# General-KS local induced history

`GeneralKSLocalInducedHistory` is the local-action component required by the
extended tilted/frequency-mixing history gate.  For every supplied homogeneous
KS history it returns separate nodal gradients and action forces for

$$
g^A(\tau)=(N,\beta,q_{\rm ADM},r).
$$

It imports the locked four-dimensional ledger

$$
S_{\rm local}=-\int\sqrt{|g|}\,[A R+C_FF^2+C_WC^2+C_EE_4+C_\Box\Box R]
$$

with $V_{\rm full}=0$.  The Einstein channel uses the existing GHY-completed
spherical action and the monopole Maxwell channel remains a bulk term.  The
single locked $C_W$ is the compact/local allocation; the stored unit-radius
stress is not added again.

On this homogeneous time-history domain the Euler and box-$R$ densities are
endpoint functionals,

$$
4\pi\sqrt{|g|}E_4=32\pi\,\partial_\tau
\left[\frac{\dot q_{\rm ADM}}N
\left(1+\frac{\dot r^2}{N^2}\right)\right],
$$

$$
4\pi\sqrt{|g|}\Box R=4\pi\,\partial_\tau
\left[\frac{q_{\rm ADM}r^2}{N}\dot R\right].
$$

They are therefore retained at both temporal endpoints and assigned no
invented bulk source.  The homogeneous scalar action is independent of the
shift; the returned beta gradient is zero at every node.

The owner differentiates the discretized action at every history node.  A
compactly supported directional finite difference checks those gradients,
rather than replacing them by constant endpoint forces.  The authenticated
control history is only a verifier probe and does not select a physical
duration or profile.

The existing compact boundary owner is scoped to reflecting internal
$Y=\pm L_\star$ endpoints, explicitly not the parent-child transmitting
interface.  It is not reused as a temporal Weyl boundary completion.  The
record consequently exposes the nonzero free-endpoint gradient of the Weyl
bulk action as an OPEN residual.  This owner is composable, but the extended
existence gate remains closed until the transmitting tilted-interface owner
supplies its declared endpoint variation and the fourth-order owner is also
present.

Reproduce only this record with

```sh
python3 scripts/derive_nsc_general_ks_local_history.py --check
```

No source generator, optimizer, metric timestep, stress/null reconstruction,
or coefficient refit is run.
