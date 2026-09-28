# Rectangular Newton-trial finite-source evaluation

This owner evaluates the selected rectangular Gauss-Newton candidate from the
[first full-source proposal](nsc-ks-coupled-newton-proposal.md) on every
original nonzero-angular family. It restores the same frozen upstream columns
through `RetainedUpstreamArchive` and does not repeat source preparation.

The candidate is a different history `g` from the committed coupled seed. Saved
group14 operators and the previous coupled retained-control payloads are
therefore not reused. All sixty operator channels are evolved on this `g`,
grid64, degree32, and the same twenty-three solve/verification nodes. The
source-cutoff edge uses the 192-node phase rule. Previous phase-error bounds
are not transferred.

The output is a **measured** finite-source residual and retarded Jacobian.
The proposal's linear prediction is retained only for comparison. The trial
is not accepted by this producer, and the physical local incoming gate remains
OPEN: source, field, finite-tail, phase-quadrature and between-node errors
stay `None`.

Per-family operator arrays resume only for the same candidate identity, hashes
and solver configuration. A budget stop leaves orphan payloads for recovery
and never issues a NON-EXISTENCE claim.

```sh
python scripts/derive_nsc_ks_coupled_newton_trial.py --run
python scripts/derive_nsc_ks_coupled_newton_trial.py --check
```
