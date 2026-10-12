# Local incoming-gate status after the n=16 leftover stall

This is a status register, not a physical certificate. The switched law

\[
C_\Sigma[g]=U_g C_{\rm up}U_g^\dagger=F[g]C_{\rm src}F[g]^\dagger,
\qquad \delta C_{\rm src}=0
\]

is evaluated on the declared class \(\delta r=\chi(s)[s w(z)+s^3 U(z)/6]\)
and interval \(I=S(1)+[0.12,0.18]\). The fifth clipped n=16 history has
measured residual maxima

\[
(0.004107193193038788,\ 0.005790788097514318).
\]

The Jacobian at that same \(g\) is already rebuilt from all sixty families.
Its unclipped linear leftover is

\[
(0.0001318944269235291,\ 0.0006335944911131271),
\]

with rectangular condition about \(7.05\times 10^8\). EXISTENCE requires
both residuals plus error bounds \(\le 3\times 10^{-11}\). The leftover
cannot reach that tolerance, so the n=16 walk is a named stall
`n16_damped_frozen_jacobian_leftover_cannot_reach_existence`.

The residual Chebyshev tail is about \(10^{-8}\), so a 32-coefficient
basis is not opened. Geometry between-node remainder is \(10^{-13}\).
UV tail, full between-node remainder, field accuracy and low/subgap stay
`None`. Failed optimization is not scoped NON-EXISTENCE.

```sh
python scripts/derive_nsc_local_incoming_gate_status.py --check
```
