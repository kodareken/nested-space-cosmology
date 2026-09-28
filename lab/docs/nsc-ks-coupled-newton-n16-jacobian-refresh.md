# n=16 Jacobian refresh at the fifth clipped history

The [fifth clipped trial](nsc-ks-coupled-newton-n16-damped-iterate4.md)
already evolved all sixty families at the accepted `g` and stored the full
retarded Jacobian. Re-evolving that same `g` would only repeat those
operators. This owner reuses that Jacobian together with the
[next5 leftover](nsc-ks-coupled-newton-n16-damped-next5.md).

The unclipped linear leftover stays about `1.3e-4` and `6.3e-4`, and the
rectangular condition stays about `7.05e8`. That leftover is far above
`3e-11`, so the frozen n=16 walk cannot decide EXISTENCE. The residual
Chebyshev tail on the owned 47 nodes is many orders smaller than the
leftover, so a 32-coefficient basis is not motivated. This names the
search stall. It is not scoped NON-EXISTENCE.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_jacobian_refresh.py --check
```
