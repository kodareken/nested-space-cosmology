# Accepted fifth clipped n=16 trial and next unit step

The [fifth clipped n=16 trial](nsc-ks-coupled-newton-n16-damped-iterate4.md)
is accepted as a solver step only. This owner forms the next clipped unit-step
proposal from that saved Jacobian. The unclipped linear leftover of the same
Jacobian is recorded so a frozen-Jacobian walk can be stopped without treating
a linear prediction as a residual. The physical local incoming gate remains
OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_next5.py --check
```
