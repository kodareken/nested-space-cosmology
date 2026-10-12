# Accepted fourth clipped n=16 trial and next unit step

The [fourth clipped n=16 trial](nsc-ks-coupled-newton-n16-damped-iterate3.md)
is accepted as a solver step only. This owner forms the next clipped unit-step
proposal from that saved Jacobian. The physical local incoming gate remains
OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_next4.py --check
```
