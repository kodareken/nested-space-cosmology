# Accepted second clipped n=16 trial and next unit step

The [second clipped n=16 trial](nsc-ks-coupled-newton-n16-damped-iterate.md)
is accepted as a solver step only. Its measured residual matched the linear
prediction to about `1e-9` and improved both components. This owner forms the
next clipped unit-step proposal from that saved Jacobian. The physical local
incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_next2.py --check
```
