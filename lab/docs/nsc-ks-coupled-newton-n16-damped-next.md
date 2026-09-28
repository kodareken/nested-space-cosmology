# Accepted clipped n=16 trial and next clipped proposal

The [clipped n=16 rectangular trial](nsc-ks-coupled-newton-n16-damped-trial.md)
is accepted as a solver step only. Its measured N,beta maxima on the complete
47-node union are below the previous n=16 zero-pad history, the radius bound
stays positive, all sixty operator channels were evolved for that `g`, and the
linear prediction matched the measured residual to about `1e-9`. This is not a
physical local-gate certificate.

This owner then forms the next clipped Newton/Gauss-Newton proposal from that
measured residual and retarded Jacobian, using the existing `_clip_step` bound
and line-search scales. The selected rectangular candidate still requires a
fresh evolution. The physical local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped_next.py --check
```
