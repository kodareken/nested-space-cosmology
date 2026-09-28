# Clipped n=16 Newton proposal from the saved residual

The unclipped n=16 rectangular Gauss-Newton step is rank 32/32 with condition
about `7.05e8` and step norm about `19.85`. That candidate is not evolved.

This owner reuses the already recorded n=16 residual and retarded Jacobian. It
applies the existing `_clip_step` bound `max_step_norm=1.0` and the existing
line-search scales `1, 1/2, ..., 1/64`. Geometric preconditioning stays on.
No new Dirac evolution is performed. A linear prediction is not a residual.

The rectangular candidate is authorized for a later nonlinear trial only if a
clipped scale keeps a positive radius, avoids the forbidden full-scale and
stalled n=8 identities, and predicts both all-node maxima strictly below the
measured n=16 values. Otherwise the search stall is named and no trial is
started.

The physical local incoming gate remains OPEN.

```sh
python scripts/derive_nsc_ks_coupled_newton_n16_damped.py --check
```
