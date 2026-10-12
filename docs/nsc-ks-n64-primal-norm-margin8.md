# Quadratic-margin step from the accepted beta-primary history

Parent is `6534e55b…`, residual `(2.61985695e-09, 7.30697217e-08)`.
Its Jacobian is the primal-controlled tangent at that history and is not
recomputed. The parent N prediction error `4.21773240e-12` at step norm
`1.31636197e-08` fixes the quadratic coefficient used here.

In the first 48 equilibrated modes, singular ratio about `1.3e3`, the
minimum-L1 step is drawn into the Euclidean ball of radius `6.5e-9`. The
requested N gain is eight times that coefficient times the square of the
ball radius, and the requested beta gain is `3.5e-9`. The realized step
must keep the predicted N gain at least eight times the coefficient times
the square of its own norm. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_margin8.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_margin8.py --check
```
