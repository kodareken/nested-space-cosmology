# Beta-primary step from the factor-200 history

Parent is `3ba9f4c2…`. Its Jacobian is the primal-controlled tangent at
that history. The proportional stall only required beta to move by about
`2.7e-13`. In the first 48 equilibrated modes, singular ratio about
`1.3e3`, a minimum-L1 step of norm about `1.3e-8` lowers the predicted N
maximum by three times the last N error and the predicted beta maximum
by `1e-8`. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_beta48.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_beta48.py --check
```
