# Five-times step from the accepted shrink

Parent is `ad759424…`, the accepted shrink of norm `2.58599619e-10`.
Its N prediction error is `5.79501185e-16`. The rejected ball step
`aa7e1982…` at norm `6.34883663e-09` has N error `6.05944142e-11`.
The log ratio of those errors is the power used here.

The step stays in the first 48 equilibrated modes, singular ratio about
`1.3e3`, and inside five times the accepted shrink norm. The requested
N gain is five times the power-law error at the cap, and the requested
beta gain is `1.2e-9`. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_trust05.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_trust05.py --check
```
