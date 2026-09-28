# Second beta-primary step

Parent is `6534e55b…`, residual about `(2.620e-9, 7.307e-8)`, prediction
error about `(4.2e-12, 7.9e-13)`. The step keeps 48 equilibrated modes and
asks for an N drop of five times that N error and a beta drop of `3e-9`.
The step norm is about `1e-7`. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_beta48b.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_beta48b.py --check
```
