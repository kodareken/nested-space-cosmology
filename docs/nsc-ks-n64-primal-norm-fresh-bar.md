# Fresh directions below an N gain of 1e-10

Best measured point remains `ad759424…`, residual
`(2.61952127e-09, 7.29271610e-08)`. The gap to `3e-11` is `2.58952127e-09`
in N and `7.28971610e-08` in beta. The `7.29e-8` figure is the beta gap.
Directions are pruned when their useful N gain stays below `1e-10`.

Spent truncations 56, 64, 70, 74, 80, 88 and 94 were not re-probed. One
new direction was measured: equilibrated modes 65–73, identity
`860b60b7…`, step norm `9.60610457e-07`, subspace ratio `47.384`.
Predicted gain `(1.50000000e-13, 1.00001133e-16)`, measured residual
`(2.75068623e-09, 7.29265250e-08)`. N rose. The N error
`1.31314961e-10` exceeds the predicted N gain. This direction's own
quadratic coefficient is `142.30481286`. Its slope `1.56150705e-07`
meets that coefficient at norm `1.09729743e-09`, where the N gain is
`1.71343767e-16`. That does not reach `3e-11`. The step is not accepted.

The same block has a larger linear N gain `1.44e-10` at norm `0.01`
with beta down by `1e-12`. Under the measured coefficient the quadratic
term at that norm is about `0.014`, so that step was not evolved.
Well-conditioned blocks plateau near `1e-11` in N. No screened ray has
a useful N gain of `1e-10`. One probe, sixty families. This is not
NON-EXISTENCE. The owned Chebyshev width stops at 64 coefficients.
UV, field, low/subgap and full between-node remain `None`.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_fresh_bar.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_fresh_bar.py --check
```
