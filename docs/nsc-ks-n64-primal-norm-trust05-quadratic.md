# N rise of the five-times shrink step

Best measured point is `ad759424…`, residual
`(2.61952127e-09, 7.29271610e-08)`, prediction minus measurement
`(5.79501185e-16, -6.83770152e-16)`.

The next 48-mode step, norm `1.16188004e-09`, singular ratio `1278.111`,
identity `baf90351…`, predicted `(2.61855315e-09, 7.17271610e-08)` and
measured `(2.74144727e-09, 7.17272537e-08)`. Beta followed. N rose by
`1.21926004e-10`. The N error `1.22894124e-10` is larger than the
predicted N gain `9.68119727e-13`. On this direction the quadratic
coefficient is `9.10350293e7`. Linear N gain meets that term at norm
`9.15291107e-12`, where the N gain is `7.62653068e-15`. That does not
move the residual to `3e-11`. The step is not accepted. Rank 74 and
above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_trust05_quadratic.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_trust05_quadratic.py --check
```
