# N quadratic of the 6.5e-9 ball step

Parent `6534e55b…` remains the best measured point, residual
`(2.61985695e-09, 7.30697217e-08)`. The 48-mode step inside the Euclidean
ball of radius `6.5e-9`, identity `aa7e1982…`, norm `6.34883663e-09`,
singular ratio `1278.111`, predicted `(2.61162988e-09, 6.95697217e-08)`
and measured `(2.67222430e-09, 6.95703600e-08)`.

Beta followed the Jacobian. N rose. The N prediction error
`6.05944142e-11` is larger than the predicted N gain `8.22706759e-12`.
The quadratic coefficient on this direction is `1.50329526e6`, about 62
times the coefficient `2.43404367e4` inherited from the parent step. The
norm where this direction's linear N gain would meet its own quadratic
term is `8.61998730e-10`, and the N gain there is `1.11701123e-12`. That
does not move the residual to `3e-11`. The step is not accepted. Rank 74
and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_margin8_quadratic.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_margin8_quadratic.py --check
```
