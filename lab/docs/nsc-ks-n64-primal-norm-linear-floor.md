# Linear floor after the primal-norm rank-29 history

Parent is `c2404e36…`, residual `(3.35476406e-09, 8.31024251e-08)`.
The rank-29 prediction error was `(−2.05637711e-11, −1.22887853e-12)`.

Well-conditioned truncations no longer move both components by more than
that error. Rank 27 is the first that barely clears it, by about `2.7e-14`
on N. Rank 74 and above clear it with step norms `0.30–0.79` and singular
ratios from `6.7e8`. A rank-74 step damped to norm `0.01` predicts an N gain
of `8.8e-11`, inside the nonlinearity already seen at that norm. Those steps
are not evolved.

This is a same-class search stall, not a NON-EXISTENCE certificate.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_linear_floor.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_linear_floor.py --check
```
