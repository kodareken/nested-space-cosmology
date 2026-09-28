# Rank-29 step from the primal-norm rank-70 history

Parent is `5e2fafb4…`, residual `(8.03677454e-09, 9.09560482e-08)`.
Equilibrated truncation keeps 29 modes. The raw step norm is about
`6.8e-8` and the singular ratio is a few hundred. Ranks below 27 miss
beta. Ranks above 70 raise the step into the ill-conditioned tail and are
not evolved. The linear prediction must beat the parent on both
components. The evolution uses primal-block DOP853 control.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_rank29.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_rank29.py --check
```
