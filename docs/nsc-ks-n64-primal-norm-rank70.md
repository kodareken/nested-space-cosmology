# Rank-70 step from the primal-norm rank-64 history

Parent is `544e2e6a…`, residual `(9.09431179e-09, 8.67931722e-07)`.
Equilibrated truncation keeps 70 modes. The raw step norm is about `0.0106`,
inside the unit ball and below the `0.03` step that stayed linear. Modes past
70 raise the step to `0.30` and then `0.79` without a singular gap, so rank 80
is not evolved. The linear prediction must beat the parent on both components.
The evolution uses primal-block DOP853 control. A prediction is not a residual.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_rank70.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_rank70.py --check
```
