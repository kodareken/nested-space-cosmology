# Truncated assembled Chebyshev n=32 trial

The raw n=32 assembled map has condition about `1.05e14`, so the full-rank
step is refused. A column-equilibrated truncated SVD at condition cutoff
`1e6` keeps a step of norm below 1 whose linear prediction beats the
iterate4 residual on both components. The prediction is not the residual.
This owner evolves that one candidate from the same upstream source.

```sh
python3 scripts/derive_nsc_ks_coupled_newton_n32_truncated_trial.py --run --cpu-budget 7200
python3 scripts/derive_nsc_ks_coupled_newton_n32_truncated_trial.py --check
```
