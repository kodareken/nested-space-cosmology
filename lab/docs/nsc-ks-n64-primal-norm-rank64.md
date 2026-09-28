# One rank-64 step under primal-block step control

Parent is the primal-block remeasure of the iterate6 radius. The equilibrated
truncation keeps 64 modes. The step norm is about `0.0296` and stays inside
the unit ball. The linear prediction beats both that parent and iterate6 on
N and beta. The evolution uses the same primal-block DOP853 controller. The
prediction is not a residual.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_rank64.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_rank64.py --check
```
