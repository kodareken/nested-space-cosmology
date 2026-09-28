# Sup-norm step from the primal-norm rank-29 history

Parent is `c2404e36…`. Its Jacobian is the primal-controlled tangent
assembled on that history. Least-squares truncations of it move both
maxima by about the last prediction error. The step evolved here is the
minimum-L1 coefficient vector in the first 32 equilibrated right singular
modes that lowers both predicted maxima by ten times that error. The
singular ratio of that block is a few hundred and the step norm is a few
times `1e-9`. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_sup32.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_sup32.py --check
```
