# Sup-norm step from the primal-norm sup-32 history

Parent is `0e97bba5…`. Its prediction error is about `(1.43e-12, 1.08e-13)`.
The step is the minimum-L1 vector in the first 16 equilibrated modes that
lowers both predicted maxima by two hundred times that error. The singular
ratio of the block is below one hundred. Rank 74 and above stay unevolved.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_sup16.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n64_primal_norm_sup16.py --check
```
