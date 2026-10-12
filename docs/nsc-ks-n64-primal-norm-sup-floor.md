# Sup-norm floor after the factor-200 step

Parent is `3ba9f4c2…`, residual about `(2.989e-9, 8.307e-8)`.
Its N prediction error is `1.244e-10`. In the first 48 equilibrated modes
the largest joint step that stays under norm `1e-6` clears that error by
a factor of about three. The previous step missed 43 percent of its
predicted N gain, so this margin is not evolved.

A `1e-8` between-node, field, UV, or low/subgap bound cannot put either
maximum under `3e-11`. Those terms stay `None`. The next owner is the
same declared `(w,U)` class on `I`.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_sup_floor.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_sup_floor.py --check
```
