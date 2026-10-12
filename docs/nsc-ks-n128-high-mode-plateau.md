# Stable-rank ceiling of modes 64 through 127

The column owner assembles the n=128 Jacobian at the zero-padded
`ad759424…` history. Its width-128 linear programs do not return a
feasible point. This register screens the singular vectors whose ratio
to the leading value stays below `1e6`.

The N gain reaches `9.899343e-11` at step norm `1e-4` and stays there at
`2e-4`. Beta does not fall. That is below the `1e-10` prune bar, so no
quadratic probe is evolved. n=256 is the same subclass pattern and is
not authorized by this plateau. UV, field, low/subgap and full
between-node stay `None`.

```sh
python3 scripts/derive_nsc_ks_n128_high_mode_plateau.py --write
python3 scripts/derive_nsc_ks_n128_high_mode_plateau.py --check
```
