# Probe of equilibrated modes 65–73

Best measured point remains `ad759424…`. This probe is a new direction:
the min-L1 step in equilibrated modes 65–73, inside the Euclidean ball of
radius `1e-6`, asking the linear N residual to fall by `1.5e-13` and beta
by `1e-16`. Spent truncations 56, 64, 70, 74, 80, 88 and 94 are not reused,
and the five-times coefficient is not used to authorize the step.

The run measures this direction's own quadratic N coefficient. A predicted
gain is not acceptance. The gate stays OPEN unless the measured residual
later meets `3e-11` with an enclosed error.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_probe65.py --propose
python3 scripts/derive_nsc_ks_n64_primal_norm_probe65.py --run
python3 scripts/derive_nsc_ks_n64_primal_norm_probe65.py --check
```
