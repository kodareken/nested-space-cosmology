# Modes 64 through 127 of the same (w, U)

`ad759424…` is the best measured n=64 history. Its 128-direction linear
image does not clear an N gain of `1e-10`. This owner zero-pads that
history to 128 Chebyshev coefficients per function and evolves the matter
response of modes 64 through 127 for both functions.

Modes 0 through 63 stay the primal-controlled Jacobian already measured at
`ad759424…`. The new columns use the same primal-block DOP853 control.
The frozen history must reproduce each parent family matter change within
`1e-12`. A Newton step is not authorized from the linear screen alone.

`LocalIncomingFamily128` is a subclass. The base allowlist remains
8, 16, 32 and 64 coefficients so earlier source hashes stay valid.

```sh
python3 scripts/derive_nsc_ks_n128_high_modes.py --preflight
python3 scripts/derive_nsc_ks_n128_high_modes.py --run --cpu-budget 14400
python3 scripts/derive_nsc_ks_n128_high_modes.py --check
```

`--check` rebuilds one family from its saved operator and compares the
screen, the assembled Jacobian and the drift bound. It does not evolve.
UV, field, low/subgap and full between-node stay `None`.
