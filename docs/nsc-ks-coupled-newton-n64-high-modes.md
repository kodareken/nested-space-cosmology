# Next Chebyshev modes of the same (w, U)

Iterate6 is the best measured n=32 history. Its 64-mode linear image
cannot reach `3e-11`. This owner freezes that history and evolves the
matter response of Chebyshev modes 32 through 63 for both functions.

```sh
python3 scripts/derive_nsc_ks_coupled_newton_n64_high_modes.py --run --cpu-budget 7200
python3 scripts/derive_nsc_ks_coupled_newton_n64_high_modes.py --check
```

`--check` hashes every family payload, rebuilds one family from its saved
operator, and compares the rank table with the register. It does not evolve.
