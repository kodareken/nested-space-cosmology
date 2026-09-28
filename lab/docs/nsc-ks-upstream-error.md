# Prepared upstream numerical-continuation error

The fixed upstream columns are obtained by homogeneous KS continuation from
the stored physical columns at `rho=1` to the archived binary64
`rho_up=1159676904047903/1125899906842624`. This owner supplies a directed
bound for the numerical continuation only. The physical reconstruction error
of the `rho=1` columns remains a separate missing input.

The exact KS generator is anti-Hermitian. For a continuous numerical
interpolant `Y`, Duhamel's identity therefore gives

\[
\lVert e(\rho_{up})\rVert_2
\le \lVert e(1)\rVert_2+
\int_1^{\rho_{up}}\lVert Y'-GY\rVert_2\,d\rho.
\]

The DOP853 dense polynomial and an order-eight background polynomial are
combined before norms. Bernstein convexity bounds their residual. A separate
Taylor remainder covers the exact background coefficients. A refinement
difference is recorded only as an indicator.

The family-14 pilot validates the four original energy fibers and their three
coherent columns. It also binds the result to the archived endpoint columns.
Because the same upstream error enters both reference and changed histories,
its eventual N/beta contribution must be propagated through the correlated
reference/difference equations on the final history. A full-reference stress
bound is deliberately not substituted.

```sh
python3 scripts/derive_nsc_ks_upstream_error_pilot.py --check
PYTHONPATH=src .venv/validation/bin/python -m pytest -q tests/test_nsc_ks_upstream_error.py
```
