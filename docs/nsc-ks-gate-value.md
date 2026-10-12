# Value-only residual

A history is evolved with `tangents='zero'`. The residual adds the baseline,
the geometry value, the matter change and the one-direction edge, each once.
No tangent block is evolved. The merit is

    m(g) = max over nodes of max(|E_N|, |E_beta|).

Node counts are 47, 129 or 257. The result is not an EXISTENCE certificate.

The first full iterate6 evaluation, 47 nodes and 9 workers, took 105.3
seconds of wall time and 525.3 CPU seconds. Merit `8.02778773e-4`.
N maximum `8.02778773e-4`, beta maximum `3.61431642e-6`. This replaces
the tangent-mode timing estimate.

```sh
python3 scripts/derive_nsc_ks_gate_value.py --run --nodes 47 --workers 9
python3 scripts/derive_nsc_ks_gate_value.py --check --nodes 47
```
