# High-energy source error through infinity

For all32 retained non-LLL groups, join the finite high interval to its
infinite tail at the exact owned endpoint:

\[
\epsilon_{\rm high\to\infty}
\le\epsilon_{\rm finite\ high}+\epsilon_{\rm tail}.
\]

The [finite-high receipt](nsc-incoming-source-update-v3.md) already includes
physical mode, positive thermal, numerical integration and source-sum
rounding errors. The [all-tail receipt](nsc-incoming-tail-quadrature-batch.md)
now bounds the numerical integral of the unchanged archived tail, its
physical mode error and outward thermal remainder. Add that combined tail
bound once; do not add the old physical-tail-only budget again.

The [joined record](../results/development/nsc-incoming-source-update-v4.json)
authenticates both inputs and all payloads, checks exact endpoint coverage
and binds the tail certificate to the source values already assembled.
Groups13/14 are covered from16 to infinity; the other groups from32 to
infinity. Their internal join lies at160 or320 according to the original
channel inventory. No source value, local action or reference response changes.

The combined lapse source-error bound is below `8.917e-12`, against `3e-11`.
The positive tail shift bound is retained with lossless interval endpoints
despite being smaller than binary floating-point range. The source-sum
rounding term dominates the covered shift budget, below `3.89e-18`.
These bounds do not certify floating arithmetic of the entire joint action
or a physical constraint residual.

Unbounded source contributions remain at energies below the declared
thresholds, with subgap, preparation, propagation and quadrature errors.
The matched horizon
initializer is a separate prerequisite for new validated low-band fields;
its endpoint certificate does not itself replace a spectral source integral.
Physical pressure errors in the finite high band also remain outside the
energy/current certificate. Full-source error and extended stationarity
therefore remain OPEN, with no physical initial data or metric evolution.

Decision: reuse the completed finite-high and tail certificates, close only
their missing disjoint-union error budget, and stop after authenticated
coverage/replay. This join runs no scientific producer.

```sh
python3 scripts/derive_nsc_incoming_source_update_v4.py --check
python3 -m pytest -q tests/test_nsc_incoming_source_update_v4.py
```
