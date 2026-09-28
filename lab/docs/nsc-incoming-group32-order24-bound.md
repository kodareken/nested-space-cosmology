# One group32 order24 certificate for two existing energy windows

Reuse the frozen retained-channel order24 wrapper, centered degree4 enclosure,
128-cell collar map, and rank-one energy projection. The missing connection
is a rigorous group32 error bound for the explicit order24 projector on the
actual LOW cell `[32,40]` and middle interval `[40,320]`. Their separate source
corrections remain independently owned. No archived order16 value acquires
this certificate until its matching explicit source correction is composed.

Calculate radial coefficients once for group32, at order24, 128 cells and
Taylor depth4. Contract the identical saved coefficients over each window
and the union `[32,320]`. Compare the lapse/shift component bounds with the
unchanged `3e-11` tolerance. Stop after this single certificate, replay and
focused verification; do not refine the partition or order a second time.

## Authenticate windows and reuse energy-independent coefficients

The existing middle-input owner authenticates the selected signed middle
panels. The LOW cell is verified against the selected `low_ref/32_+/-1`
metadata, edges, positive weights, and exact saved source nodes. The common
boundary40 and original upper endpoint320 are retained.

Radial integrals `I24,...,I48` and local cross-product coefficients contain
no energy cutoff. The legacy radial field `lower` originates as a tail
endpoint and is also checked by the energy-contraction interface. Each
window therefore receives an explicit metadata adapter: preserve
`original_radial_endpoint=320`, set that interface field to the window's
right endpoint, and record the exact contraction interval. The coefficient
arrays, order, horizon state, geometry, and multiplicities remain identical.
This is a contraction-domain change, not a radial recomputation or source fit.

The existing thermal bound is evaluated on the same windows and added with
directed rounding. The result is conditional on matching order24 source
corrections for LOW and middle separately. Their union is reported as an
alternative aggregate, never an additional contribution to the two pieces.
The [record](../results/development/nsc-incoming-group32-order24-bound.json)
retains the artifact, individual/combined bounds and interval-additivity checks.
Replay performs no recurrence or local coefficient generation.

```sh
python3 scripts/derive_nsc_incoming_group32_order24_bound.py --prepare
python3 scripts/derive_nsc_incoming_group32_order24_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_group32_order24_bound.py
```
