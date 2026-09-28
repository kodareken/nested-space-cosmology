# Explicit group12 numerical-order source correction

The unchanged incoming source admits the explicit finite-band correction

$$
\Delta K_A=\operatorname{tr}[(P_{24}-P_{16})V_A],
\qquad 40\le E\le320.
$$

The same fourth-order reference cancels in this difference. The existing
horizon source law, geometric parameters and signed multiplicity remain
fixed. This is a new numerical approximation to the same source, not a
selected state or stress term. The archived order16 values and thermal
insertion remain intact; the small physical thermal/scattering uncertainty
retains its separately owned bound.

## Decision and stopping condition

Reuse the arbitrary-precision Riccati coefficient owner, stable projector
difference from the group14 correction, and authenticated group12 middle
panels. Its centered order16 energy enclosure remains `1.50e-10`, above
`3e-11`. The exact missing quantities are an explicit higher-order source
increment and its independent physical error bound. Record this increment
on the existing interval using old quadrature and a positive48-node refinement,
at50 and70 digits. Stop after the correction and focused checks; its physical
order24 defect certificate is a separate owner. An order difference alone
does not certify the improved approximation.

Risks are asymptotic-series growth, convention/multiplicity mistakes and
conflating numerical integration indicators with physical mode error.
The first16 coefficients and an independent projector/vertex contraction
check the implementation. The old and corrected approximants remain
separate fields. Replay authenticates and contracts stored arrays without
regenerating coefficients, quadrature or physical modes.

The [record](../results/development/nsc-incoming-retained-order-correction.json)
contains all four stress increments and their action-gradient projection.
The [order24 remainder owner](nsc-incoming-retained-order24-bound.md)
must be composed with this correction before using the improved source.
No previous source/constraint record, scalar coefficient, seed, physical
metric history, action term or PDF is changed.

```sh
python3 scripts/derive_nsc_incoming_retained_order_correction.py --prepare
python3 scripts/derive_nsc_incoming_retained_order_correction.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_retained_order_correction.py
```
