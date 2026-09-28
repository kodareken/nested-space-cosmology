# Explicit group32 LOW and middle source corrections

Reuse the frozen source-window helper and the single group32 order24 radial
certificate in `490f202`. The missing connection is matching source values
for its actual LOW `[32,40]` and middle `[40,320]` windows. Their conditional
combined lapse bound is approximately `2.1336305733e-13`, below `3e-11`.

Prepare exactly these two source windows at numerical order24, using the
existing positive48/64 quadrature and precision50/70 controls. LOW receives
an explicit vacuum-source replacement relative to its actual selected
archived values. Middle receives only the stable `P24-P16` correction;
its archived thermal insertion remains present. Preserve all archived arrays.
Stop after the two source preparations, coefficient replay, focused checks
and local commit. No new radial recurrence, mode/scattering solve, geometry,
scale, state law, fitting parameter or action term belongs to this work.

## Two disjoint before/after changes

For LOW the additive correction is `new_order24_vacuum - archived_LOW`.
Physical thermal effects remain nonzero and bounded separately. For middle
the additive correction is the same-vertex contraction of `P24-P16`;
reference subtraction cancels in that difference, and archived thermal
values are retained. The existing source normalization counts copies,
angular signs and negative-frequency folding once.

Both source windows bind to the same immutable radial coefficient artifact
and its local order24 cross-product coefficients. The `[32,320]` union
certificate verifies the sum of the two physical component bounds. It is
not a third source contribution. Independent matrix/Pauli contractions
verify the LOW ad4 subtraction and the middle projector difference.

The [record](../results/development/nsc-incoming-group32-source.json)
contains each archived source, new value, explicit delta, updated approximant,
action-gradient correction, thermal policy, and numerical controls. Precision
and quadrature differences are indicators, not rigorous quadrature bounds.
The bound certifies the vacuum/thermal energy component of the order24
approximant; complete source accuracy and all constraint equations remain OPEN.

```sh
python3 scripts/derive_nsc_incoming_group32_source.py --prepare
python3 scripts/derive_nsc_incoming_group32_source.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_group32_source.py
```
