# Centered enclosure of the existing group14 order24 defect

Reuse the explicit order24 source correction in `67b4e1b`, the rank-one
energy estimate in `eb8c8eb`, and the existing directed Riccati recurrence.
The missing connection is a sufficiently tight radial enclosure for that
same 24-term projector. The existing 128-cell natural enclosure gives a
lapse-action error upper bound of approximately `2.51408525e-10` on
`16 <= E <= 160`, above the unchanged `3e-11` comparison tolerance.

The bounded experiment retains these 128 cells and uses a fourth-degree
whole-cell Taylor remainder. It changes neither the physical state nor the
24-term approximation. Stop after this one coefficient certificate, replay,
and focused tests; if the bound remains above tolerance, report OPEN without
another refinement. A passing component is not a complete source or constraint
certificate. Source correction composition remains independently owned.

## Enclosure and derivative depth

For each complex defect coefficient, intersect the natural rectangle with

$$
f(I)\subset\sum_{j=0}^{3}\frac{f^{(j)}(c)}{j!}(I-c)^j
+\frac{f^{(4)}(I)}{4!}(I-c)^4.
$$

The recurrence stores normalized Taylor coefficients. Intermediate `c_j`
needs only degrees through `24+4-j+1`; the final defect needs degrees zero
through four. Geometry through degree 28 suffices, including the derivative
of `c24` in `f24`. Trimming only inaccessible higher coefficients preserves
every required recurrence operation and never inserts a `c25` term into
the physical projector.

The existing directed horizon weight and positive radial upper sums remain
unchanged. The saved order24 local cross-product coefficients are imported
from the authenticated energy-budget artifact. They and the new radial
bounds enter the existing `projected_energy_bound` owner. The finite-offset
initialization, thermal terms, source quadrature, other groups and physical
stationarity remain separate.

The [record](../results/development/nsc-incoming-centered-order24.json)
contains the coefficient artifact, actual density/lapse bounds, comparison
with the prior natural bound, and explicit verification residuals. Replay
authenticates and contracts saved coefficients without repeating recurrence.

```sh
python3 scripts/derive_nsc_incoming_centered_order24.py --prepare
python3 scripts/derive_nsc_incoming_centered_order24.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_centered_order24.py
```
