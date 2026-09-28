# Conditional order24 energy bound for retained group12

Reuse the certified trimmed Taylor recurrence and rank-one energy bound,
with the fixed 128-cell collar partition and fourth-degree remainder. The
group12 order16 certificate in `44148c9` remains too broad, with lapse error
upper bound `1.5027532081e-10`. The missing connection is the corresponding
rigorous bound for an explicit 24-term approximation on the same authenticated
middle interval `[40,320]`.

Calculate group12 only, with numerical order24, 128 cells and Taylor depth4.
Compare its energy bound with the unchanged `3e-11` tolerance and stop after
this one certificate, replay and focused checks. No second refinement or
other-group calculation belongs to this record. A passing result is
conditional on separately applying the explicit order24-minus-order16 source
correction; it never certifies the archived order16 values.

## Reused enclosure and source conventions

The owner imports the frozen `trimmed_defect_jets`, the existing complex
Taylor enclosure/intersection, and the original horizon cell map and weight.
Geometry through degree28 supplies the normalized fourth derivatives of
`f24,...,f48`. Directed parameter intervals enclose the inherited compact mass
and magnetic angular labels. The actual angular signs and multiplicities
are unchanged. No mode, scattering equation, state or metric history is solved.

The local cross-product coefficients of the 24-term projector are computed
once by the existing `local_cross_product_coefficients` owner and frozen with
the radial bounds. The existing `projected_energy_bound` contracts both on
the metadata-authenticated interval. The separately inherited thermal bound
is added in lapse-action units. No source correction, quadrature, low/subgap,
other-group or angular/compact completeness claim is included.

The [record](../results/development/nsc-incoming-retained-order24-bound.json)
contains the actual error bound, conditional status, coefficient artifact,
and replay/intersection residuals. `--check` authenticates and contracts that
artifact without regenerating local or radial coefficients.

```sh
python3 scripts/derive_nsc_incoming_retained_order24_bound.py --prepare
python3 scripts/derive_nsc_incoming_retained_order24_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_retained_order24_bound.py
```
