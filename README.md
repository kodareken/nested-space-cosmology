# Nested-Space Cosmology: The Geometry of a Universe

**Douglas Ek**

NSC develops a geometric atlas of the world seen from inside a region.
Its central hypothesis is spherical inheritance through resolution:
approximately spherical, sourced regions form a hierarchy in which a shared
law is inherited while scale, state and resolved structure change.

The book follows one mathematical argument:

**Geometry -> inherited law and state -> physical boundaries -> expansion
and structure -> interior measurements.**

The dark-sector programme asks how inherited fields and unresolved
surrounding geometry appear gravitationally inside a region. The central
worked construction uses Einstein gravity, a scalar coupled to a four-form,
and charged membranes. The same action supplies regional vacuum relations,
a constrained three-region slice, an evolving protected parent, density
growth, conserved energy supply and finite examples of structure formation.
Its evolved sources then determine clocks, lensing and freely falling
motion.
An exact inherited-source comparison controls the matter-background
normalization and propagates it to a luminosity-distance interval;
NSC-P005 passes its independent unaveraged check. The same geometry
constrains vacuum, scalar and curvature fractions independently of observer
age and SI ruler. Absolute scale, preparation and population are physical inputs to
close and test within that connected account.

## Read the atlas

The main argument now follows the offset-to-horizon-to-vacuum test in seven
pages. The full technical atlas is preserved in Appendix A. The latest
[registered flow comparison](docs/prediction-register.md) calculates the
receiving stress; it does not yet identify today's vacuum fraction.
The chapter route below describes the preserved technical atlas.

- [The Geometry of a Universe: illustrated PDF](paper/nested-space-cosmology.pdf)
- [Full manuscript and primary-source links](paper/nested-space-cosmology.md)
- [Atlas programme and scientific owners](docs/nsc-cosmological-atlas.md)
- [Current calculations and their domains](docs/current-result.md)
- [Prior mathematics and the arXiv library](docs/arxiv-studies.md)

| Part | Question carried through the book |
| --- | --- |
| I. Spherical geometry and inherited law | What selects a sphere, and what does a shared law preserve? |
| II. Regions, interfaces and causal access | Which sourced geometries connect, and which observers can communicate? |
| III. Parent vacuum and evolving matter | What does one scalar/four-form action transmit across regions? |
| IV. Concentration and its reservoir | How does inherited energy supply, form and surround local structure? |
| V. Interior observation | What do motion, light and clocks measure in that geometry? |
| VI. Gravitational, fluid and quantum exchange | How do material interfaces and unresolved modes change the response? |
| VII. Effective gravity and scale | Which reductions produce vacuum stress, and which produce other sources? |
| VIII. The connected cosmological account | Which physical closures and independent tests complete the explanation? |

The sphere's surface-energy minimum is an established theorem for positive
isotropic tension at fixed volume. The book states repeated spherical
inheritance as NSC's physical hypothesis and develops its mathematical
realizations. Published constructions retain attribution beside the
equations used to connect them. Alternative actions and preparations occupy
identified branches within this single atlas.

## Calculated connections

The [same-parent astronomical comparison](docs/nsc-atlas-parent-supernova.md)
carries the inherited source into a controlled earlier Cauchy domain and
measured distance channels. The original SN and six AP comparisons are
preserved. The stronger [full shared-ruler scoreboard](docs/nsc-atlas-shared-ruler-scoreboard.md)
propagates the complete SN shape posterior through all 13 DESI entries:
joint chi-square is 1203.22 versus 1191.18 for a flat control with equal
fitted dimension. Its conditional Planck geometry diagnostic keeps the
missing acoustic/thermal sector explicit. A microscopic SI coupling length
is fixed by the supplied inherited ratios, but the hierarchy to the parent
curvature radius and the SI cosmological constant remain unselected.

The manuscript contains operator inheritance proofs, sourced geometric
attachments, causal maps, scalar formation and conservation calculations,
and joint observer relations. The motion-calibrated optical prediction
NSC-P004 passes its frozen synthetic test. The separate quiet-cap
supernova-distance prediction disagrees with the selected Pantheon+
contrast. Both results constrain the next connected cosmological calculation.

The [prediction register](docs/prediction-register.md) separates model
forecasts, calibration and comparison. The [research index](docs/research-index.md)
and [code map](docs/active-code-map.md) locate the derivations and implementation.
The detailed characteristic-field development and computational evidence
remain in the book's appendices and scientific owners.

## Reproduce the work

There is one active package, `src/recursive_horizons`, with root scripts,
tests and evidence. See [reproduction](docs/reproducing.md) for the isolated
Windows runtime, dependencies, numerical producers and historical replay.

```text
python -m pip install -e '.[paper,dev]'
python scripts/check_publication.py
python scripts/check_lab_snapshot.py
python scripts/check_claim_bindings.py
python scripts/check_prediction_register.py
python -m pytest tests/
python scripts/build_paper.py --check
```

Compact summaries and source hashes are retained in Git. Raw numerical
histories are regenerated when needed. [Preservation and recovery](docs/repository-consolidation.md)
describe the retained history and archive. The [paper guide](docs/papers.md)
records editions of the one living manuscript; current author-review and
submission metadata are in [the catalog](paper/catalog.json).
