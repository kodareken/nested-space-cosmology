# Nested-Space Cosmology: The Geometry of a Universe

Start with regions, boundaries and fields; assemble possible nested geometries;
specify inheritance and transmission; calculate what an interior observer measures.

The living research paper is [the manuscript](paper/nested-space-cosmology.md)
and [its illustrated PDF](paper/nested-space-cosmology.pdf).
The [overview map](paper/figures/geometry-overview.png) connects geometry,
observer access and the calculated early scalar response. Manuscript revisions
are independent of the software distribution, which remains version 0.27.0.

Spherical symmetry supplies a definite areal radius and tractable matching laws.
It does not establish that every particle is a spherical interior or that our
universe has a parent. The map distinguishes descent, geometric attachment,
causal influence and correlation. Missing physical connections remain visible.

The current geometric route reaches a sourced interior connection, a causal map,
a calibration waveform and a [computed response to independent exterior data](docs/nsc-atlas-prepared-response.md).
The eighth atlas resolves a central scalar waveform and a positive tail within
its finite arrival window, supported by an analytic bound and numerical refinement.
The [analytic inverse](docs/nsc-atlas-observer-inference.md) now determines
the source from a complete compatible signed field history in that fixed
geometry. It also excludes histories that violate a positive compact
preparation. A [quadratic stress check](docs/nsc-atlas-scalar-wall-consistency.md)
shows that the same conformal probe needs no additional singular surface
stress or tangential flux imbalance at its fixed wall. The
[characteristic constraint](docs/nsc-atlas-characteristic-initialization.md)
now fixes the pulse's initial radius correction when the corner geometry and
null normalization are held fixed. The radial correction and its affine slope
are negative after the compact trace; this is not a mass shift or an evolved
wall response. A [separate mass calculation](docs/nsc-atlas-source-mass.md)
now supplies the transverse constraint and geometric mass on that same
initial null leg. Its post-input quadratic coefficient is negative and
constant; this is a finite initial-data response. The
[normal-force check](docs/nsc-atlas-scalar-normal-force.md)
also finds zero additional normal-stress jump at the fixed conformal wall.
It constrains the later combined metric and wall response; it does not
require zero wall displacement. The [wall-mass closure](docs/nsc-atlas-wall-mass-closure.md)
now expresses the exterior geometric mass coefficient as an integral and a
frame correction built from the existing wall trace. This is an analytic
functional. Its [two-sided successor](docs/nsc-atlas-cap-wall-response.md)
now derives the cap mass and a regular wall-displacement equation, with
compatible angular, temporal and normal junction conditions. The common
model specification records the action, state, parameters and clocks.
The [smooth exterior calculation](docs/nsc-atlas-bulk-response.md) propagates
radius and lapse on the initial characteristic rectangle. The
[moving-wall calculation](docs/nsc-atlas-wall-evolution.md) now computes
the formal wall displacement and exterior clock correction from the same
field on the retained wedge. The [cap calculation](docs/nsc-atlas-cap-response.md)
now reconstructs the interior metric with its regular centre, checks its
mass against independent wall transport, and computes the central null-arrival
correction. The [next scalar calculation](docs/nsc-atlas-scalar-metric-response.md)
now evaluates the cap's outgoing backscatter and its displaced-wall boundary
input, with separate field and quadrature refinements. A clock-resolved
boundary and cached forcing now supply the [propagated exterior scalar
correction](docs/nsc-atlas-scalar-evolution.md), with separate cache, grid
and ODE controls. Its central reconstruction and finite-amplitude control
remain open.
Comparing geometries still needs an independently
justified source restriction or a source-linked additional channel. Finite transfer,
feedback and memory results supply reusable mathematics in their own domains.
Periodic collars are not established child universes. Constants, antimatter,
expansion and recycling remain research questions.

## Read and reproduce

- [Atlas programme and calculation dependencies](docs/nsc-cosmological-atlas.md)
- [Current prepared response and evidence](docs/nsc-atlas-prepared-response.md)
- [Observer inverse and preparation restrictions](docs/nsc-atlas-observer-inference.md)
- [Scalar stress and shell compatibility](docs/nsc-atlas-scalar-wall-consistency.md)
- [Normal-force compatibility and the remaining wall condition](docs/nsc-atlas-scalar-normal-force.md)
- [Gravitational initial data for the prepared pulse](docs/nsc-atlas-characteristic-initialization.md)
- [Geometric mass response on the prepared null surface](docs/nsc-atlas-source-mass.md)
- [Analytic exterior mass from the known wall field](docs/nsc-atlas-wall-mass-closure.md)
- [Cap mass and shared wall displacement](docs/nsc-atlas-cap-wall-response.md)
- [Smooth exterior metric on the initial rectangle](docs/nsc-atlas-bulk-response.md)
- [Retained-wedge field, moving wall and exterior clock](docs/nsc-atlas-wall-evolution.md)
- [Smooth cap and corrected central null clock](docs/nsc-atlas-cap-response.md)
- [Scalar backscatter and moving-wall boundary input](docs/nsc-atlas-scalar-metric-response.md)
- [Preparation certificate and retained geometry](docs/nsc-atlas-field-preparation.md)
- [Paper guide and historical editions](docs/papers.md)
- [Scientific strands and owners](docs/research-index.md)
- [Research inventory and reading boundaries](docs/research-inventory.md)
- [Prior art and calculation admission](docs/prior-art-and-open-claim.md)
- [arXiv mathematical library](docs/arxiv-studies.md)
- [Commands and reproduction](docs/reproducing.md)
- [Preservation and source pins](docs/repository-consolidation.md)

There is one active package at `src/recursive_horizons`, with root `scripts`,
`tests`, `docs` and `results`. Install with `python -m pip install -e '.[paper,dev]'`.
Run from this root:

```sh
python scripts/check_publication.py
python scripts/check_lab_snapshot.py
python scripts/check_claim_bindings.py
python -m pytest tests/
python scripts/build_paper.py --check
```

`make verify` runs current checks, the full retained test suite and the paper
check. CI runs that suite on Linux and Windows with Matplotlib 3.11.2.
The v2 [retained-file receipts](docs/lab-snapshot.json) authenticate current
bytes; unchanged [source pins](archive/source-history/manifest.json) preserve
original scientific contexts. Frozen records retain their exact hashes and
declared source contexts. Deleted carriers are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

Historical producers replay only in explicitly selected contexts under
ignored `.build/replay`; they do not replace current-source checks. The
[reproduction guide](docs/reproducing.md#preserved-editions-and-source-contexts)
distinguishes authentication from scientific recomputation and records each
original renderer requirement. The current renderer does not promise
byte-identical historical PNGs. Historical replay is separate from `make verify`.
New calculations require a named uncertainty and a finite budget.

This is a research manuscript under development. No arXiv submission or
submission-readiness claim is made. See [THEORY.md](THEORY.md) for claim language
and [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidance.
