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
stress or tangential flux imbalance at its fixed wall. Its bulk gravitational
evolution remains open. The [characteristic constraint](docs/nsc-atlas-characteristic-initialization.md)
now fixes the pulse's initial radius correction when the corner geometry and
null normalization are held fixed. The radial correction and its affine slope
are negative after the compact trace; this is not a mass shift or an evolved
wall response. A [separate mass calculation](docs/nsc-atlas-source-mass.md)
now supplies the transverse constraint and geometric mass on that same
initial null leg. Its post-input quadratic coefficient is negative and
constant; this is a finite initial-data response, with later wall motion
and observer corrections still open. The [normal-force check](docs/nsc-atlas-scalar-normal-force.md)
also finds zero additional normal-stress jump at the fixed conformal wall.
It constrains the later combined metric and wall response; it does not
require zero wall displacement. The [wall-mass closure](docs/nsc-atlas-wall-mass-closure.md)
now expresses the exterior geometric mass coefficient as an integral and a
frame correction built from the existing wall trace. This is an analytic
functional; its numerical evaluation and coupled cap/wall matching remain
open. Comparing geometries still needs an independently
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
- [Preparation certificate and retained geometry](docs/nsc-atlas-field-preparation.md)
- [Paper guide and historical editions](docs/papers.md)
- [Scientific strands and owners](docs/research-index.md)
- [Research inventory and reading boundaries](docs/research-inventory.md)
- [Prior art and calculation admission](docs/prior-art-and-open-claim.md)
- [arXiv mathematical library](docs/arxiv-studies.md)
- [Commands and reproduction](docs/reproducing.md)
- [Reassembly and recovery](docs/repository-consolidation.md)

There is one active package at `src/recursive_horizons`, with root `scripts`,
`tests`, `docs` and `results`. Install with `python -m pip install -e '.[paper,dev]'`.
Run from this root:

```sh
python scripts/check_publication.py
python scripts/check_lab_snapshot.py
python scripts/check_claim_bindings.py
python scripts/build_paper.py --check
```

Historical atlas v1–v6 and supporting articles replay in explicitly selected
original contexts under ignored `.build/replay`. Changed current evidence fails
authentication; current-source mismatches never silently select historical code.
Exact historical image replay uses a separate Windows Python 3.14.8 / Matplotlib
3.11.2 profile; the current paper uses Matplotlib 3.11.2. Follow the
[historical environment setup](docs/reproducing.md#preserved-editions-and-source-contexts)
before running `scripts/replay_history.py all`, or select it through `HISTORY_PYTHON`
when running `make verify`. GitHub checks the two profiles in separate jobs.
The v7 preparation and v8 response have explicit `atlas-v7` and `atlas-v8` historical replays;
the legacy `all` command retains the v1-v6 profile.
New calculations require a named uncertainty and a finite budget.

This is a research manuscript under development. No arXiv submission or
submission-readiness claim is made. See [THEORY.md](THEORY.md) for claim language
and [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidance.
