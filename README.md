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
a calibration waveform and a [finite exterior preparation operator](docs/nsc-atlas-field-preparation.md).
Independent data beside the early wall define a unique transmission problem;
its numerical received waveform, remote preparation and backreaction remain open. Finite transfer,
feedback and memory results supply reusable mathematics in their own domains.
Periodic collars are not established child universes. Constants, antimatter,
expansion and recycling remain research questions.

## Read and reproduce

- [Atlas programme and calculation dependencies](docs/nsc-cosmological-atlas.md)
- [Current preparation map and evidence](docs/nsc-cosmological-atlas-v7.md)
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
3.11.2 profile; the current paper uses Matplotlib 3.11.1. Follow the
[historical environment setup](docs/reproducing.md#preserved-editions-and-source-contexts)
before running `scripts/replay_history.py all`, or select it through `HISTORY_PYTHON`
when running `make verify`. GitHub checks the two profiles in separate jobs.
New calculations require a named uncertainty and a finite budget.

This is a research manuscript under development. No arXiv submission or
submission-readiness claim is made. See [THEORY.md](THEORY.md) for claim language
and [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidance.
