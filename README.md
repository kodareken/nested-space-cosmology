# Nested-Space Cosmology: The Geometry of a Universe

**Douglas Ek**

Build a geometric account of collapse inside an evolving enclosing region:
parent deformation drives the child's proposed offset, drawing from one
spatial seed into circle, sphere and nested continuation. Each parent is
itself enclosed. Derive the relationships before assigning familiar physical
meanings to them.

The active work is the geometric foundation. The previous action-based
atlas remains a library of conditional constructions, not the foundation's
assumptions. No new manuscript edition is being produced during this reset.

## Active desk

- [Foundation](THEORY.md): premises, interpretation and open connections.
- [Offset and winding construction](docs/nsc-atlas-offset-generator.md).
- [Current result](docs/current-result.md).
- [Code and verification](docs/active-code-map.md).
- [Working method](docs/instructions.md).

Run the focused checks and the complete non-local-to-local atlas readout with the isolated project runtime:

```text
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_offset_generator.py
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --atlas-only
```

## Retained work

The [research library](docs/research-index-history-20261102.md#retained-library-accessed-deliberately)
links the old results, implementation map, literature and manuscript.
Nothing there becomes a foundation premise merely because it was previously
calculated. Useful mathematics keeps its assumptions and attribution.

Default `rg` searches show the active foundation files. Use
`rg --no-ignore 'pattern' <specific-path>` to inspect retained work.
The exclusions change search discovery only; source, tests and evidence
remain in place and executable.
