# Historical archive — non-canonical

Everything below this directory is retained for provenance, not as the current scientific position of Recursive Horizons. Files may contain incorrect equations, superseded values, invalidated analyses, unsupported physical claims, combative rhetoric, raw assistant citation markers, broken links, or direct contradictions.

Do not cite an archived claim without checking the active [claim ledger](../docs/claim-ledger.md) and the current [paper](../paper/recursive-horizons.md).

## Layout

### `legacy-corpus-2026-08-16/current-documents/`

The complete fifteen-document public-facing corpus immediately before version 0.2.0, including its former `README.md`. It is preserved byte-for-byte through Git history. Notable retired framing includes “proof,” “finishing blow,” “zero free parameters,” an entropy equality attributed to unitarity, a variable-`c` interpretation of Schwarzschild coordinates, and discovery-level CMB anomaly language.

### `legacy-corpus-2026-08-16/exploratory-tests/`

The complete former thirteen-script suite and README. These are frozen exploratory prototypes. They are not imported by the active package, are not run by CI, and must not be described as a verification suite. See [the review summary](../docs/review-summary.md) for failure details.

### `pre-git-drafts/markdown/`

Earlier foundations, reviews, books, geometry notes, unresolved-problem lists, bounce/spin/SUSY drafts, and the “six proofs” synthesis. Some of the most epistemically careful language is in `review-notes.md`, `analysis-proof-black-hole.md`, `deep-research-report.md`, and `part2-synthesis.md`; later drafts often overstated what those reviews had already identified as assumptions.

### `pre-git-drafts/docx/`

Two historical Word playbooks from July 2026. They are distinct snapshots, not duplicates. The earlier file contains the longer general playbook; the `(1)` file is a shorter “bringing it home” continuation. Their original references rely partly on secondary sources and are superseded by the primary-source bibliography in the current paper.

### `pre-git-drafts/scripts/`

Retired raw FITS/HEALPix and Kantowski–Sachs utilities. They contain hard-coded local paths and unvalidated numerical methods.

### `pre-git-drafts/presentations/`

Two historical HTML visualizations. `presentation.html` is notable for correctly labeling the parent→child map and entropy equality as conjectures and explicitly listing the missing transition, information, dark-matter, and observational work. Its visual cycle remains an analogy, not a derivation.

### `historical-tests/`

Transient or live-fixture tests that no longer describe the current store tip
are preserved here with source-hash notes. They are not active tests and are
excluded from ordinary search and test discovery. See
[`historical-tests/README.md`](historical-tests/README.md).

## Archive policy

- Add corrections to canonical documents and the claim ledger.
- Keep historical snapshots immutable unless required to remove secrets, personal data, malware, or material that cannot lawfully be distributed.
- When a retired idea is revived, create a new active derivation with new tests; do not rewrite the archive to erase the previous failure.
