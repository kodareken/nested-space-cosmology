# Focused nested-qualities companion

[Read the PDF](../paper/local-incoming-gate-draft.pdf) ·
[Download TeX source](../paper/local-incoming-gate-draft-source.tar.gz) ·
[Inspect build hashes](../paper/local-gate-draft-manifest.json) ·
[Inspect theorem snapshot](../paper/nested-quality-evidence/snapshot.json) ·
[Inspect application snapshot](../paper/local-gate-evidence/snapshot.json)

The current article gives a connected proof of four nested qualities in a
specified finite coupled-operator family. It permits different regional states
under a shared law. The earlier source-fixed incoming calculation remains an
OPEN application checkpoint in an appendix. Its source records and numerical
claims retain their original meanings. The 43-page foundation is unchanged.

The prior nine-page companion remains in Git history at `134962d`. The revised
PDF is the same paper location with a new title and result focus, not a third
manuscript. See [the paper guide](papers.md) and [proof note](nsc-nested-qualities.md).

## What is verifiable here

The selected records and two mathematical owners in `paper/local-gate-evidence/`
were read from immutable laboratory Git blobs. Their original bytes are
preserved, with SHA-256, Git blob IDs, and the laboratory commit in the snapshot.
The snapshot is deliberately labelled incomplete: it is enough to inspect the
quoted checkpoint, but not to rerun all active laboratory calculations. The
existing historical result manifests retain their complete public dependencies.

The new artifact manifest binds the draft sources, presentation generator,
snapshot, PDF, generated bibliography, and deterministic source archive.
Budget-table numbers and the candidate pair are generated from the original
records. Hash checks establish identity; they do not constitute a scientific
certificate.

From the repository root, using Python 3.12 or later:

```console
python scripts/verify_local_gate_draft.py
python -m unittest discover -s tests -p test_local_gate_draft.py -v
```

These checks need only the Python standard library. Ordinary CI also checks
the historical manifests, publication tests, and preserved notebook. Expensive
scientific regeneration runs only through the explicit `make reproduce` command
or the manual GitHub workflow option.

The additional `paper/nested-quality-evidence/` snapshot contains immutable
Git blobs for the proof and its exact controls. The verifier checks every
listed byte/hash and the proof-owner links in the record. Original Markdown
source inside that snapshot retains its laboratory-relative citation paths;
the reader-facing proof under `docs/` has checked public links. Reused earlier
records are references, not a new complete gravitational certificate.

## Rebuild the PDF

Use the exact compiler identified by `paper/local-gate-draft-manifest.json`.
The earlier CI builder used Ubuntu 24.04 TeX packages, `pdflatex`, and BibTeX,
with standard `amsmath`, `amssymb`, `geometry`,
`booktabs`, `microtype`, and `hyperref` packages. It is not yet a pinned arXiv
container. A different distribution can legitimately produce different bytes.

```console
python scripts/build_local_gate_draft.py --engine tectonic --check
```

The current Mac build uses Tectonic 0.17.0 (official release binary, separately downloaded into ignored build storage). The builder also retains its pdflatex/BibTeX mode. It uses two clean temporary directories, a fixed source-date epoch,
suppressed volatile PDF metadata, and disabled shell escape. The pdflatex
mode runs BibTeX and three TeX passes; Tectonic manages its own bibliography
and convergence passes in untrusted mode. It rejects unresolved references, overfull boxes, and differing output
bytes. The source archive contains `main.tex` at its root, all input macros,
the `.bib`, and the generated `main.bbl`; it excludes logs, hidden files, and
absolute paths. Both clean PDF builds must match the tracked artifact.

For an intentional editorial update, regenerate `evidence_values.tex` using
`evidence_values` in `scripts/local_gate_draft.py`, run the builder without
`--check`, render every page for review, then refresh public provenance with
`python scripts/build_publication_provenance.py`. Review the changed manifest
and Git diff before publication. Historical scientific records must not change
to accommodate an editorial update.

## Path to submission

An arXiv submission requires a reviewed claim matching the actual theorem,
verified source/PDF, references, disclosure, account/category eligibility and
Douglas's final metadata/license/preview decision. The old physical local-gate
certificate is required only if that result is claimed; it is not a prerequisite
for this finite theorem. This update makes no submission or release claim.
