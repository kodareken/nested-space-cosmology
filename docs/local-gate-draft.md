# Focused local incoming-gate draft

[Read the PDF](../paper/local-incoming-gate-draft.pdf) Â·
[Download TeX source](../paper/local-incoming-gate-draft-source.tar.gz) Â·
[Inspect build hashes](../paper/local-gate-draft-manifest.json) Â·
[Inspect source snapshot](../paper/local-gate-evidence/snapshot.json)

The draft develops the question in English and then gives the corresponding
Schur reduction, retarded state law, lapse/shift constraints, and local error
criterion. A notation table, nine-component budget table, and scoped primary
references make the mathematical argument inspectable. Douglas Ek is the author;
significant Codex and Grok assistance is disclosed in the article.

The scientific verdict is **OPEN**. The candidate is not a physical root, the
one-cell field pilot is not a whole-cone certificate, and three known v4 error
terms do not complete the budget. The successor 257-node budget needs newly
authenticated applicability for all nine terms. Release `0.27.0` is reserved for
a closed local result and its full dependency graph.

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

## Rebuild the PDF

Use the exact compiler identified by `paper/local-gate-draft-manifest.json`.
The CI builder uses the Ubuntu 24.04 TeX packages, `pdflatex`, and BibTeX,
with standard `amsmath`, `amssymb`, `geometry`,
`booktabs`, `microtype`, and `hyperref` packages. It is not yet a pinned arXiv
container. A different distribution can legitimately produce different bytes.

```console
python scripts/build_local_gate_draft.py --check
```

The builder uses two clean temporary directories, a fixed source-date epoch,
suppressed volatile PDF metadata, disabled shell escape, BibTeX, and three TeX
passes. It rejects unresolved references, overfull boxes, and differing output
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

An arXiv submission still requires a closed local certificate, the complete
scientific dependency closure, independent verifier replay, a pinned compatible
TeX environment, account/category eligibility, and Douglas's final metadata,
license, preview, and submission actions. This draft performs none of them.
The focused draft and notebook have different purposes and version histories.
