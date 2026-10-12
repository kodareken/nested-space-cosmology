# Repository reassembly, 12 October 2026

The sole active root is `nested-space-cosmology`. Code, scripts, tests, notes and
evidence live directly in `src`, `scripts`, `tests`, `docs` and `results`.
Distribution identity remains `nested-space-cosmology` version 0.27.0.

The [migration manifest](../results/reassembly/manifest.json) records original
paths, SHA-256 hashes, byte counts, dispositions and destinations. It covers
laboratory material, the second worktree snapshot, source proofs, previews and
helpers. Identical overlaps have one owner; divergent implementations remain
under `archive/reassembly-20261012`. Deliberate current edits have separate
successor records; original hashes are never rewritten. Missing or changed
undeclared destinations fail the checker.

## Preservation and reading boundaries

The original scientific baseline is commit
`04a0cb9378977772243ebf7c9be3c74182b030f2`. The complete
[foundation edition](../paper/archive/foundation-v0.26.0/archive-manifest.json)
was archived before replacement. The second worktree's 400-file snapshot remains
under `archive/reassembly-20261012/w2`; its 294 unique research files enter the
root as exploratory candidates, not newly verified results. The differing RN
implementation remains a named variant. Temporary fixtures are historical context.

The reset backup `.research-reset-preservation` is retained. Unrelated agent
configuration is outside consolidation. Previews and helpers are preserved under
`archive/reassembly-20261012/previews-and-helpers`; obsolete containers enter the
Recycle Bin only after recovery checks. The [inventory](research-inventory.md),
[strand index](research-index.md) and [migration dispositions](../results/reassembly/README.md)
retain owners and reading boundaries. Hash authentication proves preservation,
not mathematical correctness. Unread derivations remain explicitly unverified.

## Explicit historical replay

`python scripts/replay_history.py all` materializes original sources under ignored
`.build/replay/04a0cb937897`, authenticates retained evidence against the baseline,
then runs all six atlas checks and both supporting-paper checks. Original source
proofs remain in `archive/source-history`; the historical Git object store is
process-local, never a permanent alternate.

A changed context fails rather than repairing itself. Current evidence is checked
before materialization. Current producers use root paths and reject old bindings;
they never fall back to old source. Choose a historical context deliberately.
Missing Git objects, corrupt archives and changed preserved evidence are errors.

## Recovery

Run `python scripts/reassemble_repository.py --check` before retirement or after
transfer. Each row supplies an original path and retained destination. Original
tracked sources are also recoverable from the pinned baseline; the w2 snapshot
covers its uncommitted work. Do not restore a permanent competing nested project.
The [original lab import receipt](lab-snapshot.json) remains historical provenance.

## Completed retirement

The two obsolete containers and six empty preview folders were moved to the
Windows Recycle Bin after a complete resident-file audit. The
[retirement receipt](../results/reassembly/retirement.json) accounts for 16,422
files through pinned Git blobs, preserved copies and authenticated ZIP members.
Its verifier supports `--check-recovery` after the original containers are gone.
The stale secondary-worktree metadata was archived before Git pruned its obsolete
registration; its branch and pinned commit remain recoverable. The root is now
the only registered worktree for ordinary work.
