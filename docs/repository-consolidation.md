# Repository consolidation — 28 September 2026

## Decision and resulting layout

Douglas requested one local main matching GitHub, one `codex/work-branch`, and
removal of the large old search history after confirming his backup. The existing
GitHub publication history is retained. The current laboratory is imported as
`lab/` without merging its 32 GB historical object database into GitHub.
This is an organizational change, not a new scientific result.

The canonical Mac directory is `Documents/BlackHoles-Infinity`. The old
`Documents/nested-space-cosmology` path becomes a compatibility symlink to it,
so both established paper paths continue to resolve to the same files. The
former lab Git repository and manual worktrees are retired together in recoverable
Trash, with a local receipt under the user's local state directory. They are not
active sources of instructions. Existing work outside Documents is not removed.

The old 38 laboratory branches are preserved in the retired repository; GitHub's
publication ancestry remains in the new main. The active repository has only
`main` and `codex/work-branch`. An archive is not a third working repository.

## Preserved data and evidence

[Lab import receipt](lab-snapshot.json) records every imported file, its original
Git blob and SHA-256, the original lab/public SHAs, and the 1,860 intentionally
omitted NPZ paths. About 3 GiB of current source, retained trajectories,
calibrations, records and code remain. Old search arrays (about 31 GiB expanded)
are not part of the new Git ancestry. Four untracked experiment entries and ignored
local run material are preserved in the retired original checkout, not promoted
as reviewed public results. No numerical search is launched by this migration.

The three component evidence closures require four historical source files.
`lab/.source-history/` retains their exact commit/tree/blob Merkle paths (18 Git
objects), checked both by SHA-1 object identity and SHA-256. The lab launcher
provides them as process-local alternate objects for unchanged legacy `git show`
checks. This is a limited source proof cache, not complete old Git ancestry.
`git fsck` on the live repository does not depend on that partial history.

Other historical replays may require additional old sources and omitted search
payloads. The old restore script under `lab/scripts/` refers to the OLD lab Git
history: use it only with the archived repository or the user's backup. Do not
expect old search objects to exist in the new compact repository. Do not fetch
all old history merely to satisfy a historical test unrelated to active work.

The public paper manifests and scientific lab originals retain their bytes.
Only navigation/routing files and consolidation tools change. The curated public
provenance excludes the lab subtree, which has its own complete import receipt
and credential scan. A published OPEN record remains OPEN.

## Verification and maintenance

Run `python scripts/check_lab_snapshot.py` to check imported bytes, source proof
objects and credential patterns. Run lab replays through `scripts/lab.py`; run
publication checks at the root. Deliberate new lab changes need explicit receipt
updates retaining original hashes. Do not add large transient search arrays to
Git by default. Use external content-addressed payload storage with a recorded
recovery location if future campaigns need it.

Both PDF manifests authenticate the same PDFs as before consolidation. Current
workflow changes are ordinary commits, not a v0.27.0 release. Research follows
the latest user instruction. An accepted conversation plan applies where it
stays compatible with that instruction and does not override later steering.
Douglas later removed the handover and retired repository plan copies as
authorities; see
[instructions](instructions.md#scientific-framing-and-authority).

## Completed local checks

The import includes 11,544 original resident lab files; 11,541 are byte-identical
to their original Git blobs and three navigation files carry explicit new routing.
The receipt also includes 19 source-cache files (manifest plus 18 exact objects).
No scientific record or payload was rewritten.

- Full Git object verification passed on the new repository, without the old
  object database or permanent Git alternates.
- 356 focused publication/consolidation tests passed (one pre-existing skip),
  including isolated source-object replay and missing/corrupt-object controls.
- 12 focused radius, reference-defect and endpoint-contraction tests passed.
- Component-budget, whole-cone binding, current UV-v3 and endpoint-pilot replay
  passed in the consolidated lab; the physical gate remains OPEN.
- All 60 source families loaded: 6,092 signed batches, 96,336 signed energy rows,
  64 history coordinates. No evolution or nonlinear search was launched.
- The foundation PDF rebuilt to its existing hash; the companion's source,
  selected evidence and PDF manifest passed authentication.
- After the directory swap, lab integrity and component-budget replay passed again.

The two unmerged historical FGC worktree commits have separate patches and all
four former worktrees have compressed snapshots in the local consolidation
state directory. The 24 untracked files (four experiment entries) also have a
separate checked copy there, outside Trash. The swap receipt records the exact
retired paths and successful repair of their worktree connections. Retiring the
old Git database therefore does not discard those uncommitted or unmerged files.

A bounded Grok review attempt timed out with malformed output; it is not counted
as independent approval. The checks above were executed directly by the integrator.

## Later explicit receipt successors

`scripts/update_lab_snapshot.py` registers deliberate lab successors and keeps
each file's original import hashes. The completed family `(14, 1)` source
insertion, including the direct, subgap and threshold witness trees, belongs
in that receipt. `physical_upstream_budget_component` stays null. The cubic
UV inventory and leading residual do not certify the ultraviolet tail. The
local incoming gate remains OPEN and its campaign is paused
([instructions](instructions.md#scientific-framing-and-authority)). Public provenance remains
[PUBLICATION-PROVENANCE.json](../PUBLICATION-PROVENANCE.json), rebuilt by
`scripts/build_publication_provenance.py`. That file does not copy the lab
subtree.
