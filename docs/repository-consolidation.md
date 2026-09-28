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
workflow changes are ordinary commits, not a v0.27.0 release. Research continues
from the [handover](handover.md), and the accepted [plan](PLAN.md) is unchanged.
