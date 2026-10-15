# Repository reassembly and audited cleanup

The sole active root is `nested-space-cosmology`. Code, scripts, tests, notes
and retained evidence live directly in `src`, `scripts`, `tests`, `docs` and
`results`. Distribution identity remains `nested-space-cosmology` version 0.27.0.

The 12 October 2026 migration and retirement records, superseded variants and
deleted carriers are historical evidence in
[Zenodo DOI 10.5281/zenodo.23090015](https://doi.org/10.5281/zenodo.23090015).
The retired containers are no longer runtime inputs to current-root verification.
Some original bytes still authenticate retained evidence through the sparse
object archive described below. The
[current retained-file receipt](lab-snapshot.json) is the v2 successor: each
checked row names its current repository-relative path, exact byte count
and SHA256, while retaining the original import path, hashes and Git-blob
provenance where present. Audited trash and merged sources are excluded
from its live inventory. Explicit later edits require an updated current
receipt; original provenance identities are preserved.

## Preservation and reading boundaries

The original scientific baseline is commit
`04a0cb9378977772243ebf7c9be3c74182b030f2`. The complete
[foundation edition](../paper/archive/foundation-v0.26.0/archive-manifest.json)
remains frozen. The second worktree's 400-file snapshot and 294 unique
research candidates were preserved during reassembly; their retired
carriers and variant context are cited through the DOI above. They do not
become newly verified results because their hashes were preserved.

Current scientific owners and unread boundaries remain in the
[inventory](research-inventory.md) and [strand index](research-index.md).
Hash authentication establishes byte integrity, not mathematical correctness.
Unrelated agent configuration remains outside this cleanup.

## Explicit historical replay

`python scripts/replay_history.py all` deliberately materializes historical
source contexts under ignored `.build/replay` for the preserved editions.
The exact source proofs remain in `archive/source-history`. Its original
25 Merkle objects and six source pins remain unchanged. The historical Git
object store is process-local, never a competing installed package.

The approved 15 October 2026 history compaction adds 1,654 authentication
objects (93,389,210 compressed bytes) for 2,866 explicitly inventoried paths
across 17 original contexts. The [object manifest](../archive/source-history/manifest.json)
records their identities. Byte-identical objects already retained in the
baseline tree are reused instead of copied. Original receipt arrays are
authenticated as bytes; they do not replace current numerical inputs.

Only the necessary commits, Merkle trees and blobs are transported, without
their parent ancestry. These objects are tracked as archival files and exposed
to readers through a process-local `GIT_ALTERNATE_OBJECT_DIRECTORIES`; they
are not live history refs or permanent repository alternates. The keep list
for history filtering includes the approved archive additions and runtime
setup. Current scientific records and paper files retain their baseline bytes.
Full pre-cleanup checkouts and campaigns require the Zenodo archive; the
sparse authentication store does not preserve an entire historical checkout.

Current checks use retained root evidence and reject changed undeclared
bytes. Historical records keep their original producer identity; a current
source mismatch does not silently select an older implementation. Missing
Git objects, corrupt retained artifacts and changed pinned evidence fail.

Explicit `legacy-replay` source lookup can consult the manifest's bounded
list of original commit/path/blob identities after searching the retained
Git log. It requires both the recorded SHA256 and the declared Git blob
identity. Production reads, explicit commit pins, numerical inputs and
custom lookup callbacks retain their original contracts.

## Verification and recovery

Run `python scripts/check_lab_snapshot.py` after transfer or an authorized
receipt update. It checks the v2 retained-file hashes directly and then
authenticates the declared `archive/source-history` inventory and unchanged
original source pins.
It does not require the retired reassembly manifest or retired containers.
After a later paper or code edit, refresh only the declared current byte/hash
fields before checking; retain the original import metadata and source pins.

Deleted historical content is addressed by the DOI above or its original
pinned Git context. Restore individual required evidence into this one root
when explicitly needed; do not recreate a permanent nested checkout.
The obsolete containers were previously retired after resident-file and
recovery checks. Their retirement receipt is historical DOI evidence, while
the active root remains the only registered worktree for ordinary work.
