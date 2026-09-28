# NSC storage: measured necessity and lean working tree

## Current state, 2026-09-28

**31.03 GiB of historical search NPZ payloads have left the active working tree.**
The 1,860 files are in recoverable Trash and also remain byte-identical in Git
at `8f06a65aaae6051095161f46217dc57a1bf2751f`. No Git history was rewritten.
The current checkout uses sparse patterns to omit these tracked files normally;
missing archived payloads are intentional, not an instruction to rerun science.

- [Exact payload index, sizes and hashes](storage/historical-search-payloads.json)
- [Measured timeline, verification and restore control](storage/necessity-audit.json)
- [Portable sparse rules](storage/historical-search.sparse)
- [Explicit restore command](../scripts/restore_nsc_search_payloads.py)

The Trash copy is `/Users/admin/.Trash/nsc-historical-search-20260928-123024/`.
It retains the original relative paths, a complete `restore-receipt.json` and
an append-only `moves.jsonl`. Each file was checked against its Git blob and
SHA-256 before relocation; the atomic same-volume move preserved its inode,
size and modification time. A genuine 1,702,678-byte payload was restored from
pinned Git, hash-checked and returned to Trash as a separate restore control.

## Why the working tree grew

These are sums of tracked file bytes in historical Git trees; they exclude
that checkout's historical `.git` disk allocation:

| Date / commit | Tracked files |
|---|---:|
| 20 September, `fd3c184c` | 2.01 GiB |
| 21 September, `6a4b6697` | 4.49 GiB |
| 22 September, `8f5ca8d0` | 10.02 GiB |
| 23 September, `65c824c4` | 33.28 GiB |
| Before this storage change, `8f06a65a` | 34.04 GiB |

The large increase was saved Newton/high-mode search data. Thirty directories
under `nsc-ks-coupled-newton*` and one `nsc-ks-n128-high-modes` directory account
for about 31.08 GiB including their small JSON records. Their NPZ files occupy
31.03 GiB. Of the NPZ member storage, 30.30 GiB is `operator/tangent.npy` and
`operator/tangent_z.npy`; only 0.73 GiB is other arrays. These NPZ files used
ZIP_STORED. They are distinct candidate derivatives, not literal duplicate blobs.

## What remains necessary

The unchanged source inventory and retained upstream families, current history
coefficients, corrected value measurements, current field/UV pilots, baseline
certificates and numerical methods remain resident. All small historical search
JSON registers, scripts and documents also remain: they preserve coefficients,
assembled derivatives, measured residuals, rejected steps and scientific scope.
The 1.3 GiB older validated fine trajectory is retained, as are unclassified
raw archived runs. The broad 43-page paper and focused 5-page paper are unchanged.

Grok's independent read-only audit `f0d3b18b-4824-45cc-8ffd-0cbde6707b84`
confirmed the giant n64 caches belong to older 47-node histories. The current
32+32-coordinate, 129-node path uses the fixed source and creates fresh operators.
Its starting-radius calculation reads the small iterate5/iterate6 JSON records,
not their old operator NPZ files.

## Verification of necessity

A read interceptor first denied all opens under the 31 historical directories,
including Git subprocess fallbacks. With that denial active:

- all 60 positive source families / 6,092 signed batches / 96,336 signed energy
  rows loaded; current (2,32) coefficients and the (64,64) metric initialized;
- current field binding, physical-source inventory replay, UV-v3 replay,
  reference and endpoint pilots, and the continuous-method recorder passed;
- the v4 component budget replay passed: N known error sum
  `8.999152980170087e-12`, still six missing components;
- the actual numerical-search admission still refused the incomplete budget.
  No scientific solve was launched and no admission criterion was bypassed.

After the real move, the source/setup load, component-budget replay, field
binding, UV-v3 replay and endpoint pilot passed again. The public repository
check authenticated 662 files and 100 historical steps; PDF integrity passed.
Three restore-tool tests cover exact-byte recovery, rejection of an existing
changed destination, corrupt digests, traversal and duplicate index entries.

One pre-existing check is separately recorded: reanchor-v2 `--check` refuses
its live `nsc_ks_history_evaluator.py` hash before accessing any denied directory.
That check is not reported as passing or blamed on the move.

A broad literal-reference crawl reached 21.29 GiB because it followed old
candidate ancestry and their replay payloads. That is historical provenance,
not proof that all those arrays are operands of the current computation.
The already-declared formal closure for the three enclosed components has
1,020 nodes and no reference to the removed directories (about 0.40 GiB of
current file bytes). Existing scientific records were not rewritten to trim
that closure. A future final certificate must still authenticate its actual
proof dependencies.

## Historical replay and restoration

Old full-search `--check` commands and tests such as
`test_nsc_gate_stage0_replay.py` require their archived payloads. Restore the
specific data and use the historical code version where its bindings require
it. Their absence must not be presented as a passed historical replay.

```sh
python3 scripts/restore_nsc_search_payloads.py --status
python3 scripts/restore_nsc_search_payloads.py --restore results/development/artifacts/nsc-ks-coupled-newton-trial/family-24-+1.npz
python3 scripts/restore_nsc_search_payloads.py --restore-directory results/development/artifacts/nsc-ks-coupled-newton-n64-primal-norm-trust03
```

The command retrieves original blobs from the pinned commit, verifies size and
SHA-256 and refuses to overwrite a different file. It requires no numerical
regeneration. Keep the Git history: it is the durable recovery copy even if
Trash is later emptied. Sparse rules are local to this checkout; the four older
worktrees were checked and did not acquire sparse mode. On a new lean checkout,
apply the supplied sparse rules intentionally; a normal full checkout otherwise
materializes the historical payloads again. `git sparse-checkout disable` also
restores the full tracked tree and its size.

## Remaining disk use and prior cleanup

Resident results fell from about 34 GiB to about 2.97 GiB. The whole repository
is now about 34.7 GiB including roughly 30.4 GiB of Git and 1.2 GiB of raw runs.
This does not make the whole repository 2 GiB: Git still preserves the original
history. Changing that requires a separate archive/history design which retains
pinned scientific commits and the other worktrees; no such rewrite was done.

Earlier cleanup moved 28 disposable cache/Finder/PDF-preview targets (22.55 MiB)
to `/Users/admin/.Trash/nsc-disposable-20260928-115059`, with checked hashes and
restore receipts. The compiler, environments and PDF originals were retained.

Moving files to Trash reduces the active folders. Disk space becomes available
only after Trash is emptied; this operation did not empty it. The local physical
gate remains OPEN. Storage reduction is not scientific certification.
