# Future-only evidence I/O

`src/recursive_horizons/evidence_io.py` is a Phase −1 **infrastructure**
helper for later artifacts. It is not scientific authority, not compact
evidence, not a binder, and not a runner. It does not classify a terminal,
admit a step, or decide whether a result may be promoted.

Historical binders and runners must not import it. Do not treat a successful
read, Git query, or exclusive publication as a physical or numerical claim.
The no-retrofit test freezes the source/script path set at `5ec2530`; newly
created prospective artifacts may use the utility without changing that
historical set.

## Public API

| Function | Role |
|---|---|
| `canonical_json_bytes(value)` | Encode a JSON value as sorted compact ASCII bytes. |
| `load_canonical_json(raw)` | Parse those bytes and require an exact canonical round-trip. |
| `read_regular_file(root, relative, *, max_bytes=DEFAULT_MAX_BYTES)` | Bounded no-follow regular-file read. |
| `git_read(repository, operation)` | One closed, sanitized, read-only Git query. |
| `publish_exclusive_file(root, relative, payload)` | Publish new file bytes without replacement. |
| `publish_exclusive_directory(root, relative, files)` | Publish a new directory of regular files without replacement. |

Git operations are typed and closed:

- `ResolveCommit(revision)` — `HEAD` or a full lowercase hex object name
- `ReadBlob(commit, path)` — one regular `100644`/`100755` blob
- `InspectTree(commit, path=None)` — recursive tree listing
- `InspectDelta(commit, against=None)` — name-status delta (`against=None` is `--root`)
- `CommitParents(commit)` — ordered parent identities, without assuming a single parent
- `InspectWorktree()` — changed paths and hidden index flags; `clean` is true
  only if neither exists. Rename detection is disabled, and assume-unchanged /
  skip-worktree flags cannot masquerade as an authenticated clean image.

There is no generic Git exec, no mutation command, and no shared scientific
reduction with campaign stores.

## Failure semantics

| Exception | Meaning |
|---|---|
| `CanonicalJSONError` | Duplicate keys, non-finite numbers, subclasses, cycles, over-deep nesting, unsupported types, or non-canonical bytes. |
| `UnsafePathError` | Traversal, symlink, hard link, special file, missing `O_NOFOLLOW`, root ancestor walk, or identity race. |
| `GitQueryError` | Operation outside the closed set, unsafe revision/path, missing trusted Git, cap/timeout, or sanitized Git failure. |
| `PrepublicationError` | Destination was not published. |
| `UnsupportedPublication` | Exclusive no-replace rename is unavailable; fail closed. |
| `PostpublicationUncertainty` | Destination **was** published; a later fsync, identity, or content check failed. The published result is not deleted and not retried. |

`PrepublicationError.published` is `False`.
`PostpublicationUncertainty.published` is `True`.

## Retained staging

Staging is not deleted merely because a random staging name was selected.

This call may remove a staging path only when **all** of the following hold:

1. this call exclusively created that inode (`O_EXCL` / exclusive `mkdir`);
2. the staging name still names that same `(st_dev, st_ino)`;
3. every remaining entry and byte matches the known staging snapshot; and
4. publication has not already made the destination visible.

A pre-existing collision, substituted inode, unknown entry (including an
extra empty directory), changed byte, or staging object whose ownership this
call cannot prove is **retained**.
`PrepublicationError.staging_retained` is `True` in that case. Recursive
deletion is never applied to a foreign or swapped tree.

After an exclusive rename, leftover or uncertain staging is also retained.
A success receipt is returned only after the destination inode and
bytes/inventory are independently revalidated against the snapshot.

## Path and Git hygiene

The caller-supplied root must be an existing absolute canonical path with no
symlinked ancestors. A caller may explicitly resolve its chosen root before
calling the utility; the utility does not silently reinterpret a relative or
symlinked root. Reads and publication keep the live parent-chain
`(st_dev, st_ino)` identities and re-walk the caller-visible path before
returning success. Replacement of an ancestor during an open-fd read is a
race, not a successful stale read.

Reads refuse traversal, leaf and ancestor symlinks, `st_nlink != 1` regular
files, fifos and other special files. The default byte bound is 16 MiB.
Canonical JSON is also bounded to 16 MiB and depth 64; malformed deep input
returns a typed error rather than escaping as an interpreter recursion error.

Git queries:

- exec a trusted absolute Git from `{/usr/bin,/opt/homebrew/bin,/usr/local/bin}/git`
- use a minimal environment (`PATH=/usr/bin:/bin`, no inherited `LD_*`/`DYLD_*`/`GIT_*`)
- force `GIT_CONFIG_NOSYSTEM`, empty global/system config, and `GIT_CONFIG_COUNT=0`
- pass `--no-replace-objects`, `--no-optional-locks`, `--no-ext-diff`, empty
  `diff.external` / `diff.textconv`, and alias-emptying `-c` overrides
- capture stdout/stderr with a 16 MiB cap and kill the process on overflow or timeout
- never accept a user string that is not a closed revision or a safe relative path

## Publication

File and directory payloads are snapshotted into plain `bytes` (and a frozen
path map from a built-in `dict`) **before** any `bytes(payload)` conversion of a non-bytes object
and before staging. Entry count and total size are bounded during that
snapshot so a mutable mapping or an integer/custom payload cannot change the
receipt or allocate unbounded memory.
Directory publication remains limited to 1,024 entries. TDG11's new authority
must inspect this repository's more than 1,024 tracked paths, so its prospective
Git-read contract raises only the Git listing ceiling to 16,384 entries.
Publication limits and historical consumers are unchanged.

Staging is created exclusively, written, and fsynced, including every nested
directory, before the parent publication. After
every pre-publish hook this call revalidates the parent chain, the owned
staging inode, and the exact file bytes or directory inventory. Extra,
mutated, or replaced leaves and unexpected empty directories fail closed.
Exclusive no-replace rename
(`renameatx_np(RENAME_EXCL)` on Darwin, `renameat2(RENAME_NOREPLACE)` on
Linux) then moves the owned inode onto the destination. If that primitive is
missing, publication fails closed **before** the destination becomes visible.

The destination is independently re-read after rename. A success receipt
hashes those verified bytes or that verified inventory. A mismatch after
rename is postpublication uncertainty, never a success receipt.

An existing, partial, or racing destination is a prepublication failure. The
utility does not overwrite, resume, or treat identical bytes as success.

Directory publication exists so a later diagnostic terminal can drop a small
file tree without inventing a second I/O stack. It still performs no
classification.

## Non-claims

This module does not:

- bind or verify a compact scientific result
- authorize state advance, recovery, or candidate execution
- inspect or mutate `runs/` campaign stores
- share classification or publication decisions with historical runners

Owner documents, the claim ledger, and compact hashes remain the scientific
record. See the [active code map](active-code-map.md) for current extension
points.
