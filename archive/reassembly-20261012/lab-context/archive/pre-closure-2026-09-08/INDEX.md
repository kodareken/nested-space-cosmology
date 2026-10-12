# Frozen laboratory checkpoint

Original checkpoint: `5f38712ca01ddd71e715fd265088925a73369aba`. The [manifest](manifest.json) records all 1,731 original tracked paths, modes, blobs and SHA-256 hashes, plus protected local evidence and the worktree inventory before cleanup.

Files preserve original bytes and scientific status. Relative links/imports refer to the original layout. From the laboratory root:

```sh
python3 scripts/reproduce_history.py scale-closure
python3 scripts/reproduce_history.py core-identities
python3 scripts/reproduce_history.py repository-check
```

The launcher materializes the original checkpoint. The historical repository check reports old PLAN/AGENTS contract failures without repair; it is not a current gate. The machine-specific PLAN link is preserved as text, never followed. Its resolved former content is supplemental context. Existing earlier archives and raw stores remain unchanged. Relocation and reproduction do not promote any physical claim.

## Observed historical checker outcome

The launcher was exercised at the original Git checkpoint with Git ancestry available and no raw stores. `repository-check` returned exit 1: the machine-specific PLAN link is unavailable in a portable checkout; original AGENTS/PLAN phrase contracts disagree; several raw-dependent records cannot bind without local stores; and the FGC catalog omits the later NSC records. These failures are preserved, not repaired or reclassified as physics results. The historical `scale-closure` calculation reproduces successfully.
