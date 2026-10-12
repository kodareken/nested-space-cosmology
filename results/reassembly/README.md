# Root reassembly: preservation and research boundaries

The canonical project is the repository root. Current code uses one
`src/recursive_horizons` package and root `scripts/`, `tests/`, `docs/` and
`results/`. [The migration manifest](manifest.json) records each original
laboratory file, its original bytes and its new destination. It is a recovery
and layout receipt, not a scientific claim.

The baseline is commit `04a0cb9378977772243ebf7c9be3c74182b030f2`.
All 15,451 entries of the pre-migration laboratory receipt were checked before
moving: 3,866,941,958 bytes, with no missing file or hash mismatch. The broader
resident-file manifest also preserves ignored runtime material. Existing
different root files were not overwritten; laboratory variants have their
own archived destinations. Deliberate current edits have separate successor
hashes, while their original identities remain intact.

## What came from the second worktree

The w2 base is `d273ce6c86cacdf35d57c97292b270fac5fa6713`. Its 400 dirty
files match the 10 October preservation receipt exactly; 92 are temporary test
fixtures and 308 are research or navigation files. Of those 308, 294 are unique
research files now present in the root tree, five force-oscillation files were
already identical to root owners, and nine paths have distinct variants.
The complete 400-file snapshot remains under
[`archive/reassembly-20261012/w2`](../../archive/reassembly-20261012/w2/receipt.json).

The 294 imports remain **exploratory candidates**. Importing code or records
does not adopt a cosmological interpretation. In particular, the archived
`nsc_rn_pg.py` variant does not replace the canonical implementation. The eight
other differing paths are navigation and claim documents; their assertions
must be reconciled at current owners rather than becoming competing policy.

The prior [research inventory](../../docs/research-inventory.md) and
[strand dispositions](../research-reset/strand-dispositions.json) retain their
scientific reading boundary. This migration checked bytes, paths and ownership
collisions. It did not independently rederive the charged, influence-chain,
constants, clock, or winding calculations. Universal dark-sector elimination,
preset matter fractions and cross-universe age interpretations are not admitted
by this receipt. The corresponding original notes and corrections remain
available for future investigation.

## Preserved historical contexts

- Exact source proof objects: `archive/source-history/`.
- Former laboratory archive, run material and metadata:
  `archive/reassembly-20261012/lab-context/`.
- Distinct laboratory copies of existing root owners:
  `archive/reassembly-20261012/lab-variants/`.
- Six preview folders and six obsolete helper scripts:
  `archive/reassembly-20261012/previews-and-helpers/`.

The preview collection contains intermediate variants, not just byte-identical
duplicates. It is retained for recovery, while official atlas records remain
the evidence to cite. Historical producers and their original relative layout
are replayed in an isolated ignored build directory. Current source checks
must not pretend the old bundle was produced by relocated or edited code.

`python scripts/reassemble_repository.py --check` authenticates every current
destination and separately checks original bytes for declared successors.
`--record-successors` accepts an explicit list of intentionally edited migrated
paths; it does not refresh the original hashes. Missing or changed current
files fail verification even when historical Git bytes still exist. The
script does not delete trees, stage files, or retire worktrees.

Ignored caches and saved test fixtures have a deterministic ZIP carrier,
`archive/reassembly-20261012/local-runtime-and-ignored.zip`. Each packed row
names its member and retains its original hash. The checker authenticates the
whole ZIP and every member even when local copies still exist. A fresh checkout
can recover an absent ignored copy from that explicit carrier; corrupt current
copies are still errors. No scientific file becomes optional because Git once
ignored it.
