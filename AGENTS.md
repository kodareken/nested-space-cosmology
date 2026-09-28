# NSC — one repository, active index

Read [instructions](docs/instructions.md), [the accepted plan](docs/PLAN.md),
and [the current handover](docs/handover.md). Keep this index in view during
work and update its owners when relevant files, behavior or scope change.
Do not rely on chat memory or historical continuation instructions.

```text
AGENTS.md                 entry point and maintained index
 docs/                    active instructions, PLAN.md, handover.md, paper guides
 paper/                   both manuscripts, PDFs and frozen paper evidence
 src/, scripts/, tests/   curated publication/replay owners
 results/                 published evidence and release manifests
 lab/                     active scientific workspace in this SAME Git repo
   src/, scripts/, tests/ current numerical and proof owners
   docs/                  derivations, claim ledger and scientific code map
   results/               retained source, calibration and numerical evidence
   archive/               historical code within its recorded domain
   .source-history/       exact historical source objects for current replay
```

| Need | Owner |
|---|---|
| Location, branches, commands and recovery | [Instructions](docs/instructions.md), [consolidation](docs/repository-consolidation.md) |
| Approved research goal and next step | [Plan](docs/PLAN.md), [handover](docs/handover.md) |
| Existing mathematics, assumptions and gaps | [Lab claim ledger](lab/docs/claim-ledger.md), [code map](lab/docs/active-code-map.md) |
| Active field methods | [Whole cone](lab/docs/nsc-ks-whole-cone-field-v1.md), [radius coupling](lab/docs/nsc-ks-radius-coupling-bounds.md), [endpoint contraction](lab/docs/nsc-ks-endpoint-contraction.md) |
| Both papers and public status | [Paper guide](docs/papers.md), [current result](docs/current-result.md), [catalog](paper/catalog.json) |
| Finite nested-quality theorem | [Proof](docs/nsc-nested-qualities.md), [frozen evidence](paper/nested-quality-evidence/snapshot.json) |
| Original source snapshot and large-data exclusions | [Import receipt](docs/lab-snapshot.json), [consolidation](docs/repository-consolidation.md) |
| Public verification and builds | [Reproducing](docs/reproducing.md), [draft](docs/local-gate-draft.md) |

Work on `codex/work-branch`; integrate verified work into `main` and push both
to `kodareken/nested-space-cosmology`. The local incoming gate remains OPEN.
The finite theorem does not complete the accepted local-gate plan. Different
regions can have different states under the same law; eternity is not a proof
target. Scientific work uses `lab/` owners, not the older root publication code.
