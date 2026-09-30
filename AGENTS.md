# NSC — one repository, active index

Read [instructions](docs/instructions.md) and the latest user instruction in
the conversation. That instruction overrides a stale plan. Repository
plan/handover copies are historical, not task authority. Do not recreate the
deleted handover or another plan file. Keep this index in view during work
and update its owners when relevant files, behavior or scope change. Use
explicit current user instructions and live evidence rather than remembered
or historical continuation points.

```text
AGENTS.md                 entry point and maintained index
 docs/                    active instructions, scientific status and paper guides
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
| Task authority and scientific status | Latest user instruction, which overrides a stale plan; [scope](docs/instructions.md#scientific-framing-and-authority), [status](docs/current-result.md), [code map](lab/docs/active-code-map.md) |
| Existing mathematics, assumptions and gaps | [Lab claim ledger](lab/docs/claim-ledger.md), [code map](lab/docs/active-code-map.md) |
| Source and UV certification methods | [Continuous source-error insertion](lab/docs/nsc-ks-source-operator-majorant.md), [family 14_1 completion](lab/docs/nsc-ks-source-operator-majorant-v7.md), [preparation correction](lab/docs/nsc-vacuum-source-correction.md), [high-energy remainder](lab/docs/nsc-vacuum-source-remainder.md), [full method map](lab/docs/active-code-map.md) |
| Field method owners | [Whole cone](lab/docs/nsc-ks-whole-cone-field-v1.md), [radius coupling](lab/docs/nsc-ks-radius-coupling-bounds.md), [endpoint contraction](lab/docs/nsc-ks-endpoint-contraction.md) |
| Both papers and public status | [Paper guide](docs/papers.md), [current result](docs/current-result.md), [catalog](paper/catalog.json) |
| Finite nested-quality theorem | [Proof](docs/nsc-nested-qualities.md), [frozen evidence](paper/nested-quality-evidence/snapshot.json) |
| Finite circulation and local response | [Note](lab/docs/nsc-finite-turnover.md), [record](lab/results/development/nsc-finite-turnover-v1.json) |
| Incoherent imbalance mean circulation | [Note](lab/docs/nsc-imbalance-turnover.md), [record](lab/results/development/nsc-imbalance-turnover-v1.json) |
| Original source snapshot and large-data exclusions | [Import receipt](docs/lab-snapshot.json), [consolidation](docs/repository-consolidation.md) |
| Public verification and builds | [Reproducing](docs/reproducing.md), [draft](docs/local-gate-draft.md) |

Work on `codex/work-branch`; integrate verified work into `main` and push both
to `kodareken/nested-space-cosmology`. Scientific work uses `lab/` owners, not
the older root publication code.

Current question: Douglas's proposed mechanism in the
[claim ledger](lab/docs/claim-ledger.md). Regeneration stays first. The
incoherent-imbalance calculation is the
[note](lab/docs/nsc-imbalance-turnover.md) and
[record](lab/results/development/nsc-imbalance-turnover-v1.json), not a new
roadmap. The numerical incoming-gate campaign is paused; the gate remains OPEN.
[Scope](docs/instructions.md#scientific-framing-and-authority).
