# Active laboratory — routing index

This is the `lab/` subtree of the unified NSC repository, not a separate Git repo.
Read the [root index](../AGENTS.md), [instructions](../docs/instructions.md),
and the latest user instruction in the conversation. An accepted plan in the
conversation applies where it stays compatible with that instruction and does
not override later steering. Do not use the historical lab PLAN/handover as
current task instructions or recreate the deleted handover. These pages
supersede old machine, branch and sibling-repository routes below.

| Area | Owner |
|---|---|
| Derivations and assumptions | [Claim ledger](docs/claim-ledger.md), [code map](docs/active-code-map.md) |
| Implementation and run entry points | `src/`, `scripts/`, `tests/` |
| Scientific records and payloads | [Results](results/README.md), `results/development/` |
| Current method evidence | [Code map](docs/active-code-map.md), [claim ledger](docs/claim-ledger.md) |
| Historical search exclusions | [Consolidation](../docs/repository-consolidation.md) |

Run Python through root `scripts/lab.py` so the correct source owners and pinned
historical evidence are selected. Update the owning scientific docs and root index when relevant
work changes. Current status is the [gap map](docs/claim-ledger.md#scientific-gap-map).
v5, `scripts/derive_nsc_spherical_galerkin_refinement_v5.py`, is the saved
\(T=0.005\) diagnostic and preparation, not a permanent production-driver
requirement. The episode driver
`scripts/derive_nsc_spherical_feedback_episode.py` reuses it. Its saved
\(T=0.05\) record is
`results/development/nsc-spherical-feedback-episode-v1.json`. The weak
residual is `docs/nsc-spherical-cauchy-weak.md`. The streamed reducer is
the `backend="streamed"` path in
`src/recursive_horizons/nsc_evolving_reduction.py`.
[Archived candidate notes](docs/active-code-map.md#archived-candidate-notes)
are not the queue. The local incoming gate remains OPEN and its campaign is paused;
file integrity is not scientific closure. Scope:
[root instructions](../docs/instructions.md#scientific-framing-and-authority).
