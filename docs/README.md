# Docs route

Use this page to find an owner. Claim status stays in the owning note.
This page does not merge independent strands, and it does not turn a
constrained calculation into a settled theory.

| Need | Owner |
|---|---|
| Commands, branches, and scientific framing | [instructions.md](instructions.md) |
| Current calculations and verdicts | [current-result.md](current-result.md) |
| Regeneration controls on the saved episode | [Note](../lab/docs/nsc-regeneration-controls.md), [record](../lab/results/development/nsc-regeneration-controls-v1.json) |
| Regeneration-controls successor | [Record](../lab/results/development/nsc-regeneration-controls-v2.json) |
| Conditional local response | [Note](../lab/docs/nsc-coupled-local-response.md), [record](../lab/results/development/nsc-coupled-local-response-v1.json) |
| Spherical null expansion | [Note](../lab/docs/nsc-spherical-null-expansion.md), [record](../lab/results/development/nsc-spherical-null-expansion-v1.json) |
| Local-boundary review | [Record](../lab/results/development/nsc-local-boundary-review-v1.json) |
| Initial geometry-graded window | [Note](../lab/docs/nsc-geometry-graded-window.md), [record](../lab/results/development/nsc-geometry-graded-window-v1.json) |
| Both papers | [papers.md](papers.md) |
| Reproduction | [reproducing.md](reproducing.md) |
| How the laboratory was imported | [repository-consolidation.md](repository-consolidation.md) |
| Focused OPEN draft notes | [local-gate-draft.md](local-gate-draft.md) |
| Finite nested-quality construction | [nsc-nested-qualities.md](nsc-nested-qualities.md) |
| Prior art and the open claim | [prior-art-and-open-claim.md](prior-art-and-open-claim.md) |
| 10 September 2026 development record | [development-update-2026-09-10.md](development-update-2026-09-10.md) |
| Historical road-to-arXiv plan | [PLAN.md](PLAN.md) |

`PLAN.md` is retained historical input, not a work queue. `lab/PLAN.md` is
the same plan except one claim-ledger link: this copy is repository-relative,
and the laboratory copy uses an absolute path plus a handover link.
`lab/scripts/check_current.py` still raises if `lab/PLAN.md` exists. That
checker and the retention sentence in [instructions](instructions.md)
disagree. Both copies stayed, because the plan still carries its historical
claim contract and removing the laboratory copy would change a hash-locked
import for a reason the current instructions do not ask.

`lab/handover.md` stays as a historical checkpoint. The active code map cites
it for dated measurements. It does not assign the next task.

## Copies and separate strands

Most `nsc-*.md` names also exist under `lab/docs/` with identical bytes.
Those pairs are one imported strand stored twice, so publication links and
the laboratory import receipt each have a path. They are not two results.

Same name, different bytes:

- `nsc-closure-verification-2026-09-07.md` — the `docs/` copy adds that the note is the 59th public record and that the original 58 JSON records stay byte-for-byte. The laboratory copy omits those two sentences.
- `nsc-prior-art-reuse.md` — the laboratory copy is the longer procedure, including literature boundaries and incoming-gate delivery wording. The `docs/` copy is the shorter public procedure.
- `nsc-quantum-measure-reuse.md` — the laboratory copy adds a paragraph on the canonical spherical action map. The `docs/` copy does not.

Notes that exist only under `lab/docs/` are later or laboratory-only strands,
including the claim ledger, the active code map, the spherical-coupling
series, and the saved consumers
[nsc-regeneration-controls.md](../lab/docs/nsc-regeneration-controls.md),
[nsc-coupled-local-response.md](../lab/docs/nsc-coupled-local-response.md),
[nsc-spherical-null-expansion.md](../lab/docs/nsc-spherical-null-expansion.md),
and
[nsc-geometry-graded-window.md](../lab/docs/nsc-geometry-graded-window.md).
Those notes have no `docs/` copy. The controls successor and the
local-boundary review are JSON records under `lab/results/development/`.
`lab/archive/` is historical material inside its recorded domain.

The public repository `nsc-public-polish` is a separate Git history, not a
worktree of this one. It has no `lab/` tree. Its shared `docs/` names matched
these copies at this cleanup except `current-result.md` and `reproducing.md`.
It does not carry `instructions.md`, `papers.md`, `PLAN.md`,
`local-gate-draft.md`, `repository-consolidation.md`,
`nsc-nested-qualities.md`, or `nsc-local-observer-correspondence.md`.
