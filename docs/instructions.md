# Instructions for the unified NSC repository

The laboratory and publication now belong to **one checkout and one Git history**.
The maintained branches are `main` and `codex/work-branch` at
`kodareken/nested-space-cosmology`. Work on `codex/work-branch`, test the relevant
changes, fast-forward `main`, and push both. No force-push is needed.

## Establish context after compaction

Read the root AGENTS index, [PLAN.md](PLAN.md), and [handover](handover.md).
The active numerical owners are in `lab/`; root `src/`, `scripts/`, `tests/`
and `results/` retain the curated publication chain. Do not mix the two Python
package versions. Launch scientific Python commands from the repository root:

```sh
python scripts/lab.py scripts/derive_nsc_ks_gate_budget_v4.py --check
python scripts/lab.py -m pytest tests/test_nsc_ks_endpoint_contraction.py -q
```

The launcher uses the caller's interpreter, sets the lab source path and working
directory, and provides the four pinned historical sources needed by the active
evidence closure. Use the existing validation environment for scientific work;
its requirements remain in `lab/pyproject.toml`. The root publication environment
has its own pinned requirements. Neither launcher starts a campaign by itself.

[Consolidation and recovery](repository-consolidation.md) explains the old paths,
removed search arrays and historical replay limits. Older lab README, PLAN and
handover copies are retained as historical input; this root index and the active
pages under root `docs/` own current navigation. Never resume an old Windows or
archived campaign solely because a search found its instructions.

`make check test paper-check draft-check` verifies the public surface, imported
lab integrity and both papers. It does not close the scientific gate. Keep
scientific original bytes intact; new results use successors with explicit
provenance. When intentionally updating imported lab code, update the import
receipt's current hashes and record the new source revision without altering
its original import hashes. The receipt is an integrity check, not a proof.

For an intentional lab edit, use explicit paths (no blanket re-import):

```sh
python scripts/update_lab_snapshot.py lab/src/recursive_horizons/changed_owner.py
python scripts/build_publication_provenance.py
make check
```

The helper preserves original import hashes and rejects historical excluded
payloads and the pinned source cache. Review its diff and run the checks relevant
to the actual change before committing. Removing evidence needs a separate
dependency decision; this helper does not silently drop missing files.

## Current nested-quality result

The focused companion now leads with the finite construction in
[nsc-nested-qualities.md](nsc-nested-qualities.md). The previous numerical
application remains OPEN in its appendix. That finite result does not require
eternity, heat-death or a full incoming-gate campaign. The laboratory section has
resumed its approved source-fixed local-gate programme; this public
checkpoint does not mark that programme complete. Preserve different
regional states under common laws. Before claiming two descriptions conflict,
identify their region, scale, observer, state and compared quantity.

## Keep the papers connected

The organizing inheritance postulate remains

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The broad manuscript contains the physical interpretation, operator/Schur/
recursion calculations, source history and references. Its v0.26.0 source and
PDF are preserved. The companion retains one source-fixed local incoming test as an application
under `C_Sigma[g]=U_g C_up U_g†` on `I=S(1)+[0.12,0.18]`. That application is
OPEN, separately from the finite nested-operator theorem. A later edition can incorporate verified new results with an explicit
change record; do not silently substitute a new PDF for a frozen version.

Keep assumptions, imported identities, repository calculations, numerical
diagnostics and interpretations distinct. Spatial inheritance, charge
conjugation and interference retain distinct operators. The regeneration
picture remains the motivation; global recurrence and cosmological matching
are later research. Count the declared Gaussian action and induced terms once.

## Reuse and verification

Consult the prior-art reuse map and existing evidence before compute. Name the
missing connection and the smallest check that can resolve it. Stop repeated
verification that leaves the physical gap unchanged. Historical derivations
are reference material, not a queue to repeat after moving computers.

Preserve scientific original bytes and their replay inputs. Routine edits run
relevant tests and publication checks, not old numerical campaigns. Changed PDF
source requires an updated bound build, deterministic-output check and visual
review. An unchanged PDF does not need rebuilding simply because a task resumes.

Present both papers from the same front page with dates, roles, source links and
status. Release 0.27.0 remains reserved for a closed local result and complete
scientific dependency chain. Ordinary documentation of OPEN work does not claim
that release or authorize arXiv submission.

Use relative repository links and GitHub-compatible dollar-delimited math.
Credentials and private runtime state stay outside this repository. Historical
lab records retain their original machine paths as provenance; those paths are
not active routing instructions. New instructions use repository-relative links. Douglas Ek is the accountable author; significant AI
assistance is disclosed without listing AI tools as authors.
