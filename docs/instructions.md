# Instructions for the unified NSC repository

The laboratory and publication now belong to **one checkout and one Git history**.
The maintained branches are `main` and `codex/work-branch` at
`kodareken/nested-space-cosmology`. Work on `codex/work-branch`, test the relevant
changes, fast-forward `main`, and push both. No force-push is needed.

## Establish context after compaction

Read the root AGENTS index and Douglas's latest instruction in the
conversation. That instruction overrides an older approved plan when they
conflict. On 28 September 2026 Douglas removed the active handover and
explicitly retired repository PLAN.md copies as task authorities. Do not restore
the handover, duplicate a plan, or infer a new task from those old documents.
The 28 September road-to-arXiv cursor, which made one incoming-gate campaign
the compulsory next step, does not override the later pause of that campaign.
If the current instruction is missing, ask Douglas instead of resuming a
historical file or an expensive numerical run. Current evidence is owned by the
[claim ledger](../lab/docs/claim-ledger.md), [code map](../lab/docs/active-code-map.md)
and versioned result records.
The active numerical owners are in `lab/`; root `src/`, `scripts/`, `tests/`
and `results/` retain the curated publication chain. Do not mix the two Python
package versions. Launch scientific Python commands from the repository root:

```sh
python scripts/lab.py scripts/derive_nsc_ks_gate_budget_v4.py --check
python scripts/lab.py -m pytest tests/test_nsc_ks_endpoint_contraction.py -q
```

These examples show the launcher. They are not an order to continue the paused
gate campaign. The launcher uses the caller's interpreter, sets the lab source path and working
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

## Scientific framing and authority

The question is whether specified nested surroundings produce an effective
local response under familiar local laws. Nested surroundings means structure
that is not resolved in the local description. An effective local response is
what remains in the local equations after that structure is reduced, usually
as an extra source. Familiar local laws means the local equations a local
observer already uses, not a demand for new local physics.

LambdaCDM is the standard cosmological model: general relativity with ordinary
matter, cold dark matter, and a cosmological-constant-like term. It is one
familiar description of that local response. Refuting it, or treating an open
comparison as a conflict, is not the task. A completed match to every LambdaCDM
observable is also not established. Same reduced equation form and full
physical equivalence are separate claims. Do not remove local terms or change
measured local laws. Do not claim dark matter, dark energy, antimatter, or
eternity as results already derived.

Keep three layers distinct:

- Reusable mathematics, including block reduction and the inheritance relation,
  stands only inside its stated domain. Do not recalculate that relation, and
  do not replace it with an infinity of nests.
- A proposed interpretation, including regeneration or a dark-sector reading,
  is not a theorem. Schur elimination does not by itself identify cosmological
  dark sectors.
- One chosen numerical application, the source-fixed local incoming gate, tests
  one declared geometry, preparation, source inventory, and history. Its
  interval, residual threshold, source, and history apply only there.

The gate remains OPEN: neither local existence nor a controlled exclusion is
established. Its expensive error-budget campaign is paused. The gate is a
retained application, not the centre of the next task. Its relevance can be
reconsidered when useful. That reconsideration is not a compulsory essay,
check, or permission barrier before other work. Pausing does not ban later
error-controlled numerical work on an appropriate question, and compaction is
not a reason to restart this campaign. Pausing does not discard that
application's error bounds or locked parameters. A failure or gap in that
class is not a verdict on LambdaCDM or on every nested interpretation, and it
is not a prerequisite for other nested claims or for arXiv discussion in
general. Release 0.27.0 remains reserved for a closed result of this
application.

Next research reuses existing inherited-law, regional-state, transfer and
conservation results to investigate sustained balanced turnover at stable
overall regional size and content. Relate that question to time as the
registration of spatial change, and to the matter–antimatter hypothesis.
Those readings are not proved. Choose tests that serve this mechanism. Do not
make a written assessment of the paused gate the required next act.

Pictures of an interior and its surroundings, a room, or two sides of a coin
organize the question. They are not literal equations, and they are not claims
the author must defend as physical identities. The present spherical
calculation is a controlled symmetry class. Regional states may differ under a
common law. Reuse the established operator and transfer results. A mathematical
identification must keep the domain and assumptions actually derived.

Charge conjugation is the standard map between particle and antiparticle
sectors, including reversal of gauge charge. Using that name does not derive
antimatter from the nested geometry, and it does not forbid investigating a
deeper geometric origin. It also does not make sheet exchange, an interior, or
a metaphor into antimatter. Those remain distinct questions, not new
completion requirements.

## Current nested-quality result

The focused companion now leads with the finite construction in
[nsc-nested-qualities.md](nsc-nested-qualities.md). The previous numerical
application remains OPEN in its appendix, and its campaign is paused. That
finite result does not require eternity, heat death, or the incoming-gate
campaign. This checkpoint does not mark the gate complete and does not make
the gate a prerequisite for the finite result. Preserve different regional
states under common laws. Before claiming two descriptions conflict, identify
their region, scale, observer, state and compared quantity.

## Keep the papers connected

The organizing inheritance postulate remains

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The broad manuscript contains the physical interpretation, operator/Schur/
recursion calculations, source history and references. Its v0.26.0 source and
PDF are preserved. The companion retains one source-fixed local incoming test as an application
under `C_Sigma[g]=U_g C_up U_g†` on `I=S(1)+[0.12,0.18]`. That application is
OPEN and its campaign is paused. It is separate from the finite nested-operator
theorem, and its threshold is not a prerequisite for that theorem or for nested
claims in general. A later edition can incorporate verified new results with an explicit
change record; do not silently substitute a new PDF for a frozen version.

Keep assumptions, imported identities, repository calculations, numerical
diagnostics and interpretations distinct. Spatial inheritance, charge
conjugation and interference retain distinct operators. Charge-conjugation
terminology does not close a geometric question. The regeneration picture
remains motivation, not a literal proof; global recurrence and cosmological
matching are not established. Count the declared Gaussian action and induced terms once.

## Reuse and verification

Consult the prior-art reuse map and existing evidence before compute. Name the
missing connection and the smallest check that can resolve it. Stop repeated
verification that leaves the physical gap unchanged. Historical derivations
are reference material, not a queue to repeat after moving computers.

Preserve scientific original bytes and their replay inputs. Routine edits run
relevant tests and publication checks, not the paused incoming-gate campaign
or another expensive generator. Changed PDF
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
