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

The question is the nested mechanism. Structure that one region does not
resolve can still leave a local response, and the same relations can recur
across regional gradients. Familiar local laws are the equations a local
observer already uses. Comparing or fitting LambdaCDM is not this task.

The inheritance relation below is reused inside its stated domain. Further
constructions stay finite. Readings in terms of time, antimatter, a dark
sector, or a continuing gradient stay open research. Rooms, coins, fountains,
and clocks are motivation for that research. A finite candidate can be chosen
creatively. What gets reported is the computation and an error statement that
another person can check, in a named domain.

The source-fixed local incoming gate is one retained application. It remains
OPEN, and its error-budget campaign is paused. Its interval, threshold, and
source apply only there. The pause does not create a new compulsory gate, and
the gate is not a prerequisite for the finite nested-quality theorem or for
the spherical evolution below. The next release version follows a
demonstrated result. This instruction does not reserve a version for
closure of the paused gate. Frozen v0.26.0 remains the recorded
foundation. Before treating two descriptions as the same
claim, name the region, scale, observer, state, and compared quantity.

## Present regeneration programme

The current programme is finite nested regeneration under those local laws.
Starting checkpoint `9a9090a`. The executed loop, and the code owner, is
`lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py`.
The v5 replay,
`lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py`,
remains the self-contained historical seed. The saved \(T=0.05\)
measurement is `lab/scripts/derive_nsc_spherical_feedback_episode.py`,
which reuses that preparation. Its four runs are the recorded episode.

`compose_fine_hamiltonian` recomputes the column source and the geometric
rates on the prolonged state. With the matter force included it feeds
`force_Q` into the conjugate momentum,

$$
\dot p_Q \leftarrow \dot p_Q - F_Q/\Delta x_Q.
$$

\(Q\), \(r\), \(\chi\), the momenta and the spinor columns evolve. \(L\) and
\(\beta\) stay gauge controls. `solve_initial_radius`, and the v5 replay's
full-source Newton, set the initial radius from that source
(\(\rho = F_L/\Delta x_Q\)).

The inheritance postulate used by this work is

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The active action is \(\Gamma_{\mathrm{one}}\): the canonical Gaussian plus
the same-spectrum local induced term, counted once. The vacuum-matched CTP
branch remains a separate historical choice.

The independently replayed \(n_f=512\), \(T=0.005\) window in the v5 record
is a geometric diagnostic. Its stored maxima are a full Hamilton residual of
about \(9.98\times 10^{-6}\), a full momentum residual of about
\(1.38\times 10^{-6}\), and a maximum absolute proper radial velocity of
about \(0.0079505\). That velocity figure is a max-abs value. Signed motion
is not read from it. The saved \(T=0.05\) episode reaches the requested
time. All 62 physical refinement rows meet one percent. Proper velocity
on the fine primary run ends near \([-0.061956,0.088159]\), and the field
and gravitational energies exchange about \(0.4123\). The fine full
Hamilton residual ends near \(2.83\times 10^{-5}\) and the fine momentum
residual near \(4.30\times 10^{-6}\); the coarse residual is higher.
Continuum constraint control and renewal are not proved. The weak helper
records a conditional coercivity constant and a Fourier-product
conditioning comparison. It does not certify a total initial-error bound.
The historical initial tolerance \(10^{-8}\) is not an
automatic veto of this diagnostic. A weak residual alone is not sufficient;
the full residuals stay reported. The toy \(H\) and \(B\) fractions are an
optional witness. Geometry-derived links are the target.

This loop is the present finite realization because its geometry and its
state forces already come from the same action. An extra invented rule for
the radius as a function of energy is unnecessary. The historical hand-set
\(H\) and \(B\) fractions, and the historical strong tolerance \(10^{-8}\),
do not set the objective.

Construction stays flexible. Identify the shared structure, propose one
finite realization, and check feasibility, CPU, and the scale of the effect
before any larger run. The recorded outcome guides the next model step.
Progress is an effect, a derived relation, a removed assumption, or a named
blocker. A passing label, a test count, or a longer document is not that
progress.

Recorded numbers and owners:
[current result](current-result.md#present-spherical-loop).

## Papers

The [focused companion](nsc-nested-qualities.md) leads with the finite
nested-operator construction: coupled regions, normalized inheritance,
corresponding spectra and an exact reduced local response on a finite window.
Regional states may differ under the same law. The construction is not a
proof of infinite physical time, and it does not depend on the incoming-gate
campaign. The gate remains an OPEN appendix of that companion.

The broad manuscript, v0.26.0 source and PDF, stays frozen. It holds the
physical motivation, operator and Schur calculations, source history and
references. Statements in that dated file about the next step are history.
A later edition can take a verified result only with an explicit change
record. Do not replace a frozen PDF in place.

Assumptions, imported identities, repository calculations, numerical
diagnostics and interpretations stay distinct. Spatial inheritance, charge
conjugation and interference keep distinct operators. Global recurrence and
cosmological matching are not established.

## Reuse and verification

Consult the prior-art reuse map and existing evidence before compute. Name the
missing connection and the smallest check that can resolve it. Stop repeated
verification that leaves the physical gap unchanged. Historical derivations
are reference material.

Preserve scientific original bytes and their replay inputs. Routine edits run
the relevant tests and publication checks. The paused incoming-gate campaign
and any new expensive generator stay off that path. Changed PDF source
requires an updated bound build, a deterministic-output check and visual
review. An unchanged PDF does not need rebuilding because a task resumes.

Present both papers from the same front page with dates, roles, source links
and status. Ordinary documentation of OPEN work does not assign the next
release version or authorize arXiv submission. Frozen v0.26.0 stays the
recorded foundation. Historical release plans remain dated history.

Use relative repository links and GitHub-compatible dollar-delimited math.
Credentials and private runtime state stay outside this repository. Historical
lab records retain their original machine paths as provenance; those paths are
not active routing instructions. New instructions use repository-relative
links. Douglas Ek is the accountable author; significant AI assistance is
disclosed without listing AI tools as authors.
