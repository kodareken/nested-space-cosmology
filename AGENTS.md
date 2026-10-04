# NSC — one repository, active index

Read [instructions](docs/instructions.md), the [docs route](docs/README.md),
and the latest user instruction in the conversation. That instruction
overrides a stale plan. Repository plan and handover copies are historical,
not task authority. Do not recreate the deleted handover or another plan
file. Keep this index current when ownership changes. Prefer the owning
note and live evidence over a remembered continuation point.

Scientific notes are constrained research. The method, domain, and evidence
stay with the owning file. A proposal is not an observation, and an
observation is not a settled theory. This index routes those strands. It
does not merge them, and it does not restate their numbers.

```text
AGENTS.md                 this index
 docs/                    framing, status, paper guides, and one note per public strand
   README.md              docs route, including which copies are duplicates
 paper/                   both manuscripts, PDFs, and frozen paper evidence
 src/, scripts/, tests/   curated publication and replay owners
 results/                 published evidence and release manifests
 lab/                     scientific workspace in this same Git repository
   src/, scripts/, tests/ current numerical and proof owners
   docs/                  derivations, claim ledger, and scientific code map
   results/               retained source, calibration, and numerical evidence
   archive/               historical material inside its recorded domain
   .source-history/       exact historical source objects, plus one content carrier for two recovered local-boundary sources
```

| Need | Owner |
|---|---|
| Location, branches, commands, and recovery | [Instructions](docs/instructions.md), [consolidation](docs/repository-consolidation.md) |
| Where a document lives, and which copies match | [Docs route](docs/README.md) |
| Task authority and scientific status | Latest user instruction, which overrides a stale plan; [scope](docs/instructions.md#rn-radial-programme), [status](docs/current-result.md#rn-radial-parent-4-october-2026), [code map](lab/docs/active-code-map.md) |
| Existing mathematics, assumptions, and gaps | [Lab claim ledger](lab/docs/claim-ledger.md), [code map](lab/docs/active-code-map.md) |
| RN radial parent | Scientific start `5cb775a7`. Current reproducible assessment is v2. Manual v1 remains the historical writer failure. Finite continuation to $T=20$. The numerical family is still running and is not claimed. [Status](docs/current-result.md#rn-radial-parent-4-october-2026), [rules](docs/instructions.md#rn-radial-programme), [v2 assessment](lab/results/development/nsc-rn-pilot-assessment-v2.json), [v1 record](lab/results/development/nsc-rn-pilot-assessment-v1.json), [continuation](lab/results/development/nsc-rn-neutral-continuation-v1/summary.json) |
| Historical periodic parent | NF128 exterior-to-cut response measured; the bounded NF256 confirmation is recorded. Further periodic autonomous exchange runs stay off the active path. [Status](docs/current-result.md#responsive-parent-3-october-2026), [rules](docs/instructions.md#present-regeneration-programme) |
| Strong source, kinetic-offset and frozen-column comparisons | [Method](lab/docs/nsc-discovery-parent-strong-controls.md), [consumer](lab/scripts/assess_nsc_discovery_parent_strong_controls.py), [record](lab/results/development/nsc-discovery-parent-strong-controls-v1.json) |
| Coupled exterior-to-cut response | [Method](lab/docs/nsc-discovery-parent-cut-response.md), [adapter](lab/src/recursive_horizons/nsc_discovery_parent_cut_response.py), [driver](lab/scripts/derive_nsc_discovery_parent_cut_response.py), [completed batch](lab/results/development/nsc-discovery-parent-cut-response-v1/prepare.json), [assessment](lab/docs/nsc-discovery-parent-cut-assessment.md), [record](lab/results/development/nsc-discovery-parent-cut-assessment-v1.json), [consumer](lab/scripts/assess_nsc_discovery_parent_cut_response.py) |
| Recorded NF256 cut confirmation | [Note](lab/docs/nsc-discovery-parent-cut-confirmation.md), [record](lab/results/development/nsc-discovery-parent-cut-confirmation-v1/confirmation.json), [payload](lab/results/development/nsc-discovery-parent-cut-confirmation-v1/confirmation.npz), [driver](lab/scripts/derive_nsc_discovery_parent_cut_confirmation.py). Producer `5cb775a7`. Comparison and CPU stay on the status page. |
| Parent preparation, episode, observer, response, step control and Dirac contact | [Preparation](lab/docs/nsc-discovery-parent.md), [episode](lab/docs/nsc-discovery-parent-episode.md), [observer](lab/docs/nsc-discovery-regions.md#leading-parent-observer), [response](lab/docs/nsc-discovery-parent-response.md), [step control](lab/docs/nsc-discovery-parent-step-control.md), [Dirac contact](lab/docs/nsc-discovery-dirac-contact.md) |
| Earlier discovery programme | Preserved from `ea6c699`; [record](docs/current-result.md#active-discovery-programme-2-october-2026) |
| Earlier discovery implementation | [Episode runner](lab/docs/nsc-discovery-episode.md), [observer and atlas](lab/docs/nsc-discovery-observables.md), [coupled response](lab/docs/nsc-discovery-response.md), [FFT adapter](lab/src/recursive_horizons/nsc_discovery_backend.py) |
| Discovery prediction and source-selected regions | [Matched-clock prediction](lab/docs/nsc-discovery-prediction.md), [gradient cuts and moving ledger](lab/docs/nsc-discovery-regions.md), [initial source scale](lab/docs/nsc-discovery-scale.md) |
| Discovery physical assessment and coupled reduction | [Batch assessment](lab/docs/nsc-discovery-episode-assessment.md), [four-metric tides](lab/docs/nsc-discovery-tidal.md), [causal coupled memory](lab/docs/nsc-discovery-coupled-memory.md) |
| Discovery regime family and energy accounting | [Source family](lab/docs/nsc-discovery-family.md), [disjoint RK-stage balance](lab/docs/nsc-discovery-balance.md) |
| Discovery holding balance and geometric controls | [Static spectral source](lab/docs/nsc-discovery-stationary.md), [same-action source-free solution](lab/docs/nsc-discovery-vacuum-control.md), [whole-state translation](lab/docs/nsc-discovery-translation.md) |
| Deeper inherited hand-off and future prediction | [Grandchild response](lab/docs/nsc-discovery-grandchild.md), [producer](lab/src/recursive_horizons/nsc_discovery_grandchild.py), [driver](lab/scripts/derive_nsc_discovery_grandchild.py) |
| Physical source-width response | [Width prediction](lab/docs/nsc-discovery-width-response.md), [producer](lab/src/recursive_horizons/nsc_discovery_width_response.py), [driver](lab/scripts/derive_nsc_discovery_width_response.py) |
| Ambient extent and source-consistent boundary controls | [Extent intervention](lab/docs/nsc-discovery-extent.md), [producer](lab/src/recursive_horizons/nsc_discovery_extent.py), [driver](lab/scripts/derive_nsc_discovery_extent.py) |
| General positive initial radius on a declared nonconstant geometry | [Dynamic preparation](lab/docs/nsc-discovery-dynamic-preparation.md), [constructor](lab/src/recursive_horizons/nsc_discovery_dynamic_preparation.py), [driver](lab/scripts/derive_nsc_discovery_dynamic_preparation.py) |
| Exact free-transport comparison of the saved return | [Free return](lab/docs/nsc-discovery-free-return.md), [consumer](lab/src/recursive_horizons/nsc_discovery_free_return.py), [driver](lab/scripts/derive_nsc_discovery_free_return.py) |
| Connected geometric influence of the Gaussian field | [Geometric vertices](lab/docs/nsc-discovery-geometric-vertices.md), [consumer](lab/src/recursive_horizons/nsc_discovery_geometric_vertices.py), [driver](lab/scripts/derive_nsc_discovery_geometric_vertices.py) |
| Analytic noncompact limit and conditioned static-query assessment | [Limit note](lab/docs/nsc-discovery-stationary-limit.md), [independent consumer](lab/scripts/assess_nsc_discovery_stationary_limit.py) |
| Initial geometric rate and curvature controls | [Initial sectors](lab/docs/nsc-discovery-initial-sector.md), [producer](lab/src/recursive_horizons/nsc_discovery_initial_sector.py), [driver](lab/scripts/derive_nsc_discovery_initial_sector.py) |
| Constrained curvature mode of the local action | [Assessment](lab/docs/nsc-discovery-curvature-sector.md), [symbolic owner](lab/src/recursive_horizons/nsc_discovery_curvature_sector.py), [consumer](lab/scripts/assess_nsc_discovery_curvature_sector.py) |
| Physical domain of the local curvature approximation | [Units and scope](lab/docs/nsc-discovery-eft-domain.md), [saved-state assessment](lab/scripts/assess_nsc_discovery_eft_domain.py) |
| Prepared-state geometric work and influence reference | [Pulse control](lab/docs/nsc-discovery-spectral-work.md), [producer](lab/src/recursive_horizons/nsc_discovery_spectral_work.py), [driver](lab/scripts/derive_nsc_discovery_spectral_work.py) |
| Leading common-action field–geometry comparison | [Method](lab/docs/nsc-discovery-leading-einstein.md), [model](lab/src/recursive_horizons/nsc_discovery_leading_einstein.py), [completed v2](lab/results/development/nsc-discovery-leading-einstein-v2/manifest.json) |
| Physically scaled numerical continuation | [Step control](lab/docs/nsc-discovery-leading-step-control.md), [numerical owner](lab/src/recursive_horizons/nsc_discovery_leading_step_control.py) |
| Static periodic obstruction and required exterior traction | [Virial](lab/docs/nsc-discovery-leading-virial.md), [assessment](lab/scripts/assess_nsc_discovery_leading_virial.py) |
| Source-responsive nested collars and parent continuation criterion | [Collar](lab/docs/nsc-discovery-parent-traction.md), [producer](lab/src/recursive_horizons/nsc_discovery_parent_traction.py), [withheld measurement](lab/results/development/nsc-discovery-parent-traction-v1/measurement.json), [continuation criterion](lab/docs/nsc-discovery-parent-continuation.md) |
| Source-history-selected turn and withheld full-field prediction | [Turn](lab/docs/nsc-discovery-homogeneous-turn.md), [reduction](lab/src/recursive_horizons/nsc_discovery_homogeneous_turn.py), [independent field backend](lab/src/recursive_horizons/nsc_discovery_turn_fullfield.py) |
| Completed finite parent–child pair | [Owning note](lab/docs/nsc-nested-parent-child.md), [model](lab/src/recursive_horizons/nsc_nested_parent_child.py), [response](lab/src/recursive_horizons/nsc_nested_parent_child_response.py), [confirmation driver](lab/scripts/derive_nsc_nested_parent_child_confirmation.py), [v2 evidence](lab/results/development/nsc-nested-parent-child-confirmation-v2.json), [preserved v1](lab/results/development/nsc-nested-parent-child-v1.json) |
| Parent–child replay and figure | [Portable basis](lab/results/development/nsc-nested-parent-child-replay-basis-v1.json), [independent consumer](lab/scripts/check_nsc_nested_parent_child_confirmation.py), [figure](lab/results/development/nsc-nested-parent-child-figure.png) |
| Present regeneration loop | [Status](docs/current-result.md#present-spherical-loop), [code owner](lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py), [saved v5 diagnostic](lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py) |
| Same-action conformal gauge and episode | [Status](docs/current-result.md#same-action-conformal-episode-1-october-2026), [proof](lab/docs/nsc-spherical-conformal-gauge.md), [driver](lab/scripts/derive_nsc_spherical_conformal_episode.py), [record](lab/results/development/nsc-spherical-conformal-episode-v1.json), [payload](lab/results/development/nsc-spherical-conformal-episode-v1.npz), [gauge tests](lab/tests/test_nsc_spherical_conformal_gauge.py), [episode tests](lab/tests/test_nsc_spherical_conformal_episode.py) |
| Finite conformal chain through $T=0.3$ | [Owning note](lab/docs/nsc-spherical-conformal-continuation.md), [continuation](lab/results/development/nsc-spherical-conformal-episode-v2.json), [surface flow](lab/results/development/nsc-spherical-conformal-transport-v2.json), [clock](lab/results/development/nsc-spherical-conformal-clock-v2.json), [tests](lab/tests/test_nsc_spherical_conformal_continuation.py), [preservation tests](lab/tests/test_nsc_spherical_conformal_preservation.py) |
| Same-realization conformal local response | [Note](lab/docs/nsc-conformal-local-response-successor.md), [module](lab/src/recursive_horizons/nsc_conformal_local_response_successor.py), [driver](lab/scripts/derive_nsc_conformal_local_response_successor.py), [record](lab/results/development/nsc-conformal-local-response-v2.json), [tests](lab/tests/test_nsc_conformal_local_response_successor.py) |
| Independent conformal memory control | [Owning local note](lab/docs/nsc-conformal-local-response-successor.md), [module](lab/src/recursive_horizons/nsc_conformal_memory_control.py), [record](lab/results/development/nsc-conformal-memory-control-v1.json), [tests](lab/tests/test_nsc_conformal_memory_control.py) |
| Conformal source controls and actual metric curvature | [Source controls](lab/docs/nsc-spherical-conformal-source-controls.md), [source record](lab/results/development/nsc-spherical-conformal-source-controls-v1.json), [curvature](lab/docs/nsc-spherical-conformal-curvature.md), [curvature record](lab/results/development/nsc-spherical-conformal-curvature-v1.json), [curvature tests](lab/tests/test_nsc_spherical_conformal_curvature.py) |
| Spherical feedback episode, $T=0.05$ | [Status](docs/current-result.md#recorded-t005-feedback-episode), [driver](lab/scripts/derive_nsc_spherical_feedback_episode.py), [record](lab/results/development/nsc-spherical-feedback-episode-v1.json), [payload](lab/results/development/nsc-spherical-feedback-episode-v1.npz), [tests](lab/tests/test_nsc_spherical_feedback_episode.py), [independent tests](lab/tests/test_nsc_spherical_feedback_episode_independent.py) |
| Maintained continuation, $T=0.05$ to $0.085$ | [Status](docs/current-result.md#regeneration-episode-1-october-2026), [note](lab/docs/nsc-regeneration-episode.md), [driver](lab/scripts/derive_nsc_regeneration_episode.py), [record](lab/results/development/nsc-regeneration-episode-v1.json), [payload](lab/results/development/nsc-regeneration-episode-v1.npz), [tests](lab/tests/test_nsc_regeneration_episode.py) |
| Saved episode assessment and independent audit | [Assessment](lab/docs/nsc-spherical-episode-assessment.md), [module](lab/src/recursive_horizons/nsc_spherical_episode_assessment.py), [tests](lab/tests/test_nsc_spherical_episode_assessment.py), [audit](lab/docs/nsc-regeneration-realization-audit.md), [independent tests](lab/tests/test_nsc_regeneration_realization_independent.py) |
| Regeneration controls | [Status](docs/current-result.md#regeneration-controls-1-october-2026), [note](lab/docs/nsc-regeneration-controls.md), [module](lab/src/recursive_horizons/nsc_regeneration_controls.py), [driver](lab/scripts/derive_nsc_regeneration_controls.py), [record](lab/results/development/nsc-regeneration-controls-v1.json), [payload](lab/results/development/nsc-regeneration-controls-v1.npz) |
| Regeneration-controls successor | [Status](docs/current-result.md#regeneration-controls-successor-1-october-2026), [record](lab/results/development/nsc-regeneration-controls-v2.json), [independent tests](lab/tests/test_nsc_regeneration_uniform_independent.py) |
| Conditional local response | [Status](docs/current-result.md#conditional-local-response-1-october-2026), [note](lab/docs/nsc-coupled-local-response.md), [module](lab/src/recursive_horizons/nsc_coupled_local_response.py), [driver](lab/scripts/derive_nsc_coupled_local_response.py), [record](lab/results/development/nsc-coupled-local-response-v1.json), [series](lab/results/development/nsc-coupled-local-response-v1.npz) |
| Same-trajectory local-response successor | [Status](docs/current-result.md#same-trajectory-local-response-successor-1-october-2026), [note](lab/docs/nsc-coupled-local-response.md), [record](lab/results/development/nsc-coupled-local-response-v2.json), [series](lab/results/development/nsc-coupled-local-response-v2.npz) |
| Independent local-response controls | [Note](lab/docs/nsc-coupled-local-response.md#independent-control-replay), [driver](lab/scripts/derive_nsc_coupled_local_response_audit.py), [record](lab/results/development/nsc-coupled-local-response-audit-v1.json), [tests](lab/tests/test_nsc_coupled_local_response_independent.py) |
| Spherical null expansion | [Status](docs/current-result.md#spherical-null-expansion-1-october-2026), [note](lab/docs/nsc-spherical-null-expansion.md), [module](lab/src/recursive_horizons/nsc_spherical_null_expansion.py), [driver](lab/scripts/derive_nsc_spherical_null_expansion.py), [record](lab/results/development/nsc-spherical-null-expansion-v1.json) |
| Local-boundary review | [Status](docs/current-result.md#local-boundary-review-1-october-2026), [v1 record](lab/results/development/nsc-local-boundary-review-v1.json), [v2 binding](lab/results/development/nsc-local-boundary-review-v2.json), [independent tests](lab/tests/test_nsc_local_boundary_independent.py) |
| Initial geometry-graded window | [Status](docs/current-result.md#initial-geometry-graded-window-1-october-2026), [note](lab/docs/nsc-geometry-graded-window.md), [module](lab/src/recursive_horizons/nsc_geometry_graded_window.py), [driver](lab/scripts/derive_nsc_geometry_graded_window.py), [record](lab/results/development/nsc-geometry-graded-window-v1.json), [tests](lab/tests/test_nsc_geometry_graded_window.py) |
| Weak spherical Cauchy residual | [Note](lab/docs/nsc-spherical-cauchy-weak.md), [module](lab/src/recursive_horizons/nsc_spherical_cauchy_weak.py), [driver](lab/scripts/derive_nsc_spherical_cauchy_weak.py), [record](lab/results/development/nsc-spherical-cauchy-weak-v1.json), [tests](lab/tests/test_nsc_spherical_cauchy_weak.py) |
| Initial lapse bound and nearby initial-state proof | [Note](lab/docs/nsc-spherical-cauchy-weak.md#initial-lapse-component), [certificate](lab/results/development/nsc-spherical-cauchy-error-v1.json), [tests](lab/tests/test_nsc_spherical_cauchy_weak_bound.py) |
| Source and UV certification methods | [Continuous source-error insertion](lab/docs/nsc-ks-source-operator-majorant.md), [family 14_1 completion](lab/docs/nsc-ks-source-operator-majorant-v7.md), [preparation correction](lab/docs/nsc-vacuum-source-correction.md), [high-energy remainder](lab/docs/nsc-vacuum-source-remainder.md), [full method map](lab/docs/active-code-map.md) |
| Field method owners | [Whole cone](lab/docs/nsc-ks-whole-cone-field-v1.md), [radius coupling](lab/docs/nsc-ks-radius-coupling-bounds.md), [endpoint contraction](lab/docs/nsc-ks-endpoint-contraction.md) |
| Current finite article and preserved papers | [Paper guide](docs/papers.md), [current PDF](paper/finite-regeneration.pdf), [source](paper/finite-regeneration/main.tex), [build/evidence](paper/finite-regeneration-manifest.json), [catalog](paper/catalog.json) |
| Finite nested-quality theorem | [Proof](docs/nsc-nested-qualities.md), [frozen evidence](paper/nested-quality-evidence/snapshot.json) |
| Finite circulation and local response | [Note](lab/docs/nsc-finite-turnover.md), [record](lab/results/development/nsc-finite-turnover-v1.json) |
| Incoherent imbalance mean circulation | [Note](lab/docs/nsc-imbalance-turnover.md), [record](lab/results/development/nsc-imbalance-turnover-v1.json) |
| Spherical feedback action | [Note](lab/docs/nsc-spherical-feedback-action.md), [module](lab/src/recursive_horizons/nsc_spherical_feedback_action.py) |
| Direct conformal ADM source | [Note](lab/docs/nsc-conformal-adm-source.md), [module](lab/src/recursive_horizons/nsc_conformal_adm_source.py) |
| Provisional spherical coupling v1 | [Note](lab/docs/nsc-spherical-coupling.md), [module](lab/src/recursive_horizons/nsc_spherical_coupling.py), [record](lab/results/development/nsc-spherical-coupling-control-v1.json) |
| Variational spherical Galerkin coupling v2 | [Note](lab/docs/nsc-spherical-galerkin-coupling.md), [module](lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py), [record](lab/results/development/nsc-spherical-coupling-control-v2.json) |
| Matched Galerkin refinement v3 | [Note](lab/docs/nsc-spherical-galerkin-coupling.md#matched-refinement-v3), [driver](lab/scripts/derive_nsc_spherical_galerkin_refinement.py), [record](lab/results/development/nsc-spherical-coupling-refinement-v3.json) |
| Galerkin refinement v5 | [Driver](lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py), [record](lab/results/development/nsc-spherical-coupling-refinement-v5.json), [payload](lab/results/development/nsc-spherical-coupling-refinement-v5.npz), [independent tests](lab/tests/test_nsc_spherical_galerkin_independent.py) |
| Finite-window embedding v1 | [Note](lab/docs/nsc-finite-window-embedding.md), [module](lab/src/recursive_horizons/nsc_finite_window_embedding.py), [record](lab/results/development/nsc-finite-window-embedding-v1.json) |
| Seam-regular embedding v2 | [Note](lab/docs/nsc-finite-window-embedding-smooth.md), [module](lab/src/recursive_horizons/nsc_finite_window_embedding_smooth.py), [record](lab/results/development/nsc-finite-window-embedding-v2.json) |
| Spherical Cauchy data | [Note](lab/docs/nsc-spherical-cauchy-data.md), [module](lab/src/recursive_horizons/nsc_spherical_cauchy_data.py), [record](lab/results/development/nsc-spherical-cauchy-data-v1.json) |
| Regional coordinate and normal-observer ledger | [Note](lab/docs/nsc-regional-energy-exchange.md), [module](lab/src/recursive_horizons/nsc_regional_energy_exchange.py), [record](lab/results/development/nsc-regional-energy-exchange-v1.json) |
| Archived nf256 probes | [Receipt](lab/archive/probes/README.md) |
| Archived external two-room comparison | Fixed Bronnikov–Fabris phantom black-universe bytes, 4 October 2026; [owning note](lab/docs/nsc-two-room-comparison.md). Separate from the magnetic RN programme and from the papers. |
| Evolving retained-region reduction | [Note](lab/docs/nsc-evolving-reduction.md), [module](lab/src/recursive_horizons/nsc_evolving_reduction.py) |
| Original source snapshot and large-data exclusions | [Import receipt](docs/lab-snapshot.json), [consolidation](docs/repository-consolidation.md) |
| Public verification and builds | [Reproducing](docs/reproducing.md), [draft](docs/local-gate-draft.md) |
| Windows calculation runtime | `.venv/validation` on Python 3.14.8. Numba 0.68 uses all CPU threads. HIP 7.2 reaches the Radeon. Section below. |

## Windows calculation runtime

This checkout's scientific interpreter is `.venv/validation`, Python 3.14.8
from `%LOCALAPPDATA%\Programs\Python\Python314`. Numba 0.68 is the
accelerated runtime and uses every logical core (`NUMBA_NUM_THREADS`,
16 on this Ryzen 9 6900HX). Existing producers still import SciPy. A new
numerical kernel uses Numba.

The AMD GPU path is the installed HIP SDK 7.2. `HIP_PATH` is
`C:\Program Files\AMD\ROCm\7.2\`. hipInfo names the device AMD Radeon
RX 6600M, architecture `gfx1032`. Python loads `amdhip64_7.dll` and
`hiprtc0702.dll` from that SDK's `bin`. The venv startup hook puts
those DLLs on the process path. Recheck CPU threads and a HIP saxpy
with `.venv/validation/Scripts/python.exe .venv/validation/Scripts/check_numba_hip.py`.

Do not look for Homebrew, `/opt/homebrew`, or `/Users/admin`. A macOS
virtual environment is not in this tree. Older notes may still name
those paths; they are the other machine's provenance.

The Windows branches are `win11/blackholes-infinity/main` and
`agent/blackholes-infinity/work`. Work happens on the second. The
remote branches `main` and `codex/work-branch` stay as published.

Live branch commands are in the
[instructions](docs/instructions.md). The two-room import commit
`d95175c25a9989aea681967d6089aa7b2bf27706` is historical. The RN
scientific start remains `5cb775a7`. Root owns Git
integration and publication. Preserve other owners' files and existing
scientific branches. Scientific Python uses `lab/`
through root `scripts/lab.py`;
root `src/`, `scripts/` and `tests/` own the curated publication chain.

Verdicts and the working diagnostic are in
[current result](docs/current-result.md).
Framing: [instructions](docs/instructions.md#scientific-framing-and-authority).

## References

Every explicit arXiv study cited in this repository, including archived notes and the preserved pre-closure paper. The line says what this repository referenced. A citation is not an NSC result. The same study in `docs/` and `lab/docs/` is listed once. Version suffixes are folded into the base identifier.

- [0711.1088](https://arxiv.org/abs/0711.1088) — Nemec, Tománek and Cuniberti, Appendix A.3: surface Green-function recursion.
- [0711.1178](https://arxiv.org/abs/0711.1178) — Dunne, functional determinants: equations 18–20 and the periodic examples.
- [0805.2909](https://arxiv.org/abs/0805.2909) — Codello, Percacci and Rahmede: Appendix A trace and Mellin identities used in source matching.
- [0806.2683](https://arxiv.org/abs/0806.2683) — Dabrowski, Garecki and Blaschke: five-dimensional conformal-curvature transformation reviewed for the warp term.
- [0901.0944](https://arxiv.org/abs/0901.0944) — Boonserm and Visser: analytic bounds on transmission probabilities, applied to omitted potential tails.
- [0904.0612](https://arxiv.org/abs/0904.0612) — Dappiaggi, Hack and Pinamonti: Definition 4.1 and Theorem 4.1, plus the covariant Dirac stress, trace anomaly and renormalization.
- [0907.5424](https://arxiv.org/abs/0907.5424) — Archived: Baumann, TASI lectures on inflation.
- [0910.4930](https://arxiv.org/abs/0910.4930) — Bronnikov and Donskoy: brane-world black universes, effective brane equations and the projected bulk-Weyl stress.
- [0911.5074](https://arxiv.org/abs/0911.5074) — Imported candidate for a Standard Model spectral action with skew-symmetric torsion, a different torsion class from the operator used here.
- [1001.2036](https://arxiv.org/abs/1001.2036) — Andrianov and Lizzi: bosonic spectral action induced by cancelling the scale anomaly of the regularized fermionic action.
- [1008.3630](https://arxiv.org/abs/1008.3630) — Iochum, Levy and Vassilevich: bulk and boundary spectral-action terms with skew torsion and chiral bag conditions, including the connection-convention conversion.
- [1008.3980](https://arxiv.org/abs/1008.3980) — Chamseddine and Connes: spectral action on a space with boundary, used for the extrinsic-curvature coefficient required by the gravitational Hamiltonian.
- [1009.5179](https://arxiv.org/abs/1009.5179) — Dappiaggi, Hack and Pinamonti: spinor boundary construction for the child state.
- [1101.1424](https://arxiv.org/abs/1101.1424) — Pfäffle and Stephan: torsionful spectral action, including Proposition 5.4, Theorem 5.7 and the Bochner and heat-kernel formulae.
- [1106.3263](https://arxiv.org/abs/1106.3263) — Andrianov, Kurkov and Lizzi: spectral action, Weyl anomaly and Higgs-dilaton potential, including the anomaly factorization and scale-versus-cutoff distinction.
- [1109.2524](https://arxiv.org/abs/1109.2524) — Chung et al.: gravitational Dirac production in inflation, cited as different from the anisotropic radius-pulse benchmark.
- [1109.6234](https://arxiv.org/abs/1109.6234) — Clemson et al.: primary analysis that different momentum-transfer completions can alter growth.
- [1111.3328](https://arxiv.org/abs/1111.3328) — Pusey, Barrett and Rudolph: on the reality of the quantum state.
- [1112.0475](https://arxiv.org/abs/1112.0475) — Archived: Olmo and Rubiera-Garcia, a regular-core comparator for Ricci-based gravity, not evidence for the proposed phase transition.
- [1210.2663](https://arxiv.org/abs/1210.2663) — Kurkov and Lizzi: Higgs-dilaton Lagrangian from spectral regularization, including the common local Weyl anomaly.
- [1301.2430](https://arxiv.org/abs/1301.2430) — Archived: Olmo and Rubiera-Garcia, a further regular-core comparator, not evidence for the proposed phase transition.
- [1305.3085](https://arxiv.org/abs/1305.3085) — Archived: Ziaie, Moniz, Ranjbar and Sepangi, Einstein–Cartan collapse of a homogeneous Weyssenhoff fluid.
- [1305.5178](https://arxiv.org/abs/1305.5178) — Archived: Bromm, formation of the first stars.
- [1307.3823](https://arxiv.org/abs/1307.3823) — Feng Rong: Proposition 2.1 on holomorphic Briot–Bouquet germs, imported for the incoming surface-crossing existence argument.
- [1308.4552](https://arxiv.org/abs/1308.4552) — Obukhov, Silenko and Teryaev: Dirac Hamiltonian, equations 2.14–2.16, used for the ADM source constraints.
- [1309.0820](https://arxiv.org/abs/1309.0820) — Adam, Naya, Sanchez-Guillen and Wereszczynski: BPS Skyrme model and nuclear binding energies.
- [1311.4182](https://arxiv.org/abs/1311.4182) — Fischmann: conformal powers of the Dirac operator, used for the free compact-carrier identification.
- [1312.2235](https://arxiv.org/abs/1312.2235) — Kurkov, Lizzi and Vassilevich: high-energy bosons do not propagate, used for the spectral-cutoff propagation wall.
- [1312.2960](https://arxiv.org/abs/1312.2960) — Adam, Naya, Sanchez-Guillen and Wereszczynski: BPS Skyrme energy linear in topological charge, including the static-energy equations.
- [1404.0602](https://arxiv.org/abs/1404.0602) — Bianchi and Smerlak: equations 6–9 for the peeling power and relative radiation entropy of a specified incoming vacuum.
- [1404.5635](https://arxiv.org/abs/1404.5635) — Archived: Bousso, Casini, Fisher and Maldacena, proof of a quantum Bousso bound.
- [1405.0397](https://arxiv.org/abs/1405.0397) — Castillo-Felisola, Corral, Kovalenko and Schmidt: equations 13 and 18–20, the higher-dimensional Einstein–Cartan contact reused on compact profiles.
- [1408.3203](https://arxiv.org/abs/1408.3203) — Archived: Engelhardt and Wall, quantum extremal surfaces; island and Page language is not an established parent-to-child process.
- [1408.6112](https://arxiv.org/abs/1408.6112) — Barvinsky: Euclidean-to-in-in rule, equations 2.25–2.29.
- [1409.2471](https://arxiv.org/abs/1409.2471) — Chamseddine, Connes and Mukhanov: equations 4–16, spectral volume quantization, not adopted as an NSC volume constraint.
- [1503.03095](https://arxiv.org/abs/1503.03095) — Adam et al.: self-gravitating BPS Skyrme neutron-star calculation, mean-field limit versus full field theory.
- [1504.00351](https://arxiv.org/abs/1504.00351) — Dai, Pajer and Schmidt: local density and curvature response to a long perturbation, with an isotropic and tidal split.
- [1506.02882](https://arxiv.org/abs/1506.02882) — Pagani and Percacci: quantum gravity with torsion and non-metricity, used as a Palatini-flow compatibility benchmark.
- [1508.03272](https://arxiv.org/abs/1508.03272) — Archived: Olmo, Rubiera-Garcia and Sanchez-Puente, a regular-core comparator, not evidence for the proposed phase transition.
- [1509.02878](https://arxiv.org/abs/1509.02878) — Tóth: a classical Dirac test field on a fixed extremal dyonic Kerr–Newman background can violate the null energy condition.
- [1510.06930](https://arxiv.org/abs/1510.06930) — Archived: Langlois and Noui, degenerate higher-derivative theories beyond Horndeski.
- [1512.00761](https://arxiv.org/abs/1512.00761) — Finster and Röken: loss of ellipticity and self-adjointness of the Dirac Hamiltonian for a class of horizon boundary problems.
- [1512.06074](https://arxiv.org/abs/1512.06074) — Lucat and Prokopec: equation 9 at xi=1, the axial contact coefficient, and their microscopic in-in source.
- [1601.04744](https://arxiv.org/abs/1601.04744) — Archived: Donnelly and Freidel, local subsystems in gauge theory and gravity.
- [1601.06809](https://arxiv.org/abs/1601.06809) — Natário, Queimada and Vicente: conditional theorem that a null-energy-condition test field cannot destroy an extremal Kerr–Newman horizon.
- [1605.01603](https://arxiv.org/abs/1605.01603) — Krishnan and Raju: equations 4–5, the Dirichlet metric variation of the Einstein and Gibbons–Hawking–York terms.
- [1606.02728](https://arxiv.org/abs/1606.02728) — Fan, Chen and Lü: criticality in Einstein–Gauss–Bonnet gravity, where the ordinary bilinear graviton is absent on the coalesced vacuum.
- [1606.05831](https://arxiv.org/abs/1606.05831) — Archived: Kobayashi, generic instabilities of nonsingular cosmologies in Horndeski theory.
- [1607.07791](https://arxiv.org/abs/1607.07791) — Imported candidate: an Einstein–Cartan black universe with regular torsion and a positive-kinetic scalar, not yet the same spectral torsion operator.
- [1609.01731](https://arxiv.org/abs/1609.01731) — Bondi–Sachs formulation used for the distinction between asymptotic Bondi mass and energy carried by outgoing radiation.
- [1612.07281](https://arxiv.org/abs/1612.07281) — Archived: Bonga, Gupt and Yokomizo, tensor perturbations during inflation in a spatially closed universe, including the three-sphere tensor eigenvalue.
- [1704.03143](https://arxiv.org/abs/1704.03143) — Contreras et al.: physical models for dipolar asymmetry tested with CMB polarization.
- [1704.04807](https://arxiv.org/abs/1704.04807) — Archived: exact self-dual Skyrmion identities, equations 11–16, embedded as a scale sector.
- [1708.02951](https://arxiv.org/abs/1708.02951) — Archived: Langlois et al., an extended-mimetic limiting-curvature comparator, not the healthy default.
- [1710.05834](https://arxiv.org/abs/1710.05834) — LIGO, Virgo, Fermi-GBM and INTEGRAL: GW170817 and GRB 170817A, used as the published multimessenger timing source, not a raw reanalysis.
- [1710.10155](https://arxiv.org/abs/1710.10155) — Archived: Papallo, hyperbolicity of the most general Horndeski theory.
- [1711.02080](https://arxiv.org/abs/1711.02080) — Archived: Silva et al., curvature-triggered spontaneous scalarization, cited as a comparator that does not replace an independent collapsing source.
- [1712.03876](https://arxiv.org/abs/1712.03876) — Archived: Ben Achour et al., an extended-mimetic limiting-curvature comparator, not the healthy default.
- [1803.08827](https://arxiv.org/abs/1803.08827) — Briscese and Modesto: introductory assumptions of their Euclidean-continuation unitarity result, which is not imported for the decreasing kinetic profile.
- [1805.05107](https://arxiv.org/abs/1805.05107) — Barbero, Ferreiro, Navarro-Salas and Villaseñor: section IV, the Dirac adiabatic expansion as a unitary spinor expansion.
- [1807.04726](https://arxiv.org/abs/1807.04726) — Maldacena, Milekhin and Popov: traversable wormholes, including Appendix A, section 2.3 and the charged massless-Dirac vacuum throat.
- [1807.06209](https://arxiv.org/abs/1807.06209) — Archived: Planck 2018 results VI, cosmological parameters.
- [1807.06211](https://arxiv.org/abs/1807.06211) — Planck 2018 results X: constraints on inflation, including the polarization test of a direction-dependent tilt.
- [1807.08361](https://arxiv.org/abs/1807.08361) — Archived: Mironov, Rubakov and Volkova, a bounce beyond Horndeski with general-relativity asymptotics.
- [1812.05590](https://arxiv.org/abs/1812.05590) — Archived: Silva et al., nonlinear coupling controls quenching and radial stability of scalarized compact objects.
- [1902.02222](https://arxiv.org/abs/1902.02222) — Archived: Cabral, Lobo and Rubiera-Garcia, Einstein–Cartan–Dirac identities, especially equations 8–11, 18, 20, 23 and 83.
- [1903.07543](https://arxiv.org/abs/1903.07543) — Archived: Ripley and Pretorius, scalar-Gauss–Bonnet collapse can change from hyperbolic to elliptic character.
- [1905.08255](https://arxiv.org/abs/1905.08255) — Archived: Penington, entanglement-wedge reconstruction and the information paradox; cited as motivation, not an established parent-to-child process.
- [1906.02552](https://arxiv.org/abs/1906.02552) — Planck 2018 results VII: isotropy and statistics of the CMB.
- [1908.08050](https://arxiv.org/abs/1908.08050) — Ruhdorfer, Serra and Weiler: section III.1, removing redundant curvature operators by metric redefinition.
- [1910.01648](https://arxiv.org/abs/1910.01648) — Montero, Van Riet and Venken: section 2, charged de Sitter geometry used before importing a near-extremal seed.
- [2003.08398](https://arxiv.org/abs/2003.08398) — Archived: Kovács and Reall, well-posed formulation of Lovelock and Horndeski theories.
- [2008.02136](https://arxiv.org/abs/2008.02136) — Popławski: equations 1, 29 and 32–34, a spin-fluid Einstein–Cartan collapse and bounce, not the present vacuum control.
- [2008.10995](https://arxiv.org/abs/2008.10995) — Gérard, Häfner and Wrochna: Theorem 1.1, the Unruh state for massless fermions on Kerr and its Hadamard property, on its stated domain.
- [2011.03547](https://arxiv.org/abs/2011.03547) — Archived: East and Ripley, evolution of Einstein-scalar-Gauss–Bonnet gravity in a modified harmonic formulation.
- [2012.06401](https://arxiv.org/abs/2012.06401) — Bochniak and Sitarz: stability of doubled FLRW solutions, including equation II.33.
- [2108.11630](https://arxiv.org/abs/2108.11630) — Gérard and Stoskopf: Lemma 6.4 and Proposition 6.8, Sobolev propagation and the time-dependent approximate projector.
- [2112.05819](https://arxiv.org/abs/2112.05819) — Archived: Murphy et al., ESPRESSO limit on variation of the fine-structure constant.
- [2201.03839](https://arxiv.org/abs/2201.03839) — Bochniak and Sitarz: section 2, spectral interaction between universes, the prior geometric setting of a two-sheet field connection.
- [2206.00014](https://arxiv.org/abs/2206.00014) — Archived: fixing the dynamical evolution in scalar-Gauss–Bonnet gravity, cited at the action gate.
- [2211.08477](https://arxiv.org/abs/2211.08477) — Archived: Hegade, Ripley and Yunes, where and why Einstein-scalar-Gauss–Bonnet theory breaks down.
- [2306.01695](https://arxiv.org/abs/2306.01695) — Archived: Thaalba et al., scalar-Gauss–Bonnet gravity with a Ricci coupling, equation 1.
- [2306.03780](https://arxiv.org/abs/2306.03780) — Barvinsky and Wachowski: equations 4.17, 4.32 and 4.34, the conformal-stress relation beyond conformally flat metrics.
- [2306.14880](https://arxiv.org/abs/2306.14880) — Gimeno-Amo et al.: hemispherical power asymmetry in intensity and polarization for Planck PR4.
- [2307.12190](https://arxiv.org/abs/2307.12190) — Archived: Popławski, a universe in a black hole, used as a source reduction rather than a measured speed of light.
- [2307.12677](https://arxiv.org/abs/2307.12677) — Archived: Ranocha and Giesselmann, a residual becomes a solution-error bound only through an applicable stability estimate.
- [2309.00873](https://arxiv.org/abs/2309.00873) — Nakas, Pappas and Stuchlík: a single areal-radius chart requires single-valued metric functions, so it is not applied directly to the two-sheet branches.
- [2309.10034](https://arxiv.org/abs/2309.10034) — Tristram et al.: cosmological parameters from the final Planck PR4 release.
- [2403.08705](https://arxiv.org/abs/2403.08705) — Archived: Lara et al., a fixing-the-equations approach whose equation 6 is recorded as not used for coefficient authority.
- [2405.14056](https://arxiv.org/abs/2405.14056) — Briscese et al.: sections 3.2 and 4.2, the nonlocal diagrammatic inverse is not the time-ordered correlator of physical states.
- [2406.19108](https://arxiv.org/abs/2406.19108) — Agüera y Arcas and collaborators, Computational Life: self-replicating programs in random computational systems, cited as related motivation only.
- [2407.06790](https://arxiv.org/abs/2407.06790) — Archived comparison source for a quoted NICER chiral-effective-field-theory maximum neutron-star mass.
- [2409.11398](https://arxiv.org/abs/2409.11398) — Archived: Thaalba et al., horizon-penetrating Painlevé–Gullstrand comparator and characteristic diagnostic for gravity with a Ricci coupling.
- [2410.16264](https://arxiv.org/abs/2410.16264) — Archived: gauge-invariant hyperbolicity study of spherical evolution in scalar-Gauss–Bonnet gravity.
- [2503.14738](https://arxiv.org/abs/2503.14738) — DESI Collaboration: DR2 baryon-acoustic-oscillation measurements and cosmological constraints, not a DR2 full-shape analysis.
- [2503.14744](https://arxiv.org/abs/2503.14744) — DESI Collaboration: neutrino constraints from DR2 BAO and DR1 full shape, including the stated absence of a joint 2025 full-shape plus DR2 analysis.
- [2503.19441](https://arxiv.org/abs/2503.19441) — Wright et al.: KiDS-Legacy cosmological constraints from cosmic shear.
- [2503.24263](https://arxiv.org/abs/2503.24263) — Archived: Meissner and Penrose, the physics of conformal cyclic cosmology.
- [2504.05597](https://arxiv.org/abs/2504.05597) — Gimeno-Amo et al.: cosmological-parameter variations across the sky in Planck PR4.
- [2510.23823](https://arxiv.org/abs/2510.23823) — H0DN Collaboration: the local distance network.
- [2601.14559](https://arxiv.org/abs/2601.14559) — DES Collaboration: Year 6 cosmological constraints from galaxy clustering and weak lensing.
- [2602.18761](https://arxiv.org/abs/2602.18761) — Forero-Sánchez et al.: joint DESI DR1 full-shape and DR2 BAO constraints, a 2026 preprint rather than the 2025 DESI DR2 paper.
- [2608.23861](https://arxiv.org/abs/2608.23861) — Heidari et al.: equation 6, an exterior formula matched at a stated parameter point, not a global continuation across the throat.
- [astro-ph/0305457](https://arxiv.org/abs/astro-ph/0305457) — Archived: Magueijo, review of new varying-speed-of-light theories.
- [astro-ph/0607657](https://arxiv.org/abs/astro-ph/0607657) — Bolton, Rappaport and Burles: galactic-scale constraint on the post-Newtonian parameter gamma from lensing and dynamical masses.
- [astro-ph/0608407](https://arxiv.org/abs/astro-ph/0608407) — Archived: Clowe et al., the Bullet Cluster as empirical evidence for dark matter.
- [astro-ph/9805201](https://arxiv.org/abs/astro-ph/9805201) — Archived: Riess et al., supernova evidence for an accelerating universe.
- [astro-ph/9811018](https://arxiv.org/abs/astro-ph/9811018) — Archived: Albrecht and Magueijo, a time-varying speed of light as a proposed cosmological solution.
- [astro-ph/9812133](https://arxiv.org/abs/astro-ph/9812133) — Archived: Perlmutter et al., supernova measurements of Omega and Lambda.
- [cond-mat/0209642](https://arxiv.org/abs/cond-mat/0209642) — Klich: the Fock-space trace identity behind Levitov's formula, implemented for the geometry vertex.
- [cond-mat/9404027](https://arxiv.org/abs/cond-mat/9404027) — Jauho, Wingreen and Meir: state-dependent quantum transport.
- [gr-qc/0001069](https://arxiv.org/abs/gr-qc/0001069) — Archived: Martel and Poisson, regular coordinate systems for Schwarzschild and other spherical spacetimes.
- [gr-qc/0007036](https://arxiv.org/abs/gr-qc/0007036) — Archived: Magueijo, covariant and locally Lorentz-invariant varying-speed-of-light theories.
- [gr-qc/0207080](https://arxiv.org/abs/gr-qc/0207080) — Archived: Ashtekar and Krishnan, dynamical horizons, energy, angular momentum, fluxes and balance laws.
- [gr-qc/0212028](https://arxiv.org/abs/gr-qc/0212028) — Hollands: the curved-spacetime PCT theorem, cited for an observable invariant.
- [gr-qc/0304007](https://arxiv.org/abs/gr-qc/0304007) — Archived: Ventrella and Choptuik, critical phenomena in the Einstein–massless-Dirac system.
- [gr-qc/0404074](https://arxiv.org/abs/gr-qc/0404074) — Hollands and Wald: local covariant renormalization, used for relational vacuum normalization.
- [gr-qc/0505065](https://arxiv.org/abs/gr-qc/0505065) — Archived: Barceló, Liberati and Visser, analogue gravity.
- [gr-qc/0506037](https://arxiv.org/abs/gr-qc/0506037) — Archived: hyperbolicity of second-order-in-space systems of evolution equations, section III.
- [gr-qc/0506126](https://arxiv.org/abs/gr-qc/0506126) — Archived: Hayward, formation and evaporation of nonsingular black holes.
- [gr-qc/0511109](https://arxiv.org/abs/gr-qc/0511109) — Bronnikov and Fabris: regular phantom black holes, the fixed external two-room comparison, not a replacement for their source domain.
- [gr-qc/0605031](https://arxiv.org/abs/gr-qc/0605031) — Dolan, Doran and Lasenby: section IV.2 and equation 72, the Schwarzschild Painlevé–Gullstrand connection and neutral flux control.
- [gr-qc/0607039](https://arxiv.org/abs/gr-qc/0607039) — Archived: Ashtekar, Pawlowski and Singh, improved dynamics for the quantum nature of the big bang.
- [gr-qc/0611022](https://arxiv.org/abs/gr-qc/0611022) — Bronnikov, Dehnen and Melnikov: sections 2–3, regular black holes and black universes; their phantom source is not the Dirac source used here.
- [gr-qc/9211002](https://arxiv.org/abs/gr-qc/9211002) — Parker and Simon: perturbative order reduction of local curvature-squared terms in the Einstein equation.
- [gr-qc/9211017](https://arxiv.org/abs/gr-qc/9211017) — Archived: Bekenstein, relation between physical and gravitational geometry.
- [gr-qc/9403008](https://arxiv.org/abs/gr-qc/9403008) — Archived: Garay, quantum gravity and minimum length.
- [gr-qc/9404011](https://arxiv.org/abs/gr-qc/9404011) — Smolin: the fate of black-hole singularities and Standard Model parameters, cited as a prior natural-selection proposal.
- [gr-qc/9406053](https://arxiv.org/abs/gr-qc/9406053) — Archived: Borde, topology change in classical general relativity.
- [gr-qc/9504004](https://arxiv.org/abs/gr-qc/9504004) — Archived: Jacobson, thermodynamics of spacetime and the Einstein equation of state.
- [gr-qc/9505009](https://arxiv.org/abs/gr-qc/9505009) — Camporesi and Higuchi: equations 3.51–3.54, the Dirac spectrum on spheres and real hyperbolic spaces.
- [gr-qc/9612057](https://arxiv.org/abs/gr-qc/9612057) — Archived: Borde, regular black holes and topology change.
- [gr-qc/9710089](https://arxiv.org/abs/gr-qc/9710089) — Archived: Hayward, the unified first law of black-hole dynamics.
- [gr-qc/9801079](https://arxiv.org/abs/gr-qc/9801079) — Archived: Finster, Smoller and Yau, particlelike solutions of the Einstein–Dirac equations.
- [gr-qc/9810032](https://arxiv.org/abs/gr-qc/9810032) — Fabris, Pelinson and Shapiro: section 2, the anomaly-induced-action construction.
- [gr-qc/9901069](https://arxiv.org/abs/gr-qc/9901069) — Hollands: adiabatic Hadamard states for Dirac fields, used to propagate a reference through a smooth auxiliary metric.
- [gr-qc/9904021](https://arxiv.org/abs/gr-qc/9904021) — Martín and Verdaguer: stochastic semiclassical gravity and the influence-functional approach to gravitational response.
- [gr-qc/9906015](https://arxiv.org/abs/gr-qc/9906015) — Buchert: averaged expansion, curvature and backreaction for irrotational dust, which does not supply this repository's quantum state.
- [hep-ph/0105021](https://arxiv.org/abs/hep-ph/0105021) — Ponton and Poppitz: sections 2.1–2.2, the Einstein-frame radius and distinct local counterterms; their fitted coefficients are not imported.
- [hep-ph/9912408](https://arxiv.org/abs/hep-ph/9912408) — Grossman and Neubert: section 2, five-dimensional spinor reduction and the chiral endpoint choice; their Randall–Sundrum phenomenology is not imported.
- [hep-th/0011232](https://arxiv.org/abs/hep-th/0011232) — Alimohammadi and Mohseni Sadjadi: section 2, curved-space bosonization of the massive Schwinger model.
- [hep-th/0103138](https://arxiv.org/abs/hep-th/0103138) — Flachi, Moss and Toms: equations 6–14, massless interval determinants, conformal reduction and the Wilson-line mechanism.
- [hep-th/0103186](https://arxiv.org/abs/hep-th/0103186) — Dasgupta and Loll: sections 2–4, a conditional conformal-factor cancellation in a Lorentzian proper-time gauge, not the heat-kernel regulator used here.
- [hep-th/0110007](https://arxiv.org/abs/hep-th/0110007) — Archived: Spradlin, Strominger and Volovich, Les Houches lectures on de Sitter space.
- [hep-th/0204253](https://arxiv.org/abs/hep-th/0204253) — Grumiller, Kummer and Vassilevich: dilaton gravity in two dimensions, reused for spherical reduction and first-order canonical gravity.
- [hep-th/0206054](https://arxiv.org/abs/hep-th/0206054) — Borokhov, Kapustin and Wu: section 4.1 and the appendix, the magnetic angular spectrum and degeneracies.
- [hep-th/0207060](https://arxiv.org/abs/hep-th/0207060) — Kiritsis et al.: bulk-boundary energy exchange that modifies cosmological evolution, cited as prior art still requiring the same NSC quantum state.
- [hep-th/0306138](https://arxiv.org/abs/hep-th/0306138) — Vassilevich: heat-kernel expansion, including sections 3.3 and 7.1 and equations 3.34, 4.26–4.28 and 5.30–5.32.
- [hep-th/0602146](https://arxiv.org/abs/hep-th/0602146) — Iso, Umetsu and Wilczek: equations 25–28, Hawking radiation from charged black holes and the free fermionic flux normalization.
- [hep-th/9502109](https://arxiv.org/abs/hep-th/9502109) — Mottola: section III, a configuration-space measure built from metric deformations and divided by diffeomorphism orbits.
- [hep-th/9606001](https://arxiv.org/abs/hep-th/9606001) — Chamseddine and Connes: the spectral action principle, imported as the reduction of gravity and gauge fields from one Dirac operator.
- [hep-th/9807038](https://arxiv.org/abs/hep-th/9807038) — Gusev and Zelnikov: equations 83–89, fermionic thermal images and the canonical free energy on a finite proper-time endpoint.
- [hep-th/9905177](https://arxiv.org/abs/hep-th/9905177) — Archived: Bousso, a covariant entropy conjecture.
- [hep-th/9911055](https://arxiv.org/abs/hep-th/9911055) — Garriga and Tanaka: gravity in the Randall–Sundrum brane world, including the zero mode, continuum and two-wall linearized response.
- [math-ph/0201055](https://arxiv.org/abs/math-ph/0201055) — Panati, Spohn and Teufel: sections 2 and 4.4, the space-adiabatic projector and almost-invariant subspace used for the reference band.
- [math-ph/0403050](https://arxiv.org/abs/math-ph/0403050) — Kirsten and McKane: equation 49 and the general boundary-condition construction for Sturm–Liouville determinants.
- [math/0304347](https://arxiv.org/abs/math/0304347) — Lee: Theorems 1.1 and 1.5, the Burghelea–Friedlander–Kappeler gluing formula for the zeta determinant of a Dirac Laplacian.
- [physics/0503066](https://arxiv.org/abs/physics/0503066) — Archived: Noether, invariant variation problems, in the Tavel translation.
- [quant-ph/0102125](https://arxiv.org/abs/quant-ph/0102125) — Archived: Simon, Bužek and Gisin, the no-signaling condition and quantum dynamics.
- [quant-ph/9609002](https://arxiv.org/abs/quant-ph/9609002) — Archived: Rovelli, relational quantum mechanics.
