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
| Task authority and scientific status | Latest user instruction, which overrides a stale plan; [scope](docs/instructions.md#scientific-framing-and-authority), [status](docs/current-result.md), [code map](lab/docs/active-code-map.md) |
| Existing mathematics, assumptions, and gaps | [Lab claim ledger](lab/docs/claim-ledger.md), [code map](lab/docs/active-code-map.md) |
| Active discovery programme | Accepted conversation plan; [current execution status](docs/current-result.md#active-discovery-programme-2-october-2026), [working rules](docs/instructions.md#present-regeneration-programme) |
| Discovery implementation | [Episode runner](lab/docs/nsc-discovery-episode.md), [observer and atlas](lab/docs/nsc-discovery-observables.md), [coupled response](lab/docs/nsc-discovery-response.md), [FFT adapter](lab/src/recursive_horizons/nsc_discovery_backend.py) |
| Discovery prediction and source-selected regions | [Matched-clock prediction](lab/docs/nsc-discovery-prediction.md), [gradient cuts and moving ledger](lab/docs/nsc-discovery-regions.md), [initial source scale](lab/docs/nsc-discovery-scale.md) |
| Discovery physical assessment and coupled reduction | [Batch assessment](lab/docs/nsc-discovery-episode-assessment.md), [four-metric tides](lab/docs/nsc-discovery-tidal.md), [causal coupled memory](lab/docs/nsc-discovery-coupled-memory.md) |
| Discovery regime family and energy accounting | [Source family](lab/docs/nsc-discovery-family.md), [disjoint RK-stage balance](lab/docs/nsc-discovery-balance.md) |
| Discovery holding balance and geometric controls | [Static spectral source](lab/docs/nsc-discovery-stationary.md), [same-action source-free solution](lab/docs/nsc-discovery-vacuum-control.md), [whole-state translation](lab/docs/nsc-discovery-translation.md) |
| Deeper inherited hand-off and future prediction | [Grandchild response](lab/docs/nsc-discovery-grandchild.md), [producer](lab/src/recursive_horizons/nsc_discovery_grandchild.py), [driver](lab/scripts/derive_nsc_discovery_grandchild.py) |
| Physical source-width response | [Width prediction](lab/docs/nsc-discovery-width-response.md), [producer](lab/src/recursive_horizons/nsc_discovery_width_response.py), [driver](lab/scripts/derive_nsc_discovery_width_response.py) |
| Ambient extent and source-consistent boundary controls | [Extent intervention](lab/docs/nsc-discovery-extent.md), [producer](lab/src/recursive_horizons/nsc_discovery_extent.py), [driver](lab/scripts/derive_nsc_discovery_extent.py) |
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
| Evolving retained-region reduction | [Note](lab/docs/nsc-evolving-reduction.md), [module](lab/src/recursive_horizons/nsc_evolving_reduction.py) |
| Original source snapshot and large-data exclusions | [Import receipt](docs/lab-snapshot.json), [consolidation](docs/repository-consolidation.md) |
| Public verification and builds | [Reproducing](docs/reproducing.md), [draft](docs/local-gate-draft.md) |

Work on `codex/work-branch` and integrate verified changes into `main`;
preserve existing scientific branches and publication history. Follow the
active request for integration and publication. Scientific Python uses `lab/`
through root `scripts/lab.py`;
root `src/`, `scripts/` and `tests/` own the curated publication chain.

Verdicts and the working diagnostic are in
[current result](docs/current-result.md).
Framing: [instructions](docs/instructions.md#scientific-framing-and-authority).
