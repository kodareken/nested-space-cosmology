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
| Present regeneration loop | [Status](docs/current-result.md#present-spherical-loop), [code owner](lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py), [saved v5 diagnostic](lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py) |
| Existing mathematics, assumptions and gaps | [Lab claim ledger](lab/docs/claim-ledger.md), [code map](lab/docs/active-code-map.md) |
| Source and UV certification methods | [Continuous source-error insertion](lab/docs/nsc-ks-source-operator-majorant.md), [family 14_1 completion](lab/docs/nsc-ks-source-operator-majorant-v7.md), [preparation correction](lab/docs/nsc-vacuum-source-correction.md), [high-energy remainder](lab/docs/nsc-vacuum-source-remainder.md), [full method map](lab/docs/active-code-map.md) |
| Field method owners | [Whole cone](lab/docs/nsc-ks-whole-cone-field-v1.md), [radius coupling](lab/docs/nsc-ks-radius-coupling-bounds.md), [endpoint contraction](lab/docs/nsc-ks-endpoint-contraction.md) |
| Both papers and public status | [Paper guide](docs/papers.md), [current result](docs/current-result.md), [catalog](paper/catalog.json) |
| Finite nested-quality theorem | [Proof](docs/nsc-nested-qualities.md), [frozen evidence](paper/nested-quality-evidence/snapshot.json) |
| Finite circulation and local response | [Note](lab/docs/nsc-finite-turnover.md), [record](lab/results/development/nsc-finite-turnover-v1.json) |
| Incoherent imbalance mean circulation | [Note](lab/docs/nsc-imbalance-turnover.md), [record](lab/results/development/nsc-imbalance-turnover-v1.json) |
| Spherical feedback action | [Note](lab/docs/nsc-spherical-feedback-action.md), [module](lab/src/recursive_horizons/nsc_spherical_feedback_action.py) |
| Direct conformal ADM source | [Note](lab/docs/nsc-conformal-adm-source.md), [module](lab/src/recursive_horizons/nsc_conformal_adm_source.py) |
| Provisional spherical coupling v1 | [Note](lab/docs/nsc-spherical-coupling.md), [module](lab/src/recursive_horizons/nsc_spherical_coupling.py), [record](lab/results/development/nsc-spherical-coupling-control-v1.json) |
| Variational spherical Galerkin coupling v2 | [Note](lab/docs/nsc-spherical-galerkin-coupling.md), [module](lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py), [record](lab/results/development/nsc-spherical-coupling-control-v2.json) |
| Matched Galerkin refinement v3 | [Note](lab/docs/nsc-spherical-galerkin-coupling.md#matched-refinement-v3), [driver](lab/scripts/derive_nsc_spherical_galerkin_refinement.py), [record](lab/results/development/nsc-spherical-coupling-refinement-v3.json) |
| Galerkin refinement v5 | [Driver](lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py), [record](lab/results/development/nsc-spherical-coupling-refinement-v5.json), [payload](lab/results/development/nsc-spherical-coupling-refinement-v5.npz), [independent tests](lab/tests/test_nsc_spherical_galerkin_independent.py) |
| Spherical feedback episode, $T=0.05$ | [Driver](lab/scripts/derive_nsc_spherical_feedback_episode.py), [record](lab/results/development/nsc-spherical-feedback-episode-v1.json), [payload](lab/results/development/nsc-spherical-feedback-episode-v1.npz), [tests](lab/tests/test_nsc_spherical_feedback_episode.py), [independent tests](lab/tests/test_nsc_spherical_feedback_episode_independent.py) |
| Weak spherical Cauchy residual | [Note](lab/docs/nsc-spherical-cauchy-weak.md), [module](lab/src/recursive_horizons/nsc_spherical_cauchy_weak.py), [driver](lab/scripts/derive_nsc_spherical_cauchy_weak.py), [record](lab/results/development/nsc-spherical-cauchy-weak-v1.json) |
| Finite-window embedding v1 | [Note](lab/docs/nsc-finite-window-embedding.md), [module](lab/src/recursive_horizons/nsc_finite_window_embedding.py), [record](lab/results/development/nsc-finite-window-embedding-v1.json) |
| Seam-regular embedding v2 | [Note](lab/docs/nsc-finite-window-embedding-smooth.md), [module](lab/src/recursive_horizons/nsc_finite_window_embedding_smooth.py), [record](lab/results/development/nsc-finite-window-embedding-v2.json) |
| Spherical Cauchy data | [Note](lab/docs/nsc-spherical-cauchy-data.md), [module](lab/src/recursive_horizons/nsc_spherical_cauchy_data.py), [record](lab/results/development/nsc-spherical-cauchy-data-v1.json) |
| Regional coordinate and normal-observer ledger | [Note](lab/docs/nsc-regional-energy-exchange.md), [module](lab/src/recursive_horizons/nsc_regional_energy_exchange.py), [record](lab/results/development/nsc-regional-energy-exchange-v1.json) |
| Archived nf256 probes | [Receipt](lab/archive/probes/README.md) |
| Evolving retained-region reduction | [Note](lab/docs/nsc-evolving-reduction.md), [module](lab/src/recursive_horizons/nsc_evolving_reduction.py), [streamed tests](lab/tests/test_nsc_evolving_reduction_streamed.py) |
| Original source snapshot and large-data exclusions | [Import receipt](docs/lab-snapshot.json), [consolidation](docs/repository-consolidation.md) |
| Public verification and builds | [Reproducing](docs/reproducing.md), [draft](docs/local-gate-draft.md) |

Work on `codex/work-branch`. Integrate verified work into `main` and push both
to `kodareken/nested-space-cosmology`. Scientific Python uses `lab/` through
root `scripts/lab.py`.

Starting checkpoint `9a9090a`. Essential scope: the nested mechanism and the
structure common across regional gradients. The inheritance relation is
reused. The finite realization is the spherical Galerkin loop. Its code
owner is `lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py`.
The v5 replay remains the self-contained $T=0.005$ seed. The saved
$T=0.05$ measurement is
`lab/scripts/derive_nsc_spherical_feedback_episode.py`: four runs reached
the requested time, the 62 physical refinement rows meet one percent, and
a constraint-consistent continuum solution is not validated. The weak
Cauchy helper and the streamed reducer are separate method owners. The
source-fixed incoming gate stays OPEN, its campaign paused, and it is not
a prerequisite. The next release version follows a demonstrated result.
Frozen v0.26.0 stays the foundation.
Numbers and owners:
[current result](docs/current-result.md#present-spherical-loop).
Framing: [instructions](docs/instructions.md#scientific-framing-and-authority).
