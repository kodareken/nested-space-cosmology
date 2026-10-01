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

Work on `codex/work-branch`; integrate verified work into `main` and push both
to `kodareken/nested-space-cosmology`. Scientific work uses `lab/` owners, not
the older root publication code.

Current question: Douglas's proposed mechanism in the
[claim ledger](lab/docs/claim-ledger.md). Regeneration stays first. The
incoherent-imbalance calculation is the
[note](lab/docs/nsc-imbalance-turnover.md) and
[record](lab/results/development/nsc-imbalance-turnover-v1.json), not a new
roadmap. The [spherical feedback action](lab/docs/nsc-spherical-feedback-action.md)
and the [direct conformal ADM source](lab/docs/nsc-conformal-adm-source.md)
are checked primitives. The
[provisional spherical coupling v1](lab/docs/nsc-spherical-coupling.md) is
`PROVISIONAL_CONSTRAINT_DRIFT`, not a self-consistent solution or a renewal.
Initial constraints pass and the work balances. The N=64 Hamilton constraint
reaches 36.98886722265645 by T=0.005 and is unchanged by timestep halving.
Spatial diagnostics reduce it to 2.330762882869135 at N=128 and
0.575679312602035 at N=256, still above 1e-3. The odd-lobe phase fix is a
successor and does not overwrite v1. The corrected source puts that phase
on the minus column. The
[variational Galerkin coupling v2](lab/docs/nsc-spherical-galerkin-coupling.md)
is `FAIL_HELD_OUT_CONSTRAINT`; a projected residual is not a pass. The
[matched refinement v3](lab/docs/nsc-spherical-galerkin-coupling.md#matched-refinement-v3)
is `FAIL_INITIAL_SOLVE`: internal Newton stopped at rho=0.5 with residual
1.287606e-10 and wrote no full-source result. The
[evolving reduction](lab/docs/nsc-evolving-reduction.md) is prescribed-control
convergence with no accepted coupled trajectory. Passing tests do not
validate these failed trajectories. The coupled state and geometry loop
already executes. Fixed `L` and `beta` are gauge; `Q`, `r`, `chi` and the
momenta evolve. The
[v5 short response](lab/results/development/nsc-spherical-coupling-refinement-v5.json)
was independently replayed. At nf=512 and T=0.005 it remains a diagnostic:
Hamilton about 9.98e-6, momentum about 1.38e-6, proper radial velocity about
0.0079505. It is not an initial 1e-8 pass and not a renewal. The original
H and B fractions are one witness, not a physical requirement for every
later realization.
[Embedding v1](lab/docs/nsc-finite-window-embedding.md) and the
[seam-regular v2](lab/docs/nsc-finite-window-embedding-smooth.md) compression
are an optional bridge, not a further fitting prerequisite. The
[Cauchy data](lab/docs/nsc-spherical-cauchy-data.md) leave the mean and
held-out current explicit. The
[regional ledger](lab/docs/nsc-regional-energy-exchange.md) uses
dx-normalized proper pressure and normal energy; `F_Q Qdot` is coordinate
metric work. Historical vacuum branch B is not a prerequisite. Next
diagnostic: the same coupled loop to T=0.05 with geometry-derived B, effect
and spacetime convergence, and full, projected, and weak diagnostics. A
small weak residual alone does not prove accuracy. Strong residuals stay
reported. The old initial 1e-8 line is not a blind veto. Renewal is not
demonstrated. The resolved active action remains `Gamma_one`: the canonical
Gaussian plus the same-spectrum local induced term, counted once. The numerical
incoming-gate campaign is paused; the gate remains OPEN.
[Scope](docs/instructions.md#scientific-framing-and-authority).
