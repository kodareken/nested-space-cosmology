# Nested-Space Cosmology — public working orientation

This file orients contributors to the one-equation thesis and its evidence.
Preserve the central thesis, scientific attribution and recorded result scope.

**Presentation:** lead with equations, results and their physical connections.
The README owns the general research-status statement. Elsewhere retain the
specific assumptions, domains and unresolved terms that matter to a result,
and remove repeated generic disclaimers. Use direct language for established
identities and computed quantities. Keep authenticated notes, numerical records
and source bytes intact when revising the public narrative.

The public front door uses two layers. Begin with the one-inherited-spectrum
picture in language a general reader can visualize, then give the exact Dirac,
Schur, state, stress, and recursion equations that make each correspondence
precise. Define specialist terms after the intuitive relation they formalize.
Use positive/negative frequency and charge conjugation for antimatter; use Hz
for physical frequency and reserve dB for logarithmic amplitude ratios.

For GitHub-facing pages, use dollar-delimited inline math and display math
blocks. Write the trace as `\mathrm{Tr}` and use `\lvert`/`\rvert` for
absolute values in tables, whose literal bars otherwise split cells.
Check GitHub's Markdown rendering when changing math presentation.

**Execution priority: reuse first.** Consult the [reuse map](docs/nsc-prior-art-reuse.md)
and current result before allocating scientific compute. Historical derivations
are reference material, not a queue to repeat.

Before a new derivation or run, identify the closest primary formulation and
existing implementation, what can be reused, and the exact missing NSC
connection. Compute only what can resolve that gap. Small convention/domain
checks may be needed for integration; rebuilding an established theory is not.

Keep one compact decision entry in the existing plan or delegation brief:
reused equation/artifact, missing connection, decision the result can change,
smallest check and stopping condition. A rerun additionally names the changed
dependency, failed check or specific unresolved concern. Use the
[shared procedure](docs/nsc-prior-art-reuse.md#before-allocating-scientific-compute).
Once the source and its applicability are established, reuse the reference
and conclusion; do not spend context repeatedly reconstructing settled work.

**Stop a repetition loop:** if further work only reconfirms a known result and
leaves the named physical gap unchanged, stop expanding it and return to the
missing connection or input. More tests, records or precision do not themselves
advance closure. Reuse completed verification unless relevant inputs/code
changed, a failure occurred, or a new concern makes another check informative.
Apply this rule to all six targets and delegation briefs; workers must not
duplicate baseline derivations or broad verification. Keep settled work in
references and compact notes rather than reconstructing it each continuation.
For documentation or navigation changes, check the changed links and commands.
Retain completed scientific verification; do not rerun generators or rebuild
an unchanged PDF merely because work resumes in another chat.

**Douglas Ek** is the accountable author. **ChatGPT** and **OpenAI Codex** assisted with formulation, code, computation, drafting, and audit. Treat every AI contribution as candidate work until it matches explicit mathematics, committed compact JSON, cited sources, and stated nonclaims.

## Central thesis

Black holes, universes, singularities, infinity, and local constants can be understood as parts of one inherited recursive architecture: **spaces forming inside other spaces through collapse, localization, and renewed expansion**.

```text
Parent space
└── Local collapse / black hole
    └── Child space
        └── ...
```

- Eternity is an unending process of finite rooms, not one infinite container.
- Singularity and infinity are two directions of the same local-chart wall.
- Local constants are formation parameters of one room.
- Unresolved rooms act locally through the same kernel that makes a boundary black.
- Parent collapse and child expansion are two sides of one finite transition.

The whole proposal is one recursive spectral equation:

$$
Z_{\Lambda}^{\mathrm{one}}(\mathbb D_\Theta)
=
\int \mathcal D\varphi\,
Z_\Lambda\!\left(e^{-\varphi/2}\mathbb D_\Theta e^{-\varphi/2}\right),
\qquad
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The current concrete regulator is the heat-spectral action plus its fermionic
quadratic form. The configuration changes from parent to child. The law,
\(\Theta\), energy meaning, and boundary kernel do not.

## Constructive rules

1. Build the strongest coherent realization of the one operator before judging the wider architecture.
2. Use established physics as a launch surface, not as a prohibition against changing interpretations outside their tested domains.
3. Treat contradictions, divergences, and singularities as instructions to inspect the local chart, scale, or missing parent/child map.
4. Distinguish exact equality from structural resemblance.
5. Do not convert failure of one truncation, coordinate chart, or numerical method into rejection of Nested-Space Cosmology.
6. Do not spend scientific compute reproving imported interference, solitons, nuclear tables, black-universe existence, spectral-action Standard Model terms, or shadow-matter geometry in isolation.
7. If a compact JSON `nonclaims` flag is `false`, public prose must not claim the corresponding statement.

## Evidence labels

Use the labels defined in [THEORY.md](THEORY.md): Postulate, Imported result, Repository derivation, Numerical diagnostic, Open prediction, Interpretive hypothesis.

## Current frontier

Use the [September 10 development update](docs/development-update-2026-09-10.md)
before the earlier preprint checkpoint below. The common-source equations,
boundary state, measure, ADM source maps and flux condition are committed.
The full causal geometric/UV source remains an explicit missing definition.
The separate development snapshot authenticates twelve later records while
preserving the 100-record preprint manifest. A push of authenticated existing
work requires publication/path checks, not replay of the old scientific chain.

The latest source includes the phase-resolved parent Dirac covariance and
its canonical massless stress on the expanding child. Both neck null
contractions are negative, while density, anisotropy and flux leave
independent metric residuals. The raw finite proper-time determinant has a
different thermal state/source relation from the canonical field; the
state-dependent conversion is explicit in its equilibrium control.

Use [the current result](docs/current-result.md) as the work cursor. The next
owner is the common causal/state completion, including the remaining
compact, gauge, boundary and recursive contributions. Reuse the computed
source and normalization evidence. Another angular refinement or a local
coefficient scan does not supply this missing functional. Full physical
couplings, self-sourcing and observational predictions remain open.

At common dimensional energy, parent normalization gives

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

The retarded map must be accompanied by quantum-state data. Induced geometric terms must be counted once. The present compact evidence is indexed by the release manifest; the physical requirements and six targets are in [docs/current-result.md](docs/current-result.md).

## What not to claim

Do not claim that Nested-Space Cosmology, a final \(\zeta\), a particle spectrum, a dark-sector fit, or our universe being inside a black hole has been proved. Do not treat \(1006/1015\), \(54/503\), or any determinant-only scale minimum as a child-matched prediction. Do not repair a failed truncation by inserting sector masses, a second \(\Theta\), or an extra dark function.

## Reading order

1. [README.md](README.md) — one equation, import/derive/invalidate split, map.
2. [THEORY.md](THEORY.md) — labeled thesis.
3. [docs/prior-art-and-open-claim.md](docs/prior-art-and-open-claim.md) — launch surfaces.
4. [docs/current-result.md](docs/current-result.md) — compact checkpoint.
5. [docs/reproducing.md](docs/reproducing.md) — compact-record inspection and regeneration.
6. [CONTRIBUTING.md](CONTRIBUTING.md) — how to change the public record.
