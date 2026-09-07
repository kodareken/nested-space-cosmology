# Contributing

Nested-Space Cosmology welcomes criticism, calculations, counterexamples, and clearer formulations. Contributions are evaluated by whether they make the one-equation claim easier to falsify or reproduce—not by whether they protect the current narrative.

**Douglas Ek** is the accountable author. **ChatGPT** and **OpenAI Codex** have been used for formulation, code, computation, drafting, comparison, and audit. Pull requests that include AI-generated text or code must still be checked against explicit mathematics, committed compact JSON, cited sources, and stated nonclaims. AI assistance does not transfer scientific responsibility.

## Before opening a change

1. Read [README.md](README.md), [THEORY.md](THEORY.md), and [AGENTS.md](AGENTS.md).
2. Read [docs/prior-art-and-open-claim.md](docs/prior-art-and-open-claim.md) and [docs/current-result.md](docs/current-result.md).
3. State whether the change is a postulate, an imported result, a repository derivation, a numerical diagnostic, an open prediction, or an interpretive hypothesis.
4. Cite primary literature for technical factual claims. Do not present imported black-universe, Skyrme, QFT, spectral-action, or shadow-matter results as novelty.
5. For a numerical result, include inputs, units, uncertainty treatment, code path, package versions, and an explicit failure condition.
6. If a compact JSON `nonclaims` flag is `false`, do not write prose that treats the corresponding statement as shown.

## Current frontier

**The current frontier is the unsolved energy-resolved recursive child tail.**

Do not restore a scale root by adding an independently weighted geometric action, retuning \(\Phi\) after seeing the target, or inserting a dark-matter or dark-energy function. The derived unsolved equation is

$$
\Gamma(x)=K(x)-B^{\dagger}\Gamma(x/\Omega)^{-1}B.
$$

## Verification

From the repository root:

```bash
make verify
```

or the pieces:

```bash
make check
make test
make reproduce
make paper
```

See [docs/reproducing.md](docs/reproducing.md) for pins, overwrite policy, and what regeneration does not imply.

A calculation may be valuable when it falsifies a preferred truncation. Negative results belong in the compact record and in [docs/current-result.md](docs/current-result.md), not in a private discard pile.

## Scientific writing rules

- Say “conditional on” when a result depends on a postulate.
- Say “inferred from” when the target observable is used to determine an upstream parameter.
- Reserve “prediction” for a quantity fixed without using the data later compared with it.
- Reserve “measurement” for a validated estimator with uncertainty and provenance.
- Separate a coordinate statement from an invariant or locally observed one.
- State the domain of a theorem or model result; a symmetry-reduced identity is not a generic collapse result.
- Do not claim that Nested-Space Cosmology, a final \(\zeta\), a particle spectrum, a dark-sector fit, or our universe being inside a black hole has been proved.

## What not to mix

Keep each repository a coherent public object. Do not combine this programme with unrelated private infrastructure, internal authorization language, or unpublished chronology. Use only relative links inside the repository, plus ordinary literature citations.

The technical paper lives at
[paper/nested-space-cosmology.md](paper/nested-space-cosmology.md).
