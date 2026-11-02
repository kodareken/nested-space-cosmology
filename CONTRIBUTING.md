# Contributing

Nested-Space Cosmology welcomes derivations, calculations, constructive criticism
and clearer formulations. Contributions should make the common-operator
construction more explicit, reproducible and physically informative.

**Douglas Ek** is the accountable author. **ChatGPT** and **OpenAI Codex** have been used for formulation, code, computation, drafting, comparison, and audit. Pull requests that include AI-generated text or code must still be checked against explicit mathematics, committed compact JSON, cited sources, and stated nonclaims. AI assistance does not transfer scientific responsibility.

## Before opening a change

Start with the [foundation](THEORY.md) and [active code map](docs/active-code-map.md).
When reusing mathematics, state the equation, its assumptions and the
specific connection left to calculate. A change that merely repeats known
mathematics or adds another passing baseline does not establish a new connection.
State the decision the result can change, the smallest useful check and its
stopping condition. Use the existing work entry; no extra tracking system is
needed. A rerun must identify a changed dependency, failure or concrete
unresolved concern. Apply the same scope to delegated work.

1. Read [AGENTS.md](AGENTS.md) and [current instructions](docs/instructions.md). Do not edit AGENTS.md without Douglas's explicit instruction.
2. Use [the current result](docs/current-result.md). Consult retained owners only when they supply a needed input.
3. State whether the change is a postulate, an imported result, a repository derivation, a numerical diagnostic, an open prediction, or an interpretive hypothesis.
4. Cite primary literature for technical factual claims. Do not present imported black-universe, Skyrme, QFT, spectral-action, or shadow-matter results as novelty.
5. For a numerical result, include inputs, units, uncertainty treatment, code path, package versions, and an explicit failure condition.
6. If a compact JSON `nonclaims` flag is `false`, do not write prose that treats the corresponding statement as shown.

## Current frontier

The current frontier is the geometric offset construction and the contact
law for opposite windings. The contraction and conditional two-output
result are derived in the [construction owner](docs/nsc-atlas-offset-generator.md).
The contact dynamics and output selection are open. Old campaigns and
their specific constraints remain in the [retained library](docs/research-index.md).
They are not a work queue after compaction.

## Verification

Match verification to the change. Navigation edits need link checks, not
fresh numerical derivations. For
scientific changes, use focused checks while developing and one integrated
verification pass at a coherent boundary. Repeat an expensive run only after
a relevant change, failure or unresolved concern makes it informative.

The active foundation check is:

```text
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_offset_generator.py
```

The broader historical checks are available when that work is explicitly selected:

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

Lead with the equation, result and physical connection. The README contains
the general research-status statement; avoid repeating it throughout the
paper. Keep the domain, assumptions, provenance and remaining terms beside
the particular result they qualify. Write established identities directly.

Write for two readers at once. Give a general reader the visual relationship
first—gradient, oscillation, spectrum, resonance, boundary, inheritance—then
give the exact operator statement a physicist can inspect. Jargon should name
a relationship already made understandable, not replace the explanation.

- Say “conditional on” when a result depends on a postulate.
- Say “inferred from” when the target observable is used to determine an upstream parameter.
- Reserve “prediction” for a quantity fixed without using the data later compared with it.
- Reserve “measurement” for a validated estimator with uncertainty and provenance.
- Separate a coordinate statement from an invariant or locally observed one.
- State the domain of a theorem or model result; a symmetry-reduced identity is not a generic collapse result.
- Do not claim that Nested-Space Cosmology, a final $\zeta$, a particle spectrum, a dark-sector fit, or our universe being inside a black hole has been proved.

## What not to mix

Keep each repository a coherent public object. Do not combine this programme with unrelated private infrastructure, internal authorization language, or unpublished chronology. Use only relative links inside the repository, plus ordinary literature citations.

The [paper guide](docs/papers.md) identifies both sources, PDFs and manifests.
The foundational manuscript is preserved; the focused companion does not
replace it.
