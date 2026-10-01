# The NSC papers and their relationship

These documents belong to one research programme organized by

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The relation is the inheritance postulate. Its recorded finite operator,
Schur-complement, recursion and transfer calculations remain reusable within
their stated domains. Continuing the local calculation does not erase them.

The latest user instruction directs the active task. An accepted conversation
plan applies where it stays compatible with that instruction and does not
override later steering. NSC seeks the shared structure of nested regions and
the local response that structure leaves. Dated manuscript wording and
archived execution notes do not override that scope. Cosmological
background/perturbation matching remains subsequent work.

## Foundational manuscript

[Nested-Space Cosmology, v0.26.0 (43-page PDF)](../paper/nested-space-cosmology.pdf) ·
[Markdown source](../paper/nested-space-cosmology.md) ·
[Metadata](../paper/metadata.json) · [Build manifest](../paper/build-manifest.json)

This is the broad working preprint dated 15 September 2026. It contains the
physical motivation, one-spectrum/inheritance formulation, operator and source
calculations, developments up to that date, limitations and references. Its
source and PDF are preserved byte-for-byte. Statements about the next step in
that dated version are historical; current work is described in
[the current-result page](current-result.md).

## Current focused companion

[Regional Transfer, Geometric Feedback, and Local Memory (14-page PDF)](../paper/finite-regeneration.pdf) ·
[LaTeX source](../paper/finite-regeneration/main.tex) ·
[arXiv source package](../paper/finite-regeneration-arxiv.tar.gz) ·
[Build manifest](../paper/finite-regeneration-manifest.json)

The 1 October 2026 edition presents the measured finite spherical realization:
regional source differences, surface exchange, actual metric response,
maintained localization, and a local response on the same generated trajectory.
The finite operator inheritance construction is reused in its stated domain.
Source controls, matched proper clocks, independent curvature jets, and causal
omission errors are reported with their own numerical comparisons. This is
the current edition of the focused companion, alongside the preserved foundation.

The pinned science is commit 8a52257fc2829c73d485c7b4c6e310baadfe3c38.
The [evidence snapshot](../paper/finite-regeneration/evidence/snapshot.json)
authenticates its complete declared source/input graph, including historical
source contexts. Tectonic 0.17.0 produced two byte-identical clean cached builds;
fonts are embedded and text is searchable. The source package contains the
manuscript, bibliography, generated bibliography, and four used figure PDFs.
The actual arXiv TeX Live processor and human author review remain pending.

## Historical focused edition

[Nested Qualities and Local Responses (5-page PDF)](../paper/local-incoming-gate-draft.pdf) ·
[Preserved source](../paper/local-gate-draft/main.tex) ·
[Source archive](../paper/local-incoming-gate-draft-source.tar.gz) ·
[Build manifest](../paper/local-gate-draft-manifest.json)

The 27 September 2026 edition retains its finite nested-operator theorem and
historical incoming-gate appendix. Its source, PDF, archive and manifest stay
unchanged. The [proof note](nsc-nested-qualities.md) and
[frozen exact-control snapshot](../paper/nested-quality-evidence/snapshot.json)
remain available. Earlier editions are history of this same companion;
they do not create another active paper role.

## Status and verification owners

[Current result](current-result.md) owns scientific status, and
[the finite realization note](../lab/docs/nsc-spherical-conformal-continuation.md)
routes the numerical owners. [The catalog](../paper/catalog.json) has two
active document roles: the foundation and the current focused companion,
with historical companion links. The measured finite result and the OPEN,
paused incoming-gate application have separate status fields.

Use the read-only command:

~~~sh
make finite-check
~~~

It authenticates the declared evidence graph, source package and PDF without
repeating a scientific campaign or compiling TeX. An intentional new article
build uses the pinned project builder through make finite-paper. The historical
[OPEN draft guide](local-gate-draft.md) remains its old edition's verification
owner. [Reproduction](reproducing.md) explains both paths. The measured result
supports the new v0.27.0 publication; the frozen foundation remains v0.26.0.

The [v0.27.0 research release](https://github.com/kodareken/nested-space-cosmology/releases/tag/v0.27.0)
is published. Its [verification receipt](../results/releases/v0.27.0.json)
records the science and publication commits, annotated tag, successful CI run,
and the hashes of all five downloaded assets.
