# Paper

## Canonical manuscript

- Source: [`recursive-horizons.md`](recursive-horizons.md)
- PDF: [`recursive-horizons.pdf`](recursive-horizons.pdf)
- Bibliography data: [`references.bib`](references.bib)

The Markdown manuscript is the semantic source of truth. The PDF is reproducibly rendered with pinned ReportLab by `scripts/build_paper.py`; byte-identical output is expected for the same ReportLab version and selected font set. It contains the same text, numbered headings, tables, links, and page metadata.

The cover, manuscript, and PDF metadata credit the work as **authored and
produced by Douglas Ek & ChatGPT 5.6 Sol**. The manuscript's contribution
disclosure records the human research direction and publication
accountability, the AI collaborator's material role, and the verification
boundary applied to AI-generated work.

The [FGC local-defocusing workspace](fgc-local-defocusing/README.md) is a
separate pre-result gate for the future focused technical paper. It is not
included in the prospectus PDF and does not yet contain a result manuscript.

## Build

From the repository root:

```bash
python3 -m pip install -r requirements-paper.txt
python3 scripts/build_paper.py
```

The command writes `paper/recursive-horizons.pdf` atomically and reports its page count and SHA-256 hash. It does not require LaTeX or network access.

## Title choice

The technical title leads because the work is a conjectural research programme.
“Can Einstein and Hawking Finally Rest Easy Together? Here Is What Must Be
Proven” is retained as the popular subtitle: it expresses the intention to
preserve Einstein’s local causal geometry and Hawking/Bekenstein horizon
thermodynamics while making clear that the transition joining them is still an
open, testable construction.

## Citation status

The bibliography favors primary papers. A citation supports only the established ingredient or published observation described in the sentence that cites it. It does not transfer authority to the project’s new postulates.
