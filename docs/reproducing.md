# Reproducing the compact results

This repository publishes 58 compact JSON records and the scripts that generated them. Regeneration checks those records. It does not create a final \(\zeta\), a particle spectrum, a dark-sector fit, or a proof that our universe is inside a black hole.

The live scientific checkpoint is [docs/current-result.md](current-result.md). Claim labels are in [THEORY.md](../THEORY.md). Contribution rules are in [CONTRIBUTING.md](../CONTRIBUTING.md).

## Requirements

- Python 3.12 or newer
- the pinned dependencies in `requirements.txt`:
  - `mpmath==1.3.0`
  - `numpy==2.5.1`
  - `scipy==1.17.1`
  - `sympy==1.14.0`

Optional paper rendering uses `matplotlib==3.11.1` and `reportlab==4.4.9`
from `requirements-paper.txt`.

From the repository root:

```bash
python3 -m pip install -e ".[paper]"
```

or:

```bash
make install
```

## Public commands

The Makefile is the declared public gate:

```bash
make check
make test
make reproduce
make paper
make verify
```

- `make check` validates manifests, claims, paths, and links.
- `make test` runs the focused publication tests.
- `make reproduce` recomputes the 58-artifact chain in portable mode.
- `make reproduce-exact` requires byte-identical recomputation.
- `make paper` rebuilds the tracked paper PDF.
- `make verify` runs check, test, reproduce, and the paper byte check.

Portable mode requires the generator hashes, graph, schemas, classifications,
symbolic expressions, gates, nonclaims, and declared headline observables to
agree. Headline numerics use only the tolerances recorded in the manifest.
Incidental coordinates of a numerically flat diagnostic argmin are not public
observables and may vary between LAPACK implementations. Exact mode remains
the stronger same-environment check and compares every output byte.

Committed JSON can be inspected without regeneration. The 58 individual generators in `scripts/run_*.py` are the present per-artifact writers. Use them only to create a *new* compact record; they refuse to overwrite a file that already exists.

## Committed JSON is the public result

Each compact record under `results/` is the scientific object. Typical runners refuse to overwrite an existing file. Do not edit committed JSON by hand. Do not delete a committed record to manufacture a “fresh” success.

To inspect the live frontier record:

```bash
python3 - <<'PY'
import json
from pathlib import Path
p = Path("results/nsc-2-zeta1-recursion-map.json")
d = json.loads(p.read_text())
print(d["artifact_id"])
print(d["classification"])
print(d["nonclaims"])
print(d["gate"]["mode_resolved_tail_solved"], d["gate"]["zeta_derived"])
print(d["next_result"])
PY
```

The expected nonclaims are all `false`: \(\Omega\) is not selected, the functional tail is not computed, the scale root is not restored, and physical \(\zeta\) is not promoted.

## Individual runners

The 58 generators live in `scripts/` and write a single named file in `results/`. They authenticate declared input JSON by SHA-256, enforce their gates, and write canonical JSON. A representative frontier runner:

```bash
python3 scripts/run_nsc_zeta1_recursion_map.py
```

If the output file already exists, that command exits with an error rather than rewriting history. Use `make reproduce` for the portable public chain. Use an individual runner only when adding a *new* compact record that does not yet exist.

## What regeneration does not mean

Passing `make reproduce` means the committed chain can be recomputed from the published scripts and pins. It does not promote a numerical diagnostic, close the energy-resolved child tail, or convert an imported black-universe, Skyrme, spectral-action, or shadow-matter result into Nested-Space novelty.

The manuscript source and generated preprint are
[paper/nested-space-cosmology.md](../paper/nested-space-cosmology.md) and
[paper/nested-space-cosmology.pdf](../paper/nested-space-cosmology.pdf).
