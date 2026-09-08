# Reproducing the compact results

This repository publishes the original 58 v0.1.0 compact JSON records, preserved byte-for-byte, plus scoped follow-ups 59–81. Regeneration checks those records. It does not create a final \(\zeta\), a particle spectrum, a dark-sector fit, or a proof that our universe is inside a black hole.

The live scientific checkpoint is [docs/current-result.md](current-result.md). Claim labels are in [THEORY.md](../THEORY.md). Contribution rules are in [CONTRIBUTING.md](../CONTRIBUTING.md).

Subsequent [compact-interaction](nsc-compact-interaction.md),
[five-dimensional UV matching](nsc-torsion-uv-map.md), and
[published-flow compatibility](nsc-flow-compatibility.md) development evidence
is stored separately in `results/development/`. It does not change the
immutable v0.3.0 release or its 81-record count. Run
`make reproduce-development` to authenticate its inputs and compare every
recorded field. The ordinary `make reproduce` and `make verify` routes include
this focused check; `make reproduce-exact` refers to the released collection.

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
- `make reproduce` recomputes the historical 58-artifact chain portably and follow-ups 59–81 by comparing every recorded field.
- `make reproduce-exact` requires byte-identical recomputation.
- `make paper` rebuilds the tracked paper PDF.
- `make verify` runs check, test, reproduce, and the paper byte check.

Portable mode requires the generator hashes, graph, schemas, classifications,
symbolic expressions, gates, nonclaims, and declared headline observables to
agree. Headline numerics use only the tolerances recorded in the manifest.
Incidental coordinates of a numerically flat diagnostic argmin are not public
observables and may vary between LAPACK implementations. Exact mode remains
the stronger same-environment check and compares every output byte.

For the smooth-geometry convergence diagnostic, the release specification
also declares a conditioned comparison of `observed_orders`. These values
are derived from small differences between independently computed band edges:

\[
e_j=|\lambda_j^{\mathrm{lattice}}-\lambda^{\mathrm{continuum}}|,
\qquad p_j=\log_2(e_j/e_{j+1}).
\]

The checker recomputes both identities within **each** expected and regenerated
record, allowing at most 16 binary64 ULPs for arithmetic/library rounding.
Every error must be finite and larger than the absolute comparison budget
\(b=2\times10^{-12}\), and each regenerated error must differ from its frozen
value by at most \(b\), with no additional relative tolerance. This tightens
the ordinary \(10^{-8}\) absolute error comparison by a factor of 5,000.
The admissible interval for each regenerated order is then calculated from
the frozen errors, rather than chosen as a blanket order tolerance:

\[
\log_2\frac{e_j-b}{e_{j+1}+b}
\;\le p_j^{\mathrm{regenerated}}\le\;
\log_2\frac{e_j+b}{e_{j+1}-b}.
\]

Only a 32-ULP cushion is added when evaluating these interval endpoints; it
covers the permitted identity rounding and endpoint arithmetic. An order
inconsistent with its own errors still fails even when it lies inside the
interval. The policy affects exactly nine `absolute_errors` and six
`observed_orders` values. Raw continuum and lattice spectra, every other
numeric field, keys, conventions, gates, nonclaims, and source hashes retain
their existing comparisons. Exact mode still compares every byte.

The \(2\times10^{-12}\) budget is an explicit portability acceptance threshold,
not a certified total solver uncertainty. The immutable
[smooth-geometry implementation](../src/recursive_horizons/nsc_smooth_geometry.py)
uses a band-root `xtol=1e-12`, an ODE discriminant with `rtol=2e-11` and
`atol=2e-12`, and floating-point eigensolutions. These distinct numerical
errors cannot be bounded from root `xtol` alone. The ratio interval accounts
for how accepted small raw-error differences are amplified in the derived
order; it does not weaken or replace the underlying calculations.

Committed JSON can be inspected without regeneration. The original 58 generators in `scripts/run_*.py` are the historical per-artifact writers; their bytes are frozen. Scoped follow-up checkers are `scripts/check_nsc_scale_closure.py`, `scripts/check_nsc_regulated_recursion.py`, `scripts/check_nsc_radial_spectrum.py`, and `scripts/check_nsc_geometric_chain.py`. Use an individual writer only to create a *new* compact record; they refuse to overwrite a file that already exists.

The unit-closure follow-up authenticates selected v0.1.0 scripts and JSON by SHA-256. Isolated reproduction copies those committed files as auxiliary provenance inputs; it does not add them as extra steps of the historical 58-record chain and does not repin the original laboratory import commit. Direct all-field verification is:

```bash
python3 scripts/check_nsc_scale_closure.py --check
python3 scripts/check_nsc_regulated_recursion.py --check
python3 scripts/check_nsc_radial_spectrum.py --check
python3 scripts/check_nsc_geometric_chain.py --check
python3 scripts/check_nsc_covariant_source.py --check
python3 scripts/check_nsc_measure_normalization.py --check
python3 scripts/check_nsc_influence.py --check
python3 scripts/check_nsc_response_matching.py --check
```

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

The historical 58 generators live in `scripts/` and write a single named file in `results/`. They authenticate declared input JSON by SHA-256, enforce their gates, and write canonical JSON. A representative frontier runner:

```bash
python3 scripts/run_nsc_zeta1_recursion_map.py
```

If the output file already exists, that command exits with an error rather than rewriting history. Use `make reproduce` for the portable public chain. Use an individual runner only when adding a *new* compact record that does not yet exist.

## What regeneration does not mean

Passing `make reproduce` means the committed chain can be recomputed from the published scripts and pins. It does not promote a numerical diagnostic, close the energy-resolved child tail, or convert an imported black-universe, Skyrme, spectral-action, or shadow-matter result into Nested-Space novelty.

The manuscript source and generated preprint are
[paper/nested-space-cosmology.md](../paper/nested-space-cosmology.md) and
[paper/nested-space-cosmology.pdf](../paper/nested-space-cosmology.pdf).

## Short demonstration

Run `make demonstrate` for the computed finite boundary response, prescribed geometry-pulse vacuum work, and invariant sheet/chirality sector. The command checks committed evidence and states the domain of each result. It does not fit observations or generate a new theory claim.
