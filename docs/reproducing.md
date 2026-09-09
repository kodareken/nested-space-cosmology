# Reproducing the compact results

This repository publishes the original 58 v0.1.0 compact JSON records, preserved byte-for-byte, plus scoped follow-ups 59–81 and nineteen nested development records 82–100. Regeneration checks those records. It does not create a final \(\zeta\), a particle spectrum, a dark-sector fit, or a proof that our universe is inside a black hole.

The live scientific checkpoint is [docs/current-result.md](current-result.md). Claim labels are in [THEORY.md](../THEORY.md). Contribution rules are in [CONTRIBUTING.md](../CONTRIBUTING.md).

The nested records live under `results/development/`:
[compact interaction](nsc-compact-interaction.md),
[five-dimensional UV matching](nsc-torsion-uv-map.md),
[published-flow compatibility](nsc-flow-compatibility.md),
[charged-sector parity](nsc-charged-self-sourcing-route.md),
[vacuum/charge coefficients](nsc-vacuum-charge-matching.md),
[compact boundary action](nsc-compact-boundary-action.md), and
[curved compact source](nsc-compact-casimir.md).
Records 82–85 remain the v0.4.0 nested set. Records 86–88 are the v0.5.0
follow-ups. Records 89–91 add the [horizon source](nsc-horizon-source.md),
[warped quantum source](nsc-warped-source.md), and
[compact/light matching](nsc-compact-matching.md) in v0.6.0. Their JSON still lacks `artifact_id`; manifest identities are
assigned outside those files. Internal status strings that mention v0.3.0
remain historically true. The compact-mass check is an explanatory application
note, not a separate record. The later generators recursively authenticate
input records; isolated reproduction retains those raw dependency files and
normalizes only result input digests after authentication. Run
`make reproduce-development` to recompute the nineteen nested records with
their recorded all-field policies. `make reproduce` and `make verify` cover
the complete 100-record graph; `make reproduce-exact` requires byte identity.

Records 92–100 add the gauge source, spherical action, curvature EFT,
full spectral endpoint, child-state bound, absolute massless reference,
finite angular reference, parent-matched source and state/regulator conversion.
All nine retain their laboratory JSON, generator bytes and all-field tolerances.
The original 91 manifest step objects remain unchanged.

The angular-stress record authenticates the laboratory `pyproject.toml`.
Its immutable copy is kept under `tests/fixtures/source-checkpoints/` and
restored to its original path only inside an isolated reproduction workspace.
The public installation configuration is separate. The C++17 angular
transport source is included in both the source closure and the package.
Its large angular integration and the parent-state calculation dominate the
new runtime; the Linux reproduction job has a 120-minute ceiling. This is
an execution limit, not a claimed measured runtime.

## Requirements

- Python 3.12 or newer
- A C++17 compiler with thread support for angular-stress
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
- `make reproduce` recomputes the historical 58-artifact chain portably and follow-ups 59–100 by comparing every recorded field. Nested development paths are copied and compared explicitly.
- `make reproduce-exact` requires byte-identical recomputation.
- `make paper` rebuilds the tracked paper PDF.
- `make verify` runs check, test, reproduce, and the paper byte check.

Portable mode requires the generator hashes, graph, schemas, classifications,
symbolic expressions, gates, nonclaims, and declared headline observables to
agree. Headline numerics use only the tolerances recorded in the manifest.
Incidental coordinates of a numerically flat diagnostic argmin are not public
observables and may vary between LAPACK implementations. Exact mode remains
the stronger same-environment check and compares every output byte.

The spectral-endpoint record also contains the numerical residual
`independent_derivative_error = extrapolated_log_derivative - analytic_derivative`.
At radius4, it subtracts derivatives near4754.59 to obtain a value near1e-7.
Linux CI exposed a change from -1.40e-7 to +1.56e-7, although both estimates
satisfy the generator's original accuracy requirement. Portable mode now
checks the Richardson extrapolation and subtraction identities within eight
ULPs, enforces the original `3e-8*max(1,abs(analytic_derivative))` accuracy,
and propagates the two operands' original comparison budgets to their
difference. Every raw derivative retains its original tolerance. Only the
four named residual fields receive this declared derived-quantity policy;
the scientific JSON, generator, comparator and exact mode remain unchanged.
This is a portable acceptance rule, not a new physical error estimate.
The generator's direct `--check` retains its original comparison; use the
public portable route for this cross-runtime comparison. If CI reproduction
fails, it retains generated public result JSON as a short-lived diagnostic
artifact. Private laboratory raw stores and credentials are not part of
that isolated result directory.

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
python3 scripts/check_nsc_compact_interaction.py --check
python3 scripts/check_nsc_torsion_uv_map.py --check
python3 scripts/check_nsc_flow_compatibility.py --check
python3 scripts/check_nsc_charged_sector.py --check
python3 scripts/check_nsc_vacuum_charge_matching.py --check
python3 scripts/check_nsc_compact_boundary_action.py --check
python3 scripts/check_nsc_compact_casimir.py --check
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

Run `make demonstrate` to inspect authenticated stored results and their physical scope. It checks the release specification, record hashes and source/input closure, then displays concise results. It does not launch scientific generators by default. To rerun the selected checks, use `make demonstrate-recompute` or `python3 scripts/demonstrate.py --recompute`. Neither route fits observations or creates a new theory claim.
