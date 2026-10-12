# Actual18 direct-vacuum rows, pi/2 <= E < 2

This owner records resumable direct-vacuum evidence for nine original
positive rows of family 14, angular sign +1, and each row's actual negative
partner. The rows are `group14/low32_1` 32 through 40. Their archived
energies run from `0x1.92f9814c33825p+0` to `0x1.fae56f00f6a85p+0`, inside
`pi/2 <= E < 2`.

The frozen direct window `2 <= E < 16` starts at `group14/low32_1` row 41.
It does not own these nine rows. This recorder does not call that window,
does not rewrite its schema, and does not write its campaign directory.

## Transport

The homogeneous law is the reviewed direct-vacuum owner:

`n_y = 2 h(y) cross n`

from log-distance `-18`, DOP853 tolerances `2e-13` and `2e-15`, 192 bits,
frame order 16, 48 metric terms, analytic radius `0.1`, maximum step
`0.0125`, and a degree-12 defect integral. Mass and angular labels are the
archived `pi/2` and `sqrt(5)`, each enclosed by union with the exact value.
The archived columns and covariances are not replaced. No new source solver
is used.

The negative partner keeps its own columns and covariance. Those columns are
the S3 conjugate of the positive columns. The covariance satisfies
`Cminus = I - conj(Cplus)`. The Bloch endpoint compared with that partner is
`(nx, -ny, -nz)`, the image of `I - S3 conj(Q) S3`. The opposite full-vector
flip is not used. Adjudication `c82b600b` is recorded and is not inverted.

Occupation is the reviewed `max(sqrt(f), n)`. It is stored separately and
added once:

`[L, U] + [0, Occ] = [L, U + Occ]`.

Errors stay unweighted. The archived quadrature weight and the stored column
weight are recorded so a later aggregate can apply them.

## Checkpoints

Each completed positive row is one exclusively published directory

`rows/row-group14__low32_1-NNNNN/`

with canonical `record.json` and `witness.npz`. Resume proves every saved
row before another capture. Replay validates the saved trace and does not
call the solver. A missing or changed binding is rejected before a new
solve. An incomplete or extra-file row directory is rejected rather than
repaired. A timeout keeps completed directories and does not publish a
partial row. Missing rows contribute an identity and null bounds, not zero.

The authoritative source is commit `eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9`
of `kodareken/nested-space-cosmology`. The reviewed snapshot copies of the
direct-vacuum and occupation owners match that commit. Baseline family
source hashes and archive input digests are stored with each row.

## Commands

From `threshold18`, with the project-local validation environment:

```sh
.venv/bin/python -m pytest -q tests/test_nsc_threshold18_source.py
.venv/bin/python scripts/derive_nsc_threshold18_source.py --check --output results/development/nsc-threshold18-actual18-v1
.venv/bin/python scripts/derive_nsc_threshold18_source.py --resume --row-budget 9 --cpu-budget 180 --output results/development/nsc-threshold18-actual18-v1 --report evidence/actual18-report.json
```

`--check` writes nothing when the directory is absent and never solves.
`--resume` will not run unless both budgets are present. The default output
is `results/development/nsc-threshold18-actual18-v1`.

## Mac integration

Do not apply these paths on the Mac from this host. The later drop, under
`/Users/admin/Documents/BlackHoles-Infinity/lab`, is:

- `src/recursive_horizons/nsc_threshold18_source.py`
- `scripts/derive_nsc_threshold18_source.py`
- `tests/test_nsc_threshold18_source.py`
- `docs/nsc-threshold18-actual-source.md`
- `results/development/nsc-threshold18-actual18-v1/`

The pinned lab here is `threshold18/pin/lab` at the commit above. The
checkpoint rows are the evidence to copy. This window does not edit that
Mac tree and does not publish Git.

## Limits

This is a local unweighted comparison of 18 signed identities. It is not an
inheritance identity, not the 108-row subgap batch, and not local-gate
closure. Quadrature weights are stored and not applied. Completing the nine
pairs leaves `physical_upstream_budget_component` null and the physical
local gate OPEN.
