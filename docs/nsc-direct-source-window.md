# Direct-vacuum source window, 2 <= E < 16

This owner records resumable direct-vacuum evidence for one original slice of
family 14, angular sign +1. It calls the reviewed homogeneous Bloch transport.
It does not call the correction-forcing equation, and it does not rewrite the
source law or the archived arrays.

The slice is the original positive rows with \(2\le E<16\), together with each
row's actual negative partner. The live archive contains 87 such positive
rows: `group14/low16_1` rows 16 through 47, and `group14/low32_1` rows 41
through 95. Their energies run from 2.071557525252559 to 15.9788018699833.
Those 174 signed identities are disjoint from the current 868-row source
coverage. Other energies, the other positive angular family, quadrature error,
and the local gate stay outside this window.

## What each completed row contains

The positive fiber is transported by

\[n_y=2\,h(y)\times n\]

from log-distance \(-18\), with the direct owner's fixed DOP853 tolerances
\(2\times10^{-13}\) and \(2\times10^{-15}\). The settings bound here are 192
bits, frame order 16, 48 metric terms, analytic radius 0.1, maximum step
0.0125, and a degree-12 defect integral. Mass and angular labels are the
archived \(\pi/2\) and \(\sqrt5\), each enclosed by union with the exact value.
The archived columns and covariances are not replaced.

The archived distance and the analytic occupation allowance are stored
separately. Occupation is the reviewed \(\max(\sqrt f, n)\) and is not inserted
into the archived interval. Each signed total is the interval sum of that
distance with the occupation allowance, formed once:

\[
[L,U]+[0,\mathrm{Occ}]=[L,U+\mathrm{Occ}].
\]

The negative partner uses the complemented endpoint \((n_x,-n_y,-n_z)\) on that
partner's own columns and covariance. The positive scalar bound is not copied.
Errors stay unweighted. The archived quadrature weight and the stored column
weight are recorded so a later aggregate can apply them; this window does not.

## Checkpoints, resume and refusal

Each completed positive row is one exclusively published directory

`rows/row-group14__low16_1-NNNNN/` or `rows/row-group14__low32_1-NNNNN/`

containing canonical `record.json` and `witness.npz`. The witness holds the
exact positive and negative columns, both source covariances, and the saved
Bloch trace. The record binds the archive digests, both source digests, both
preparation digests, group, panel, row, all three energy-fiber entries on each
sign, angular sign, mass, \(\rho_\mathrm{up}\), horizon, \(\kappa\), \(\omega\),
the settings above, and the implementation hashes of this owner and the direct
vacuum method.

Resume proves every saved row before it may capture another. Replay validates
the saved trace and does not call the solver. A missing or changed binding is
rejected before any new solve. An incomplete, non-canonical, or extra-file row
directory is rejected rather than repaired. A timeout keeps completed row
directories and does not publish a partial row. `--row-budget` is the maximum
number of new positive rows, and the CPU budget is process time for those new
solves only. Zero rows asks for no new solve.

The JSON status stays `OPEN` when this finite window is incomplete and when it
is complete. Missing rows contribute an identity, not a zero bound. Completing
these 87 pairs does not fill `physical_upstream_budget_component`, does not
cover every family, and does not close the local gate.

## Commands

From the repository root, with the validation environment:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_direct_source_window.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_direct_source_window.py --check --output <checkpoint-directory>
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_direct_source_window.py --resume --row-budget 1 --cpu-budget 120 --output <checkpoint-directory>
```

`--check` writes nothing when the directory is absent and never solves.
`--resume` will not run unless both budgets are present. The default output
path is `results/development/nsc-direct-source-window-v1`. No production
campaign is started by the tests. The frozen middle-source campaign directory
is refused as an output.

The tests solve at most two original positive rows, in one temporary directory,
with their actual negative partners. Resume and replay of those saved traces
do not capture again. Changed settings, a changed source, a changed or missing
payload, a duplicate identity, and a missing or changed binding are rejected.
The census is the exact 87-row window above. Existing production directories
are compared by bytes and are not required to be absent.

## Limitations

This is a local unweighted source comparison on one declared energy window. It
is not an inheritance identity, not an infinite nest, and not an N or beta
residual. Quadrature weights are stored and not applied. The distances are not
contracted through the response kernel. Group 14's other positive angular
family, energies outside \(2\le E<16\), and the remaining covered or uncovered
source rows are not supplied here. The physical gate remains OPEN.

## Observed bounded pilot

The measurement below is filled from the two-row temporary-directory test.
Those two rows are the lowest energies in the window, not a uniform sample of
all 87. The 87-row campaign has not been run.

Pilot figures are recorded after that test. Until then, no full-window CPU
cost is claimed.
