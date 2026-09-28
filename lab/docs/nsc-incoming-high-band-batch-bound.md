# Remaining retained high-band energy certificates

Reuse the frozen order24 trimmed Riccati recurrence, rank-one energy estimate,
128-cell collar map and centered degree4 interval remainder. The missing
connection is matched order24 vacuum/thermal error control for the remaining
28 retained channels, after excluding completed groups12,14,22,32. The batch
contains 55 disjoint source windows: LOW `[32,40]` plus each actual middle
interval for 27 groups, and only the existing `[16,160]` middle for group13.

First certify group6 once. If its component budget passes, continue the other
27 channels at the same resolution with at most three CPU worker processes.
Each worker caches common geometry; each channel gets one radial pass and a
content-addressed receipt. Completed compatible receipts are resumed instead
of recalculated. Stop after this one batch at this resolution, whether PASS
or OPEN. No source, mode, scattering, metric, parameter or action generation
belongs to this work. Passing bounds require separately matched order24 source
corrections and do not certify archived order16 values or source quadrature.

## Exact sign symmetry, without an occupation assumption

The existing recurrence starts with `c1=L+i*m`. For real geometric jets its
angular-sign transformation is

$$
c_n(-\lambda)=(-1)^n\overline{c_n(\lambda)},\qquad
f_p(-\lambda)=(-1)^{p+1}\overline{f_p(\lambda)}.
$$

The derivative and nonlinear terms have these same parities separately.
An exact symbolic induction/coefficient certificate checks the base case,
every recurrence degree through24, and every defect degree24 through48.
Conjugation and sign preserve rectangular complex norm bounds. A small actual
interval recurrence control checks both signs independently, and saved
completed-channel receipts confirm equal signed norm arrays. The batch can
therefore compute one positive-angular norm enclosure and explicitly copy it
to both angular signs with their existing per-sign multiplicities. Physical
source signs remain distinct; no occupation symmetry or state is asserted.

Each saved energy-independent radial coefficient set is contracted over its
actual windows. A combined interval is an alternative aggregate, not another
source contribution. The record sums each group's union once, preserving the
55-window accounting. Thermal terms retain their positive inherited bound.

`--prepare` saves individual receipts as workers finish. `--check` validates
all producer/input hashes, replays coefficient contractions and verifies
aggregation without any recurrence. A running process is not restarted when
an observation call times out.

```sh
python3 scripts/derive_nsc_incoming_high_band_batch_bound.py --prepare
python3 scripts/derive_nsc_incoming_high_band_batch_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_high_band_batch_bound.py
```
