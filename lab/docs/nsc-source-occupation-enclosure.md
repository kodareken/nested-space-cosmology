# Vacuum distance of the unchanged affine source

The three-channel source already has horizon occupation f, coherence s and
incoming occupation n, with s squared = f(1-f). Subtracting the vacuum projector
P=diag(0,1,0) gives the horizon block [[f,-i s],[i s,-f]] and the third entry n.
The square of the horizon block is exactly f times the identity. Consequently

    ||C_source - P||_operator = max(sqrt(f), n).

This keeps coherence. Replacing sqrt(f) by f would discard its contribution.
The original source law sets n=0 for positive energy E<omega*mass. The helper
uses that branch only when the entire input enclosure satisfies the strict
inequality. An interval crossing the threshold retains the nonzero upper
bound; equality is on the nonzero branch, as in the existing source law.

A normalized 2x3 sewing map has operator norm one. Its homogeneous unitary
continuation therefore cannot increase this source-to-vacuum distance. This
is the finite-occupation term that can be added to a separately proved vacuum
preparation error. It neither sets the reflection to zero nor changes source
columns, temperatures, occupations, or the upstream preparation.

For original group14 E=2.071557525252559, the incoming occupation is already
zero, and the analytic bound is below 1.383e-12. The older generic exponential
allowance was about 1.072e-6 because it did not use the incoming mass threshold.
That older bound remains valid but is unnecessarily large for this row.
This improves a proof estimate, not the underlying state.

The identity concerns the analytic source. Binary storage can round 1-f to 1;
its literal matrix distance can then exceed the analytic bound by a tiny
rounding term. The tests demonstrate that distinction explicitly. When the
helper is used in preparation comparison, the unchanged stored covariance is
already compared separately against the exact vacuum; its arithmetic must not
be silently identified with the analytic source matrix.

Tests check the exact block identity, normalized sewings, signed partners,
coherence and incoming-channel omission controls, and threshold-crossing
intervals. Five tests pass. Independent review passed all five tests and accepted the scoped analytic
bound. Integration into successor source records remains pending. No existing record or gate budget is changed.

Owners: `src/recursive_horizons/nsc_source_occupation_enclosure.py` and
`tests/test_nsc_source_occupation_enclosure.py`.

## Moderate-energy diagnostic

On original group14/low32_1 row41, E=2.071557525252559, the sharper occupation
allowance was combined with the existing validated vacuum correction:

| Expansion order | Initial log distance | Maximum step | Cells | Positive covariance-error upper (rounded up) | CPU seconds |
|---|---:|---:|---:|---:|---:|
| 4 | -40 | 0.025 | 1553 | 3.071e-11 | 6.69 |
| 8 | -18 | 0.00625 | 2682 | 1.969e-12 | 22.50 |

The eighth-order Bloch error was 7.5051e-13; the retained analytic occupation
allowance was 1.38252e-12. This is a positive-row diagnostic, not a saved-witness
source record. No negative partner or neighboring energy is counted by it.
A separate near-threshold test at E=1.5594631247710442 failed the current phase
tube check (cell1138 of1544); that setting did not certify the row. A later [wider-tube control](nsc-near-threshold-phase-pilot.md) succeeds with explicit bounds.

Reproduce the two positive-row diagnostics with
`python scripts/lab.py scripts/pilot_nsc_moderate_source_occupation.py --run`
using the validation interpreter. No witness or certificate is written.
