# High-angular preparation-method pilot

Three original positive-angular, positive-energy source rows were tested on
29 September 2026: the first E>=32 row of groups12,22,32. These sample the largest
angular level in each retained compact-mass sector. The source arrays and
occupation law remain unchanged.

| Group | Energy | Direct M16 covariance-error upper | Corrected M8 covariance-error upper | Correction cells | CPU seconds |
|---|---:|---:|---:|---:|---:|
| 12 | 32.019251120011916 | 1.634e-9 | 4.470e-13 | 952 | 6.11 |
| 22 | 32.010944552602076 | 1.215e-10 | 3.130e-13 | 886 | 7.94 |
| 32 | 32.010944552602076 | 3.288e-10 | 2.894e-13 | 932 | 8.12 |

The displayed error bounds are rounded upward. Timing includes both methods.
The direct method uses 128 cells for its coefficient remainder; the correction
uses its existing nonzero horizon initializer and continuously validated driven
Bloch trajectory, with max_step=0.1 and a 30-second capture CPU cap per row.
Both retain the finite thermal/coherent allowance. These are original-source
covariance comparisons, not N/beta residuals or integrated source-budget terms.

The direct expansion alone is too loose to serve as a universal replacement
for the correction method near E=32 in these high-angular families. This is a
proof-method limitation, not evidence excluding the physical class. The smaller
corrected bounds support using witnesses where the direct bound is excessive.
Actual source-weighted aggregation still decides the contribution to the gate.

No negative partners or neighboring rows were included by this pilot. It does
not expand the authenticated coverage in the v2 source-insertion record.
The script prints original input hashes, measured bounds and timing; it does
not write scientific checkpoint files.

```sh
.venv/validation/bin/python scripts/lab.py scripts/pilot_nsc_high_angular_preparation.py --run
```

Owners reused: [vacuum remainder](nsc-vacuum-source-remainder.md),
[correction proof](nsc-vacuum-source-correction.md). The original three-row
control was executed; the repository wrapper adds explicit opt-in and compact
output formatting to that same calculation.
