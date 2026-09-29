# Preparation-correction range pilot

This diagnostic selects the first retained group-14, positive-angular,
positive-energy source row at or above E=8, 4 and 2. It keeps the original
preparation and uses the existing order-eight vacuum correction with a nonzero
horizon remainder. These results choose a method; they do not certify a source
panel, fill the upstream budget, or change the local gate's OPEN status.

Run from the repository root:

```sh
.venv/validation/bin/python scripts/lab.py scripts/pilot_nsc_vacuum_correction_range.py --run
```

Each capture has a 30-second CPU limit. The script prints source input hashes,
then the bounds and cost. No source or certificate files are written. The
continuous defect validation runs after capture. Total CPU cost observed on
29 September 2026 was about 15.4 seconds for the four cases.

| Energy | Maximum step | Cells | Vacuum Bloch error | Thermal/coherent allowance | Total positive covariance error upper |
|---|---:|---:|---:|---:|---:|
| 8.0211981300 | 0.1 | 499 | 1.74739e-12 | 7.66255e-24 | 1.46988e-12 |
| 4.0054722763 | 0.1 | 364 | 1.61108e-9 | 2.86392e-12 | 8.09024e-10 |
| 2.0715575253 | 0.1 | 259 | 1.24454e-6 | 1.07189e-6 | 1.69419e-6 |
| 4.0054722763 | 0.025 | 673 | 3.80238e-11 | 2.86392e-12 | 2.24883e-11 |

The refined E=4 case reduces the total bound by about 36 times. Its coarse
bound therefore should not be treated as a measured source discrepancy.
The E=2 vacuum surrogate has a substantial finite-occupation allowance that
step refinement cannot remove. A useful bound there requires the full prepared
state, including thermal/coherent information, or a sharper bound retaining
that information. None of these observations excludes the physical class.

Only selected positive rows were examined. Other rows, negative partners,
quadrature and the N/beta aggregation have not been calculated by this pilot.
The [validated correction owner](nsc-vacuum-source-correction.md) and
[continuous insertion](nsc-ks-source-operator-majorant.md) define their separate
proof scopes. The resumable middle-energy coverage driver is being developed
separately; this pilot does not claim its completion.

A subsequent [source-occupation calculation](nsc-source-occupation-enclosure.md)
uses the incoming mass threshold already present in the original source law.
It removes the unnecessarily large incoming allowance at the tested E~2 row;
a smaller-step correction then gives a positive covariance bound below
1.969e-12. The table above remains the historical coarse-method measurement.
