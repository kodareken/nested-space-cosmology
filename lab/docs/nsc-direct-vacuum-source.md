# Direct homogeneous vacuum Bloch transport

This is a reusable per-row method. It transports one positive-energy vacuum
branch with the homogeneous Bloch equation and compares that branch with one
archived source fiber. It is not full source coverage, not a tail-remainder
certificate, and not local-gate closure. The source occupation law is unchanged.

Owner: `src/recursive_horizons/nsc_direct_vacuum_source.py`.

## Same channel for the frame and the generator

`DirectVacuumBloch` takes a positive energy, a nonnegative mass, a signed
angular label and a horizon hint. It builds `metric_horizon_frame` from those
four values and uses that frame's horizon root and energy in
`BlochSource.hamiltonian_series` and `BlochSource.rhs_numeric`. There is no
frame argument, so an unrelated frame cannot be paired with the generator.
Any positive energy inside the frame and metric domains is allowed, including
energies above the mass. `BlochSource`'s own constructor remains restricted to
the positive subgap channel and is not used to build this model.

`VacuumCorrection` is not a base class. Its expansion forcing is not evaluated
and is not subtracted. The numeric right-hand side is the homogeneous cross
product `2 h × n`.

## Vacuum branch, tail and norm

On the interior chart the first frame column is the positive-energy vacuum
branch, including when `E > m`. Its component balls keep the Cauchy tail of the
frame series. The initial Bloch vector is the projector of that column:

```text
n = (2 Re(v0 conj(v1)), -2 Im(v0 conj(v1)), |v0|^2 - |v1|^2).
```

The integrator stores binary midpoints. Validation passes the balls, not those
midpoints, back through `validate_bloch`. The resulting initial error is the
ball-to-float distance. It is nonzero and is kept inside the Bloch error.
The interior generator is anti-Hermitian, so the exact evolution preserves the
column norm. The validated bridge uses the enclosed initial norm rather than
replacing it by one.

The declared start must match the first stored time exactly. Adjacent cells
must join in time and in the three Bloch components. The target is a radial
coordinate `rho >= 1` in the same horizon chart; the metric collar `|u| < 0.4`
is required, and the log-distance round trip is the existing endpoint bridge.

## Archived matrices and occupation, kept separate

The archived comparison calls `original_covariance_distance_bounds` on the
supplied 2×3 columns and 3×3 covariance. It does not rewrite those arrays.

The negative partner is not given the positive bound. Its endpoint is the
complement

```text
(n_x, -n_y, -n_z),
```

which is the Bloch image of `I - S3 conj(Q_+) S3`, and that endpoint is
compared with the actual archived negative columns and covariance. The lower
and upper bounds are recomputed. Copying the positive upper is not a negative
result.

The finite-occupation allowance is the reviewed `max(sqrt(f), n)` from
`occupation_vacuum_distance`. `compare_signed_archives` returns it beside the
two archived distances and does not add it. A caller that wants one positive-row
source upper adds that allowance once. Adding it inside the matrix distance, or
adding it twice, is a different claim. Missing, boolean, or omitted archive
entries are rejected rather than replaced by zero.

Capture takes an explicit positive CPU budget and step. If the budget is
exhausted, the integrator raises and no trajectory is returned.

## One demonstrated row

The reviewed pilot, and the module test, use the original group14 positive
family, the smallest archived energy at least 2, and that fiber's actual
negative partner. The mass and angular label are the archived `π/2` and
`√5`, enclosed by union with the exact values. The start is log-distance `-18`,
the step is `0.00625`, and the precision is 192 bits. On that positive row the
pilot reported vacuum error `1.59008e-12` and, after adding the occupation
allowance once, positive source error `2.38908e-12`, with 2680 cells in about
2.3 CPU seconds. Those two figures are the demonstrated positive-row
comparison only. They do not certify the negative partner, neighboring
energies, quadrature, the tail remainder, or the local gate. The gate stays
OPEN.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_direct_vacuum_source.py
```

## Parent integration controls

After correcting pytest patch calls and comparing interval endpoints rather
than nonexact Arb objects by equality, all eight tests passed. The method now
also rejects complex Bloch traces and missing/boolean Bloch errors explicitly.
The tests compare both actual archived signs; the earlier prototype figures
above refer only to its positive-row report.

A same-row step control used the unchanged frame and source, and added the
occupation allowance once to each signed archived distance:

| Maximum step | Cells | Positive upper (rounded up) | Negative upper (rounded up) | Total CPU seconds |
|---|---:|---:|---:|---:|
| 0.025 | 753 | 2.242e-11 | 2.233e-11 | 0.63 |
| 0.0125 | 1340 | 2.410e-12 | 2.313e-12 | 1.13 |

Total time includes capture, validation and both comparisons. These are
selected-row calibration measurements, not uniform settings for every source
family. No production coverage or gate budget is filled by this table.

Independent Grok review of the integrated method passed all eight tests and
found no bound, phase, domain or input defect within this per-row contract.
The occupation allowance must still be added once by the caller. Only the
selected moderate-energy fiber has been demonstrated; a source-window recorder
and per-row bindings are separate work.
