# Active-ramp continuous-field pilot

The Stage-2 pilot in `scripts/validate_nsc_ks_field_ramp_pilot.py` selects cell
150 of the unchanged four-energy family-14_1 stream. This is in the active
normal-cutoff transition, rather than the initial undeformed part of the
trajectory. It is a new cell, not another replay of the existing cell 122.

The adapter reuses the original trajectory producer's source, 1024-node grid,
129 target nodes, zero-tangent DOP853 settings, lock and 120-second capture CPU
cap. It changes only the cell selection and successor output paths. The
capture binds the adapter and the original producer. The enclosure reuses the
authenticated v4 Fourier-profile payload and current whole-cone method.

Before enclosure, the saved weights, energies, coherent covariance, initial
columns, grid and target nodes must equal the current declared inputs. The
saved cell keeps its exact segment identity. The normalized residual integral
equals the bound on the normalized-cell residual; a second cell-width factor
is rejected. A scientific digest excludes runtime and detects altered results.

Recorded integral uppers for derivative orders 0 and 1 are approximately
`2.73421533437e-16` and `1.24093214881e-12`. Exact dyadic endpoints are in
`results/development/nsc-ks-field-ramp-cell150-v1.json`. The initial enclosure
took about 56 CPU seconds on this Mac. The polynomial contribution dominates
the time/radius remainder in this cell. These are field-residual integrals,
not N/beta constraint error bounds and not directly comparable to the total
constraint budget without the propagation and contraction owners.

Only one cell and four selected energy rows are enclosed. The upstream
physical-source uncertainty, whole-trajectory coverage and full family coverage
remain missing. No 60-family campaign or nonlinear search was launched. A
production scaling decision needs the integrated propagated and contracted
budget, including the remaining source terms.

```sh
.venv/validation/bin/python scripts/lab.py scripts/validate_nsc_ks_field_ramp_pilot.py --check
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_field_ramp_pilot.py
```

`--check` authenticates the saved record and its dependencies; it does not rerun
the enclosure. `--capture` and `--enclose` refuse to replace existing records.
The captured payload is retained when a resource stop prevents an enclosure.
