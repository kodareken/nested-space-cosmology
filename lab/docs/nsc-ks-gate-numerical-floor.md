# One-knob numerical floor

The base is the value-only iterate6 residual at 47 nodes. Each later row
changes one solver setting and records the change in the nodal maxima.
Those changes are indicators. They are not error bounds, and no drift is
subtracted.

Loop numerics are the cheapest setting whose indicator floor is at or below
`1e-10` per component. Certificate numerics need an indicator floor at or
below `3e-12`. Neither is declared until the ladder has been measured.

```sh
python3 scripts/derive_nsc_ks_gate_numerical_floor.py --record
python3 scripts/derive_nsc_ks_gate_numerical_floor.py --next
python3 scripts/derive_nsc_ks_gate_numerical_floor.py --check
```

`--next` evolves one missing label, and only when no other `derive_nsc_`
process is running with `--run`.
