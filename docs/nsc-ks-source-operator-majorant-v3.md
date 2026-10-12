# Source insertion including the replayed low row15 pair

The v3 successor replays the complete middle window and the new low-row witness,
then adds that pair to the existing high-energy and first-low-row coverage.
Original group14/low16_1 row15 is distinct from both the old low row0 and the
middle-panel row15. Duplicate identities, missing signed bounds, a changed
upstream slice and already-weighted error inputs are rejected.

The same continuous exact-response insertion and original weights are used.
The result covers868 signed rows of family14_1, with300 remaining in this
family. It retains the complete JSON/NPZ dependency bindings and leaves the
full upstream-budget component null. No other family, quadrature or numerical
field error is supplied, and the local gate stays OPEN. Earlier records retain
their bytes.

Owners: `scripts/derive_nsc_ks_source_operator_majorant_v3.py` and
`tests/test_nsc_ks_source_operator_majorant_v3.py`. The driver supports exclusive
`--record` and full saved-witness `--check`; replay launches no new source solve.
