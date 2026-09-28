# Resolved local envelope indicators

The independent KS comparison has a128/256-node matter change up to8.57e-11
on I, above1e-11. Its source preparation is fixed and its temporal tolerance
comparison already passes. This separate two-case check changes only the
remaining Fourier resolution and then the integration tolerance.

Run512 nodes with the existing tight settings2e-13/2e-15/.0005, then512
with5e-14/5e-16/.00025. Retain the same actual rho_up, canonical initial
columns, source covariance/weights, alpha=.001 history, .4 computational
length and21 target nodes. The budget is two solves and30 CPU seconds.
Do not repeat old PG or KS cases and do not automatically add a third case.

The256/512 and512/tighter differences are tested against1e-11 as numerical
indicators. They do not certify continuum error, between-node residuals or
the full spectrum. Current matter and its tangent use actual envelope
derivatives; source energies are not substituted as outgoing momenta.
No drift is removed, no source is changed, and the physical gate stays OPEN.

```sh
python3 scripts/derive_nsc_ks_envelope_resolution.py --run
python3 scripts/derive_nsc_ks_envelope_resolution.py --check
```
