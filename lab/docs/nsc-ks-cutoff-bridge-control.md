# Saved all-family control with the derived source-cutoff edge

This control joins the [analytic kernel bridge](nsc-source-cutoff-bridge.md)
to the saved raw response of all60 nonzero-angular source families and their
120 explicit signed contributions. It preserves each original family,
source-panel record, upstream slice and raw result. No Dirac or source
preparation run is repeated.

The phase integral and its retarded derivative use the same saved analytic
radius history, incoming nodes and upstream coordinate. The common
unit-angular coefficient is multiplied by the original per-family ledger
sum of `multiplicity*ell²`, separately for each energy sign. Both angular
signs are retained and the zero-angular baseline allocations remain once.

Two geometry quadratures, with96 and192 Gauss nodes on each normal-window
panel, measure a convergence indicator. A centered change of the original
history amplitude checks the weighted edge derivative. Neither comparison
is a rigorous quadrature or finite-source-tail bound.

The record retains the raw and joined gradients and Jacobians separately.
The joined value still omits an explicitly unbounded finite-cutoff raw tail,
as well as unresolved source, field and between-node errors. The control
history is not a root and this is not a local EXISTENCE or NON-EXISTENCE
certificate. The edge follows from the same source/action prescription;
no local/reference coefficient, source law or `Gamma_rest` is changed.

```sh
python scripts/derive_nsc_ks_cutoff_bridge_control.py --check
```
