# First validated lower-band source replacement

The [validated group32 window](nsc-incoming-validated-energy-window.md)
replaces its original LOW contribution on `[24,32]` once:

\[
T_{\rm v5}=T_{\rm v4}-T_{\rm old,32,[24,32]}
                         +T_{\rm validated\ vacuum,32,[24,32]}.
\]

Its density correction is `-3.709381090210057e-11`; the complete window
lapse-error bound is `1.8628405254514626e-16`. The physical thermal remainder
is positive and bounded separately. The true vacuum current is zero by
the exact projector identity; raw polynomial current and norm defects remain
in the original window receipt. No endpoint column is normalized.

The [new composition record](../results/development/nsc-incoming-source-update-v5.json)
authenticates the window and previous source, verifies that the new cell
meets the old coverage at32, and propagates its explicit source increment
into the existing joint-constraint approximants. Group32 is now covered
from24 through infinity; other groups retain their previous16/32 lower
thresholds. Earlier high regions, tails, LLL and group13 subgap data are
preserved. Incremental source-sum rounding is added separately.

The combined covered lapse-source error remains below `8.917e-12` against
`3e-11`. Lower-energy/subgap errors and full numerical source accuracy are
still unresolved. This partial source budget does not certify the complete
joint action arithmetic or its constraint residual. Physical initial data
and extended stationarity remain OPEN.

Reuse decision: compose the completed field/integral certificate and its
stored delta; no successful preflight, point-energy propagation, full window,
high source or tail calculation is repeated. Stop at authenticated
disjoint-region accounting and source-to-action replay.

```sh
python3 scripts/derive_nsc_incoming_source_update_v5.py --check
python3 -m pytest -q tests/test_nsc_incoming_source_update_v5.py
```
