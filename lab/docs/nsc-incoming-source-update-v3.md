# Complete finite high-band source replacement

The high source is replaced once by the directly integrated order24 vacuum:

\[
T_{\rm new}=T_{\rm original}-T_{\rm original,high}
             +T_{24,\rm high}^{\rm vacuum}.
\]

The original source includes the existing group13 subgap replacement, paired
tail approximation and one physical LLL state. Their values are retained;
previous high-band corrections are not added again. Physical thermal stress
is a nonzero bounded remainder, including its current. The zero current of
the high vacuum approximation does not set the physical current to zero.

The [record](../results/development/nsc-incoming-source-update-v3.json) joins
the [direct numerical integrals](nsc-incoming-high-band-quadrature-batch.md)
with the [remaining28 physical certificates](nsc-incoming-high-band-batch-bound.md)
and existing groups12/14/22/32 certificates. It authenticates every referenced
producer and payload. Groups13/14 cover `[16,160]`; the other groups cover
`[32,160]`, extended to320 for groups10/11/12/31/32. Adjacent physical pieces
partition the same direct integration interval; their union is counted once.

The sum of physical projector, thermal, source quadrature and arithmetic
bounds closes the finite-high lapse/shift source budget against `3e-11`.
Its lapse bound is `8.90304e-12`: `5.35558e-12` physical mode/thermal plus
`3.54746e-12` numerical integration. Rounding when embedding this replacement
in the original recorded source sum is bounded separately; the other fixed
pieces retain their unresolved physical errors. The baseline raw action
coefficients become approximately `(-0.69593057250,-0.05712098699)`.
The accompanying pressure integrals have numerical enclosures, but the
energy-specific physical certificate does not certify their physical errors.
No radial recurrence, quadrature producer or old control is rerun by this join.

The infinite-tail physical projector certificate remains valid. Its source
integral uses stored finite-precision Taylor coefficients and numerical
convergence indicators; those do not constitute a rigorous numerical integral
bound. The record therefore keeps tail numerical error separate and null.
Lower-energy and subgap physical/numerical errors are also unresolved.
The full source and constraint residual remain OPEN. The local constraint
plane and positively enclosed lapse coefficient are unchanged; no initial
data, physical root, endpoint solution or metric timestep is selected.

Reproduce only the saved-data composition:

```sh
python3 scripts/derive_nsc_incoming_source_update_v3.py --check
python3 -m pytest -q tests/test_nsc_incoming_source_update_v3.py
```
