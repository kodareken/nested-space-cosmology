# Original-source control of the coupled history evaluator

The bounded control restores the first16 original positive-energy rows of
`group14/low16_1` and their explicit negative-energy partner. Their saved
upstream columns and preparation digests are reused. No scattering or source
preparation is repeated and no occupation, phase, weight or coupling is fit.

The two local functions have eight Chebyshev coefficients each. The central
history keeps the saved nonzero `w_1=0.001` control, with all other
coefficients zero. A finite difference varies a declared combination of
`w_0,w_1,U_0,U_1`. Every candidate propagates all16 retarded directions
through the actual KS operator before assembling both constraints.

Eight solve nodes and17 independent verification nodes share one fixed
target union. Exact subsets and signed partners reuse each candidate's
operator. Baseline and geometry enter once. Three candidate evolutions on
the64-point numerical cell, with an eighth-degree energy interpolant on
`[0,0.25]`, are limited to30 CPU seconds in total.

The finite-difference comparison is a derivative **control**, with tolerance
`3e-8`. It is not a field-accuracy bound and this coarse control is not a
candidate physical root. Source/low-energy accuracy, interpolation/field
errors, the source-cutoff/coincidence bridge and between-node residuals are
not certified here. All missing physical error bounds remain explicit.

Saved operators reconstruct the original columns, derivatives and retarded
responses during replay. The verifier repeats their assembly and the
derivative comparison without rerunning Dirac evolution.

```sh
python scripts/derive_nsc_ks_history_evaluator_control.py --check
```
