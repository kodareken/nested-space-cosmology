# Directed Taylor enclosure of the source-phase derivative

The [phase evaluator](nsc-dirac-source-phase-evaluation.md) integrates

\[
\frac{f_{s,z}}{\ell^2}
=-\int_1^{\rho_u}
\frac{r_z}{r^3}\bigl(\rho,\,z-s D(\rho)\bigr)\,d\rho,
\qquad
D(\rho)=\int_1^\rho a(\rho')^{-2}\,d\rho',
\]

with Gauss--Legendre panels. Those node differences are indicators. This
owner supplies a directed value ball for the same unit-angular integrand
at one incoming control point, for both original source signs.

The history is the saved one-direction family
`w=(0,1)`, `U=(0)`, amplitude `.001`, axial windows `.03/.06` and normal
windows `.007/.03`. Its analytic identity is checked against the
[cutoff-bridge control](nsc-ks-cutoff-bridge-control.md). The upstream
coordinate is that record's exact binary64 `rho_up`. The original ledger
weight per energy sign is `107952`. No source or field solve is repeated.

## Reused jets

`D` uses `acb.integral` of `inv_a_integrand**2` at the cell midpoint,
the Lipschitz bound `|D'|<=25/16` on a rho ball, and higher Taylor
coefficients from `BackgroundSeries.inv_a2`. The characteristic
`z_char=z-s D` is composed with an increment whose constant term is
exact `0`. Axial `w,U` and `w_z,U_z` are the owned
`AnalyticRadiusFamily` jets along that increment, with derivative order
at most 16. The actual radius and `r_z` use the same shift, plateau and
profiles as the geometry owner. Smooth-flat cutoffs raise
`SubdivisionNeeded` until a cell is classifiable.

On a dyadic cell of half-width `h`, degrees `0..13` of the center jet
are integrated (odd powers vanish by symmetry). The remainder is at most

\[
\frac{2 h^{15}}{15}\sup\bigl|a_{14}\bigr|,
\]

where `a_14` is coefficient 14 of the full-cell jet. Adaptive bisection
stops on that remainder, on a wall, or on the declared cell/depth/CPU
limits. Exact dyadic endpoints keep the left of the first cell at
`rho=1`; `rho_up-1` is decomposed into power-of-two lengths so balls do
not inflate below the incoming slice. Interval inverses retain lower
denominators.

The [continuum local embedding](nsc-ks-local-embedding.md) is reused. No
new global matching or metric step is introduced.

## One-point pilot

The recorded control is the leftmost incoming `z` only. The target is a
unit-angular radius `<=2e-16` per source sign, which is enough for an
aggregate `N,beta` phase-value error below `~1e-11` after the original
angular-square weight. Focused tests check vanishing constant/zero
profiles, containment of an independent high-precision sample integral,
chart identity, both characteristics, and context restoration. The
Taylor remainder is the proof; Gauss-node agreement is not.

Source-tail, field, between-node and other missing errors remain `None`.
The physical local gate stays `OPEN`.

```sh
PYTHONPATH=src /tmp/nsc-validated-runtime.8pg_7t1z/bin/python -m pytest -q tests/test_nsc_dirac_source_phase_bound.py
/tmp/nsc-validated-runtime.8pg_7t1z/bin/python scripts/derive_nsc_dirac_source_phase_bound_control.py --check
```
