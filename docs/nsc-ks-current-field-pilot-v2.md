# Current-history one-cell field-proof pilot

This adapter bounds the continuous KS residual of **one saved DOP853
cell** of the live Chebyshev-32 `(w,U)` history `0b0e4ced…`. It does not
transfer the old fine-history proof. That proof belongs to amplitude
`ALPHA=0.001`, `U=0`, 246 cells and four source energies; its profile
table owner hardcodes the four keys `(n,0)`. The current history has
both `w` and `U` nonzero. `AnalyticRadiusFamily.time_series` already
emits the fourteen monomials `w^p U^q` with `p+q<=4`.

Parent Codex owns the capture
`results/development/nsc-ks-current-trajectory-pilot-v2.json` and the
conditional module `nsc_ks_difference_error.py`. This adapter owns only
the difference residual polynomial, the one-cell validator, its tests
and this note. There is no Dirac or source evolution here.

## Declared cell

The capture is family `14_1`, four original energy rows times three
original coherent columns, grid 1024, numerical period 0.4, dense
segment index 122. The reconstruction is the production split
`D=X-A` with `A'=G_ref A` and `D'=L_g D+M A`. `A` is z-independent.

## Radius and reciprocal remainder

Power-coefficient absolute sums of the Chebyshev-32 `U` profile are of
order `3.9e8` and destroy the radius margin. The adapter uses
`nsc_ks_current_history_bounds.chebyshev_profile_bounds`, which retains
the represented basis. On this history that gives `radius_lower≈1.39976`
and `eta≈1.70e-4`. The reciprocal-radius remainder of order four is
taken from those Chebyshev norms, not from the power-sum owner.

## Polynomial split

`DifferenceResidualPolynomial` builds `R_X` with the full potential and
`R_A` with matching zero potential keys, then subtracts Bernstein
coefficients. That is `R_D=R_X-R_A` before any norm. Independent
controls with zero `D` and nonzero `A` check that the background
residual cancels and that only `M A` remains. Subtracting already
bounded `||R_X||` and `||R_A||` would hide that cancellation: both
norms stay positive while `R_D` is identically zero when `D=0` and the
potential is empty.

Time-Taylor and omitted-radius remainders of `L_g` on `D` and of `M`
on `A` are added. Background `G_ref` remainders belong to `R_A` and are
not included in `R_D`.

## Fourier profiles

The fourteen spatial keys are enclosed with the owned alias/tail
theorem. On Chebyshev-32 `U` the eighth-derivative `L^1` norms are
enormous, so omitted-band tails dominate any global Fourier remainder
at a retained index that fits a 180 CPU-second one-cell budget. That
is a named analytic obstruction,

`cheb32_mixed_profile_fourier_omitted_band_from_global_derivative_l1`,

not a physical NON-EXISTENCE result. The polynomial finite band and
the time/radius remainder are still recorded.

## Conditional characteristic and matter majorants

`propagate_difference_error` matches the production error equations
`e_A'=G_ref e_A+R_A` and `e_D'=L_g e_D+M e_A+R_D` with `e_{A,z}=0`.
The `F_z` split matches the envelope: reference axial error is
`E|e_A|`; difference axial error is `|e_{D,z}|+E|e_D|` for
`D_z-i E D`. `difference_matter_error` bounds the difference
contraction `D C A*+A C D*+D C D*` with the same vertex norms as
`finite_matter_error`, and it does not treat `A` and `D` as independent
fields.

It does **not** match a whole-field certificate on this cell:

- the integrals are conservative whole-interval Gronwall majorants
- `M e_A` uses a uniform `|e_A|` times `∫||M||`
- residuals, `M_integral` and source error are caller obligations
- one saved cell is not the backward cone of `I`

This pilot therefore does not call those majorants as a field proof.
Sampled DOP853 defects are not promoted to bounds.

## Recorded one-cell result

Cell 122 on `rho ∈ [1.0149853515625002, 1.0151074218750002]`, 12 source
columns, 14 reciprocal keys, Chebyshev `radius_lower≈1.3997614174467379`,
`eta≈1.7041610947297358e-4`. The rejected power-sum `U` bound is
`3.873e8`.

| Contribution on the whole `xi,z` cell | order 0 | order 1 |
|---|---:|---:|
| Polynomial finite Fourier band | `9.121e-6` | `4.620e-3` |
| Time-Taylor and reciprocal-radius remainder | `6.914e-14` | `8.109e-11` |
| Polynomial plus omitted profile band | `312.1` | `3.661e5` |
| Total continuous normalized residual | `312.1` | `3.661e5` |

The finite band and time/radius remainder are a numerical continuous
residual bound for this cell. The total is dominated by the global
eighth-derivative Fourier tail of `U^4` (`B_8≈3.12e52`, omitted-band
`≈1.37e30` at retained index 64). Named obstruction:
`cheb32_mixed_profile_fourier_omitted_band_from_global_derivative_l1`.

Measured cost: `97.46` CPU seconds, `103.5` wall seconds, cap 180.
Storage: capture payload `3163524` bytes, field record `52562` bytes,
cell proof `1343` bytes. Capture payload hash
`c51701d2304c8567…`, history `0b0e4ced…`, source digest `d48a07ada2749c6f…`.

The physical local gate remains OPEN. This is not EXISTENCE, not scoped
NON-EXISTENCE, and not a 246-cell or full-source certificate.

The versioned record is
`results/development/nsc-ks-current-field-pilot-v2.json`.

```sh
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v2.py --run
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v2.py --check
```
