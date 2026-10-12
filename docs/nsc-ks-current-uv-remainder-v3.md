# Current finite-history UV remainder v3

Mac integration 2026-09-28 separates v3 report functions into
`nsc_ks_current_uv_report.py`, preserving the v2 remainder owner byte-for-byte.
The recorder writes `nsc-ks-current-uv-remainder-v3-rebound.json`; the original
worker v3 JSON/NPZ remain unchanged. No new physical coefficient is supplied by
this provenance repair. Both v2 and the rebound v3 must replay successfully.


The v2 remainder successor named A2, A3, n4, C4, C_M, current-history
`L0 A4` H2 integrals, the upstream higher-order H2 remainder, transported
majors and the vacuum N/beta tail, and bound those names to the switched
state law, live profile `0b0e4ced…`, horizon preparation, source inventory
and KS subtraction. Its production A2/A3/n4 transports were not constructed.
The massless vanishing of the local `Im(T_s` major `A2)` right-hand side was
not allowed to set the upstream datum or the full-envelope coefficient to
zero. Historical v1/v2 records remain immutable.

## Owned characteristic transports

The v3 owner is `src/recursive_horizons/nsc_ks_current_uv_transport.py`,
recorded by `scripts/derive_nsc_ks_current_uv_remainder_v3.py` as
`results/development/nsc-ks-current-uv-remainder-v3.json` with payload
`results/development/artifacts/nsc-ks-current-uv-remainder-v3.npz`.

On the actual Dirac operator `L0`, with `T_s = ∂_ρ - s a^{-2} ∂_z`,

\[
\operatorname{Im}(T_s A_{2,s})
=-\frac{m\ell}{4}\frac{s a^2 r_ρ + r_z}{r^2},
\]

\[
\operatorname{Re}(T_s A_{2,s})
=\frac{-a_ρ m^2 a r^3
+\ell^2(-a_ρ a r + r_ρ a^2 + s r_z + 2 q r)
+2 m^2 q r^3}{4 r^3}.
\]

`Re(T_s` major `A2)` retains the transported major-A1 phase `q`. Dummy
`L0 A_j` symbols are not used. Minor A3 retains major A2. For `m=0`,
`Re(T_s` major `A3)=\ell^2 h/(2 r^2)` and `T_s n4` is independent of the
unowned upstream imaginary datum `h`.

The two histories share the unowned upstream imaginary major-A2 datum and,
when `m=0`, the local imaginary transport right-hand side vanishes. The
history-minus-reference field is therefore identically zero. The
full-envelope coefficient and the upstream datum stay `None`.

Interval majorants, not samples or fitted decay, bound

- `|q|` from `|T_s q|=(m^2+\ell^2/r^2)/2` and the owned radius/chart slab,
- `|δq|` from `|r_g^{-2}-r_{\rm ref}^{-2}|`,
- `|δ \operatorname{Re} A_{2,s}|` from the closed real parent formula and
  duration times `|F_g|+|F_{\rm ref}|`.

## History-minus-reference C4

Production `C4` is the equal-`μ` original-pair contraction of the `e^{-3}`
N and beta action densities of `g` minus the same densities of the
reference, on the declared interval `I`, with action minus once and measure
`de/(2π)`. It is not the generic full-envelope cubic coefficient of `g`.
Massless `J_3` drops the unowned `h`. The numerical value remains `None`
until minor A3 `= a^2/(2i) Π_{-s} L0 A2` is bounded on the slab.

## C_M and the vacuum tail

The owned `B_z` and `B_{zz}` commutator integrals are unchanged. Connecting
them to a finite `C_M` still requires production current-history `L0 A4` H2
integrals and the same-column upstream H2 remainder above order four.
Manufactured periodic and small-matrix H2 fixtures cannot supply those
inputs. No frozen-background tail certificate is transferred. `C_M` and the
vacuum N/beta tail remain `None`. Gate 2 stays OPEN.

## Next primitive

Interval or analytic bound of minor A3 `= a^2/(2i) Π_{-s} L0 A2` on the
current-history slab (second geometry jets of `r` and of the transported
phase `q`, together with `Re(major A2)`); then `T_s` major A3 and `T_s n4`;
then evaluation of the history-minus-reference `e^{-3}` contraction on `I`.
After that, production `L0 A4` H2 integrals and the same-column upstream H2
remainder for a finite `C_M` and vacuum tail.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 /Users/admin/Projects/BlackHoles-Infinity/.venv/validation/bin/python -m pytest -q tests/test_nsc_ks_finite_history_uv_remainder.py tests/test_nsc_ks_current_uv_remainder_v2.py tests/test_nsc_ks_current_uv_remainder_v3.py
python3 scripts/derive_nsc_ks_current_uv_remainder_v3.py --check
```
