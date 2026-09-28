# Current-history midpoint-jet profile bound, v3

This owner replaces the v2 global interval-Clenshaw L1 of the fourteen
compact profiles `w^p U^q` by a midpoint Chebyshev jet with endpoint
remainder. Production evaluation, the hash-bound module
`nsc_ks_current_history_bounds.py`, and the v2 adapter are unchanged.
There is no Dirac or source evolution. A one-cell bound is not a field
certificate and does not close the local incoming gate.

## Declared proof parameters

These were printed before both runs.

| Parameter | Value |
|---|---|
| Comparison CPU cap | 60 s |
| One-cell CPU cap | 180 s |
| Bits | 120 |
| Derivative order | 8 |
| Transition panels | 16 |
| Interior panels | 4 |
| Fourier `M`, `K` | 1024, 64 |
| Spatial period | `1759218604441/4398046511104` (actual grid spacing, not nominal `0.4`) |
| Profiles | actual 32-coefficient `(w,U)`, same analytic geometry as `one_direction_metric` |
| `free_contains_zero` as proof | not used |

The v2 record `695f2eb3` enclosed history `0b0e4ced` family `14_1` cell 122
with the same Fourier counts. Old `B_8` values are taken from that
immutable record. New `B_8` is computed live.

## Method

Chebyshev polynomials of the first kind obey DLMF 18.9.1,

\[
T_0=1,\quad T_1=x,\quad T_{n+1}(x)=2x T_n(x)-T_{n-1}(x).
\]

Differentiating \(j\) times gives

\[
T_{n+1}^{(j)}(x)=2j\,T_n^{(j-1)}(x)+2x\,T_n^{(j)}(x)-T_{n-1}^{(j)}(x).
\]

For \(M\ge 1\) and \(|x|\le M\), \(|T_n^{(j)}(x)|\le T_n^{(j)}(M)\). On
\([-1,1]\) this is the Chebyshev endpoint bound (DLMF 18.9.21:
\(T_n'=n U_{n-1}\), \(U_{n-1}(1)=n\), so \(T_n'(1)=n^2\)). For \(x\ge 1\),
\(T_n^{(j)}\) is positive and increasing. Ultraspherical derivatives
(DLMF 18.9.19) with \(U_n=C_n^{(1)}\) recover the first-derivative
identity. The axial map is the exact binary `mapparms()` dyadic. Then
\(P^{(j)}(z)=\mathrm{scale}^j\sum a_n T_n^{(j)}(x(z))\). Expanding a
high-degree profile into monomials and summing absolute power
coefficients is a different, looser bound; that is why the power-sum
owner is not used here.

On an interval of midpoint \(c\) and radius \(r\), two valid enclosures
of \(P^{(j)}\) are intersected:

1. Mean-value remainder: \(P^{(j)}(c)+[-r,r]B_{j+1}\), with \(B_{j+1}\)
   the global endpoint bound of \(|P^{(j+1)}|\).
2. Exact finite Taylor of a degree-\(N\) polynomial,
   \(P^{(j)}(c+h)=\sum_{m=0}^{N-j} P^{(j+m)}(c) h^m/m!\), together with
   truncated Taylor plus Lagrange remainder \(r^k/k!\,B_{j+k}\).

Midpoint Clenshaw has no interval wrapping. Broad-interval Clenshaw of
\(T_{31}\) is a valid natural interval extension, but dependency
explosion made v2's \(U^4\) eighth-derivative L1 \(3.12\times 10^{52}\).
That looseness is not an analytic non-existence obstruction.

The owned `plateau_series` cutoff jet is multiplied in the
\(f^{(n)}/n!\) convention. Interior cells have \(\chi\equiv 1\). Derivative
factorials and the declared inner/outer domains are preserved. Tests
check independent numpy/mpmath polynomial derivatives, interval
containment of samples, mixed Leibniz products, and that the new jet is
strictly tighter than interval Clenshaw on a broad interior ball.

## Comparison at unchanged Fourier counts

Measured comparison CPU \(0.639\) s, cap \(60\), \(56\) derivative cells.

| Key | old \(B_8\) | new \(B_8\) | ratio | old tail \(T_0\) | new tail \(T_0\) |
|---|---:|---:|---:|---:|---:|
| \(w\) | \(6.156\times 10^{29}\) | \(1.135\times 10^{25}\) | \(1.84\times 10^{-5}\) | \(2.70\times 10^{7}\) | \(4.97\times 10^{2}\) |
| \(U\) | \(1.683\times 10^{35}\) | \(2.489\times 10^{29}\) | \(1.48\times 10^{-6}\) | \(7.38\times 10^{12}\) | \(1.09\times 10^{7}\) |
| \(w^2\) | \(2.090\times 10^{30}\) | \(9.894\times 10^{21}\) | \(4.73\times 10^{-9}\) | \(9.16\times 10^{7}\) | \(4.34\times 10^{-1}\) |
| \(wU\) | \(5.699\times 10^{35}\) | \(1.361\times 10^{26}\) | \(2.39\times 10^{-10}\) | \(2.50\times 10^{13}\) | \(5.96\times 10^{3}\) |
| \(U^2\) | \(1.554\times 10^{41}\) | \(1.913\times 10^{30}\) | \(1.23\times 10^{-11}\) | \(6.81\times 10^{18}\) | \(8.38\times 10^{7}\) |
| \(w^3\) | \(4.508\times 10^{30}\) | \(1.349\times 10^{19}\) | \(2.99\times 10^{-12}\) | \(1.98\times 10^{8}\) | \(5.91\times 10^{-4}\) |
| \(w^2U\) | \(1.229\times 10^{36}\) | \(1.819\times 10^{23}\) | \(1.48\times 10^{-13}\) | \(5.38\times 10^{13}\) | \(7.97\) |
| \(wU^2\) | \(3.348\times 10^{41}\) | \(2.468\times 10^{27}\) | \(7.37\times 10^{-15}\) | \(1.47\times 10^{19}\) | \(1.08\times 10^{5}\) |
| \(U^3\) | \(9.124\times 10^{46}\) | \(3.370\times 10^{31}\) | \(3.69\times 10^{-16}\) | \(4.00\times 10^{24}\) | \(1.48\times 10^{9}\) |
| \(w^4\) | \(5.708\times 10^{30}\) | \(1.663\times 10^{16}\) | \(2.91\times 10^{-15}\) | \(2.50\times 10^{8}\) | \(7.29\times 10^{-7}\) |
| \(w^3U\) | \(1.553\times 10^{36}\) | \(2.227\times 10^{20}\) | \(1.44\times 10^{-16}\) | \(6.80\times 10^{13}\) | \(9.76\times 10^{-3}\) |
| \(w^2U^2\) | \(4.223\times 10^{41}\) | \(3.000\times 10^{24}\) | \(7.10\times 10^{-18}\) | \(1.85\times 10^{19}\) | \(1.31\times 10^{2}\) |
| \(wU^3\) | \(1.149\times 10^{47}\) | \(4.063\times 10^{28}\) | \(3.54\times 10^{-19}\) | \(5.03\times 10^{24}\) | \(1.78\times 10^{6}\) |
| \(U^4\) | \(3.124\times 10^{52}\) | \(5.537\times 10^{32}\) | \(1.77\times 10^{-20}\) | \(1.37\times 10^{30}\) | \(2.43\times 10^{10}\) |

The \(U^4\) eighth-derivative L1 dropped twenty orders of magnitude. This
is a usable named continuous-proof method. It is still a coarse bound:
interior-panel refinement does not move \(B_8\); remaining mass sits on
the cutoff ramps.

## One-cell residual

The comparison improved \(B_8\) by more than \(10^{-6}\), so cell 122 was
bounded with `DifferenceResidualPolynomial` and the new jet L1. Measured
cost \(121.85\) CPU seconds, \(126.93\) wall seconds, cap \(180\).
Radius lower bound \(1.3997614174467379\), \(\eta\approx 1.70\times 10^{-4}\).
Period \(0.3999999999998636=1759218604441/4398046511104\).

| Contribution on the whole \(\xi,z\) cell | order 0 | order 1 |
|---|---:|---:|
| Polynomial finite Fourier band | \(1.678\times 10^{-11}\) | \(9.053\times 10^{-9}\) |
| Time-Taylor and reciprocal-radius remainder | \(3.319\times 10^{-20}\) | \(3.894\times 10^{-17}\) |
| Polynomial plus omitted profile band | \(5.481\times 10^{-4}\) | \(6.429\times 10^{-1}\) |
| Total continuous normalized residual | \(5.481\times 10^{-4}\) | \(6.429\times 10^{-1}\) |

v2 totals were \(312.1\) and \(3.661\times 10^{5}\). The finite band and
time/radius remainder also shrank, because v2's alias balls and profile
suprema inherited the huge global L1. Dominating contribution:
omitted profile Fourier band of Cheb32 \(U\)-powers. Named obstruction
`cheb32_mixed_profile_fourier_omitted_band_from_midpoint_jet_l1`.

`free_contains_zero` is not recorded as a proof. Interval containment of
zero is not an identity certificate.

## Scope

One saved cell. Unknown global field error and unknown source error.
Not a 246-cell certificate, not EXISTENCE, not scoped NON-EXISTENCE.
The physical local gate remains OPEN.

Named next gap:
`transition_cutoff_jet_L1_still_dominates_omitted_Fourier_band;
one_cell_unknown_global_field_and_source_error`.

A local (not global-L1) Fourier remainder on the cutoff ramps, or a
sharper cutoff jet, is the next bound improvement. Raising `K` at this
\(B_8\) does not put the omitted band under the finite residual.

## Reproduction

```sh
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v3.py --pilot
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v3.py --run
PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v3.py --check
PYTHONPATH=src .venv/validation/bin/python -m pytest -q tests/test_nsc_ks_chebyshev_jet_bound.py
```

Versioned records:
`results/development/nsc-ks-current-field-pilot-v3-profile-bounds.json`,
`results/development/nsc-ks-current-field-pilot-v3.json`.
