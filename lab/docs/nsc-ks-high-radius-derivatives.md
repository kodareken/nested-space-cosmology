# Spatial supremum and Fourier-l1 bounds for the reciprocal-radius difference

On the existing pure-radius nonzero control,

\[
q=\frac1{r_{\rm ref}+\delta r}-\frac1{r_{\rm ref}}.
\]

This owner supplies directed **spatial/axial** derivative uppers of \(q\)
through order \(9\), and the stronger **Fourier-\(\ell^1\) moments**
\(A_j(q)\) used by exact Weyl-shift Taylor remainders. It reuses the
committed \(p=16\) fine Fourier profile table and the two-derivative radius
enclosure. It does not evolve a field or source, fit a profile, or change
the action. Mixed \(T/z\) jets through order \(9\) are a separate root owner.
This file remains spatial supremum and Fourier-\(\ell^1\) only.

A covariance UV bound may use auxiliary \(P\) in \(C=P+D\) without a
band-action change. Generating \(\Gamma_{\rm sub}\) is a different use.
Normal-time integrals, high derivatives of the fourth-order symbol, and the
exact Moyal remainder remain missing. This is not a UV or physical-gate
certificate.

## Directed denominators

Norm uppers (\(W_j\), \(|c_k|\), \(\sigma\), support width) stay directed
upward. Period length and radius lowers are **positive enclosing intervals**,
never collapsed `.upper()` records used as denominators. Collapsed
upper-only dyadic records are rejected.

For exact \(L\), \(\omega=2\pi/L\) is formed from the enclosing interval, then
the resulting moment is taken upward. For a geometry lower \(r_\star\),
\(1/r_\star\) uses the enclosure (equivalently a directed lower) so that a
low-bit evaluation cannot undercover the exact reciprocal. Strict tests at
working precision \(8,12,30\) on \(L=7/6\) show that \(1/L^{\mathrm{upper}}\)
fails to cover \(6/7\), while the enclosing reciprocal covers it as an exact
dyadic versus the exact rational. Float tolerances are not used.

## Reused inputs

| Input | Owner |
|---|---|
| Alias and omitted-band tails, \(j<p-1\) | `nsc_ks_profile_fourier_bound.alias_and_tail_bounds` |
| Certified coefficient balls of physical \(W\) | [fine profile table](nsc-ks-fine-profile-table.md), payload `fc20111d0c72f1ac81d56ee5ea64cdd1bb94a9bfc6e86edfa19d1f8992f0ce33` |
| Analytic profile fingerprint | `a4c061f297027bcb9934fa2338bb81fb799ca4ec26a025f6bc563384c70cd69d` |
| Radius lower, \(r_{\rm ref}\) lower, and \(W_j\) through order \(2\) | [radius enclosure](nsc-ks-radius-enclosure.md) |
| Directed dyadic packing | `nsc_ks_ball_trajectory.exact_upper` |

Physical \(W\) already contains \(\alpha=0.001\). That factor is not applied
again. \(U\equiv 0\). The period remains \(205/512\). Compact axial support
remains width \(2\times\texttt{outer}\) (binary image of \(0.12\)).

## Fourier-\(\ell^1\) moments and spatial suprema

\[
A_j(f)
=\sum_{|k|\le K}|\omega_k|^j\,|c_k|
+T_j,\qquad
\omega_k=\frac{2\pi k}{L}.
\]

\(A_j\) is the Wiener-algebra moment of \(\partial_z^j f\) and also a spatial
supremum majorant. The same number is reported as both. It is the input to
sums \(\sum_n|\omega/2|^n\|V_n\|\) in an exact finite-history shift remainder.

For the reciprocal, \(A_0(1/r)\) is **not** inferred from a pointwise
\(r_{\min}\) alone. With \(d_0=A_0(\delta r)\) and \(d_0<r_{\mathrm{ref,min}}\),
Neumann in the Wiener algebra gives

\[
A_0(1/r)\le\frac1{r_{\mathrm{ref,min}}-d_0}.
\]

Higher moments obey the same Leibniz recursion with this inverse \(A_0\)
bound. The difference satisfies

\[
A_0(q)\le\frac{d_0}{r_{\mathrm{ref,min}}(r_{\mathrm{ref,min}}-d_0)},
\qquad
A_n(q)=A_n(1/r)\quad(n\ge 1).
\]

Pointwise bounds remain separate: \(q_0\le\delta r_{\max}/(r_{\mathrm{ref,min}} r_{\min})\)
and \(q_{{\rm full},0}\le 1/r_{\min}\), with \(L^1\le w\sup\) and
\(L^2{}^2\le w\sup^2\) on the owned axial support. Those are not Fourier
moments.

On this control the recorded pointwise \(r_{\min}\) is a smaller lower bound
than the Neumann margin \(r_{\mathrm{ref}}-d_0\), so it happens to be
conservative for \(A_0(1/r)\). That comparison is checked, not assumed in
general. The Fourier moments use Neumann.

## Recorded constants

| Constant | Exact | Display |
|---|---|---|
| Period \(L\) | \(205/512\) | \(0.400390625\) |
| Normal support \(\sigma\) | \(1080863910568919/36028797018963968\) | \(0.03\) |
| Axial support width \(w\) | \(1080863910568919/9007199254740992\) | \(0.12\) |
| \(r_{\rm ref,min}\) | \(7/5\) | \(1.4\) |
| Pointwise \(r_{\min}\) | radius-record rational | \(1.39997\) |
| Neumann margin \(r_{\rm ref}-d_0\) | dyadic lower | \(1.399975605153736\) |
| \(d_0=A_0(\delta r)\) | dyadic upper | \(2.439484626404717\times10^{-5}\) |
| Amplitude in \(W\) | binary64 of \(0.001\) | \(0.001\) |
| Fourier settings | \(p=16\), \(M=65536\), \(K=8192\) | table bits \(120\), sum bits \(160\) |

Uppers use `kind=dyadic` (`mantissa`,`exponent`) or `kind=rational`.
Lowers use `kind=dyadic_lower` or `kind=rational` with `binary64_lower`.
`binary64_*` fields are display only.

## Distinct pointwise and Fourier-\(\ell^1\) bounds of \(q\)

| \(j\) | pointwise \(\sup\|q^{(j)}\|\) | Fourier \(A_j(q)\) | pointwise \(L^1\) |
|---:|---:|---:|---:|
| 0 | \(1.244661684793611\times10^{-5}\) | \(1.244656701478140\times10^{-5}\) | \(1.493594021752334\times10^{-6}\) |
| 1 | \(1.176704425440322\times10^{-3}\) | \(1.176695002996133\times10^{-3}\) | \(1.412045310528386\times10^{-4}\) |
| 2 | \(2.363285495716073\times10^{-1}\) | \(2.363266571585380\times10^{-1}\) | \(2.835942594859288\times10^{-2}\) |
| 3 | \(9.753104470024994\times10^{1}\) | \(9.753026371242463\times10^{1}\) | \(1.170372536402999\times10^{1}\) |
| 4 | \(7.393962625063174\times10^{4}\) | \(7.393903417308674\times10^{4}\) | \(8.872755150075809\times10^{3}\) |
| 5 | \(9.191897545178308\times10^{7}\) | \(9.191823940441528\times10^{7}\) | \(1.103027705421397\times10^{7}\) |
| 6 | \(1.739708546190926\times10^{11}\) | \(1.739694615386337\times10^{11}\) | \(2.087650255429111\times10^{10}\) |
| 7 | \(4.727894771120511\times10^{14}\) | \(4.727856912329344\times10^{14}\) | \(5.673473725344613\times10^{13}\) |
| 8 | \(1.795063412318621\times10^{18}\) | \(1.795049038302943\times10^{18}\) | \(2.154076094782345\times10^{17}\) |
| 9 | \(1.644655262255394\times10^{22}\) | \(1.644642092673062\times10^{22}\) | \(1.973586314706472\times10^{21}\) |

\(A_0(W)=8.131615421349057\times10^{-4}\) already includes \(\alpha\);
\(\delta r\) uses \(\sigma A_j(W)\), not \(\alpha\sigma A_j(W)\). For \(j=0,1,2\)
this lies below the Chebyshev enclosure \((0.001,0.283333,125.556)\). For
\(j\ge 1\), difference moments equal the full reciprocal moments.

## What is still missing for an exact UV residual

- mixed \(T/z\) remainder integrals (spatial/Fourier only here);
- high derivatives of the actual fourth-order reference symbol;
- the exact Moyal remainder;
- the Bloch trace inequality (separate job);
- the actual fourth-order reciprocal remainder (root review).

The physical local gate on \(I=S(1)+[.12,.18]\) remains OPEN.
