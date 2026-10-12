# Uniform-in-k operator-norm bounds of the auxiliary projector

Root already owns the actual mixed coefficients of auxiliary `P_0..P_4`
through the regular `mu=1/|k|` jets of
[the scaled reference projector](nsc-scaled-reference-projector.md). Those
jets are the small-transfer owner. This owner supplies the missing
**all-real-`k`** operator-norm majorant, needed when a large axial Fourier
transfer evaluates the symbol at `k ± omega/2` and that argument crosses
`0`. The symbols remain auxiliary. They are not `C_{\mathrm{up}}` or
`C_\Sigma`, they do not replace the physical state, and they do not add a
`Gamma` term. The exact split remains `C_g=P_g+D_g`.

The bound is an operator-norm statement about Taylor coefficients of `P_j`,
not a reconstruction of the symbols themselves and not a physical gate.

## Result

On an owned [KS geometry box](nsc-ks-spacetime-geometry-jets.md) with
`a>0`, `r>0` and strictly positive rest gap

\[
c_{\min}=\inf\frac1a,\qquad
g_{\min}=\sqrt{m_{\min}^2+\bigl(\lvert\ell\rvert_{\min}/r_{\max}\bigr)^2}>0,
\]

every mixed Taylor coefficient of auxiliary `P_0,\ldots,P_4` obeys an
explicit operator-norm majorant that is **independent of the expansion
momentum** `k_0\in\mathbb R`. Zero rest gap is rejected; there is no
denominator floor. The `ell=0` pure-radius response is a different owner.
The third axis of the jet is ordinary `k`, not `mu`. Default orders are
physical `9` and `k` `5`, enough for `d_z^4 P_4`; order `j` is trimmed to
`(9-j,5-j)`.

Grid samples at finite `k` are sanity checks. Uniformity is the normalized
gap algebra below, not extrapolation from those samples.

## Hamiltonian and gap increment

Pure KS, `N=1`, `beta=0`:

\[
H=c\,k\,\Sigma_3+V,\qquad
c=\frac1{a(T)},\qquad
V=-m\Sigma_1+\frac{\ell}{r}\Sigma_2,
\]

\[
g(k)^2=c^2k^2+v^2,\qquad
v^2=m^2+\ell^2/r^2.
\]

At any expansion basepoint, `g_0=g(k_0)` is a scalar constant in the Taylor
algebra. The increment `Delta=g^2/g_0^2-1` is majorized by a nonnegative
jet `X` with `X_{000}=0`:

| `k`-order | physical multiindex | `X` majorant |
|---|---|---|
| `0` | `alpha \neq 0` | `\|(c^2)_alpha\|/c_{\min}^2 + \|(v^2)_alpha\|/g_{\min}^2` |
| `1` | any `alpha` | `\|(c^2)_alpha\|/(c_{\min} g_{\min})` |
| `2` | any `alpha` | `\|(c^2)_alpha\|/g_{\min}^2` |
| `>=3` | any | `0` |

The `k`-order-`1` line uses AM-GM `|k|/g^2\le 1/(2cv)`. Times the `2` from
`d_k(k^2)` this is `|c^2_alpha|/(c_{\min}g_{\min})`, not twice that. The
`k`-order-`2` line is already a Taylor coefficient (`2!` removed). Strict
interval denominators take **lowers** of `c` and of the rest gap; an upper
is never rounded and then inverted.

Because `X_{000}=0`, the finite-order formal compositions

\[
\frac{g_0}g\;\prec\;(1-X)^{-1/2},\qquad
\frac{g_0^2}{g^2}\;\prec\;(1-X)^{-1}
\]

bound derivatives at each point without asserting convergence on a disk.
The reused [scalar mixed jet](nsc-scaled-reference-projector.md) performs
those compositions.

## `H/g_0` and `P_0`

Pauli orthogonality gives `\|H_0\|_{\mathrm{op}}=g_0`, so the majorant `A`
of `H/g_0` has `A_{000}=1` exactly. For `alpha\neq0`,

\[
A_{alpha,0}\le\frac{|c_alpha|}{c_{\min}}+\frac{\|V_alpha\|_{\mathrm{op}}}{g_{\min}},
\qquad
A_{alpha,1}\le\frac{|c_alpha|}{g_{\min}},
\]

and higher `k` vanish. The potential uses

\[
\|V_alpha\|_{\mathrm{op}}
=\sqrt{(m\,\delta_{alpha,0})^2+\bigl(\ell\,(1/r)_alpha\bigr)^2}
\]

from coefficient absolute uppers; there is no trace-one projection on
`P_j` for `j\ge1`.

\[
P_0=\frac I2-\frac H{2g}
\quad\Rightarrow\quad
\|P_0\|\;\le\;
\tfrac12 A(1-X)^{-1/2}+\tfrac12.
\]

The constant term is `1`. `Q=I-P_0` has the same majorant (`\|Q_0\|_{\mathrm{op}}=1`,
higher derivatives `-P_0`).

## PST recurrence majorants

The recurrence is
[Panati–Spohn–Teufel, sections 2 and 4.4](https://arxiv.org/abs/math-ph/0201055)
as already executed by
[the spatial reference symbol](nsc-spatial-reference-symbol.md). Each star
coefficient is replaced by its absolute value
`1/(2^l l!)\binom{l}{h}` on the mixed `(z,k)` split; both factors receive
ordinary `k` derivatives.

\[
G_n
=\sum_{i,j<n,\,l=n-i-j\ge0}
\operatorname{star}_l(P_i,P_j).
\]

The transport defect, after dropping kinetic Weyl orders `>=2` (`a=a(T)`
has no `z` derivative), is

\[
\|F_n\|
\le
\|d_T P_{n-1}\|
+\|C\|\,\|d_z P_{n-1}\|
+\sum_{j=0}^{n-1}\frac{2}{2^l l!}\,\|d_z^l V\|\,\|d_k^l P_j\|,
\quad l=n-j.
\]

`C` is the absolute Taylor jet of `c`; its `z` derivatives vanish. Off-diagonal
and diagonal pieces of `P_n` use the matrix commutator factor `2` and the
same `P_0` majorant for `Q`:

\[
\|P_n^{\mathrm{off}}\|
\le
\frac12 A F_n(1-X)^{-1}/g_{\min},
\qquad
\|P_n^{\mathrm{diag}}\|
\le
2\,P_0\,G_n\,P_0.
\]

All mixed product coefficients and factorials are kept. The output is an
upper bound on **operator norms** of Taylor coefficients; the actual
derivative bound is that coefficient times `t!\,z!\,k!`. Frobenius and
nuclear conversions on `C^2` are the explicit factors `\sqrt{2}` and `2`.
They are not applied to the majorant.

Constant geometry makes `P_j=0` for `j\ge1` identically in `k`, while `P_0`
retains a nontrivial `k` jet.

```sh
PYTHONPATH=src /tmp/nsc-validated-runtime.8pg_7t1z/bin/python -m pytest -q tests/test_nsc_global_projector_bounds.py
```

## Remaining UV inputs

These bounds close the all-`k` coefficient envelope of auxiliary `P`. They
do not close the finite-history remainder. Still required, and not computed
here, are the exact Weyl-defect seminorms of
[the shift remainder](nsc-weyl-commutator-remainder.md), spatial radius
norms from the existing enclosure, nuclear-norm propagation of `D=C-P`
(using the explicit `sqrt(2)` or `2` conversion if a Schatten norm is
chosen), the initial `C_{\mathrm{up}}-P_{\mathrm{up}}` bound, and the
changed-history source UV remainder. The local incoming gate stays OPEN.
No physical `C` is identified with auxiliary `P`.
