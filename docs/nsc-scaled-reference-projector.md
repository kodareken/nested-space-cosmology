# Regular high-momentum jets of the auxiliary fourth-order projector

The [spatial PST owner](nsc-spatial-reference-symbol.md) already solves the
Weyl-product projector of
[Panati, Spohn and Teufel, sections 2 and 4.4](https://arxiv.org/abs/math-ph/0201055)
on the canonical half-density Hamiltonian. Its `SymbolJet` tables are
hardcoded to total physical order four, and `P_4` is trimmed to a value.
Substituting `k=1/mu` on an interval that meets `mu=0` therefore cancels
`P_0` against `P_\infty` and loses the true decay.

This owner keeps that same recurrence and the same pure KS history
`N=1`, `beta=0`, `a=a(T)`, `r=r(T,z)` from
[coordinate-time geometry jets](nsc-ks-spacetime-geometry-jets.md).
The new work is the regular representation

\[
\mu=\frac1{|k|},\qquad e=\mathrm{sign}(k),\qquad
P_0=P_\infty+\mu B_0,\qquad
P_j=\mu^{j+1}B_j\quad(j\ge1),
\]

with each `B_j` finite at `mu=0`. The symbols are auxiliary. They are not
`C_{\mathrm{up}}` or `C_\Sigma`, they do not replace the physical state, and
they do not add a `Gamma` term. The exact split remains
`C_g=P_g+D_g`. `G_rho=i H_T/a` is not formed here; the clock is coordinate
`T`.

## Regular leading term

With `C=Sigma_3/a`, `V=-m Sigma_1+(ell/r) Sigma_2`, `b=a V`,
`b_2=a^2(m^2+ell^2/r^2)` and `s=sqrt(1+mu^2 b_2)`,

\[
P_\infty=\frac{I-e\,Sigma_3}2,\qquad
B_0=-\frac{b}{2s}+e\,Sigma_3\frac{mu\,b_2}{2s(s+1)}.
\]

`P_0=P_\infty+mu B_0` is not a difference of two near-equal projectors.
The formulas were checked against the unscaled PST identities; no sign,
power or factorial correction was required.

Canonical `k` derivatives use `d_k=-e mu^2 d_mu`,

\[
d_k^h(mu^p B)=mu^{p+h}D[p,h]B,
\qquad
D[p,0]B=B,\qquad
D[p,h+1]B=-e\bigl((p+h)D[p,h]B+mu\,\partial_\mu D[p,h]B\bigr).
\]

Any positive `T`, `z` or `k` derivative kills `P_\infty`, so `P_0` uses
`p=1` and `B_0`. For `j>=1`, `p=j+1`.

## Scaled defects

Idempotence and transport defects keep the owned sums, with every star
factor written through `D[i+1,h]B_i`. Then

\[
G_n=mu^{n+2}\widetilde G_n,\qquad
F_n=mu^n\widetilde F_n.
\]

Kinetic Weyl orders `>=2` vanish because `a` has no `z` dependence. The
off-diagonal reconstruction uses `Hbar=e C+mu V` and
`gapbar^2=a^{-2}+mu^2(m^2+ell^2/r^2)`,

\[
B_n=\frac{[\bar H,\widetilde F_n]}{4\,\overline{\mathrm{gap}}^2}
+mu\bigl[-P_0\widetilde G_n P_0+(I-P_0)\widetilde G_n(I-P_0)\bigr].
\]

Scalar trace and Berry terms are retained. There is no pointwise
trace-one or CAR projection on `P_j`. Formal PST compatibility is
inherited; a numerical UV bound is not inferred from power counting.

## Jets

Default tables are large enough for a fourth axial derivative of the
defect: input physical order 9 and mu order 5, so `B_j` retains physical
`9-j` and mu `5-j`. Coefficients are mixed Taylor terms
`d_T^t d_z^z d_mu^m f/(t! z! m!)`. Denominators use strictly positive
enclosing intervals; an upper bound is never rounded and then inverted.
Explicit callables or free metric shapes are rejected.

```sh
PYTHONPATH=src /tmp/nsc-validated-runtime.8pg_7t1z/bin/python -m pytest -q tests/test_nsc_scaled_reference_projector.py
```

## Remaining UV inputs

These jets supply the missing mixed derivatives of the auxiliary symbol.
They do not close the finite-history remainder. Still required, and not
computed here, are the exact Weyl-defect seminorms of
[the shift remainder](nsc-weyl-commutator-remainder.md), spatial radius
norms from the existing enclosure, nuclear-norm propagation, the initial
`C_{\mathrm{up}}-P_{\mathrm{up}}` bound, and the changed-history source
UV remainder. The local incoming gate stays OPEN.
