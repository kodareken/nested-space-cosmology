# Actual fourth-order defect coefficients at small axial transfer

On one owned geometry box, the [regular projector](nsc-scaled-reference-projector.md)
supplies all momenta `|k|>=K`, including the shifted segments required by
`|omega|<=|k|`. Its mu interval must contain `[0,2/K]`; a pointwise energy
value is rejected. The coefficients bound the actual Weyl defect, using the
formal transport identities only as exact algebraic cancellations.

For `R_T=P_T+i[H_T,P]`, the [exact remainder identity](nsc-weyl-commutator-remainder.md)
gives the kinetic/temporal part

\[
|k|^{-5}\left[\partial_TB_4+\tfrac12\{\sigma_3/a,\partial_zB_4\}\right].
\]

The Pauli anticommutator cancels its off-diagonal entries algebraically.
No measured drift or field contribution is removed. An axial derivative of
order q uses `B4_(T z^q)` and `B4_(z^{q+1})`, so q=4 requires total retained
physical order9.

For each projector order j, put `n=5-j`. Its potential Taylor remainder
uses the nth momentum derivative. The regular derivative core obeys

\[
\partial_k^n P_j=\mu^6D[j+1,n]B_j.
\]

Both shifted momenta keep their sign and have magnitude at least `|k|/2`.
For q axial derivatives, Leibniz and the reciprocal Fourier moments give

\[
\|\partial_z^q R_{T,\rm small}\|_F
\le C_{5,q}|k|^{-5}+C_{6,q}|k|^{-6},
\]
\[
C_{6,q}=\sum_{j=0}^4\sum_{h=0}^q
{128\binom qh\,|\ell|\,A_{5-j+h}(1/r_g-1/r_{\rm ref})
\over 2^{5-j}(5-j)!}
\sup\|\partial_z^{q-h}D[j+1,5-j]B_j\|_F.
\]

The factor128 is the product of the two shifted segments and `2^6` from
their minimum momentum. All needed potential moments have positive order,
so the spatially constant reference potential contributes zero to them.
Full moments conservatively bound the small-transfer subset. The rho-clock
bound is obtained from `R_rho=-R_T/a`, since a has no z dependence.

The recorded pilot uses the original nonzero-history profile, the group14
channel parameters, both signs of canonical k, and a small geometry box in
the preparation/axial transition. It is not a new field or source run.
Large transfer, coverage of the entire preparation domain and time
integration are not supplied. In particular, canonical-k bounds are not
silently identified with the tail of the original input-energy quadrature.
The physical local incoming gate remains OPEN.
