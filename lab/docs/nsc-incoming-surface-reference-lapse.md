# Closed reference coefficient of the compatible spatial lapse response

For `B_ref=partial E_N/partial(r_Tzz)` at the inherited incoming baseline,
write `p=k/a`, `L=ell/r`, `M2=m^2+L^2`, `Ha=a_T/a`, `Hr=r_T/r`, and
`D=copy_count*degeneracy`. The per-group result is

\[
\boxed{B_{{\rm ref},g}=
\frac{D L^2\{(3H_r-2H_a)M^2+8H_r m^2\}}
 {120\pi a r M^6}.}
\]

This includes the full spatial Weyl recursion through order4. It retains
the first-normal background terms and uses the existing positive subtraction
action sign. `D` already includes both angular signs; no extra frequency
folding is inserted. Massive zero-angular channels give zero. The LLL retains
its separate constant reference and geometric allocation.

Reuse the frozen projector/idempotence equations, incoming intrinsic data,
the band-bulk identity and the standard full-line beta moments. The missing
entry is B_ref in the compatible-function principal matrix. Stop after the
closed kernel/integral and exact checks; no cu/cv calculation, reference
quadrature, baseline probe or mode/radial generator is repeated.

## Lower-order identity suffices for the complete fourth-order trace

Idempotence gives `Tr(H Pn)=omega Tr(Gn)`, where Gn excludes the unknown
`P0 Pn+Pn P0-Pn` terms. For `delta r=b*T*z^2/2`, the changes in P0,P1,P2
at the evaluation point vanish, but `P1_zz`, `P2_zz` and the induced P3
change do not. The surviving fourth-order trace is

\[
\delta\operatorname{tr}G_4/b
=2\operatorname{tr}(P_1\delta P_3/b)
-\tfrac14\operatorname{tr}(P_{0,kk}P_{2,zz}/b)
-\tfrac14\operatorname{tr}(P_{1,kk}P_{1,zz}/b).
\]

Odd Weyl orders cancel in paired traces. The order4 Weyl term involving only
P0 vanishes because all its pure intrinsic spatial derivatives vanish on
the incoming surface. The unknown's grade3 means that only Ha,Hr can multiply
it at order4; higher-normal background derivatives would exceed that grade.

The mixed P3 coefficient is extracted from the exact shifted linear projector
response with `H+/-=H(k+/-q/2)`, and formal mode `exp(-i nu T+i q z)`.
The time/spatial labels are bookkeeping variables only. Its transport and
idempotence identities are checked through the needed `nu*q^2` coefficient.
The P2 spatial derivative reuses the homogeneous lower recurrence, keeping
`delta(h_T)_zz=V*b` and `delta(h_TT)_zz=-4Hr*V*b`, with `V=-L*sigma2/r`
and fixed field `r_TT`. The undifferentiated changes still vanish at the
evaluation point. Dropping the second relation would lose background terms.

## Point kernel and full-line integral

With `omega^2=p^2+M2`, the fourth-order kernel multiplying b is

\[
K_4(p)=\frac{L^2}{32a^2r\omega^9}
\left[-7H_aM^2p^2+H_r\{5m^2M^2+(9L^2+7m^2)p^2+2p^4\}\right].
\]

The third-order kernel is nonzero before integration,
`K3=-L*m*p/(8*a^2*r*omega^5)`. Its full-line integral vanishes by odd parity
and absolute convergence; it is not set to zero pointwise.

The standard full-line moments of `1/omega^9`, `p^2/omega^9`, and
`p^4/omega^9` are respectively `32/(35M^8)`, `16/(105M^6)`, and `4/(35M^4)`.
Multiplying by `D/(2pi)` and `dk=a dp` gives the boxed coefficient.

The [helper](../src/recursive_horizons/nsc_incoming_surface_reference_lapse.py)
checks the matrix, trace, Weyl and moment identities exactly. Its closed
formula is for the principal worker's directed evaluation and integration;
this note makes no independent principal-determinant, physical-surface or
extended-stationarity claim. All sources, scales, the rho0 seam and other
declared action terms remain unchanged.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_reference_lapse.py
```
