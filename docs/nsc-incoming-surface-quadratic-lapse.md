# Closed quadratic reference lapse response

For `v=r_Tz`, the coefficient multiplying v² in the included lapse equation
has the following per-group reference contribution:

\[
\boxed{C_{{\rm ref},g}=
\frac{D L^2(L^2+3m^2)}{60\pi a r^2(m^2+L^2)^3}},\qquad L=\ell/r.
\]

`D=copy_count*degeneracy` includes both angular signs. This is the coefficient
of v², not the second derivative with respect to v. The zero-angular massive
channels give zero; the LLL retains its separate reference/geometry owners.

Reuse the full spatial Weyl projector, its idempotency energy trace identity,
and standard full-line moments. The missing coefficient C_ref is independent
of the other background derivatives because v² already has derivative grade4.
Thus a constant a,r background with `delta r=v*T*z` retains every needed term.
Stop after exact identities and comparison with the existing n24 control;
no new reference quadrature, cu/cv calculation, mode, state or metric run.

## Spatial terms are retained before taking the evaluation point

Let `h=(-m,L,p)`, `p=k/a`, `omega²=m²+L²+p²`, and `V=partial_r h` at
fixed ell. The first projector on a general nearby scalar-radius jet is

\[
P_1=\frac{r_T(h\times V)\cdot\sigma}{4\omega^3}
-\frac{mLr_z}{4ar\omega^3}I.
\]

Although P1 vanishes at T=z=0, its z derivative does not. Its scalar spatial
term has a nonzero time derivative. Keeping both contributions makes the
full second-order transport source vanish at that point, while idempotence
gives the nonzero scalar

\[
P_2(0,0)=-\frac{vL^2p}{8ar\omega^5}I.
\]

On T=0, `r_T(z)=v*z`. The homogeneous P2 dependence is quadratic in r_T
even though the field r_TT remains fixed, because differentiating ell/r
produces `h_TT=2L*r_T²/r²` in its angular component. Therefore P2_zz cannot
be omitted. The surviving fourth-order trace is

\[
\operatorname{tr}G_4=
\operatorname{tr}(P_2^2)
-\tfrac14\operatorname{tr}(P_{0,kk}P_{2,zz})
+\tfrac14\operatorname{tr}(P_{1,zk}^2).
\]

All pure intrinsic spatial derivatives of P0 vanish on T=0, removing the
order4 P0/P0 Weyl term. P1's point value removes the ordinary P1/P3 term.
Odd Weyl trace pairs cancel by the previously verified identity. The helper
checks the complete first- and second-order transport/idempotence equations
before using these simplifications.

## Kernel and inherited full-line measure

Using `Tr(H P4)=omega Tr(G4)`, the v² point kernel is

\[
K_C(p)=\frac{L^2\{L^4+L^2m^2+(12m^2-8L^2)p^2+12p^4\}}
 {32a^2r^2\omega^9}.
\]

The standard moments of `1/omega^9`, `p²/omega^9`, `p⁴/omega^9` and
the existing `D/(2pi)` measure with `dk=a dp` give the boxed result.
The archived n24 control is compared only after both momentum and actual
angular signs are combined; linear-v odd node contributions are not
misidentified as quadratic response. Original arrays remain unchanged and
no old projector or quadrature generator is executed.

The [helper](../src/recursive_horizons/nsc_incoming_surface_quadratic_lapse.py)
exposes the closed formula and exact residuals for the coefficient worker.
Local C and the other compatible-ODE coefficients belong to that independent
owner. No surface boundary, physical initial datum or extended stationarity
is selected here.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_quadratic_lapse.py
```
