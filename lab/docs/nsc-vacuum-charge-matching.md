# A vacuum–gravity–gauge matching condition for the charged throat

The charged source route has a concrete consistency condition. If the
retained positive fermionic spectral contribution is used by itself to set
the local vacuum, Einstein and Maxwell coefficients, it cannot supply the
near-extremal charged starting geometry of the imported throat construction.
Changing the positive cutoff profile, the number of identical charged copies,
or the algebraic warp moments does not remove this obstruction.

This is a result about a specified leading-coefficient realization. It does
not exclude the completed NSC action or all charged geometries. It identifies
the vacuum/gauge/gravity combination that the common functional must supply
before reusing the [charged-fermion self-sourcing solution](nsc-charged-self-sourcing-route.md).

## One normalization for all three coefficients

Use N identical unit-charge complex Dirac fields on the common five-dimensional
carrier. The four-dimensional metric and gauge field are independent of Y,
A_Y=0, and the compact size is held fixed for this coefficient projection.
Let h(z) be a nonnegative scalar spectral profile with finite nonzero moments

\[
Q_\alpha[h]=\frac1{\Gamma(\alpha)}\int_0^\infty
z^{\alpha-1}h(z)dz,\qquad
f_5=Q_{5/2},\quad f_3=Q_{3/2},\quad f_1=Q_{1/2}.
\]

The published Dirac heat coefficients and the [checked gauge normalization](nsc-charged-self-sourcing-route.md)
give the following **volume, R4 and gauge-strength terms** after integrating
the common conformal warp. Write I_p=integral exp(p sigma)dY and
c5=(4 pi)^(5/2):

\[
V_{\rm vol}=\frac{4N\Lambda^5 f_5 I_5}{c_5},\qquad
A=\frac{N\Lambda^3 f_3 I_3}{3c_5},\qquad
C=\frac{2N\Lambda f_1 I_1}{3c_5}.
\]

The additional warp-curvature contribution is retained separately below.
These coefficients come from one contribution to the same functional;
they are not three independently fitted weights. Other field-dependent
potentials of the complete operator, including a nonuniform scalar link,
have not been evaluated in this coefficient assignment.

For the corresponding local Lorentzian action in -+++ conventions,

\[
S_4=\int\sqrt{-g}\,[A R-C\mathcal F_{\mu\nu}\mathcal F^{\mu\nu}-V],
\]

the identifications are G_N=1/(16 pi A), g4²=1/(4C), and
lambda4=V/(2A). Here lambda4 is a cosmological curvature parameter,
not the spectral cutoff Lambda. With nonzero integer magnetic flux q,
the charge-radius normalization of
[Maldacena–Milekhin–Popov, equation2.3](https://arxiv.org/html/1807.04726v3)
is r_Q²=pi q² G_N/g4². Therefore

\[
\boxed{\Xi\equiv\lambda_4r_Q^2=\frac{q^2VC}{8A^2}.}
\]

For V=V_vol this becomes

\[
\Xi_{\rm vol}=3q^2\frac{f_5f_1}{f_3^2}\frac{I_5I_1}{I_3^2}.
\]

N and the overall cutoff scale cancel. This is a dimensionless connection
between three physical sectors, evaluated before assigning their total
renormalized values.

## Bounds that do not require a spectrum scan

Cauchy–Schwarz applied to the Mellin integrals gives

\[
\frac{f_5f_1}{f_3^2}\geq
\frac{\Gamma(3/2)^2}{\Gamma(5/2)\Gamma(1/2)}=\frac13.
\]

A second application to exp(5 sigma/2) and exp(sigma/2) gives
I5 I1≥I3². Consequently Xi_vol≥q² for any such nonnegative profile.

If h has a nonnegative proper-time representation, its moments are
positive averages of t^(-5/2), t^(-3/2) and t^(-1/2). Cauchy–Schwarz
then gives the stronger f5 f1≥f3², hence

\[
\boxed{\Xi_{\rm vol}\geq3q^2.}
\]

This subclass includes the exponential heat profile and the adopted
proper-time window. For u=nu_match/Lambda in (0,1), that window has
f5=(1-u^5)/5, f3=(1-u^3)/3 and f1=1-u. Its bound is also checked through

\[
9(1-u^5)(1-u)-5(1-u^3)^2
=(1-u)^4(4u^2+7u+4)\geq0.
\]

A nonnegative resolution kernel does not guarantee that every possible
completed action profile is nonnegative. For example, the previously
compared weighted-log profile changes sign. The result is not extended
to such profiles or to the complete scale-integrated quantum functional.

## The charged starting geometry imposes the opposite requirement

Use the standard Reissner–Nordstrom–de Sitter lapse, reviewed in
[Montero, Van Riet and Venken, section2](https://arxiv.org/html/1910.01648v4),
with the charge radius normalized as above:

\[
f(r)=1-\frac{2G_N M}{r}+\frac{r_Q^2}{r^2}-\frac{\lambda_4 r^2}{3}.
\]

At a double horizon, f=f'=0. The identity
f+r f'=1-r_Q²/r²-lambda4 r² gives

\[
\Xi=z(1-z)\leq\frac14,\qquad z=\lambda_4r_h^2,
\qquad \frac14-z(1-z)=(z-\tfrac12)^2.
\]

Thus the retained volume coefficients lie outside this charged extremal
range for every nonzero integer q. The approximately flat magnetic seed
used in the imported construction needs an even stronger separation from
cosmological curvature. This does not re-solve its Einstein/Casimir problem;
it checks whether our proposed coefficient assignment can be inserted into it.

## The known warp term is not silently omitted

The standard conformal-curvature transformation in five dimensions gives
R5=exp(-2 sigma)(R4-8 sigma''-12 sigma'^2), where primes are Y derivatives.
This is an application of the transformation formulas reviewed in
[Dabrowski, Garecki and Blaschke](https://arxiv.org/html/0806.2683v3).
The remaining bulk Einstein contribution to the four-dimensional constant
term is, with A5=A/I3 including the same N fermion copies,

\[
V_{\rm warp,bulk}=A_5(8B-12J),\qquad
B=[e^{3\sigma}\sigma']_{-L_\star}^{L_\star},\qquad
J=\int e^{3\sigma}(\sigma')^2dY.
\]

For the specified sigma=-18(Y/L_star)²/1015, the exact Gaussian integrals
give

\[
\frac{V_{\rm warp,bulk}}{V_{\rm vol}}
=-\frac{f_3}{f_5}\frac{0.0235029858774536\ldots}{\zeta}.
\]

In the proper-time window, f3/f5≤5/3. The magnitude is therefore at most
0.0391716431290893…/zeta. For zeta≥1, retaining this complete bulk
Einstein warp term still leaves Xi≥2.882485…q². The controlled bulk
expansion requires stronger scale separation, including nu_match L_star≫1;
the small-size region is not used to manufacture a root.

For comparison, the standard Dirichlet-metric Gibbons–Hawking–York term
would cancel 8A5 B. Its remaining warp factor is only
0.000418230858372… times (f3/f5)/zeta. This boundary completion is **not**
assumed to be the one derived by the actual Dirac transmission domain.
Independent boundary terms remain part of the open matching problem.
Higher-curvature terms and other quantum contributions are not evaluated
by this leading-order comparison. Small controlled corrections cannot
bridge an order-one, growing-with-q mismatch; corrections of that size
would require a different or untruncated realization.

## What the complete action must determine next

For a positive-vacuum Einstein–Maxwell seed, the necessary condition on the
**complete** coefficients is

\[
\boxed{\frac{V_{\rm full}C_{\rm full}}{A_{\rm full}^2}\leq\frac{2}{q^2}.}
\]

If A and C were unchanged and only V corrected, the positive proper-time
baseline would require V_full/V_vol≤1/(12q²), even before demanding an
approximately flat exterior. This is a required effect to derive, not a
counterterm value selected here.

The missing contribution may belong to the full invariant/scale measure,
the remaining quantum fields and gauge/constraint determinants, the
actual compact boundary action, further operator-field potentials, or a
nonlocal/multimetric realization.
All must come from the specified common functional. A change of regulator
or a vanishing vacuum subtraction is not inserted to obtain the desired
geometry. The positive-trace bound does not apply to a full supertrace;
neither does counting bosonic and fermionic degrees of freedom establish
the needed cancellation without their actual operators, spectra and measure.

The [development record](../results/development/vacuum-charge-matching.json)
authenticates the already released input records. The
[runner](../scripts/check_nsc_vacuum_charge_matching.py) verifies the exact
normalization and horizon identities, the window polynomial, and the
Gaussian warp integrals against independent adaptive integration. All
fields are compared; float tolerances are 3e-13 absolute and relative.
No magnetic-throat simulation or prior record reproduction was repeated.

This is a new matching result for the project using established mathematics.
No world-priority claim, full-theory exclusion, or completed self-sourcing
claim follows. The v0.4.0 release remains the published checkpoint while
this additional matching condition is developed.
