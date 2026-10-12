# Fixed coefficients of the compatible incoming surface operators

Reuse `delta r=T*w(z)+T^3*U(z)/6`, `delta a=0`, the same local/reference
action and fixed intrinsic source. The operator shape is already established:

$$
E_N=A(w)U+B(w)w''+C(w')^2+D(w),\qquad
E_\beta=dU'+ew'''+F(w)w'+S_\beta.
$$

## Decision and stopping condition

The missing quantities are the fixed polynomials F(w), D(w)-S_N and C.
Reuse the closed shift-response formulas with `Hr -> Hr+w/r`, while
`R2=r_TT/r` remains fixed, for F. Use the existing homogeneous Bloch
normalization and standard full-line moments for the lapse-reference change;
use the covariant homogeneous local Euler operator for its local part.
Extract the local v² coefficient only AFTER the full two-dimensional lapse
Euler variation. Its derivative grade4 permits zero other derivative
backgrounds for that coefficient. The reference C coefficient has a separate
independent owner.

Stop after one closed polynomial/interval record and focused replay tests,
or a precise missing term. Do not integrate w(z), choose initial/boundary
data, change a source/state/scale/action parameter, rerun old Cu/Cv or
reference/probe generators, or expand beyond this operator closure.

Likely issues are missing radius-chain-rule terms in higher time derivatives,
incorrect treatment of the full spatial lapse variation, and confusing the
source-independent coefficient proof with numerical source accuracy. Keep
raw r_TT,r_TTT fixed for D(w), apply all Euler derivatives before the v²
restriction, and preserve exact formal S_N,S_beta throughout. Stored mixed
and v responses remain validation controls.

## Closed coefficients

Write `H=Hr+w/r`, with `a,r,Ha,A2,R2,A3,R3` fixed to the original slice.
Here `A2=a_TT/a`, `R2=r_TT/r`, `A3=a_TTT/a`, `R3=r_TTT/r`;
the shifted values of `Hr` are not inserted into their baseline relations.
The exact reference and local F formulas are imported verbatim from the
[closed shift coefficient](nsc-incoming-shift-coefficient.md), then evaluated
at H. Its saved constant interval is reused unchanged. This gives the new
linear and quadratic coefficients without repeating the Cv proof.

The full local lapse Euler expression is

\[
L_N-\partial_T L_{N_T}-\partial_z L_{N_z}
 +\partial_z^2 L_{N_{zz}}.
\]

It is formed before restricting the radius jets. Setting all derivative
backgrounds except `r_Tz=v` to zero then gives

\[
C_{\rm local}=\frac{64\pi C_W}{3a}
 -\frac{h_q+\log r}{15\pi a}=-\frac{2K}{a},\qquad
K=-\frac{32\pi C_W}{3}+\frac{h_q+\log r}{30\pi}.
\]

The Einstein, constant Euler, LLL, WZ Euler and WZ BoxR contributions to
this coefficient are zero after the full Euler variation. Import the
[independently proved full-Weyl reference coefficient](nsc-incoming-surface-quadratic-lapse.md)
from `b74f5fc`:

\[
C_{{\rm ref},g}=\frac{D_g L^2(L^2+3m^2)}
 {60\pi a r^2(m^2+L^2)^3},\quad L=\ell/r.
\]

Both local and per-group reference expressions obey
`C=−(A1+d/r)/a²` exactly; the regular-w branch has its own owner.
Together with the frozen `B=−(2A+(Hr+w/r)d)/a²`, `e=−d/a²`,
this gives `Delta_w=−d*C` for `Delta=A*e−B*d`; the artifact retains
the exact algebraic residual. No profile or crossing is inferred.
`D_g=copy_count*degeneracy` includes the angular signs once.

For the homogeneous energy change, put `h=(-m,L,p)` and
`omega²=m²+L²+p²`. The existing normalized Bloch recursion gives

\[
b_0=-h/\omega,\quad b_1=(h\times h_T)/(2\omega^3),\quad
b_2=-h\times\partial_Tb_1/(2\omega^2)
     -|b_1|^2b_0/2,
\]
\[
E_2=\omega|b_1|^2/2,\quad E_3=\omega b_1\cdot b_2,\quad
E_4=(h\times b_1)\cdot\partial_T b_2/(2\omega)
       +\omega|b_2|^2/2.
\]

E3 integrates to zero by odd momentum parity. E2 and E4 have positive-gap
rational kernels. Standard beta-function full-line moments and the inherited
measure `D_g*a*dp/(2*pi)` integrate them exactly; no momentum quadrature is
used. The module retains the raw derivative chain rules for `ell/r` and
`k/a`. Differentiating the integrated E4 with respect to R3, divided by r,
reproduces the frozen Cu equation exactly. Subtracting the baseline after
`Hr -> H` leaves degree at most four. The local homogeneous lapse polynomial
uses the same full Euler expression. The zero-angular massive channels have
zero reference change; the LLL geometry remains included once locally.

The content-addressed proof artifact stores every exact per-channel
polynomial. At the unchanged ledger the coefficient intervals are:

| Coefficient | Directed interval |
|---|---|
| C | [0.11119593102288508, 0.11119593102288536] |
| F0, imported Cv | [16.575975622058497, 16.575975622058518] |
| F1 | [−0.07632764397901040, −0.07632764397900982] |
| F2 | [−0.003168250383882580, −0.003168250383882434] |
| D1 | [−23.722260206004698, −23.722260206004627] |
| D2 | [9.350373714084240, 9.350373714084265] |
| D3 | [−0.04426658181020114, −0.04426658181020083] |
| D4 | [0.011149125347319842, 0.011149125347319903] |

Here `D(w)=S_N+D1*w+D2*w²+D3*w³+D4*w⁴` and
`F(w)=F0+F1*w+F2*w²`. S_N is the exact finite baseline constraint constant,
including its original local allocation; S_beta is treated likewise.
Neither is assigned a numerical source approximation by this coefficient
record. The executable momentum primitive is
`G(w)=F0*w+F1*w²/2+F2*w³/3`, so
`E_beta=partial_z[d*U+e*w''+G(w)]+S_beta` without choosing any endpoint or
integration constant.

## Verification and scope

The 50/70-digit interval evaluations use the same exact-geometric/stored-value
hulls and unchanged binary ledger as the prior coefficients. They have
intersecting enclosures. The one new homogeneous w-direction point control
uses four signed momentum/angular points, without integration, and agrees
with the original projector energy orders within 1.74e−18. The full local
density adapter agrees with the original density partial derivatives within
8.89e−16 (both tolerance 3e−13). Existing v and mixed action responses,
including local channel allocations, are replayed against the new C and the
saved Cu/Cv; tolerance is the unchanged 3e−11 in action units. The old
coefficient estimate differs by about 1.43e−8, but its v=.01 action response
difference is about 1.43e−12; these quantities are recorded separately.

The result closes these polynomial coefficients of the included lapse/shift
operators. It supplies neither the other metric equations nor a physical
boundary/history choice. Source numerical accuracy, global endpoint balance,
metric evolution and Gamma_rest remain their separately owned questions.
The script's `--check` evaluates the saved exact expressions and verifies
hashes; it runs no symbolic coefficient, projector or old probe producer.

```sh
python3 scripts/derive_nsc_incoming_surface_coefficients.py --prepare
python3 scripts/derive_nsc_incoming_surface_coefficients.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_coefficients.py
```
