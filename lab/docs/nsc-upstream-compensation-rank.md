# A fixed upstream interval preserves the metric-response rank

\[
\boxed{\sigma_{\min}(\overline V_b)>\frac{51}{200}>0,
\qquad \operatorname{supp}b\subset
[203/200,41/40],\quad b\ge0,\ b\not\equiv0.}
\]

Every smooth common bump in this interval gives three independent real
responses of the actual raw KS lapse, axial scale and radius at the selected
diagonal energy. This closes the width/rank gap in the
[upstream compensation lemma](nsc-upstream-metric-compensation.md).
The result concerns the **normalized** response: the physical response has
smallest singular value greater than `mu_b * 51/200`, with an unspecified
positive bump-dependent scale. No bump shape, compensating amplitudes or
history solution is selected.

The [exact record](../results/development/nsc-upstream-compensation-rank.json)
authenticates the previous operator, geometry and selected-label owners.
Its [producer](../scripts/derive_nsc_upstream_compensation_rank.py) uses
symbolic identities and rational arithmetic only.

## Actual metric vertices and input guards

At `N_K=1,beta_K=0,k=-E`, the
[raw Dirac vertices](../src/recursive_horizons/nsc_incoming_cauchy_jets.py)
give

\[
V=\begin{pmatrix}-m&0&0\\ \ell/r&0&-\ell/r^2\\
-E/a&E/a^2&0\end{pmatrix},\qquad
V^{-1}=\begin{pmatrix}-1/m&0&0\\-a/m&0&a^2/E\\
-r/m&-r^2/\ell&0\end{pmatrix},\qquad
\det V=-\frac{m\ell E}{a^2r^2}.
\]

The selected block is `14_1`, with the positive energy
`E=0.5562120090641313` from the
[radial response record](../results/development/nsc-retarded-radial-response.json).
The [retained channel record](../results/development/nsc-mode-resolved-cauchy-state.json)
supplies stored binary mass `1.5707963267948966` and angular label
`2.23606797749979`. Their exact binary rationals satisfy

\[
3/2<m<8/5,\qquad 2<\ell<9/4,\qquad 1/2<E<3/5.
\]

The defining labels `m=pi/2,ell=sqrt(5)` obey the same conservative guards;
the proof does not identify their exact values with rounded stored values.
It applies to either choice within these bounds without refitting a label.

Set `rho_c=51/50,h=1/200`. The support lies strictly in `(1,103/100)`.
For the [fixed geometry](../src/recursive_horizons/nsc_lorentzian.py),

\[
r=\sqrt{1+\rho^2},\qquad
a^2=3[(1+\rho^2)\arctan(1/\rho)-\rho]-1,
\qquad (a^2)'=6[\rho\arctan(1/\rho)-1]<0.
\]

Reuse the [rational geometry bound](nsc-incoming-fourier-matching.md):
at `rho=103/100`, `atan(1/rho)>pi/4-3/203` and
`pi>314159/100000` imply

\[
a^2>\frac{547699824079}{812000000000}>\frac{16}{25}.
\]

Also `a(1)^2=3pi/2-4<1` using `pi<22/7`. Throughout the interval,
`4/5<a<1`, `2<r^2<21/10`, `r'=rho/r<1`, and `ell/r<8/5`
since `(9/4)^2/2<(8/5)^2`. The canonical generator
`G=iH/a0` has `H=-m sigma1+(ell/r)sigma2-(E/a)sigma3`, so

\[
\|G\|^2<\frac{2(8/5)^2+(3/4)^2}{(4/5)^2}
=\frac{2273}{256}<9.
\]

The clock `a0` is the fixed reference coordinate factor, equal to the
background a in these bounds. It is not differentiated under a raw axial
metric variation.

## Exact rational rank margin

The inverse gives

\[
\|V_c^{-1}\|_F^2
=\frac{1+a^2+r^2}{m^2}+\frac{a^4}{E^2}+\frac{r^4}{\ell^2}
<\frac{82}{45}+4+\frac{441}{400}
=\frac{24929}{3600}<7.
\]

Thus `sigma_min(V_c)>1/sqrt(7)>3/8`. Independently,

\[
\|V\|_F^2<2(8/5)^2+(3/4)^2+(15/16)^2+(9/8)^2
=\frac{50093}{6400}<9.
\]

For the derivative, `atan(x)>=x/(1+x^2)` implies
`rho atan(1/rho)>1/2` on this interval. Consequently
`a'=3[rho atan(1/rho)-1]/a<0` and `|a'|<15/8`. The only nonzero entries are

\[
V'_{21}=-\ell r'/r^2,\quad V'_{31}=Ea'/a^2,\quad
V'_{32}=-2Ea'/a^3,\quad V'_{23}=2\ell r'/r^3.
\]

Their absolute bounds, in that order, are
`9/8,225/128,1125/256,9/5`. The last uses `r^3>5/2`, which follows already
from `r^2>2`. Therefore

\[
\|V'\|_F^2<\frac{44085141}{1638400}<36.
\]

Let `O(rho)` be the SO(3) adjoint action of the exact canonical transport
`U_E(1,rho)`. Unitary transport gives `||O'||_2<=2||G||<6`; hence

\[
\|(OV)'\|_2\le2\|G\|\|V\|_2+\|V'\|_2<24.
\]

With the common positive reference weight,

\[
\overline V_b=
\frac{\int (b/a_0)OV\,d\rho}{\int (b/a_0)\,d\rho},\qquad
\|\overline V_b-O_cV_c\|_2<24h=\frac3{25}.
\]

Rotations preserve singular values. The singular-value perturbation bound
then yields

\[
\sigma_{\min}(\overline V_b)>
\frac38-\frac3{25}=\frac{51}{200}.
\]

No propagation values or transport phases enter this certificate.
The zero-transfer axial factor is the existing positive `w_hat(0)`;
`mu_b=w_hat(0) integral(b/a0)>0` multiplies all three columns equally.
There is no bump-independent numerical lower bound on mu_b here.

## Meaning for the next parent-history construction

For real amplitudes alpha, the actual linear response is
`K_comp=-i mu_b (Vbar_b alpha).sigma`. Its rank permits cancellation of one
selected diagonal traceless K by actual metric directions. The coefficients
and metric functions must be the **same for every energy and channel**;
independent spectral compensation is not an allowed metric history.
This one-pair condition does not establish full C0 matching, source
preservation, the parent constraints or action stationarity.

Support strictly away from Sigma leaves every incoming metric jet and its
restriction frame unchanged. The exact source and earlier preparation stay
fixed; retarded initial/inflow support guards must be checked for the union
of supports when a history is implemented. A sufficiently small overall
family amplitude is required to preserve positive metric fields and the
connected chart. Rank alone gives no finite-amplitude admissibility bound.

Direct KS Fourier/Weyl vertices at zero transfer implement this linear map.
A PG implementation must use the full
[fixed-chart pullback](nsc-pg-ks-metric-pullback.md), including principal
variations. The construction introduces no matrix counterforce, coupling
fit, state change or varying reference clock.

**Decision entry.** Reuse: the raw-vertex rank lemma and exact geometry.
Gap closed: an explicit nonempty common support interval with certified
rank. Decision: actual response-column evaluation is now justified on this
interval if selected as the next history construction. Stopping condition:
the positive exact margin `51/200`; no width search, amplitudes or field
solve are part of this gate.

Verification: `python3 scripts/derive_nsc_upstream_compensation_rank.py --check`
and `python3 -m unittest discover -s tests -p test_nsc_upstream_compensation_rank.py`.
