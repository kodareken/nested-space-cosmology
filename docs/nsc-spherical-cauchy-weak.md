# Weak initial residual from the spherical action

Diagnostic only. The production evolution remains
`nsc_spherical_galerkin_coupling`. This helper does not step a state, rerun
the radius Newton, or retune a tolerance. Record
`lab/results/development/nsc-spherical-cauchy-weak-v1.json`, status
`PARTIAL_TERMS_NO_TOTAL_BOUND`. Historical execution checkpoint `9a9090a`.
The saved v5 JSON and NPZ bytes are unchanged. Sealed assessment CPU
`6.549308` s.

`--check` reads that sealed record. It remeasures the scientific fields
from the declared v5 bytes and does not write the record, the v5 files, or
a PDF. `checkpoint_head`, `cpu_seconds`, and `wall_seconds` are historical
execution metadata. The live Git HEAD is not required to equal `9a9090a`.
A fresh assembly records the live checkpoint only through `--output` on a
new path. That path cannot be the sealed JSON.

The sealed generator hashes are the helper, driver, and owned test that
wrote the record. Those historical bytes are not stored again. The current
files differ, so replay reports that limit and does not invent a matching
hash. The coupling, Galerkin, feedback-action, and Cauchy-data hashes remain
binding. A missing declared source is an unavailable-byte limit. A changed
declared hash or a changed scientific field is rejected. No total bound is
added to that sealed diagnostic.

The CLI `--check` now names formula-source domain
`5f10ecd365843d1616e50eb16a20d7acd8377e2c` explicitly. The binding resolver
authenticates those exact historical bytes; the current helper remeasures
the saved scientific fields without invoking the original producer.
`binding_report` and `verify_saved` remain strict on current formula bytes
unless a caller supplies `source_ref`. Historical reuse of the separate
initial certificate likewise calls
`initial_lapse_binding_report(record, source_ref=...)`; its default mode
rejects a changed current formula source. Numerical payload hashes remain
binding in both modes.

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_cauchy_weak.py tests/test_nsc_spherical_cauchy_weak_independent.py -q
python scripts/lab.py scripts/derive_nsc_spherical_cauchy_weak.py --check
```

Owned and independent tests passed together in 6.34 s. `--check` on live
HEAD `e62f804` remeasured in 5.822429 s, kept historical head `9a9090a`,
reported three generator limits, and did not write. Scientific comparison
uses a relative tolerance of `1e-12` and an absolute tolerance of `0`.
That tolerance is not a geometry bound.

Measurement files at the earlier recorded replay differed from the sealed generator hashes:

| Path | sha256 |
|---|---|
| `lab/src/recursive_horizons/nsc_spherical_cauchy_weak.py` | `c993d628b5df4632a92080c2f249e84fce4fea97be481443d2108e43081acfe1` |
| `lab/tests/test_nsc_spherical_cauchy_weak.py` | `78867955ba8dbdc20682d6d779d0cb36ed818e8769a007abe56b2a8f11486d1a` |
| `lab/tests/test_nsc_spherical_cauchy_weak_independent.py` | `3599ae84d911f88138a01bd3876c28253c2afd29d294d08961152fb30ed50b76` |
| `lab/scripts/derive_nsc_spherical_cauchy_weak.py` | `3b2554a6e23fd15eb25f82f536471299e13a692a3e372466c231e9645e084fd7` |

## Reduction

For constant \(Q\), \(\chi=p_r=p_\chi=0\), \(r=y^2\) and \(\rho=\mathrm{force}_L/\mathrm{d}x\),

\[
G=Q^2 r_{\mathrm{mag}}^2+\frac{Q\rho}{8\pi A},\qquad
R_y=-4y''+Q^2 y-\frac{G}{y^3},\qquad
C= -\frac{8\pi A}{Q} y^3 R_y.
\]

SymPy simplification of `lapse_constraint + rho` against this \(C\) is 0.
The owner residual is `_radius_residual`. Its \(D(F)\) rearrangement,

\[
C=\frac{8\pi A}{Q}\big(D^2(r^2)-3(Dr)^2-Q^2 r^2+G\big),
\]

matches that owner. The \(y\) formula adds the product-rule defect
\(D^2(r^2)-3(Dr)^2-4 y^3 D^2 y\). On a smooth mode the owner gap is
`1.460520593354886e-11` and the action gap is `1.825933182431072e-09`.
On mode 40 the algebra gap stays `1.0913936421275139e-11` and the action
gap is `3143.991531266489`.

For \(G>0\), \(y>0\), and \(G\) independent of \(y\),

\[
J(y)=\int\Big(2(y')^2+\frac{Q^2 y^2}{2}+\frac{G}{2y^2}\Big)\,\mathrm{d}x
\]

has Hessian \(4\|v'\|^2+\int(Q^2+3G/y^4)v^2\). Declared \(G>0\) and
\(\rho\) independent of \(y\) give the uniform constant

\[
\mu=\min(4,Q^2)=Q^2=0.0631642220827373,
\]

so \(H(v,v)\ge 4\|v'\|^2+Q^2\|v\|^2\ge\mu\|v\|_E^2\) without \(y_*\).
An omitted independence flag does not establish that hypothesis. The
sharper value \(\min(4,c_{\mathrm{seg}})\) still uses the segment. Neither
constant is a total bound: the full residual dual, a continuum positivity
enclosure of \(G\), and a rounding enclosure are open.
`geometry_error_upper_bound` stays null.

If \(y_*\) is critical in the tested space and
\(c=Q^2+3G/z^4\ge c_{\mathrm{seg}}>0\) on the segment,

\[
\|y-y_*\|_E\le \|R_y\|_*/\min(4,c_{\mathrm{seg}}),
\]

with \(\|v\|_E^2=\|v'\|^2+\|v\|^2\) and
\(\|R_y\|_*=\sup|\int v R_y|/\|\cdot\|_E\). The supremum is only over the
space actually tested. A geometry-band number omits held-out modes. The
quadrature grid omits the continuum tail. The shift current is outside \(J\).

A supplied discrete segment has energy error `0.0025431085506270297` and
upper value `0.016958726993982685`. With \(G=-1\) the constant Hessian is
`-23.494686223338096`, and no inequality is applied. Mode 20 at amplitude
`0.02` has represented energy dual `0.00010721334128839067`, unresolved
strong maximum `19.744301113926195`, and omitted energy
`0.6295904828089242`. The represented piece is the nonlinear image of that
mode. `shift_momentum` leaves a constant current `0.3` in the residual.

## Saved v5 seed

\(G\) runs from `0.06316419930011759` to `10.951431857279351`. Pointwise
\(c_{\min}\) is `0.07126986300473899` at \(n_f=256\) and
`0.07127002523127005` at \(n_f=512\). Varying \(r\) by `1.7` changes
\(\mathrm{force}_L\) by `0`.
The radius Newton was not rerun. Column bytes were unchanged.

| | \(n_f=256\), \(n_q=1024\) | \(n_f=512\), \(n_q=2048\) |
|---|---:|---:|
| Owner strong max | `0.0017681375243085995` | `9.981455056262689e-06` |
| Represented pull max | `9.608390314276184e-10` | `8.49761520175214e-09` |
| Represented energy dual | `5.6611510229305975e-09` | `1.5869576914267803e-11` |
| Held-out strong max | `0.0017681375412189999` | `9.97865903324404e-06` |
| Energy dual of \(R_y\) | `2.3735379313702133e-07` | `1.0212964160773616e-09` |
| Product-rule defect in \(C\) | `2.746484871086371e-06` | `2.3743109576536128e-05` |
| Algebraic identity gap | `3.345977266902469e-09` | `2.5966427413993834e-08` |
| Doubled-grid strong max | `0.0017723150093473805` | `4.539313647278211e-05` |
| Current mean | `-2.2413621942646887e-14` | `1.3585588092395491e-14` |
| Held-out momentum | `2.489029399344809e-12` | `7.241782087509932e-12` |

Replay against the saved v5 maxima differs by at most
`1.368327673390013e-11`. The algebraic gap is the \(D(F)\) cancellation.
The largest \(D^2(r^2)\) term is `9.426201686802028` and
`9.43418275132038`, and one-multiply \(\varepsilon\) times that condition
is `2.399209290239975e-11` and `4.253452800603874e-09`. Nyquist amplitudes
are grid samples, `3.3964321157373294e-11` and `5.220894150511344e-11`.
No tail bound is stored. The \(y=\sqrt{r}\) product-rule defect at
\(n_f=512\) still exceeds that seed's dense owner maximum. The represented
dual does not certify the total initial error.

The strong numerator that does not pass through \(y\) is
\(P=2rr''-(r')^2-Q^2 r^2+G\), with \(C=(8\pi A/Q)P\). Antiperiodic column
products have degree at most \(nf-1\), and this odd-band radius has degree
at most \(ng/2\), so \(P\) has degree at most \(nf-1\). Both baselines have
\(nq=4nf\), and the measured `1e-8` support sits inside that
cap. On these four grids the dense owner and the Fourier product compare as
follows. The \(n_f=512\) dense maximum moves from
`9.981455056262689e-06` to `4.539313647278211e-05`, while the product
maximum stays from `5.068517579063356e-06` to `5.400936559883778e-06`.
That movement is a conditioning indicator. It is not a rounding enclosure.

## Initial lapse component

The sealed record is unchanged: `geometry_error_upper_bound` stays null
there, and `--check` still does not write. A separate enclosure is
`initial_lapse_component_bound` in the same helper. Its record is
`lab/results/development/nsc-spherical-cauchy-error-v1.json`
(`e0a37e1177b2a483eca0eb3e44243893bbdc04667496b732dd7e63467a97bdd2`),
status `FINITE_CONDITIONAL_INITIAL_LAPSE_BOUND`. It does not step the state,
rerun the radius Newton, or change the source formula. v5 JSON
`1176c5a10eab4fe61a7571ecf0a7708a6fc74b6f0ef20c54cdab53409a0e0efc` and
v5 NPZ `0c92d3d545afffb47f22c9e8d705dc4d70b2df2ca2fdf5c6abc254f37af35b0d`
were read and not written. Sealed weak bytes remain
`ce5e4345f75b187661f965f96c680a3076331e2df709845c22251edb99b3f148`.

The saved `n_f=512` columns are half-density samples of an antiperiodic
trigonometric interpolant. Fermion modes run from `-255.5` to `255.5`.
`P=-i∂_x` multiplies mode `m` by `2π m/L`, which is the symbol of the
owner momentum matrix at `η=1/2`. With the declared occupations,
`M=4κ`, and constant `Q`,

\[
\mathrm{force}_L=M(K/Q+\kappa\,\mathrm{Re}\,S_1),\qquad
\rho=\mathrm{force}_L/\mathrm{d}x.
\]

`force_L` does not read `r`, `χ`, the momenta, the lapse profile, or the
shift, so `G` is independent of `y` by that formula. Conjugate products of
the fermion modes are periodic of degree at most `511`. The odd radius
band has degree `255`, so `r r''` and `(r')^2` have degree at most `510`.
Thus `P` has degree at most `511`. The saved `n_q=2048` Nyquist is `1024`,
so that polynomial is alias-free on the owner grid. Quadrature refinement
is not the rounding certificate.

On the saved state the interval enclosure gives

| Quantity | Value |
|---|---:|
| `‖P‖_L2` upper | `7.093398495891281e-07` |
| `‖P‖_L2` directed gap | `6.947035723309529e-10` |
| Tail `‖P‖_L2` above degree 511 | `1.972067957485316e-09` |
| Continuum `r` lower | `4.210475520492436` |
| Continuum `G` lower | `0.053733032411766975` |
| `μ=min(4,Q²)` | `0.06316422208273727` |
| `‖Ry‖_L2` upper | `8.210273762221587e-08` |
| `‖y-y_*‖_E` upper | `1.2998297915340663e-06` |
| Pointwise `‖r-r_*‖_∞` upper | `4.060849694930819e-06` |

`P=2rr''-(r')^2-Q²r²+G=-y³ Ry` for `y=√r>0`. The energy dual satisfies
`‖Ry‖_*≤‖Ry‖_L2≤‖P‖_L2/y_min³`. On `{y>0}` with `G≥0`,
`H(v,v)≥4‖v'‖²+Q²‖v‖²≥μ‖v‖_E²`, so any positive critical point of `J`
obeys `‖y-y_*‖_E≤‖Ry‖_*/μ`. The circle embedding
`‖v‖_∞≤√(coth(L/2)/2) ‖v‖_E` then lifts that distance to the areal
radius. Positivity uses a `65536`-node interval evaluation plus the
coefficient `ℓ¹` bound on the first derivative. Arithmetic is
`mpmath.iv` at 25 decimal digits, with `math.nextafter` on the stored
endpoints. The recorded run used `20.465117` s of process time.

The smaller saved areal-radius effect scale is `0.0018488643658258752`.
One percent of it is `1.8488643658258752e-05`. The pointwise radius
bound is `0.21964021644804999` of that bar and
`0.0021964021644805` of the effect. This comparison is only the initial
chart. It does not propagate the error through the episode.

The sealed error record still stores `minimizer_existence_proved` as
null. That field is historical. The existence argument below discharges
it from the sealed positivity hypotheses and does not recompute `P`.
`evolution_error_bound` stays null. `total_certified` is false. A small
weak number is not by itself a continuum solution. Reuse checks the v5
bytes, the declared occupations, and the formula-source hashes. The
helper hash is a generator limit, so this later proof does not rehash
the sealed scientific record. A changed episode hash drops only the
effect comparison.

## Existence on the circle

Let \(L>0\) and \(S^1=\mathbb{R}/L\mathbb{Z}\). Let \(Q^2>0\) and let
\(G\in C^0(S^1)\) satisfy \(G\ge\delta>0\). The admissible set is the
open subset of \(H^1(S^1)\) whose continuous representative is strictly
positive. There is no boundary term.

\[
J(y)=\int_{S^1}\Big(2(y')^2+\frac{Q^2}{2}y^2+\frac{G}{2y^2}\Big)\,dx
\ge 2\|y'\|_{L^2}^2+\frac{Q^2}{2}\|y\|_{L^2}^2.
\]

The constant \(y\equiv 1\) is finite. A sublevel is bounded in \(H^1\).
Cauchy–Schwarz on the shorter arc gives
\(|y(x)-y(z)|^2\le\mathrm{arc}(x,z)\,\|y'\|_{L^2}^2\), so the sublevel
is equicontinuous. The Fourier embedding
\(\|y\|_\infty\le\sqrt{\coth(L/2)/2}\,\|y\|_E\) bounds it in \(C^0\).
Arzelà–Ascoli and weak compactness of the Hilbert ball produce a
subsequence that converges uniformly and weakly in \(H^1\).

If the limit vanishes at \(x_0\), the same Cauchy–Schwarz estimate gives
\(y(x)^2\le|x-x_0|\,\|y'\|_{L^2}^2\). The integral of \(1/y^2\) is then
at least a constant times \(\log(a/\varepsilon)\), which diverges as
\(\varepsilon\to 0\). Fatou's lemma on the nonnegative densities
\(G/(2y_n^2)\) lifts that divergence to the minimizing sequence, so the
limit cannot touch zero. Hilbert weak lower semicontinuity of
\(\|y'\|\) and uniform convergence of the potential produce a minimizer
in the positive set.

On \((0,\infty)\), \(t\mapsto t^2\) and \(t\mapsto t^{-2}\) are strictly
convex, and the second derivative of \(t^{-2}\) is \(6t^{-4}>0\). With
\(Q^2>0\) and \(G\ge\delta>0\), \(J\) is strictly convex and \(C^1\) on
the positive set. The minimizer is therefore the unique positive
critical point. The weak equation is
\(y''=(Q^2 y-G/y^3)/4\). Trigonometric \(G\) is \(C^\infty\), and the
right-hand side bootstraps to a \(C^\infty\) positive solution.

The sealed enclosure has \(G\ge 0.053733032411766975\), \(Q^2=0.0631642220827373\),
\(L=8\), and degree at most `511`. Those hypotheses hold, so this \(J\)
has a unique smooth positive critical point. The distance
`1.2998297915340663e-06` is the distance to that point. It remains an
initial-chart bound.

## Shift current

With \(\chi=p_r=p_\chi=0\) and constant \(Q\), the shift constraint is
\(-Q p_Q'+j\). A periodic \(p_Q\) has mean-zero derivative, so the mean
of \(j\) remains in the residual. The saved columns satisfy
`partner = ω conj(plus)` exactly, and the three occupation pairs are
equal. For the odd symbol \(P=-i\partial_x\),

\[
j=-M\big(1-|\omega|^2\big)\,\mathrm{Re}\big(\overline{\psi}_+ P\psi_+\big).
\]

The stored phase factor has exact binary defect

\[
1-|\omega|^2=-\frac{4585061804339445}{81129638414606681695789005144064}.
\]

That factor is not treated as zero. The enclosed mean of \(j\) is
`4.2386437564048724e-17` to `4.238643756404882e-17`, which does not
contain `0`. The same identity and coefficient \(\ell^1\) bounds give
\(\|j\|_\infty\le 1.2587292180034653\times 10^{-14}\). The geometry-band
tail, in \(L^2\), is at most `3.560223862911562e-14`. Divided by the
one-percent areal-radius bar, the uniform current is
`6.808120926930189e-10` of that bar. The source array is unchanged and
its mean is not deleted. The enclosure used `0.909995` s.

The mean is incompatible with an exact periodic solution at
\(p_r=0\). The existing action already contains \(p_r\). The fixed
profile \(p_r=\lambda r'_{\mathrm{saved}}\) with

\[
\lambda=-\frac{\int j}{\int (r'_{\mathrm{saved}})^2}
\]

uses \(\lambda\in[-1.2038809379699564\times 10^{-15},\,-1.203880937969953\times 10^{-15}]\).
Then \(j+\lambda (r')^2\) has mean zero, so a periodic \(p_Q\) cancels
the continuum shift residual. The consistent lapse update is
\(p_r^2/(4ZQ)\), and its uniform size is at most
`9.288637449732407e-32`. Moving from the saved radius to the critical
point changes that mean by at most `3.857333633418152e-21`. This
correction does not alter \(j\).

`evolution_error_bound` remains null. The initial chart is still not a
propagated episode error.

## Exact nearby initial state

The saved state still has `p_r = 0`. The correction is not installed
and no record is reset. For the fixed profile `f = r'_saved` and
`p_r = λ f`, the owned lapse constraint gains `p_r²/(4 Z Q)`. Because
`C = (8πA/Q)(P_geom + G)` and `Z = feedback_Z(A) < 0`, that addition is

\[
G_\lambda=G+\frac{p_r^2}{32\pi A Z}.
\]

On `|\lambda|\le 2\times 10^{-15}` the change in `G` is at most
`5.694312657413407e-32`, so `G_λ` stays above `0.053733032411766975` and each such
source still has a unique positive critical radius. With
`D_0=\int f r_0'\,dx` bounded below by `0.28166210833983585` from the
sealed `H^1` distance,

\[
H(2\times 10^{-15})\ge 9.024157792148394\times 10^{-16},
\qquad
H(-2\times 10^{-15})\le -2.242327781900592\times 10^{-16}.
\]

The scalar intermediate-value theorem gives a root `λ_*` in that
interval. A periodic `p_Q` is then the antiderivative of
`(λ_* f r' + j)/Q`, and both initial constraints are exact for that
nearby state. Its areal radius stays within
`4.06084969493082e-06` of the saved radius, and
`|p_r|\le 9.352990036589177\times 10^{-16}`. The source `j` is
unchanged.

The scalar margins use outward interval arithmetic. The source perturbation
is evaluated at the exact critical point, whose minimum is bounded by the
saved minimum minus the sealed Sobolev distance. It is not evaluated with
the saved minimum alone. For either radius comparison, the identity
$r_1'-r_0'=2y_1(y_1-y_0)'+2(y_1-y_0)y_0'$ bounds its derivative using
the $H^1$ distance, the Sobolev constant, and the enclosed derivative norm.
This supplies the two integral margins above without assuming a pointwise
derivative bound for the unknown critical point. The sealed initial-error
record remains unchanged.

## Constraint transport

Under the owned Hamiltonian `H=∫(L C+β D)`, with lapse and shift
prescribed and with no Dirac force, the geometric constraints satisfy

\[
\begin{aligned}
C_t&=\beta C_x+\beta_x C+(L/Q^2)D_x+(2L_x/Q^2-2LQ_x/Q^3)D,\\
D_t&=\beta D_x+2\beta_x D+L_x C
\end{aligned}
\]

identically off shell. A `p_Q` force `-F_Q` adds
`-p_χ F_Q/(2 F_χ)` to `C_t` and `Q ∂_x F_Q` to `D_t`. The undifferentiated
`L^2` product of `(C, D)`, and the same product with `D/Q`, both leave
`(L/Q^2)∂_x D`. On a trigonometric band of degree `N` that derivative is
at most `(2π N/L)\|D\|_2`, so a band-limited constraint estimate closes
with a coefficient proportional to `N`. The saved series stores maxima,
not those coefficients. That estimate would still be a constraint budget.
It is not the physical observable error.

The quadrature constraints are `Ĉ = C_geom + force_L/dx` and
`D̂ = D_geom + force_beta/dx`. The owned observer map is
`ρ = force_L/(4π r^4 Q dx)` and `j = -force_beta/(4π r^4 Q^2 dx)`, so
the physical residuals are `Ĉ/(4π r^4 Q)` and `-D̂/(4π r^4 Q^2)`.
Their absolute proper-volume integrals are `∫ |Ĉ|/r dx` and
`∫ |D̂|/(r Q) dx`. The coordinate smearing `∫ |L Ĉ + β D̂| dx` is an
energy-unit budget for the constraint Hamiltonian. None of these is a
radius error or a proper-velocity error, and none is a trajectory-error
certificate.

Using the recorded quadrature maxima and the stored positive `r_min` and
`Q_min`, the four continuations from `T=0.05` to `T=0.085` give the
following source budgets. The static lapse and shift integrals used in
the smearing are at most `5.489130825200043` and `7.280257677860477`.

| Case | `max \|Ĉ\|` | normal-energy majorant | smearing majorant | field-energy exchange | smearing/exchange |
|---|---:|---:|---:|---:|---:|
| `nf256_dtmax_0_0005` | `0.009241190100823481` | `0.017534289186211676` | `0.059419715229532` | `0.2860794429811655` | `0.2077035476940718` |
| `nf256_dtmax_0_00025` | `0.009241182567421104` | `0.017534274892288232` | `0.05942027089732437` | `0.28607942497508176` | `0.20770550312207886` |
| `nf512_dtmax_0_00025` | `0.0001395734148426783` | `0.00026483195364246294` | `0.0007892382299303408` | `0.2860794297480709` | `0.002758808036723804` |
| `nf512_dtmax_0_0005` | `0.00013953044187636887` | `0.0002647504150870611` | `0.0007889809879183477` | `0.28607944774320515` | `0.0027579086653808297` |

These majorants multiply the recorded maximum by a positive weight.
They are not continuum maxima. On `nf512` the matter chain-rule Ward
maximum remains `2.776235187920695e-17`, the source gap is `0`, and
unitarity reaches `7.098334142696672e-10`. The missing observable
primitive is the map from this source-normalized budget and a stored
Hermitian geometry mismatch `δH(t)` to the physical observables.

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_cauchy_weak_bound.py -q
```

Owned bound tests and the two existing weak-test modules passed together
in `8.34` s. The historical lapse enclosure remains the recorded
`20.465117` s run. The shift enclosure above is the separate
`0.909995` s pilot. The nearby-state proof took `1.187442` s.

The sealed record keeps the generator hashes below. They are not the
current files. Replay does not rewrite them. The record bytes are
`ce5e4345f75b187661f965f96c680a3076331e2df709845c22251edb99b3f148`.

| Path | sha256 |
|---|---|
| `lab/src/recursive_horizons/nsc_spherical_cauchy_weak.py` | `14e8c4de9028385a44b903391d07815b300a55caf82a158767141593757e3051` |
| `lab/tests/test_nsc_spherical_cauchy_weak.py` | `0b82dbad4ad5ba74c066831be350a12199b2da299f3ceea8b611f8bb344de790` |
| `lab/scripts/derive_nsc_spherical_cauchy_weak.py` | `c699f641b1a246e3336548bb0063fdc39cc2f3da5802e0c5715d4976500f7582` |
| `lab/results/development/nsc-spherical-cauchy-weak-v1.json` | `ce5e4345f75b187661f965f96c680a3076331e2df709845c22251edb99b3f148` |
