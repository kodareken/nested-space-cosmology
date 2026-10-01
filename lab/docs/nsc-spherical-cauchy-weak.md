# Weak initial residual from the spherical action

Diagnostic only. The production evolution remains
`nsc_spherical_galerkin_coupling`. This helper does not step a state, rerun
the radius Newton, or retune a tolerance. Record
`lab/results/development/nsc-spherical-cauchy-weak-v1.json`, status
`PARTIAL_TERMS_NO_TOTAL_BOUND`. Checkpoint `9a9090a`. The saved v5 JSON and
NPZ bytes are unchanged. Assessment CPU `6.787458` s.

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_cauchy_weak.py tests/test_nsc_spherical_cauchy_weak_independent.py -q
```

Owned and independent tests passed together in 7.48 s.

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

| Path | sha256 |
|---|---|
| `lab/src/recursive_horizons/nsc_spherical_cauchy_weak.py` | `14e8c4de9028385a44b903391d07815b300a55caf82a158767141593757e3051` |
| `lab/tests/test_nsc_spherical_cauchy_weak.py` | `0b82dbad4ad5ba74c066831be350a12199b2da299f3ceea8b611f8bb344de790` |
| `lab/scripts/derive_nsc_spherical_cauchy_weak.py` | `c699f641b1a246e3336548bb0063fdc39cc2f3da5802e0c5715d4976500f7582` |
| `lab/results/development/nsc-spherical-cauchy-weak-v1.json` | `974c1614468ae3dbfb1f4e13d1be026b0aa2f99489ac258287fbdfcfc097da72` |
