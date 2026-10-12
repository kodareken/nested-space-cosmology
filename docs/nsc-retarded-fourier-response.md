# Retarded Fourier-pair matching response

The full first-order covariance coefficient at a selected pair is

\[
\delta\widehat C(E_o,E_i)
=\delta A(E_o,E_i)C_i f_i^\dagger
+f_oC_o\delta A(E_i,E_o)^\dagger,
\qquad
\delta A(E_o,E_i)=\int dz\,e^{iE_o z}\delta F(z,E_i).
\]

For the tested group14 radius direction, the fine diagonal pair at
`Eo=Ei=0.5562120090641313` has spin00 coefficient
`-1.0554769134063148e-5`. Spatial and temporal pair differences are
`6.91800e-9` and `2.46272e-11`; the CAR-pair residual is `3.58011e-11`.
This is stable numerical evidence of an incoming state response, with the
rigorous matching certificate still OPEN.

The [analytic matching connection](nsc-incoming-fourier-matching.md)
performs the source-energy integral exactly at first order using the
stationary background plane waves. Thus this selected-pair calculation
has no omitted source-energy quadrature. It uses the same three source
columns and original coherent covariance at the two selected energies.
It remains one angular block and one metric variation, not the whole action
or a test of matching at every pair.

## Reuse and numerical prescription

Reuse the authenticated raw traces from the
[targeted spatial resolution](nsc-retarded-response-resolution.md). No field
or radial propagation is performed. The stationary amplitude f comes from
the initial authenticated global reference columns and the fixed surface
restriction, with its analytic `exp(i E S1)` conversion. The evolving,
drifting SBP reference is not used as a replacement for that amplitude.

Only deltaF is numerically transformed. The reference plane-wave transform
is kept analytic on the whole line. This is not a double transform of the
old four-frequency position kernel, which would introduce finite-window
sinc factors. In particular, the old energy quadrature factors
`sqrt(weight/(2*pi))` and angular multiplicities do not enter this pair.

The compact variation has exact incoming trace support inside
`(.0525,.2475)`. The full recorded PG interval `[0,.3]` therefore covers the
continuum response. All numerical trace values on that interval, including
leakage tails, are retained. Simpson integration is compared with trapezoidal
integration; phase-origin covariance is checked independently between tau
and z conventions. The toy normalization test uses the exact Gaussian
Fourier transform and a canonical phase modulation, not the model output.

The exact complete-source CAR pair identity is evaluated separately by
replacing C with the original source projector P. This identity is distinct
from the earlier nonzero, truncated position-projector kernel. Its numerical
residual is reported and never subtracted from the covariance coefficient.
The closed third input source column remains exactly zero in this band.

## Decision and scope

The [record](../results/development/nsc-retarded-fourier-response.json)
retains all sixteen complex2x2 covariance and CAR pairs for each of
1601/128,3201/128,3201/256, plus the mode-pair response. Diagnostic limits
are25% spatial,10% temporal, and1% for the CAR/response and time-quadrature
comparison ratios. Algebraic Hermiticity, phase-origin and initial-column
checks use3e-11. These are numerical consistency indicators, not rigorous
continuum error bounds or a tolerance on action stationarity.

A stable nonzero pair supplies numerical evidence that this retarded
extension does not automatically preserve the fixed C0 condition. A
rigorous exclusion would additionally need a response error bound that
excludes zero. A zero selected pair would not prove matching everywhere.
Other extensions and other metric directions remain available; this test
cannot reject the nested architecture or a different parent matching.

Full global constraints and extended stationarity remain OPEN. No source
law, scale, seed, stress, action term or physical duration is changed, and
no metric equation is evolved. The original spatially unresolved pilot
and its numerical results remain intact.

```sh
python3 scripts/derive_nsc_retarded_fourier_response.py --check
python3 -m pytest -q tests/test_nsc_retarded_fourier_response.py
```
