# Spatial diagnosis of the retarded incoming response

The saved spatial discrepancy is concentrated in the slow characteristic
component near rho1. Applying the same local generator to exactly matching
common-node forcing produces45–52% differences before any time propagation.
Together these diagnostics select **finer spatial evaluation of the same
pulse**, with appropriate reference-column accuracy and time control. They
do not select a wider pulse or a changed source.

Reuse the [three-solve pilot](nsc-retarded-compatible-response.md), its
unchanged global reference columns and the owned SBP(4,2)/SAT generator.
The gap is localization of the failed spatial indicator and discrimination
between forcing/transport resolution, baseline drift and boundary leakage.
Stop after saved-array norms and three local `F,dL,LF` evaluations. No field
propagation, radial solve, old producer, automatic refinement or metric
step is performed by this owner.

## Saved response localization

The maximum final-field difference is1.31007449e-4 at rho0.992,
characteristic spin0, E=-0.5562120091 and source column0. Spin0 carries
99.9436% of the unweighted squared difference norm;90.293% of its error
lies in rho[0.95,1.05]. Its peak moves fromrho0.968 to0.940, but a diagnostic
integer-row shift barely reduces its relative shape discrepancy
(0.9691 to0.9658). The disagreement is not explained by a simple translation.
No field or coordinate is actually shifted.

The incoming-trace difference peaks at1.03074492e-4 at PG time0.253125,
canonical spin0 and the same source column. The coarse trace norm there is
1.04190e-4; the fine trace norm is1.24743e-6. About72.13% of the time-norm
squared error lies aftertime0.21. The covariance difference peaks on the
spin0 diagonal at this same pair of times.

The spatial SBP-norm relative field difference is0.68715. The trapezoidal
time-norm trace relative difference is0.41318, and the finite sampled
covariance Frobenius ratio is0.46102. These finite-column norms describe
shape disagreement; they are not physical-state norms or continuum bounds.
The four energies contribute approximately equally to the squared error;
the closed infinity column contributes exactly zero.

## Drift and boundary diagnostics remain separate

The stationary-reference drift maxima occur at the right boundary in
characteristic spin1. At the tangent-discrepancy maximum, the baseline-drift
difference is only2.25e-8 compared with1.31e-4 for the tangent difference.
Normalized correlations are approximately5.3e-4 for the field discrepancies
and9.6e-5 for the trace discrepancies. This separates their observed patterns;
it does not prove that baseline discretization has no effect on tangents.

Final right-boundary tangents fall from2.1863e-6 to1.2978e-7 with spatial
refinement, while final left-boundary tangents are negligible. Complete
boundary histories were not saved, so their maxima cannot be localized in
time from this artifact. The existing continuum left-going argument does
not remove numerical boundary/dispersion effects.

## Conditional continuum support and numerical leakage

For this radius-only variation, both principal cones remain fixed. Only
sources upstream of rho1 can reach Sigma. The fast and slow travel times
from rho to1 are integrals of `1/(beta+1)` and `1/(beta-1)`. After subtracting
`S(rho)-S(1)=integral beta/(beta²-1)`, the earliest/latest arrivals relative
to the axial pulse edges differ by minus/plus

\[
X=\int_1^{\rho_+}\frac{d\rho}{\beta(\rho)^2-1}.
\]

Thus the continuum tangent trace is supported in `[0.09-X,0.21+X]` under
the stated unchanged-initial-past, compact-source and fixed-cone assumptions.
The cheap numerical quadrature gives X=0.035952045175436084 and the window
[0.05404795482456391,0.24595204517543606]. This edge estimate is not a directed
certificate, and the finite SBP evolution has no asserted exact causal cone.

Outside that window, the saved trace maximum and unweighted squared-norm
fraction are:

| Case | Maximum | Squared-norm fraction |
|---|---:|---:|
|401/64 |1.04190267e-4 |0.06509499 |
|801/64 |7.69596845e-6 |0.00044195 |
|801/128 |8.27776473e-6 |0.00047344 |

The record also gives trapezoidal time-norm fractions, peak coordinates and
distances of sampled nodes from the estimated edges. In particular the coarse
late lobe lies outside the continuum response window: it is numerical leakage
under this conditional support argument, not a physical parent response.
This conclusion is kept separate from a global causality proof or a numerical
error bound.

## Independent local discriminator

At times0.11,0.15,0.19, construct the exact same discrete forcing on both
saved grids:

\[
F=\delta L\,\Phi_{\rm archive}e^{-iEt},\qquad LF=L F.
\]

The pure-radius deltaL is a local potential insertion. Common-node F
comparisons are exactly zero. The relative common-node LF differences are
0.523965997,0.452112747,0.493550129, respectively. No field evolution or
incident forcing integration enters this comparison.

For the owned interior fourth-order centered derivative, the Fourier group
factor is `g(theta)=(4*cos(theta)-cos(2*theta))/3`. A DFT of the same full-grid
geometry factor delta_r/r² gives:

| PG time |401-grid power with abs(g−1)>0.25 |801-grid power |
|---|---:|---:|
|0.11 |0.337125 |0.050146 |
|0.15 |0.365063 |0.023410 |
|0.19 |0.338676 |0.044180 |

The sampled reversed-group-factor fractions are also recorded. This is a
local frozen-coefficient dispersion diagnostic, not a global variable-speed
causality theorem, a continuum error bound or a new periodic boundary
condition. It is consistent with the late coarse-grid slow-component lobe;
it does not by itself prove every contribution to the observed error.

The [record](../results/development/nsc-retarded-spatial-diagnosis.json)
authenticates the saved pilot and reused owners, preserves the prior OPEN
status, and performs local matrix/FFT algebra only. It makes no claim that
the four-energy covariance response decides full C0 or CAR matching. No
sampled projector response is subtracted as a correction; source, pulse,
Gamma allocations, seam and physical scope remain unchanged.

```sh
python3 scripts/derive_nsc_retarded_spatial_diagnosis.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_retarded_spatial_diagnosis.py
```
