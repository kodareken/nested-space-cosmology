# Classical RN reference

The owner is
[`nsc_rn_reference.py`](../src/recursive_horizons/nsc_rn_reference.py).
It fixes the action dictionary, the ingoing Reissner–Nordström
Painlevé–Gullstrand solution, an independent curvature reduction, and a
read-only assessment of sealed source-free states. It does not evolve a
shell, write a campaign record, or modify sealed NPZ or JSON.

## Action

The classical bulk term, signature \(+---\), is

\[
S=\int\sqrt{|g|}\,\bigl[-A R_L-C_F F^2\bigr].
\]

There is no Weyl-squared pole and no cosmological term. The executable
flux is the \(q\) of the spherical reduction, not the prose label \(2\pi q\):

\[
G_N=\frac1{16\pi A},\qquad
r_Q^2=\frac{C_F\,\mathrm{flux}^2}{4A}.
\]

\(A>0\) and \(C_F>0\). The Maxwell coupling recorded with this reduction is
\(g_4^2=1/(4C_F)\). Direct contraction of the resulting magnetic RN metric,
with the laboratory Riemann convention

\[
R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{cb}
+\Gamma^a{}_{ce}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{cb},
\qquad R_{bd}=R^a{}_{bad},
\]

gives \(R=0\), because the Maxwell stress is traceless,

\[
R_{ab}R^{ab}=\frac{4 r_Q^4}{r^8},
\]

\[
K=\frac{48\,\mathrm{mass}^2}{r^6}
-\frac{96\,\mathrm{mass}\,r_Q^2}{r^7}
+\frac{56 r_Q^4}{r^8}.
\]

The cross term in \(K\) is fixed by that contraction; it is not copied in
from a flipped charge sign. On the ingoing orthonormal frame

\[
R_{\hat0\hat1\hat0\hat1}=\frac{2\,\mathrm{mass}\,r-3r_Q^2}{r^4},
\qquad
R_{\hat0\hat2\hat0\hat2}=R_{\hat0\hat3\hat0\hat3}
=-\frac{\mathrm{mass}\,r-r_Q^2}{r^4}.
\]

## Chart

\[
ds^2=N^2dt^2-(dr+\beta\,dt)^2-r^2d\Omega^2.
\]

For the source-free solution the Einstein equation and unit asymptotic
Killing normalization give \(N=1\) and

\[
f=1-\frac{2\,\mathrm{mass}}r+\frac{r_Q^2}{r^2},\qquad
\beta=\sqrt{\frac{2\,\mathrm{mass}}r-\frac{r_Q^2}{r^2}}=\sqrt{1-f}.
\]

\(N=1\) is that vacuum result. A coupled shell must derive the lapse; it is
not kept at 1. The shift is real for \(r\ge r_Q^2/(2\,\mathrm{mass})\).
Horizons are the roots of \(r^2-2\,\mathrm{mass}\,r+r_Q^2=0\). Both radial
null speeds \(dr/dt=1-\beta\) and \(dr/dt=-(1+\beta)\) are negative on the
open interval between those roots, which is the pure-outflow excision.
The outer boundary is nonperiodic and absorbing. The parent interval lies
outside the outer horizon. The child collar is near that horizon and uses
a different areal interval.

\(V_4=0\). \(G_N\) and \(P^2=r_Q^2\) are the fixed-charge branch. The bare
Misner–Sharp mass on this signature is

\[
m=\frac r2\bigl(1+\nabla^ar\nabla_ar\bigr).
\]

On the vacuum chart \(\nabla^ar\nabla_ar=-f\), so
\(m=\mathrm{mass}-r_Q^2/(2r)\). The charged mass is \(m+r_Q^2/(2r)\).

The adapter entry point is `sourcefree_static_rn(radius, mass, *, A, C_F, flux)`.
It returns that static chart at the requested radius: horizons, metric jets,
invariants, tetrad components, Misner–Sharp mass and charged mass. `field_state`
and `campaign` are `None`. It does not build a Dirac packet.

## Momenta and the conformal gradient

A sealed parent checkpoint is read with `load_checkpoint`. Decode accepts only
`momentum_representation="canonical_pi"`. A missing representation is rejected.
With the stored frame \(W\) and spacing \(dx_g\),

\[
p=W(\pi/dx_g).
\]

\(dx_g\) is the stored spacing, or \(8/n\) when the record's \(n_f\) matches
\(W\) and \(n=n_f-1\). The strong empty control has \(n=127\) and
\(dx_g=8/127\). The declared anchor \(k\) is not an input, so changing it
does not change \(p\).

On the historical conformal chart \(N=rQ\) and \(a=-8\pi A\),

\[
r_t=\frac{Q p_Q}{2ar},\qquad
\nabla r\cdot\nabla r=\frac{p_Q^2}{4a^2 r^4}-\frac{r_x^2}{r^2 Q^2}.
\]

The charged mass uses that gradient on the decoded \(p_Q\). It is not fitted
to \(k\) or to a horizon. A flat sealed radius contributes \(r_x=0\). \(k=0.02\)
has no sealed checkpoint.

## Matter control

The characteristic basis still has \(\sigma_2=\mathrm{diag}(1,-1)\), and the
angular operator on the chart is \(\kappa/r\). `manufactured_quadrature_control`
only checks SBP weights and the weighted Gram. It sets
`positivity_unestablished` and is not a prepared source with \(E>0\). The
scattering-packet owner supplies that preparation.

## Curvature status

Analytic RN invariants above are the source-free reference. They are not a
coupled calibration. Historical conformal \(R_4\) is unresolved: this slice
does not call `leading.metric_jets`. An owned Fourier estimate previously
reached about \(0.757\), while prior authenticated metric jets on the same
family are about \(0.167\) to \(0.160\). Neither value is adopted. The saved
Bertotti–Robinson control is authenticated and is not used as a rescaled RN test.

## Writer

`write_assessment` can create only
`results/development/nsc-rn-reference-assessment-v1.json`, outside the sealed
input directories. It requires a frozen 40-character commit and the
assessment `input_hashes`, refuses overwrite, keeps the file read-only, and
rejects a payload above 64 MiB. This slice does not call it, and that file
is not created.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_rn_reference.py -q
```
