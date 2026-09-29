# Signed subgap Jost phase transport

This is a successor of the positive-only finite-band proof in
[nsc-massive-jost-transport-bound.md](nsc-massive-jost-transport-bound.md).
It transports an original subgap phase inward when the phase derivative is
allowed to be negative. Inward propagation may then grow by a finite factor.
The local gate stays **OPEN**. Original source columns are not replaced, and
this page does not rederive the metric, the horizon frame, or the source law.

Implementation: `src/recursive_horizons/nsc_massive_jost_mixed_transport.py`.
The predecessor module is not modified.

## What was blocking row 15

`enclose_phase_transport_cell` requires a strictly positive lower bound on

$$
f_{y,\theta}=\frac{2e^{y}}{\sqrt{A}}
\Bigl(m\cos\theta-\frac{\ell\sin\theta}{\sqrt{1+r^{2}}}\Bigr)
$$

and raises `cell metric or phase monotonicity is unresolved` otherwise.
On original `group14/low16_1` row 15,
$E=0.24867511687395624$, $m=\pi/2$, $\ell=+\sqrt{5}$, that hypothesis
is false on 23 of the 843 DOP853 cells that remain at least $1.01\times10^{-4}$
outside the coordinate horizon. The most negative enclosed rate on that band is

$$
k=-h\lambda\approx -0.03056.
$$

The growth factor of one such cell is about $e^{0.031}$, not an infinite
instability. The predecessor still rejects a negative-derivative cell. This
owner accepts it only when the metric enclosure and the phase tube both close.

## Retained comparison

The independent variable, metric series, Hairer interpolant, normalized
residual, and exterior initializer are the predecessor's:

$$
\delta(x)=\theta_x-h f_y,\qquad h=y_{\mathrm{end}}-y_{\mathrm{start}}<0.
$$

Dyadic subdivision changes only the Taylor majorant of that same residual.
The original cell width enters the rate once. Metric tails remain the
composed majorant from the predecessor helpers. The outer error is still
`subgap_phase_bound` at the stored initial node and the solver's actual angle.

Let $\lambda$ be a lower bound of $f_{y,\theta}$ on the phase tube, and
set $k=-h\lambda$. The constant-rate majorant obeys

$$
|e(1)|\le e^{-k}|e(0)|+\varepsilon\,\varphi(k),
$$

where

$$
\varphi(k)=\int_0^1 e^{-kt}\,dt.
$$

Thus $\varphi(0)=1$ and $\varphi(k)=(1-e^{-k})/k$ for $k\neq 0$. A rate
ball containing 0 uses the integrated exponential series when it is smaller
than 1, and the monotone endpoint values when it is wider. Neither branch
divides by an interval containing 0. Positive $k$ is contraction, negative
$k$ is finite expansion, and $k=0$ adds the whole defect.

The derivative of that majorant is $e^{-kx}(\varepsilon-k e_{\mathrm{in}})$,
which does not change sign on the cell. The maximum is therefore at an
endpoint. The cell is rejected if either endpoint leaves the tube, if the
normalized defect does not lie strictly inside the tube, or if $A>0$ and
$r>1$ are not enclosed. Those failures raise. They do not become a bound.

## Row 15 composition

The pilot follows `/tmp/nsc-full-low-row-pilot.py` for this single row and
both archived energy signs. It does not scan a second row and it does not
write a result record.

Settings: original `solve_jost` with order 8, collar $10^{-10}$,
`rtol=2e-13`, `atol=2e-15`; transport degree 16, 48 metric terms, 4 defect
subdivisions, tube $10^{-5}$, 192 bits; inner request
`horizon_rho+1.01e-4`. The certified inner offset is about $1.10886\times10^{-4}$.
After that node the existing horizon frame, `reflection_from_phase`, Bloch
capture, and `validate_bloch` are called unchanged. Positive and negative
covariance errors use the archived columns for row 15.

| Quantity | Observed upper |
|---|---|
| Cells (contracting / expanding / zero / straddling) | 843 (820 / 23 / 0 / 0) |
| Initializer phase error | $1.778\times10^{-15}$ |
| Maximum normalized cell defect | $4.993\times10^{-10}$ |
| Inner phase error | $1.814\times10^{-9}$ |
| Frame reflection tail | $6.460\times10^{-22}$ |
| Bloch error (671 cells, 10067 evaluations) | $2.027\times10^{-10}$ |
| Positive covariance error | $1.047\times10^{-10}$ |
| Negative covariance error | $1.047\times10^{-10}$ |
| CPU / wall time | 6.88 s / 6.90 s |

The pilot calls the row too loose only if the phase upper bound reaches
$10^{-4}$ or either covariance upper bound reaches $10^{-3}$. This row
does not. Those scales are a separation check, not a gate certificate.

The two covariance uppers are the same order as the Bloch error. They do not
show an order-one mismatch with the archived columns. They also do not change
those columns. Zero and straddling rates did not occur on this row; the tests
cover them on manufactured cells and on the analytic weight.

## Still open

One bounded row is not a source campaign. The phase owner stops at the
declared inner radius. The collar from that node is the existing frame's
evaluation, not a new sewing theorem in this file. No physical $\rho=1$
source column is replaced. The local gate remains OPEN.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_massive_jost_mixed_transport.py
.venv/validation/bin/python scripts/lab.py scripts/pilot_nsc_massive_jost_mixed_transport.py
```

Parent integration rejects missing and boolean proof inputs before Arb conversion.
In particular, None cannot become a zero error or zero rate. Numeric zero remains
valid where the analytic formula permits it. Regression tests cover this guard.

Independent Grok review passed all 15 tests and reproduced the one-row pilot
without finding an under-enclosure. This accepts the scoped method, not a full
source-family result. Saved-witness production records are a separate next step.
