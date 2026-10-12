# Directed radial subgap Jost phase transport

This is a Stage-2 source-accuracy method. It reuses the existing massive
Riccati/phase equation and the original `solve_jost` DOP853 interpolants. It
supplies a finite-band phase error from the solver-selected outer radius `R`
down to a declared inner radius lying strictly outside the horizon. The r>=8
band is the control; a second bound reaches an offset just above 1e-4. It
does not enclose the near-horizon collar, matched sewing, or a physical
`rho=1` source column.

Implementation: `src/recursive_horizons/nsc_massive_jost_transport_bound.py`.

## Equation and domain

Here `r` is the original exterior radial coordinate, with areal radius
`sqrt(1+r^2)` and horizon value `ĥ = horizon_rho`. The original subgap
integrator uses

\[
y=\log(r-\hat h),\qquad
\frac{d\theta}{dy}=f_y(y,\theta)
=\frac{2e^{y}}{A(r)}\Bigl(\sqrt{A}\Bigl(\frac{\ell\cos\theta}{\sqrt{1+r^2}}+m\sin\theta\Bigr)-E\Bigr).
\]

The metric is the exact NSC expression

\[
A(r)=1-3\bigl[(1+r^2)\arctan(1/r)-r\bigr],
\]

never a broad interval subtraction of that formula. On the validated band the convergent
expansion in `u=1/r` is used, with the same positive majorant of the omitted
odd powers as the exterior half-line bound:

\[
A=1-6\sum_{j\ge0}\frac{(-1)^j u^{2j+1}}{(2j+1)(2j+3)}.
\]

The independent variable of the interpolant is the original `y`. A declared
inner endpoint with offset at least `10^{-4}` is required, so the band never
meets the horizon. A>0, r>1, positive phase derivative and tube invariance
must be enclosed on every included cell; r>1 is not itself a horizon proof.
The remaining collar below the certified endpoint, matched sewing and interior
transport remain missing.

## Interpolant

`solve_jost` retains the original SciPy `OdeSolution` in the `run.sol` closure
named `dense`. Each local interpolant is the DOP853 degree-7 Hairer polynomial
in the cell fraction `x=(y-y_{\mathrm{old}})/h`, `h<0` inward:

\[
\theta(x)=(1-x)\theta_{\mathrm{old}}+x\theta_{\mathrm{end}}
+\sum_{j=1}^{6} F_{j}\,x^{p}(1-x)^{q},
\quad p=\lfloor(j+2)/2\rfloor,\; q=\lfloor(j+1)/2\rfloor.
\]

The stored endpoints are solver nodes: the current interpolant `y_old` and the
*next* interpolant `y_old`. Their binary64 addition check matches
`y_old+F[0]`; the proof uses the exact dyadic difference of the stored
endpoints instead of treating the rounded subtraction F[0] as exact. The six
correction vectors are `F[1:]`, the same basis as `TrajectorySegment`.

## Normalized defect and contraction

On each cell, the interpolant defect in the normalized coordinate is

\[
\delta(x)=\theta_x-h f_y\bigl(y_{\mathrm{old}}+hx,\theta(x)\bigr).
\]

Taylor expansion at `x=1/2` plus a rigorous remainder from the `N`th coefficient
on `x\in[0,1]`, with factors `2^-j` and `2^-N`, encloses `\varepsilon=\sup|\delta|`. Metric tail derivatives are
*composed* with `u=1/r(y)`: if `D_k` are the `u`-derivative uppers of the
majorant `6u^{2J+1}/(1-u)` at `1/r_{\min}`, the Faà-di-Bruno generating series
`\sum (D_k/k!)w^k` is composed with the positive series of
`|d^j(1/r)/dx^j|/j!`. Copying `D_k` into the `y` (or cell-`x`) basis is not a
proof.

The mean-value error `e=\theta_{\mathrm{true}}-\theta_{\mathrm{interp}}` satisfies
`e_x=h f_{y,\theta}e-\delta` with

\[
f_{y,\theta}=\frac{2e^{y}}{\sqrt{A}}\Bigl(m\cos\theta-\frac{\ell\sin\theta}{\sqrt{1+r^2}}\Bigr).
\]

For forward `x` and negative `h`, the contraction rate is `k=-h f_{y,\theta}`.
A positive lower bound `\lambda\le f_{y,\theta}` on the phase tube gives
`k\ge -h\lambda`. Then

\[
|e(1)|\le e^{-k}|e(0)|+\varepsilon\frac{1-e^{-k}}{k},
\]

and the cell is accepted only if this bound stays strictly inside the tube.
The degree is at least seven so differentiating the retained phase polynomial
does not omit a coefficient required by the remainder. Cell widths are exact
differences of the stored dyadic times. The initial radius is enclosed as
`horizon_rho+exp(y_start)` rather than assuming an exact log/exp round trip, and
the actual first solver angle is used, including a different series-order seed.

The outer `|e(0)|` is the existing entire-exterior initializer bound from
`nsc_massive_jost_phase_bound.py`. Failure of positivity, metric enclosure or
tube invariance raises; it does not emit a conditional wrapper in place of a
bound.

## Coverage and remaining work

The representative original group-14 label is
`E=0.0013248831260437577`, `m=\pi/2`, signed `\ell=\sqrt{5}`, with the actual
solver-selected `R` (here `60`). Analytic `π/2` and `±√5` may be enclosed
together with those binary labels. Directed bounds are exact dyadics. Global
Arb precision and series cap are restored.

This does not claim a horizon-collar bound, sewing, interior transport, source
covariance arithmetic or a physical `N/β` component. The incoming gate remains
OPEN.

## Subdivision and recorded result

A dyadic subdivision changes only the enclosure. Each subinterval uses its
midpoint and half-width in Taylor's theorem on the SAME normalized residual.
The maximum subinterval residual and minimum positive phase-derivative lower
are valid over the complete original cell. Its original step remains in the
contraction formula; no partition-width factor is applied a second time.

The source inventory supplies the original group-14 energy, mass, angular label
and background configuration. The driver saves the exact binary phase-node
and dense-correction trace in a small NPZ witness. It recomputes the original
mode and all bounds during --check, comparing the record, witness and complete
source hashes. The physical source columns are never replaced.

At the last certified node, the radial-coordinate offset is approximately
1.0218052894e-4. Degree16 with48 metric terms gives an error upper about
6.40e-11 with one bounding interval per cell; four bounding intervals reduce
that upper to about9.04e-12 for the SAME numerical trajectory (906 cells).
These are phase errors, not N/beta constraint errors. Sewing and interior
transport to rho=1 remain unbounded. The gate remains OPEN.

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_massive_jost_transport_bound.py --check
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_massive_jost_transport_bound.py
```
