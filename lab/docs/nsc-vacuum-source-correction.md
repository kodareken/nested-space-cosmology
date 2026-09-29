# An enclosed correction to the preparation expansion

The [uniform source expansion](nsc-vacuum-source-remainder.md) is effective at
high energy. Its absolute defect integral is looser near the lower end of the
retained middle panels. This successor evolves only the small correction to
that expansion and encloses its actual equation defect. It compares with the
same archived source; it does not replace the original preparation.

## Equation and initial uncertainty

Let $A_E(y)n=2\mathbf h_E(y)\times n$ be the existing homogeneous Bloch
generator and let $r_M=(n_M)_y-A_E n_M$ be the known finite-expansion defect.
For $e=n-n_M$,

$$e_y=A_E e-r_M.$$

The exact transverse defect is reused from the previous derivation. The
numerical correction starts at zero at $y_0=-18$, but its true initial value
is not assumed zero. The original uniform remainder method bounds it over
the entire interval from the exact horizon limit to $y_0$:

$$\|e(y_0)\|\le C_M^{\rm prep}(y_0)/E^M.$$

The correction's generator is skew-symmetric. Its forcing is prescribed, so
the error propagator still has norm one. No constant norm is assumed for the
driven correction itself.

## Continuous defect and coordinate composition

Each saved DOP853 cell is an anchored polynomial $p(x)$ on $0\le x\le1$,
with $y=y_0+h x$. Its normalized residual is

$$d=p_x-2h\,\mathbf h_E(y)\times p+h r_M(y).$$

The forcing is evaluated from coordinate derivatives of $P_{M+1}(\delta)$,
then composed with $\delta=e^y$. Its expansion increment has identically zero
constant coefficient; interval subtraction is not mistaken for that identity.
Midpoint Taylor coefficients and a whole-cell derivative bound enclose the
residual norm continuously. Summing the normalized cell integrals adds no
second width factor.

The final endpoint is $n_M(y_{\rm end})+p(1)$. The difference between the
binary log endpoint and the exact target rho is bounded using the exact
vacuum vector's norm one. Thus

$$\|n_{\rm true}-n_{\rm endpoint}\|
\le\epsilon_{\rm initial}+\sum\epsilon_{\rm cell}+\epsilon_{\rm bridge}.$$

The stored curve is a proof witness. Replay revalidates it without asking
the integrator to produce the same floating-point trajectory again.

## Comparison with the original source

The pilot is the first original group14/mid24_1 row, at E about 16.01925,
on the archived upstream slice. The correction provides a tight vacuum
covariance enclosure. The separately bounded finite occupations/coherence
are added to obtain bounds for the full source. Direct comparison with the
original matrices gives both lower and upper operator-distance bounds.
The negative partner is compared separately after the exact homogeneous
complement map; numerical Gram discrepancy is retained.

An archived-source error is not a gravitational constraint residual. It must
be propagated and contracted through the same effective-source kernel before
it can enter an N/beta budget. A nonzero preparation discrepancy does not
exclude the declared physical class.

The source arrays, quadrature and occupation law are unchanged. Other rows,
families, quadrature error and the changed-history UV remainder remain outside
this pilot. All full physical budget fields remain null and the local gate
remains OPEN.

Owners: `src/recursive_horizons/nsc_vacuum_source_correction.py`,
`scripts/derive_nsc_vacuum_source_correction.py`, its tests and v1 record/payload.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_vacuum_source_correction.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_vacuum_source_correction.py --check
```
