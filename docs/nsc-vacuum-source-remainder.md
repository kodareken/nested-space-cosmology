# A uniform high-energy preparation covariance remainder

This owner bounds the homogeneous source covariance at the declared upstream
slice. It reuses the exact metric and normalized horizon branch from the
[horizon-frame proof](nsc-metric-horizon-frame.md) and the Bloch equation from
the [source covariance method](nsc-subgap-source-covariance.md). It does not
replace the source, compute the inheritance identity, or bound the evolved
changed-history ultraviolet tail.

## Exact vacuum branch and finite occupations

For positive energy the zero-temperature three-port covariance is
$\operatorname{diag}(0,1,0)$. The original sewing matrix maps it to
$e_0e_0^\dagger$, independently of the reflection and transmission. The
constant phase conventions leave that projector invariant. The homogeneous
interior fundamental matrix is unitary, so its vacuum covariance has trace
one and can be written $Q=(I+\mathbf n\cdot\boldsymbol\sigma)/2$.

Set $\delta=q_h-q=t^2$, $y=\log\delta$ and
$H(\delta)=W(q_h-\delta)/\delta>0$. The radius is
$r=\csc(q_h-\delta)$. With $L=\sqrt H(mr-i\ell)$ and
$w=n_x+in_y$, the existing full-metric Bloch equation becomes

$$
w_y=-\frac{2iE}{H}w+\frac{2itL}{H}n_z,
\qquad (n_z)_y=-\frac{2t}{H}\operatorname{Im}(\overline Lw).
$$

The vacuum boundary condition is $(w,n_z)\to(0,1)$ as $t\to0$.
This is a condition on the normalized horizon branch, not on a finite-radius
approximation or a new choice of source occupation.

## Recurrence and the equation defect

Write $w\sim t\sum_{j\ge1}P_jE^{-j}$ and
$n_z\sim\sum_{j\ge0}Z_jE^{-j}$, with $P_0=0$, $Z_0=1$.
Coefficient matching and the unit-norm identity give

$$
P_{j+1}=\frac{iH}{2}\left(\frac{P_j}{2}+\delta\partial_\delta P_j\right)+LZ_j,
$$

$$
Z_n=-\frac12\sum_{j=1}^{n-1}
\left(\delta\overline{P_j}P_{n-j}+Z_jZ_{n-j}\right).
$$

The sum in the second equation is real: conjugate terms occur in pairs.
The implementation uses real and imaginary jets, retaining only derivative
orders that remain known after each recurrence. The exact metric's omitted
series and derivative tails are included by the existing `metric_H` owner.

For the truncation $\mathbf n_M$, the transverse equation defect is exactly

$$r_w=-\frac{2it}{H E^M}P_{M+1},\qquad r_z=0.$$

The second identity holds to every retained order. To see why, the recurrence
makes $|w_M|^2+(n_z)_M^2-1$ have no powers through $E^{-M}$. Its derivative
has the same property. The identity
$\frac12\partial_y|\mathbf n_M|^2=\mathbf n_M\cdot\mathbf r$
then implies that every coefficient of $r_z$ through order M vanishes:
$w_M r_w$ starts only at order M+1, while $(n_z)_M$ starts with 1.
There are no higher powers in $r_z$, which is linear in the truncation.
The finite truncation is not renormalized; renormalization would change this
equation defect.

## Continuous remainder, uniform in energy

The exact Bloch generator is skew-symmetric. Its error propagator has Euclidean
norm one. Both the exact solution and its expansion have the same horizon
limit, and the defect is integrable there. Since $dy=2dt/t$,

$$
\|\mathbf n-\mathbf n_M\|_2
\le\frac{C_M^{\rm prep}}{E^M},\qquad
C_M^{\rm prep}=4\int_0^{\sqrt{\delta_{\rm up}}}
\frac{|P_{M+1}(t^2)|}{H(t^2)}\,dt.
$$

The covariance operator error is at most half this value. The constant is
enclosed by whole-cell Arb ranges times each cell width, with a positive
lower bound for H. This is a directed integral upper bound. Refinement
differences and successive expansion orders are not used as error estimates.
The bound is uniform for every positive E; the high-energy pilot reports its
value at E=160 and its formula for all larger E. A separate comparison with
the original retained modes uses the same valid formula at their own energies.

The complete source also has finite occupations and horizon coherence. Using
the unchanged source scales, unitarity of the homogeneous frame and the
owned $|R|^2+T=1$ current identity gives

$$
\|Q_{\rm source}-Q_{\rm vac}\|_{\rm op}
\le f+n_{\rm in}+c_h
\le e^{-2\pi E/\kappa}+e^{-2\pi E/(\omega\kappa)}+e^{-\pi E/\kappa}.
$$

The declared massive threshold can only reduce the incoming occupation.
The source kappa is not replaced by the geometric horizon slope. This bound
does not require inventing an exterior reflection coefficient.

## Evidence and limits

Exact manufactured controls verify the transverse and longitudinal defects,
normalization coefficients and the factor four in the continuous integral.
A numerical integration from an independently enclosed full-metric horizon
frame checks the expansion against the original ODE. That comparison is a
diagnostic; the error certificate is the analytic defect bound. Controls also
reject absent bounds, which must not be converted to zero by Arb's constructor.

The order-8, 64-cell bound at E=160 is below 8.46e-15 for the covariance.
This is a bound relative to the declared analytic expansion, not automatically
the error of the archived numerical columns.

For the archived columns the owner explicitly adds the operator distance
between that expansion and the stored $A_{\rm up}C_{\rm src}A_{\rm up}^\dagger$.
The original source and preparation digests are retained per row. An order-16,
128-cell pilot covers 384 positive-energy rows with E>=32 and their 384
opposite-angular negative partners. For a negative partner the expansion is
complemented using $Q_- = I-\sigma_3\overline{Q_+}\sigma_3$, then compared with
the actual negative columns. This includes their numerical Gram discrepancy.
No source array is replaced, no preparation ODE is repeated, and no quadrature
weight is applied to these unweighted covariance errors.

The record concerns the original group-14 positive angular source family and
its signed partners at the upstream slice. Lower energies and the other
positive angular family are not covered by this pilot. It does not give the whole changed-history field remainder
called C_M in the UV gate records, a same-column H2 bound, or an integrated
N/beta ultraviolet tail. Those remain separate missing connections. Source
covariance error must be propagated through the actual local response and
contracted with its true spatial derivative. The physical gate stays OPEN.

Owners: `src/recursive_horizons/nsc_vacuum_source_remainder.py`,
`scripts/derive_nsc_vacuum_source_remainder.py`, and the v1 result record.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_vacuum_source_remainder.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_vacuum_source_remainder.py --check
```
