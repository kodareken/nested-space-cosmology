# Directed enclosure of the reduced massless cubic current

This method encloses the formal history-minus-reference cubic coefficient for
the original massless group-1 pair. The first saved pilot covers **one incoming
point**, for both characteristics and the equal-weight angular pair. It does
not cover all of I, massive channels, or the higher-order UV remainder. The
physical incoming gate and the physical UV budget remain OPEN.

Owner: `src/recursive_horizons/nsc_ks_cubic_uv_enclosure.py`.
The equations come from the six exact
[reduced-current identities](nsc-ks-cubic-uv-current.md#reduced-system-for-the-next-directed-quadrature).

## Fixed geometry and source

The radius remains $r=r_0+\chi(s)(s w+s^3U/6)$, with the original binary
Chebyshev coefficients and domain maps. Arb encloses their analytic evaluation.
The reference chart supplies $a$, $s=-\int_1^\rho a^{-1}$ and
$D=\int_1^\rho a^{-2}$. The characteristic is
$z(\rho)=z_{\rm in}-s_cD(\rho)$, where $s_c=\pm1$ labels the characteristic
and is distinct from the normal coordinate $s$.

The existing background, plateau and Chebyshev jet owners bound derivatives
throughout real boxes. Flat cutoff strips use their existing Cauchy bounds;
unresolved boxes request subdivision. No sampled derivative is used as a
supremum. The upstream endpoint must lie outside the normal history support,
so the preparation remains unchanged. The three initial **differences** are
zero; no unknown physical source-normalization error is assigned zero.

## Coordinate partials before characteristic composition

For $b=a\ell/(2r)$, form the bivariate Taylor coefficients of $r(\rho,z)$.
The reciprocal recurrence is

$$
v_{00}=1/r_{00},\qquad
v_{ij}=-v_{00}\sum_{(p,q)\ne(0,0)}r_{pq}v_{i-p,j-q}.
$$

Here $0\le p\le i$, $0\le q\le j$. Compute the coordinate partials of b,
then compose their jets with the zero-constant increments
$\Delta\rho=t$ and $\Delta z=-s_c\sum_{n\ge1}(a^{-2})_{n-1}t^n/n$.
The factorial ratios convert normalized bivariate coefficients to the
required coordinate derivatives before composition. This keeps $b_\rho$
distinct from the total characteristic derivative $T_{s_c}b$.

The reduced system has the triangular form

$$
q_z'=f_1(\rho),\qquad q_{zz}'=f_2(\rho),\qquad
J_3'=g(\rho)+k(\rho)q_{zz},
$$

with the geometry-only forcing from the proved cubic identities.
J3 here denotes the changed-minus-reference current; the reference has no
spatial derivatives and its common constant cancels. The N difference on
Sigma is obtained using the same proved surface contraction.

## Two directed integration controls

For a decreasing cell with $h=\rho_1-\rho_0<0$, uniform forcing intervals give

$$q_{zz}(\rho)\in Q_0+[h,0]F_2.$$

Consequently $J_1\in J_0+h(G+K[Q_0+[h,0]F_2])$. This elementary enclosure
retains the nested integral. Freezing qzz at the starting value would omit a
term and is rejected by a polynomial negative control.

The higher-order method uses Taylor coefficients $f_j=f^{(j)}/j!$ at the
starting point and throughout the cell. For either phase-gradient integral,

$$
Q_1\in Q_0+\sum_{j=0}^{n-1}f_j(\rho_0)\frac{h^{j+1}}{j+1}
+[f_n]_{\rm cell}\frac{h^{n+1}}{n+1}.
$$

For J, use the coefficient convolution of $g+kQ$. The starting Q jet is
$Q_0,f_{2,0},f_{2,1}/2,\ldots$. For the final Lagrange coefficient, its
constant is replaced by the **whole-cell tube**, while higher Q coefficients
are bounded by $[f_{2,j-1}]_{\rm cell}/j$. This encloses the required highest
derivative of J over the complete cell. All arithmetic is outward-directed.

Increasing Taylor degree with broad boxes can make bounds worse. Coarse pilot
cells exhibited that overestimation near cutoff transitions. Narrower boxes
resolve it. Every reported interval is derived independently from its own
Taylor remainder; differences between resolutions are never used as errors.

## Pilot and remaining work

The point is the original history's binary incoming center. The pilot uses
the unchanged upstream rho, both $s_c$ values, $m=0$, $|\ell|=\sqrt5$ (enclosing
both its analytic value and original binary label), and original multiplicity
6. All reduced forcing terms are even under $\ell\mapsto-\ell$, so the
angular pair contributes a factor of two. Source signs remain separate.

The per-channel N-bracket difference receives $-\mu/(2\pi)$ and J3 receives
$+\mu/(2\pi)$ for beta. These are coefficient ranges, not complete residual
or epsilon ranges. Their eventual $E^{-3}$ integral is only one term of the
UV tail. The record leaves the higher-order remainder, uniform C4 on I,
physical source error and physical UV component null.

The 1,024-cell, order-8 single-characteristic prototype took about 13 CPU
seconds and narrowed J3's interval to below $10^{-12}$. That timing is an
observed cost, not a proof quantity. The record stores exact dyadic interval
endpoints and can replay without saving large trajectories.

Next, cover the declared incoming interval and bound the remaining UV terms.
A nonzero cubic tail can be integrated as a known contribution with its own
uncertainty; its size must not be confused with that uncertainty. Any such
assembly must preserve the owned subtraction and unchanged physical action.
The record also gives the formal leading term's integral, coefficient divided
by $2E_c^2$ at the original cutoff $E_c=160$. The missing higher-order terms
prevent identifying it with the complete UV tail. No value from this pilot
has been inserted into the physical residual or its error budget.

The Mac controls also verify that the distance integrand is the same
$a^{-2}$ as the background jet owner, that the characteristic interval keeps
its radius through `_binary_arb`, and that both characteristic enclosures
contain the independent h3/n4 numerical transport. The latter is an
implementation check; the Taylor remainder supplies the directed bound.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_cubic_uv_enclosure.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_cubic_uv_enclosure.py --check
```
