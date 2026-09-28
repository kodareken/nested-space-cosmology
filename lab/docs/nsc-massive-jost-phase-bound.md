# Directed exterior subgap Jost phase bound

This is a source-accuracy method for Stage 2 of the approved conversation plan.
It reuses the existing background and massive Riccati equation. It supplies a
finite outer-boundary estimate, not the complete physical source-column error.
Implementation: `src/recursive_horizons/nsc_massive_jost_phase_bound.py`.

## Equation and domain

Here r denotes the original exterior radial coordinate rho, with areal radius
sqrt(1+r^2). On the fixed unit NSC exterior, set

\[
A(r)=1-3[(1+r^2)\arctan(1/r)-r],\qquad
V_1=\ell\sqrt A/\sqrt{1+r^2},\quad V_2=m\sqrt A.
\]

The existing ratio equation is

\[
A z'=i(V_1+iV_2)-2iEz+i(V_1-iV_2)z^2.
\]

In the subgap decaying sector the current is zero, so write \(z=e^{i\theta}\).
The exact phase equation and its outgoing limit are

\[
\theta'=f(r,\theta)=\frac2A(V_1\cos\theta+V_2\sin\theta-E),
\qquad \theta(\infty)=\arcsin(E/m),\quad 0\le E<m.
\]

The principal outgoing branch has positive cosine. Production source rows have
positive E; E=ell=0 is used solely as an exact control. A failed tube/positivity
test does not exclude a physical source; it means this bound has not certified it.

Let \(p(x)+iq(x)\), \(x=1/r\), be the *supplied binary64* ratio polynomial
from `outgoing_ratio`. No recurrence coefficient is assumed exact. For p>0,
\(\phi=\arctan(q/p)\) and

\[
d(r)=\phi'-f(r,\phi),\qquad
\phi'=-x^2\frac{pq_x-qp_x}{p^2+q^2}.
\]

## Continuous residual enclosure

The metric's convergent expansion on 0<=x<1 is

\[
A=1-6\sum_{j\ge0}\frac{(-1)^jx^{2j+1}}{(2j+1)(2j+3)}.
\]

After J terms, every normalized derivative of the omitted series on [0,1/R]
is bounded in absolute value by the corresponding derivative coefficient of
\(6x^{2J+1}/(1-x)\) evaluated at 1/R. Its coefficients are nonnegative.
This deliberately loose positive majorant covers all omitted odd powers and
their derivatives. Arb automatic differentiation then encloses the Nth
normalized derivative of d on the entire compact interval, not a sampled mesh.
With 2J+1>=N, derivatives at zero below N are unaffected by the omitted metric
tail. Taylor's theorem gives

\[
\sup_{r\ge R}|d(r)|\le
\sum_{j=0}^{N-1}|d^{(j)}(0)/j!|R^{-j}
+\sup_{[0,1/R]}|d^{(N)}(x)/N!|R^{-N}=\varepsilon.
\]

Here the Taylor derivatives refer to d as a function of x. Rounded ratio
coefficients and parameter intervals are retained in those derivatives.
There is no assumption of exact low-order cancellation or extrapolated remainder.

## Backward comparison and the boundary at infinity

Choose a phase tube |theta-phi|<=eta. The sine and cosine Lipschitz bound gives
the following uniform lower estimate, with every endpoint directed outward:

\[
f_\theta=\frac2{\sqrt A}
\left(m\cos\theta-\frac{\ell x}{\sqrt{1+x^2}}\sin\theta\right)
\ge\lambda>0.
\]

For e=theta-phi, the mean-value equation is e'=a(r,e)e-d, a>=lambda.
Integrating **backwards** from a finite L with e(L)=0 yields

\[
|e(r)|\le\frac{\varepsilon}{\lambda}
(1-e^{-\lambda(L-r)}).
\]

If epsilon/lambda<eta, the solution stays strictly inside the tube and continues
to R. Two such solutions with increasing terminal radii differ on a fixed
interval by at most a bounded endpoint difference times exp(-lambda(L-r)).
They therefore converge to a unique bounded solution on [R,infinity).

The asymptotic root arcsin(E/m) is checked to lie in the limiting tube. As
A->1 and V1->0, comparison with this constant root has a forcing tending to
zero. The same backwards integral estimate shows theta(r)->arcsin(E/m).
More explicitly, for sufficiently large r the constant root lies inside the
tube. Put b(r)=f(r,theta_infinity), which tends to zero. The bounded solution
satisfies |theta(r)-theta_infinity| <= sup_{t>=r}|b(t)|/lambda, by the same
backwards integral formula with its exponentially decaying terminal term. This
upper tends to zero. The limiting ratio also gives log-amplitude derivative
-sqrt(m^2-E^2), selecting decay rather than growth.

Thus the selected solution is the decaying outgoing branch, not an arbitrary
finite-radius terminal condition. In particular,

\[
|\theta(R)-\phi(R)|\le\varepsilon/\lambda.
\]

Finally add the directed difference between phi(R) and the actual supplied
floating initializer angle. The resulting bound includes ratio evaluation and
angle rounding for that initializer. It is also an upper bound on the complex
unit-ratio difference, since |exp(iu)-exp(iv)|<=|u-v|.

## Coverage and remaining work

The recorder selects original group-14 low-energy inventory rows, both angular
signs, and keeps their exact binary energy labels. It encloses both the binary
channel parameters and their analytic pi/2, +/-sqrt(5) values. Each pilot binds
the inventory, polynomial coefficients, angle, implementation and this derivation.
Changing the outer radius is an independently bounded control, not a convergence
estimate used as a certificate. Directed bounds are serialized as exact dyadics.

This does not yet enclose DOP853 radial transport, the matched horizon collar,
Frobenius sewing, interior transport to rho=1, source covariance arithmetic or
energy-panel interpolation. The producer's actual outer radius must be bound
when applying this to a stored column. No physical N/beta budget component is
filled by this outer-phase pilot. The incoming gate remains OPEN.
