# Conditional analytic crossing of the compatible principal root

A nonzero-velocity crossing of `Delta=0` is locally possible when the exact
crossing compatibility holds. This is a conditional statement about the
included compatible lapse/shift operators; it selects no physical source,
flux constant, crossing position, boundary data or profile.

## Reuse, missing connection and stopping condition

Reuse the fixed polynomial coefficients in `3cabfe9`, the regular-branch
certificate `0c2453f`, and the separately owned flux-form reduction. The
missing connection is whether the singular highest-derivative coefficient
prevents a finite-jet crossing even after its algebraic compatibility is met.
Verify the inverse-coordinate transformation and apply an existing
Briot–Bouquet theorem. Stop with a conditional analytic crossing germ or a
precise failed hypothesis. No coefficient/source producer, profile solver,
exceptional zero-velocity case or global extension is needed.

## Exact compatibility and inverse coordinate

The frozen identities give

\[
\Delta(w)=\Delta_1(w-w_D),\qquad \Delta_1=-dC>0,
\]
\[
(\Delta(w)w')'=R(w,z),\qquad
R(w,z)=A(w)\{\Pi_0-S_\beta z-G(w)\}+dD(w),
\]
where `G'=F`, `G(0)=0`, `D(0)=S_N`, and all coefficients are the already
proved fixed polynomials. The integration constant `Pi0` and finite exact
source constants remain real symbols. In particular, numerical source
approximations are not substituted into R.

At any real crossing position z_c with `w(z_c)=w_D` and finite
`v_c=w'(z_c) != 0`, the necessary condition is

\[
\boxed{R(w_D,z_c)=\Delta_1v_c^2>0.}
\]

Assume this exact condition. Let `x=w-w_D`, `z=Z(x)`, and
`t=Z'=1/w'`. Then `w''=-t'/t³` and the flux equation is exactly

\[
Z'=t,\qquad x t'=t-\frac{R(w_D+x,Z)}{\Delta_1}t^3.
\]

Choose either compatible real nonzero slope symbol v_c, set `t0=1/v_c`,
and center `y=Z-z_c`, `s=t-t0`. The system becomes

\[
x\binom{y'}{s'}=
\binom{x(t_0+s)}
 {t_0+s-\Delta_1^{-1}R(w_D+x,z_c+y)(t_0+s)^3}.
\]

The right side is holomorphic near `(x,y,s)=(0,0,0)` because R is a fixed
real polynomial in w,z. Compatibility makes its constant term zero. Its
Jacobian with respect to `(y,s)` is

\[
J=\begin{pmatrix}
0&0\\
-R_z(w_D,z_c)t_0^3/\Delta_1&-2
\end{pmatrix},\qquad \operatorname{spec}J=\{0,-2\}.
\]

Multiplying the ordinary equation `y'=t0+s` by x introduces no extra
analytic solution: equality for x different from zero extends to x=0 by
analyticity. In particular `Z'(0)=t0 != 0`.

## Imported existence result and reality

Feng Rong, *The Briot–Bouquet systems and the center families for holomorphic
dynamical systems*, Proposition 2.1, p.2, states that a holomorphic system
`x y'=f(x,y)` with `f(0,0)=0` has a unique holomorphic germ with `y(0)=0`
when the state Jacobian has no positive-integer eigenvalue. Zero is allowed.
This precise theorem was checked in the
[author's arXiv v1 PDF](https://arxiv.org/pdf/1307.3823v1#page=2);
the [published article](https://www.sciencedirect.com/science/article/pii/S0001870813002302)
also contains Proposition 2.1. The general convergence theorem is imported,
not rederived here.

Both eigenvalues meet that hypothesis, so there is a unique holomorphic
centered `(Z,t)` germ for the specified compatible real tuple. Complex
conjugation preserves its equation and initial value. Uniqueness therefore
makes it real on the real axis. Since `Z'(0)=t0` is nonzero, the analytic
inverse-function theorem supplies a real-analytic inverse x(z) on a
neighborhood of z_c. Thus

\[
w(z)=w_D+x(z),\qquad w(z_c)=w_D,\quad w'(z_c)=v_c\ne0
\]

is a real-analytic crossing through both sides of the principal root.
This establishes smooth existence as a consequence of analyticity; uniqueness
is asserted within the analytic germs, not as a new theorem for arbitrary
finite-regularity solutions.

The first finite curvature follows directly by differentiating the flux
identity once:

\[
w''(z_c)=\frac{R_w(w_D,z_c)+R_z(w_D,z_c)/v_c}{3\Delta_1}.
\]

All higher derivatives are finite because the germ is analytic. Recover

\[
U(z)=\frac{\Pi_0-S_\beta z-G(w(z))-e w''(z)}{d},
\]
using the already certified `d != 0`. This is analytic. Its substitution
satisfies the momentum constraint identically; `d E_N` is the negative of
the flux-equation residual, including at the root. No division by A or Delta
is required in this reconstruction. The compatible metric Taylor family
therefore has finite jets locally; no field evolution is inferred.

## Verification and limits

The small symbolic module verifies the coordinate change, centered Jacobian,
characteristic polynomial, reconstruction residual, and first-curvature
relation exactly. It imports the existing d,C and root intervals and checks
`-d*C>0`; the root is an exact named value enclosed by the receipt, never
replaced by a floating midpoint in the crossing condition. Record replay
checks input hashes and calls no scientific coefficient generator.

For any symbolic real tuple with `R(w_D,z_c)>0`, the two compatible signs
`v_c=±sqrt(R(w_D,z_c)/Delta1)` each define their own analytic crossing germ.
This does not establish that a particular earlier trajectory reaches the
root or meets compatibility, nor select one of the signs. The zero-velocity
and zero-R cases are outside this result. No global solution, boundary match,
parent preparation, pressure equation, physical source accuracy or metric
initial-value problem is claimed.

```sh
PYTHONPATH=src python3 -m recursive_horizons.nsc_incoming_surface_crossing --write
PYTHONPATH=src python3 -m recursive_horizons.nsc_incoming_surface_crossing --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_crossing.py
```
